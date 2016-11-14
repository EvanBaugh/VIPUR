#!/usr/bin/env python
# :noTabs=true:

"""
run methods specifically for submitting jobs in the SLURM queuing system
use the identical framework for PBS e.g. collecting tasks first during
preprocessing, then running them together as a batch, monitoring the queue etc.
set arbitrary limits for now...

so single_relax = False by default for this implementation
"""

################################################################################
# IMPORT

# common modules
import subprocess
import time

# bigger modules

# custom modules
from vipur_settings import SLURM_USER , SLURM_QUEUE_QUOTA , SLURM_QUEUE_MONITOR_DELAY , SLURM_BASH_SCRIPT , SLURM_JOB_OPTIONS
from helper_methods import run_local_commandline , create_executable_str

from pre_processing import *
from run_methods import determine_check_successful_function
from rosetta_feature_generation import remove_intermediate_ddg_monomer_files
from post_processing import *

################################################################################
# SLURM RUN METHODS

def run_VIPUR_SLURM( pdb_filename = '' , variants_filename = '' ,
        out_path = '' , write_numbering_map = True ,
        single_relax = False , delete_intermediate_relax_files = True ,
        demo = False , rerun_preprocessing = False ):
    # the following should probably be a separate method...

    # for the example input
    if demo:
        pdb_filename = PATH_TO_VIPUR + '/example_input/2C35.pdb'
        variants_filename = PATH_TO_VIPUR + '/example_input/2C35.txt'

        out_path = PATH_TO_VIPUR + '/example_output'

    # alternatively, run on an entire directory
    if not pdb_filename and not variants_filename:
        # current directory
        print 'no input provided, assuming you want to run on every (.pdb,.txt) file pair found in the current directory'
        pdb_filename = os.getcwd()

    if os.path.isdir( pdb_filename ) and not variants_filename:
        # assume variants_filename from pdb_filename
        variants_filename = '.txt'

    if os.path.isdir( pdb_filename ) and variants_filename[0] == '.':
        # look for file extension
        # instead, run on the directory
        if not out_path:
            out_path = os.path.abspath( pdb_filename )
        
        fa_filenames = [(out_path +'/')*bool( out_path ) + i for i in os.listdir( pdb_filename ) if get_file_extension( i ) == 'fa']
        fa_filenames = [[i , get_root_filename( i ) + variants_filename] for i in fa_filenames if os.path.isfile( get_root_filename( i ) + variants_filename ) and not os.path.isfile( get_root_filename( i ) + '.pdb' )]

        print 'running VIPUR on all (.pdb,' + variants_filename + ') file pairs found in ' + pdb_filename
        # find .pdb files
        pdb_filenames = [(out_path +'/')*bool( out_path ) + i for i in os.listdir( pdb_filename ) if get_file_extension( i ) == 'pdb']

        # look for pairs
        pdb_filenames = [[i , get_root_filename( i ) + variants_filename] for i in pdb_filenames if os.path.isfile( get_root_filename( i ) + variants_filename )]

        print str( len( pdb_filenames ) ) + ' pairs found'
        print str( len( fa_filenames ) ) + ' pairs found for sequence only mode'

        # go there...
#        os.chdir( pdb_filename )

        if not pdb_filenames:
            if not fa_filenames:
                raise IOError( '!!! no (.pdb,' + variants_filename + ') file pairs found in ' + pdb_filename + '!!?!\nAND no (.fa,' + variants_filename + ') file pairs were found...' )
            else:
                print '...only (.fa,' + variants_filename + ') file pairs were found, running in sequence only mode'

    else:
        # file extension etc.
        file_extension = get_file_extension( pdb_filename )
        root_filename = get_root_filename( pdb_filename )

        # normal execution, generalize by turning into list
        pdb_filenames = []
        fa_filenames = []
        if file_extension == 'pdb':
            pdb_filenames = [[(out_path +'/')*bool( out_path ) + pdb_filename , (out_path +'/')*bool( out_path ) + variants_filename]]
        else:
            fa_filenames = [[]]


    # combine all "filenames" to run into unified framework
    target_proteins = []#None]*(len( pdb_filenames ) + len( fa_filenames ))
    for i in pdb_filenames:
        this_out_path = get_root_filename( i[0] ) +'_VIPUR'    # directory to create
        target_proteins.append( i + [False , this_out_path] )
    for i in fa_filenames:
        this_out_path = get_root_filename( i[0] ) +'_VIPUR'    # directory to create
        target_proteins.append( i + [True , this_out_path] )

    # setup environment variables BEFORE pre processing
    # not needed with current SLURM setup...

    # pre processing
    task_summaries = []
    for i in target_proteins:
        # guess what the task summary filename 'would' be, if it exists, keep going...
        task_summary_filename = i[3]*bool( i[3] ) +'/'+ get_root_filename( i[0] ).split( '/' )[-1] + '.task_summary'
        if os.path.isfile( task_summary_filename ) and not rerun_preprocessing:
            print 'hmmm, ' + i[0] + ' seems to have run preprocessing already, skipping now'
            #continue    # skip this one, do not add to list of tasks...?
            # actually, skip running pre-processing BUT DO add it to the list of tasks
            
            # is this actually working?
        else:
            task_summary_filename = run_preprocessing( i[0] , i[1] ,
                sequence_only = i[2] , out_path = i[3] ,
                task_summary_filename = task_summary_filename ,
                write_numbering_map = write_numbering_map , single_relax = single_relax )


        # modify for SLURM script
        task_summary = load_task_summary( task_summary_filename )
        
#        task_summary['filenames']['slurm_script_filename'] = 'slurm_' + get_root_filename( i[0] ) + '.sh'
        task_summary['filenames']['slurm_script_filename'] = out_path + '/slurm_script_this_batch.sh'
        task_summary['filenames']['slurm_output_filename'] = out_path + '/slurm_output_batch.out'
        task_summary['filenames']['slurm_error_filename'] = out_path + '/slurm_error_batch.err'
        # ...awkward...they all have individual task summarization of the same script...but nowhere else to put it...
        
        for j in xrange( len( task_summary['commands'] ) ):
            slurm_options = {}

            command = task_summary['commands'][j]['command']

            # add for relax
            if task_summary['commands'][j]['feature'].replace( '_native' , '' ) == 'relax' and not 'rescore' in task_summary['commands'][j]['feature']:
                command = command.replace( '.linuxgccrelease' , '.mpi.linuxgccrelease' )
                command = 'mpiexec -n 40 ' + command
                command += ' -jd2:mpi_file_buf_job_distributor false'
                command += ' -run:multiple_processes_writing_to_one_directory'
                
            slurm_options.update( SLURM_JOB_OPTIONS )

            # special...
            if task_summary['commands'][j]['feature'] == 'psiblast' and not 'num_threads' in task_summary['commands'][j]['command']:
                task_summary['commands'][j]['command'] += ' -num_threads 40'
            
            # modify the task summary
            task_summary['commands'][j]['command'] = command
            
            # MUST still do ddg_monomer on single process...
            if 'ddg_monomer' in task_summary['commands'][j]['feature'] or 'rescore' in task_summary['commands'][j]['feature']:
                if 'rescore' in task_summary['commands'][j]['feature']:
                    # sanity check
                    if not 'variant' in task_summary['commands'][j].keys():
                        raise Exception( 'rescore command without the variant information...!?' )
                    
                    # need variant in the script, otherwise overwrite :(
                    script_filename = i[3] + '/'*bool( i[3] ) + get_root_filename( i[0] ).split( '/' )[-1] +'.'+ task_summary['commands'][j]['feature'] +'_'+ task_summary['commands'][j]['variant'] + '.slurm_script.sh'
                else:
                    # ddg monomer is "per protein", no need for more detail
                    script_filename = i[3] + '/'*bool( i[3] ) + get_root_filename( i[0] ).split( '/' )[-1] +'.'+ task_summary['commands'][j]['feature'] + '.slurm_script.sh'
                task_summary['commands'][j]['script_filename'] = script_filename

                # only write ONE submission script per batch = run of VIPUR           
                f = open( script_filename , 'w' )
                f.write( SLURM_BASH_SCRIPT( command ) )
                f.close()
            
                # use the script filename as the source for any log files
                # control the output and error paths
                for k in slurm_options.keys():
                    if '__call__' in dir( slurm_options[k] ):
                        slurm_options[k] = slurm_options[k]( script_filename )

                slurm_options['N'] = '1'
                slurm_options['n'] = '1'

                # also generate the pbs call? might as well, keep it simple...
                # srun or sbatch?
                task_summary['commands'][j]['sbatch_command'] = create_executable_str( 'sbatch' , [script_filename] , slurm_options )
            
        # rewrite the task summary
        write_task_summary( task_summary , task_summary_filename )

        task_summaries.append( task_summary )#_filename )


    # run them all
    run_VIPUR_task_summaries_SLURM( task_summaries , single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )


    # post processing
    # this look identical!!! :)
    for i in xrange( len( task_summaries ) ):
        # always okay to rerun post processing...should not make any difference
        sequence_only = target_proteins[i][2]
        task_summaries[i] = run_postprocessing( task_summaries[i] , sequence_only = sequence_only )

    return task_summaries

# different from serial method
# collect up all jobs at once, fire off serially (lol) to queue
# NOPE collect up all jobs at once - WRITE ONE SUBMISSION SCRIPT, and fire off this single job...I guess give progressive updates on status...?
# need to be cognizant of queue status, rescore must occur later
def run_VIPUR_task_summaries_SLURM( task_summaries , single_relax = False , delete_intermediate_relax_files = True , unsupported_many_jobs_version = False ):
    # queue command can be derived from task itself
    # this method will run...until all jobs are complete
    
    # reload? really how this should be done...
    for i in xrange( len( task_summaries ) ):
        if isinstance( task_summaries[i] , str ) and os.path.isfile( task_summaries[i] ):
            task_summaries[i] = load_task_summary( task_summaries[i] )
    
    # make a list of ALL jobs instead of per protein tasks - all commands have "cd" if there is an "out_path"
    # as a reference, need to track attempts etc. directly, rewrite task summaries frequently
    non_rescore_tasks = sum( [[(i , j) for j in xrange( len( task_summaries[i]['commands'] ) )] for i in xrange( len( task_summaries ) )] , [] )
    rescore_tasks = [i for i in non_rescore_tasks if 'rescore' in task_summaries[i[0]]['commands'][i[1]]['feature']]
    # only rescore once ALL are complete
    non_rescore_tasks = [i for i in non_rescore_tasks if not i in rescore_tasks]


    # actually...just give a whole different setup here
    
    if unsupported_many_jobs_version:
        run_VIPUR_tasks_SLURM( task_summaries , non_rescore_tasks )
    
        # should just be give dependency on the matching relax jobs
        # --dependency=afterok:<relax_script_filename.sh>
        # run the rescore tasks
        run_VIPUR_tasks_SLURM( task_summaries , rescore_tasks )
    else:
        # the only supported usage
        # gather up all the jobs, create a single batch script
        # could make separately, and submit a collective script to submit these...
        # but for now, going to just run crudely combined with a single sbatch
        # ...how to collect for salloc? still a mystery to me...
        
        # gather the commands
        # ...until its working, still just do them in series
        # MUST do ddg_monomer on single processor -N1 -n1 ...wtf...
        ddg_monomer_tasks = [i for i in non_rescore_tasks if 'ddg_monomer' in task_summaries[i[0]]['commands'][i[1]]['feature']]
        non_rescore_tasks = [i for i in non_rescore_tasks if not i in ddg_monomer_tasks]
        
        # do ddg_monomer in betwee...cause why not?
        # as the many-job batch variant
        run_VIPUR_tasks_SLURM( task_summaries , ddg_monomer_tasks )
        
        run_VIPUR_tasks_in_batch_SLURM( task_summaries , non_rescore_tasks )       

        # actually, do this with the ddg_monomer stuff
        run_VIPUR_tasks_SLURM( task_summaries , rescore_tasks )
    
    # return anything?
    # task summaries should be updated with all the necessary files...


################################################################################
# three obvious ways to proceed for me
# one big sbatch script, submitted requesting many resources
# tested, still need ddg_monomer and rescore as separate jobs...

# one big sbatch script - but run with the "parSlurm" approach, requires non-multithreading and no mpi

# many small jobs run separately with sbatch
# need to explicitly test, but should be identical to PDB in form

################################################################################
# ONE BIG SBATCH SCRIPT

def run_VIPUR_tasks_in_batch_SLURM( task_summaries , task_list , max_slurm_tries = 2 , ddg_monomer_cleanup = True , single_relax = True ):
    # also setup to do start-stop
    completed = [i for i in task_list if 'run' in task_summaries[i[0]]['commands'][i[1]] and 'success' in task_summaries[i[0]]['commands'][i[1]]['run']]

    attempt = 1
    while not len( completed ) == len( task_list ):
        # do not worry about the queue in this mode
        # assume you will be able to submit etc.
        # not such thing as "running or queued" either, just run one batch at a time
            
        # gather all the commands into a single script
        # choose the next job
        jobs_to_run = [i for i in task_list if
            not i in completed and
            not ( 'run' in task_summaries[i[0]]['commands'][i[1]] and
                ('success' in task_summaries[i[0]]['commands'][i[1]]['run'] or
                'failure' in task_summaries[i[0]]['commands'][i[1]]['run']) )
            ]
        print str( len( jobs_to_run ) ) + ' processes still need to finish'
    
        # write the script
        master_script_text = '\n\n'.join( [task_summaries[i[0]]['commands'][i[1]]['command'] for i in jobs_to_run] )
        # test without relax processes
#        master_script_text = '\n\n'.join( [task_summaries[i[0]]['commands'][i[1]]['command'] for i in jobs_to_run if not 'relax' in task_summaries[i[0]]['commands'][i[1]]['command']] )
        master_script_text = SLURM_BASH_SCRIPT( master_script_text )
        
        # check if they have different names!?...wait...they will...
        slurm_script_filename = task_summaries[0]['filenames']['slurm_script_filename']
        slurm_output_filename = task_summaries[0]['filenames']['slurm_output_filename']
        slurm_error_filename = task_summaries[0]['filenames']['slurm_error_filename']
        
        slurm_script_filename = slurm_script_filename.replace( '.sh' , '_'+ str( attempt ) + '.sh' )
        slurm_output_filename = slurm_output_filename.replace( '.out' , '_'+ str( attempt ) + '.out' )
        slurm_error_filename = slurm_error_filename.replace( '.err' , '_'+ str( attempt ) + '.err' )
        # can just use the first one now...

        f = open( slurm_script_filename , 'w' )
        f.write( master_script_text )
        f.close()
        
        # save a copy of this script for reference?
        # successive runs with overwrite the file...
        
        # debug
#        raw_input( 'everything okay?' )
    
        # submit sbatch
        # simple for now...
        command = 'sbatch -n 40'
        if slurm_output_filename:
            command += ' -o ' + slurm_output_filename
        if slurm_error_filename:
            command += ' -e ' + slurm_error_filename
        command += ' ' + slurm_script_filename
        batch_job_id = run_local_commandline( command , collect_stdout = True )
#        batch_job_id = run_local_commandline( 'sbatch -n 40 ' + slurm_script_filename , collect_stdout = True )
        # srun or sbatch?
        batch_job_id = batch_job_id.strip().split( ' ' )[-1]
        print 'submitted ' + batch_job_id

        
        # monitor the job until it is complete

        # pause...
        batch_complete = False
        while not batch_complete:
            queue_status = get_slurm_queue_status( only_job_status = True )
#            batch_complete = bool( [i for i in queue_status if i[0] == batch_job_id] )
            batch_complete = not batch_job_id in queue_status.keys()
            # could be an immediate failure...but don't want to linger here anyway in that case

            # debug
#            print queue_status
#            print queue_status.keys()
#            print batch_complete , batch_job_id , batch_job_id in queue_status.keys()
            for i in queue_status.keys():
                print i + '\t' + queue_status[i]

            # can be sure it doesn't need to wait if empty
            if queue_status:
                print 'waiting ' + str( SLURM_QUEUE_MONITOR_DELAY ) +'s...'
                time.sleep( SLURM_QUEUE_MONITOR_DELAY )


        # evaluate if it ran successfully
        for job_pair in jobs_to_run:
            command_dict = task_summaries[job_pair[0]]['commands'][job_pair[1]]

            check_successful = determine_check_successful_function( command_dict , single_relax = single_relax )

            success = check_successful( command_dict )

            failure_summary = ''
            if isinstance( success , bool ):
                complete = success
            elif len( success ) > 1 and isinstance( success[0] , bool ):
                complete = success[0]
                failure_summary += ' '+ ';'.join( [str( j ) for j in success[1:]] ) +' '
 
            # track the number of attempts?
            # try until failure - how many times?
            tries = 0
            if 'run' in command_dict.keys() and command_dict['run'] and not 'success' in command_dict['run'] and not 'failure' in command_dict['run']:
                tries = int( command_dict['run'] )
            tries += 1

            job_task = task_summaries[job_pair[0]]['root_filename'].split( '/' )[-1]
            job_feature = task_summaries[job_pair[0]]['commands'][job_pair[1]]['feature']
            job_variant = ''
            if 'variant' in task_summaries[job_pair[0]]['commands'][job_pair[1]].keys():
                job_variant = task_summaries[job_pair[0]]['commands'][job_pair[1]]['variant']
            this_job_description = job_task +' '+ job_feature + (' ' + job_variant)*bool( job_variant )

            if tries >= max_slurm_tries:
                # its a failure
                print this_job_description + ' completed successfully'*complete + (' failed with ' + str( tries ) + ' attempts')*(not complete)
                failure_summary = 'success'*complete + (str( tries ) +' tries;failure ' + failure_summary)*(not complete)
                completed.append( job_pair )
            elif complete:
                print this_job_description + ' simply completed successfully'
                failure_summary = 'success' #+ str( tries ) + ' tries'
                completed.append( job_pair )
            else:
                # record the number of tries
                print this_job_description + ' completed' + ' successfully'*complete
                failure_summary = str( tries )
                completed.append( job_pair )
                
            # update the record
            task_summaries[job_pair[0]]['commands'][job_pair[1]]['run'] = failure_summary
            
            # no need to be here anymore
            # optionally cleanup
#            if ddg_monomer_cleanup and command_dict['feature'] == 'ddg_monomer':#'ddg' in i['output_filename']:
#                print 'ddg_monomer writes useless output files, deleting these now...'
#                remove_intermediate_ddg_monomer_files()

        # update task_summaries e.g. write them!
        # modified: so the task summary records its own name...bah!
        for i in task_summaries:
            if not 'task_summary_filename' in i['filenames'].keys():
                raise NotImplementedError( 'should input the task summary filename (not the summary itself)...' )
            else:
                # write it out
                print 'updating: ' + i['filenames']['task_summary_filename']
                write_task_summary( i , i['filenames']['task_summary_filename'] )

        # debug
#        print attempt
#        print len( completed ) , len( task_list )
#        raw_input( 'start the next round?' )

        # need to run another batch?
        attempt += 1

    # rerun as necessary
    # implicit above...?
    
    # cleanup at the end
    # need to merge relax output...?


    # should be outside the outermost loop...which is unclear right now...        
    # write one last time?
#    for i in task_summaries:
#        if not 'task_summary_filename' in i['filenames'].keys():
#            raise NotImplementedError( 'should input the task summary filename (not the summary itself)...' )
#        else:
            # write it out
#            write_task_summary( i , i['filenames']['task_summary_filename'] )


################################################################################
# SUBMIT MANY INDIVIDUAL JOBS

# ex. need to distinguish serial vs parallel queues etc.

# submit jobs until complete or too many attempts
def run_VIPUR_tasks_SLURM( task_summaries , task_list , max_pbs_tries = 2 , ddg_monomer_cleanup = True , single_relax = False , delete_intermediate_relax_files = False ):
    # run the non_rescore tasks
    completed = [i for i in task_list if 'run' in task_summaries[i[0]]['commands'][i[1]] and 'success' in task_summaries[i[0]]['commands'][i[1]]['run']]
    # should running_or_queued be saved? written to file?
    running_or_queued = {}
    rounds = 0
    while not len( completed ) == len( task_list ):
        rounds += 1
        print '\n\nQUEUE MONITOR ROUND ' + str( rounds )
        
        # debug
#        print running_or_queued
    
        # check queue status
        queue_status = get_slurm_queue_status( only_job_status = True )

        # update "running_or_queued" list (?)
        # err, no, does not have information on which job it is...:(
        #for i in queue_status.keys():
        #    if queue_status[i] in ['R' , 'Q']:


        # used to be after submission, not occurs first

        # debug, need to know
        running_jobs = len( [i for i in queue_status.values() if i in ['R']] )
        if running_jobs:
            print str( running_jobs ) + ' are still running...'
        
        # assess outcome of completed jobs
#        still_running = 0
        # need to add the jobs that completed, removed themselves from the queue in SLURM
#        print queue_status.keys() + [j for j in running_or_queued.keys() if not j in queue_status.keys()]
        for job_id in queue_status.keys() + [j for j in running_or_queued.keys() if not j in queue_status.keys()]:
            # debug
#            if job_id in queue_status.keys():
#                print '\t'+ job_id , queue_status[job_id] , job_id in running_or_queued.keys()
#            else:
##                print '\t'+ job_id , None , job_id in running_or_queued.keys()
#                print '\t'+ job_id , job_id in running_or_queued.keys()
        
            if (not job_id in queue_status.keys()) or (queue_status[job_id] == 'C' and job_id in running_or_queued.keys()):
                task_id = running_or_queued[job_id][0]
                command_index = running_or_queued[job_id][1]
                command_dict = task_summaries[task_id]['commands'][command_index]

                check_successful = determine_check_successful_function( command_dict , single_relax = single_relax )

                success = check_successful( command_dict )

                failure_summary = ''
                if isinstance( success , bool ):
                    complete = success
                elif len( success ) > 1 and isinstance( success[0] , bool ):
                    complete = success[0]
                    failure_summary += ' '+ ';'.join( [str( j ) for j in success[1:]] ) +' '
 
                # track the number of attempts?
                # try until failure - how many times?
                tries = 0
                if 'run' in command_dict.keys() and command_dict['run'] and not 'success' in command_dict['run'] and not 'failure' in command_dict['run']:
                    tries = int( command_dict['run'] )
                tries += 1
                
                if tries >= max_pbs_tries:
                    # its a failure
                    print job_id + ' completed successfully'*complete + (' failed with ' + str( tries ) + ' attempts')*(not complete)
                    failure_summary = 'success'*complete + (str( tries ) +' tries;failure ' + failure_summary)*(not complete)
                elif complete:
                    print job_id + ' simply completed successfully'
                    failure_summary = 'success' #+ str( tries ) + ' tries'
                else:
                    # record the number of tries
                    print job_id + ' completed' + ' successfully'*complete
                    failure_summary = str( tries )
                
                # update the record
                task_summaries[task_id]['commands'][command_index]['run'] = failure_summary
            
                # optionally cleanup
                if ddg_monomer_cleanup and command_dict['feature'] == 'ddg_monomer':#'ddg' in i['output_filename']:
                    print 'ddg_monomer writes useless output files, deleting these now...'
                    remove_intermediate_ddg_monomer_files()

                # jobs that have since been completed - consider them complete?
                completed.append( running_or_queued[job_id] )
                del running_or_queued[job_id]
                
                # write out "completed"? or "running_or_queued"?

#            else:
#                still_running += 1
#        print str( still_running) + ' jobs still running (or queued)...'

        # update task_summaries e.g. write them!
        # modified: so the task summary records its own name...bah!
        for i in task_summaries:
            if not 'task_summary_filename' in i['filenames'].keys():
                raise NotImplementedError( 'should input the task summary filename (not the summary itself)...' )
            else:
                # write it out
                print 'updating: ' + i['filenames']['task_summary_filename']
                write_task_summary( i , i['filenames']['task_summary_filename'] )


        # used to be first, submit jobs them check complete
        # but slurm removes jobs from the list

        queue_space_occupied = len( [i for i in queue_status.values() if not i in ['C' , 'PD']] )    # ignore "C"ompleted jobs, "R"unning job quota are not set by us...
        # if your queue system does not have a separate "R"un quota, remove 'R' from the above!
        available_space = SLURM_QUEUE_QUOTA - queue_space_occupied

        
        # launch next jobs in available slots
        if available_space:
            print str( queue_space_occupied ) + ' jobs queued or running, could submit up to ' + str( available_space ) + ' more'
            # choose the next job
            jobs_to_run = [i for i in task_list if
                not i in completed and
                not i in running_or_queued.values() and
                not ( 'run' in task_summaries[i[0]]['commands'][i[1]] and
                    ('success' in task_summaries[i[0]]['commands'][i[1]]['run'] or
                    'failure' in task_summaries[i[0]]['commands'][i[1]]['run']) )
                ]
            print str( len( jobs_to_run ) ) + ' jobs still need to finish (after the currently running jobs complete)'
            
            # only the next few
            for i in jobs_to_run[:available_space]:
                command_dict = task_summaries[i[0]]['commands'][i[1]]
            
                # write scripts as part of pre processing?...yeah...
                # write the command to a script
                #script_filename = command_dict['out_path'] +'/'*bool( command_dict['out_path'] )+
                #script_filename = command_dict['script_filename']
            
                # submission is specific to the job
                slurm_command = ''
            
                # if its a rescore and relax jobs were separated, need to recombine them!
                if 'rescore' in command_dict['feature']:
                    # combine the individual relax runs
                    #relax_commands = [i for i in task_summary['commands'] if i['feature'].replace( '_native' , '' ) == 'relax']
                    #silent_filenames = [j['output_filename'] for j in relax_commands if j['variant'] == i['variant'] and 'run' in j.keys() and j['run'] == 'success']
                    silent_filenames = [j['output_filename'] for j in task_summaries[i[0]]['commands'] if
                        j['feature'].replace( '_native' , '' ) == 'relax' and
                        j['variant'] == command_dict['variant'] and
                        'run' in j.keys() and
                        'success' in j['run']
                        ]
                    # actually need to identify the combined_silent_filename, be sure the relax files have not already been merged
                    # which variant
                    target_variant = [j for j in task_summaries[i[0]]['variants'].keys() if j.split( '_' )[-1] == command_dict['variant'] and j.split( '_' )[0] in command_dict['command']]
                    if not target_variant:
                        # its native
                        combined_silent_filename = task_summaries[i[0]]['other']['combined_native_silent_filename']
                        combined_score_filename = task_summaries[i[0]]['other']['combined_native_score_filename']
                    elif len( target_variant ) > 1:
                        raise Exception( '??? found more than on matching variant ???\n' + ', '.join( target_variant ) )
                    else:
                        # found it
                        combined_silent_filename = task_summaries[i[0]]['variants'][target_variant[0]]['combined_silent_filename']
                        combined_score_filename = task_summaries[i[0]]['variants'][target_variant[0]]['combined_score_filename']

                    #if not single_relax:    # AND post processing has not already be run...scan for the combined silent file
                    if not single_relax and not os.path.isfile( combined_silent_filename ):
                        if not len( silent_filenames ) == ROSETTA_RELAX_OPTIONS['nstruct']:
                            raise Exception( '??? somehow the matching relax run(s) has failed ???\n' + str( i ) )
                        score_filenames = [j.replace( '.silent' , '.sc' ) for j in silent_filenames]

                        merge_rosetta_relax_output( silent_filenames , combined_silent_filename , score_filenames , combined_score_filename , delete_old_files = delete_intermediate_relax_files )
                        # rescore already knows the proper filename
                    else:
                        # just a single match for each
                        # output filename should be correct as is :)
                        None

                # do this above, in initial processing
#                elif 'ddg_monomer' in command_dict['feature']:
                    # must do it this way for now
                    # write the run script
#                    ddg_monomer_script_filename = task_summary['filenames']['slurm_script_filename'].replace( 'slurm_script_this_batch.sh' , 'run_ddg_momomer_script.sh' )
#                    f = open( ddg_monomer_script_filename , 'w' )
#                    f.write( command_dict['command'] )
#                    f.close()

            
                # submit this script using a queue command
                # srun or sbatch?
                slurm_command = command_dict['sbatch_command']    # SHOULD already have an abspath to the script
                new_job_id = run_local_commandline( slurm_command , collect_stdout = True )
#                new_job_id = new_job_id[:new_job_id.find( ' ' )]
                new_job_id = new_job_id.strip().split( ' ' )[-1]
                print 'submitted ' + new_job_id
                
                # save the job id
                # assume its queue
                running_or_queued[new_job_id] = i
#                print 'added ' + new_job_id

        else:
            print 'no new \"positions\" are available'

        # OKAY, move the "updating" to just after the status check
        # problem with ddg_monomer, runs so fast...
        # make a specific exception:
        # ...move the pause to here
        # prevents odd behaviour...um...sorta? maybe not
#        if 'ddg_monomer' in command_dict['feature']:


        # used to be where the updates were compared

        
        # pause...
        if queue_space_occupied or jobs_to_run:
            print 'waiting ' + str( SLURM_QUEUE_MONITOR_DELAY ) +'s...'
            time.sleep( SLURM_QUEUE_MONITOR_DELAY )

    # return anything?
    # write one last time?
    for i in task_summaries:
        if not 'task_summary_filename' in i['filenames'].keys():
            raise NotImplementedError( 'should input the task summary filename (not the summary itself)...' )
        else:
            # write it out
            write_task_summary( i , i['filenames']['task_summary_filename'] )

################################################################################
# PBS INTERACTION METHODS

# ASSUMES it is the ONLY application submitting jobs as this user...
# simple parse/scan the queue, use qstat -u 
def get_slurm_queue_status( user = SLURM_USER , header_lines = 1 , trailer_lines = 0 , only_job_status = True ):
    # header_lines = 2 for FULL queue, = 5 for USER queue ("-u")
    command = 'squeue'
    if user:
        command += ' -u ' + user

#    queue_info = subprocess.Popen( command.split( ' ' ) , stdout = subprocess.PIPE , stdin = subprocess.PIPE , stderr = subprocess.STDOUT ).communicate()[0]
    queue_info = run_local_commandline( command , collect_stdout = True )

    # debug
#    print queue_info

    # simple parsing
    queue_info = queue_info.split( '\n' )
    if trailer_lines:
        queue_info = queue_info[header_lines:-1*trailer_lines]
    elif header_lines:
        queue_info = queue_info[header_lines:]
    queue_info = [[j for j in i.split( ' ' ) if j.strip()] for i in queue_info if i.strip()]

    # debug
#    print queue_info

    # optionally only report the job statuses
    if only_job_status:
#        queue_info = [(i[0][:i[0].find( '.' )] , i[-2]) for i in queue_info]
        queue_info = [(i[0] , i[4]) for i in queue_info]
        # make into a dict? job ids should be unique...
        queue_info = dict( queue_info )

    # debug
#    print queue_info

    return queue_info

# primarily concerned with job status only
# need 1st and 2nd to last columns

# simple wrapper for running
# grab the output job ID
def run_slurm_job( script_filename , slurm_run_command = 'sbatch' , output_filename = 'temp_slurm.sh' ):    
    # optionally write a file
    if not os.path.isfile( script_filename ) and isinstance( script_filename , str ):
#        if os.path.isfile( output_filename ):    # should really check if it exists first...
        f = open( output_filename , 'w' )
        f.write( script_filename )
        f.close()
        
        script_filename = output_filename
    
    # submit it, grab the name :)
#    job_id = run_and_record( slurm_run_command +' '+ script_filename , )
    job_id = run_local_commandline( slurm_run_command +' '+ script_filename , collect_stdout = True )
    job_id = job_id[0].split( ' ' )[-1]
    
    return job_id


