#!/usr/bin/env python
# :noTabs=true:

"""
run methods specifically for submitting jobs in the PBS queuing system
"""

################################################################################
# IMPORT

# common modules
import subprocess
import time

# bigger modules

# custom modules
from vipur_settings import PBS_USER , PBS_QUEUE_QUOTA , PBS_QUEUE_MONITOR_DELAY , PBS_SERIAL_JOB_OPTIONS , PBS_PARALLEL_JOB_OPTIONS
from helper_methods import run_local_commandline , create_executable_str

from pre_processing import *
from post_processing import *

################################################################################
# PBS RUN METHODS

def run_VIPUR_PBS( pdb_filename = '' , variants_filename = '' ,
        out_path = '' , write_numbering_map = True ,
        single_relax = True , delete_intermediate_relax_files = True ,
        demo = False , rerun_preprocessing = False ):
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
#            print out_path
        
        fa_filenames = [(out_path +'/')*bool( out_path ) + i for i in os.listdir( pdb_filename ) if get_file_extension( i ) == 'fa']
        fa_filenames = [[i , get_root_filename( i ) + variants_filename] for i in fa_filenames if os.path.isfile( get_root_filename( i ) + variants_filename ) and not os.path.isfile( get_root_filename( i ) + '.pdb' )]

        print 'running VIPUR on all (.pdb,' + variants_filename + ') file pairs found in ' + pdb_filename
        # find .pdb files
        pdb_filenames = [(out_path +'/')*bool( out_path ) + i for i in os.listdir( pdb_filename ) if get_file_extension( i ) == 'pdb']

        # look for pairs
        pdb_filenames = [[i , get_root_filename( i ) + variants_filename] for i in pdb_filenames if os.path.isfile( get_root_filename( i ) + variants_filename )]
#        print [i for i in pdb_filenames if os.path.isfile( pdb_filename +'/'+ get_root_filename( i ) + variants_filename )]

        print str( len( pdb_filenames ) ) + ' pairs found'
        print str( len( fa_filenames ) ) + ' pairs found (for sequence only)'

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


    # pre processing
    task_summaries = []
    for i in target_proteins:
        # guess what the task summary filename 'would' be, if it exists, keep going...
        task_summary_filename = i[3]*bool( i[3] ) +'/'+ get_root_filename( i[0] ).split( '/' )[-1] + '.task_summary'
        if os.path.isfile( task_summary_filename ) and not rerun_preprocessing:
            print 'hmmm, ' + i[0] + ' seems to have run preprocessing already, skipping now'
            #continue    # skip this one, do not add to list of tasks...?
            # actually, skip running pre-processing BUT DO add it to the list of tasks
        else:
            task_summary_filename = run_preprocessing( i[0] , i[1] ,
                sequence_only = i[2] , out_path = i[3] ,
                task_summary_filename = task_summary_filename ,
                write_numbering_map = write_numbering_map , single_relax = single_relax )


        # modify for PBS script
        task_summary = load_task_summary( task_summary_filename )
        for j in xrange( task_summary['commands'] ):
            pbs_options = {}

            command = task_summary['commands'][j]

            # add for relax            
            if task_summary['commands'][j]['feature'].replace( '_native' , '' ) == 'relax' and not 'rescore' in task_summary['commands'][j]['feature']:
                command = command.replace( '.linuxgccrelease' , '.mpi.linuxgccrelease' )
                command = 'module load mvapich2/gnu/1.8.1; /share/apps/mvapich2/1.8.1/gnu/bin/mpiexec -n 36 ' + command
                command += ' -jd2:mpi_file_buf_job_distributor false'
                command += ' -run:multiple_processes_writing_to_one_directory'
                
                # also use the parallel options
                pbs_options.update( PBS_PARALLEL_JOB_OPTIONS )
            else:
                pbs_options.update( PBS_SERIAL_JOB_OPTIONS )

            # put "cd" in front
            command = ('cd '+ i[3] +';')*bool( i[3] ) + command
            
            # modify the task summary
            task_summary['commands'][j] = command
            
            
            # actually write the script...
            # don't worry about optional #PBS header info
            script_filename = i[3] + '/'*bool( i[3] ) + get_root_filename( i[0] ).split( '/' )[-1] +'.'+ task_summary['commands'][j]['feature'] + '.pbs_script.sh'
            task_summary['commands'][j]['script_filename'] = script_filename
            
            f = open( script_filename , 'w' )
            f.write( command )
            f.close()
            
            # use the script filename as the source for any log files
            # control the output and error paths
            for i in pbs_options.keys():
                if '__call__' in dir( pbs_options[i] ):
                    pbs_options[i] = pbs_options[i]( script_filename )
            
            # also generate the pbs call? might as well, keep it simple...
            task_summary['commands'][j]['qsub_command'] = create_executable_str( 'qsub' , [script_filename] , pbs_options )

        # rewrite the task summary
        write_task_summary( task_summary_filename )

        task_summaries.append( task_summary_filename )


    # run them all
#    run_VIPUR_task_summaries_serially( task_summaries , single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )
    run_VIPUR_task_summaries_PBS( task_summaries , single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )


    # post processing
    # this look identical!!! :)
    for i in xrange( len( task_summaries ) ):
        # always okay to rerun post processing...should not make any difference
        sequence_only = target_proteins[i][2]
        task_summaries[i] = run_postprocessing( task_summaries[i] , sequence_only = sequence_only )

    return task_summaries

# different from serial method
# collect up all jobs at once, fire off serially (lol) to queue
# need to be cognizant of queue status, rescore must occur later
def run_VIPUR_task_summaries_PBS( task_summaries , single_relax = False , delete_intermediate_relax_files = True ):
    # queue command can be derived from task itself
    # this method will run...until all jobs are complete
    
    # make a list of ALL jobs instead of per protein tasks - all commands have "cd" if there is an "out_path"
    # as a reference, need to track attempts etc. directly, rewrite task summaries frequently
    non_rescore_tasks = sum( [[(i , j) for j in xrange( len( task_summaries[i]['commands'] ) )] for i in task_summaries.keys()] , [] )
    rescore_tasks = [i for i in non_rescore_tasks if 'rescore' in task_summaries[i[0]]['commands'][i[1]]['feature']]
    # only rescore once ALL are complete
    non_rescore_tasks = [i for i in non_rescore_tasks if not i in rescore_tasks]


    run_VIPUR_tasks_PBS( task_summaries , non_rescore_tasks )
    
    # run the rescore tasks
    run_VIPUR_tasks_PBS( task_summaries , rescore_tasks )
    
    # return anything?
    # task summaries should be updated with all the necessary files...

#    for i in xrange( len( task_summaries ) ):
        # tasks will check if the task summaries indicates they have already be run
#        task_summaries[i] = run_task_commands_serially( task_summaries[i] , single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )


# submit jobs until complete or too many attempts
def run_VIPUR_tasks_PBS( task_summaries , task_list , max_pbs_tries = 2 , ddg_monomer_cleanup = True , single_relax = False , delete_intermediate_relax_files = False ):
    # run the non_rescore tasks
    completed = [i for i in task_list if 'run' in task_summaries[i[0]]['commands'][i[1]] and 'success' in task_summaries[i[0]]['commands'][i[1]]['run']]
    # should running_or_queued be saved? written to file?
    running_or_queued = {}
    while not len( completed ) == len( task_list ):
        # check queue status
        queue_status = get_pbs_queue_status()

        # update "running_or_queued" list (?)
        # err, no, does not have information on which job it is...:(
        #for i in queue_status.keys():
        #    if queue_status[i] in ['R' , 'Q']:

        queue_space_occupied = len( [i for i in queue_status.values() if not i in ['C' , 'R']] )    # ignore "C"ompleted jobs, "R"unning job quota are not set by us...
        # if your queue system does not have a separate "R"un quota, remove 'R' from the above!
        available_space = queue_space_occupied - PBS_QUEUE_QUOTA

        
        # launch next jobs in available slots
        if available_space:
            # choose the next job
            jobs_to_run = [i for i in task_list if
                not i in completed and
                not i in running_or_queued.values() and
                not ( 'run' in task_summaries[i[0]]['commands'][i[1]] and
                    ('success' in task_summaries[i[0]]['commands'][i[1]]['run'] or
                    'failure' in task_summaries[i[0]]['commands'][i[1]]['run']) )
                ]
            
            for i in jobs_to_run:
                command_dict = task_summaries[i[0]]['commands'][i[1]]
            
                # write scripts as part of pre processing?...yeah...
                # write the command to a script
                #script_filename = command_dict['out_path'] +'/'*bool( command_dict['out_path'] )+
                #script_filename = command_dict['script_filename']
            
            
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

            
                # submit this script using a queue command
                pbs_command = command_dict['qsub_command']    # SHOULD already have an abspath to the script
                new_job_id = run_local_commandline( pbs_command , collect_stdout = True )
                new_job_id = new_job_id[:new_job_id.find( '.' )]
                
                # save the job id
                # assume its queue
                running_or_queued[new_job_id] = i

        
        # assess outcome of completed jobs
        for job_id in queue_status.keys():
            if queue_status[job_id] == 'C' and job_id in running_or_queued.keys():
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
                    failure_summary = 'success'*complete + (str( tries ) +' tries;failure ' + failure_summary)*(not complete)
                else:
                    # record the number of tries
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


        # update task_summaries e.g. write them!
        # modified: so the task summary records its own name...bah!
        for i in task_summaries:
            if not 'task_summary_filename' in i['filenames'].keys():
                raise NotImplementedError( 'should input the task summary filename (not the summary itself)...' )
            else:
                # write it out
                write_task_summary( i , i['filenames']['task_summary_filename'] )

        
        # pause...
        time.sleep( PBS_QUEUE_MONITOR_DELAY )

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
def get_pbs_queue_status( user = PBS_USER , header_lines = 5 , trailer_lines = 1 , only_job_status = True ):
    # header_lines = 2 for FULL queue, = 5 for USER queue ("-u")
    command = 'qstat'
    if user:
        command += ' -q ' + user

#    queue_info = subprocess.Popen( command.split( ' ' ) , stdout = subprocess.PIPE , stdin = subprocess.PIPE , stderr = subprocess.STDOUT ).communicate()[0]
    queue_info = run_local_commandline( command , collect_stdout = True )

    # simple parsing
    queue_info = queue_info.split( '\n' )[header_lines:-1*trailer_lines]
    queue_info = [[j for j in i.split( ' ' ) if j.strip()] for i in queue_info]

    # optionally only report the job statuses
    if only_job_status:
        queue_info = [(i[0][:i[0].find( '.' )] , i[-2]) for i in queue_info]
        # make into a dict? job ids should be unique...
        queue_info = dict( queue_info )

    return queue_info

# primarily concerned with job status only
# need 1st and 2nd to last columns

################################################################################
# MAIN

if __name__ == '__main__':
    None


