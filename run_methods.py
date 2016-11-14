#!/usr/bin/env python
# :noTabs=true:

"""
methods for running groups of task summaries
"""

################################################################################
# IMPORT

# common modules
import os
import sys
import shutil
import subprocess
import time    # for debugging only

# bigger modules

# custom modules
from vipur_settings import *
from helper_methods import *
from pre_processing import *
from post_processing import *

from psiblast_feature_generation import *
from probe_feature_generation import *
from rosetta_feature_generation import *

from classification import *

################################################################################
# SERIAL RUN METHODS

# different options/approaches for setting up job processing

def run_serially_until_complete( command_dict , run_command , check_successful , max_tries ):
    tries = 0
    complete = False
    failure_summary = ''
    while tries < max_tries and not complete:
        # already done?
        if 'run' in command_dict.keys() and command_dict['run'] == 'success':
            complete = True
            break
    
        # run it
        run_command( command_dict['command'] )
        
        # check successful
        success = check_successful( command_dict )
        if isinstance( success , bool ):
            complete = success
        elif len( success ) > 1 and isinstance( success[0] , bool ):
            complete = success[0]
            failure_summary += ' '+ ';'.join( [str( j ) for j in success[1:]] ) +' '
                
        tries += 1
    return complete , tries , failure_summary


# simple task manager, just fire off sequentially + locally
def run_task_commands_serially( task_summary_filename ,
        ddg_monomer_cleanup = True , max_tries = 2 ,
        single_relax = False , delete_intermediate_relax_files = False ):
    # alternate input type
    if isinstance( task_summary_filename , str ):
        task_summary = load_task_summary( task_summary_filename )
    else:
        # shrug? assume its okay
        task_summary = task_summary_filename
        # modified: task summary SHOULD have this field
        if 'task_summary_filename' in task_summary['filenames'].keys():
            task_summary_filename = task_summary['filenames']['task_summary_filename']
        else:
            raise NotImplementedError( 'should input the task summary filename (not the summary itself)...' )
    
    # move into the target directory before continuing, for ddg_monomer local
    current_dir = os.getcwd()
    os.chdir( task_summary['out_path'] )    # for PBS etc. instead add "cd " to run script
        
    # skip rescore, must wait until corresponding relax job is finished
    print 'launching jobs locally...\n'
    for i in task_summary['commands']:
        # these MUST occur AFTER the relax jobs have been run and combined
        if 'rescore' in i['feature']:    # check for now...
            continue
        
        # skip those that have alreay run
        if 'run' in i.keys() and i['run'] == 'success' and os.path.isfile( i['output_filename'] ):
            print i['output_filename'] + ' appears to have been successfully generated, do not run it again'
            continue

        check_successful = determine_check_successful_function( i , single_relax = single_relax )

        # alternate method, run until complete
        completed , tries , failure_summary = run_serially_until_complete( i , run_command = run_local_commandline , check_successful = check_successful , max_tries = max_tries )
    
        # optionally cleanup
        if ddg_monomer_cleanup and i['feature'] == 'ddg_monomer':#'ddg' in i['output_filename']:
            print 'ddg_monomer writes useless output files, deleting these now...'
            remove_intermediate_ddg_monomer_files()
        
        # check for complete? failed? how many tries?
        i['run'] = 'success'*completed + (str( tries ) +' tries;failure ' + failure_summary)*(not completed)

    # not do rescore, "but only if complete" (not currently checking)
    relax_commands = [i for i in task_summary['commands'] if i['feature'].replace( '_native' , '' ) == 'relax']
    for i in task_summary['commands']:
        if 'run' in i.keys() and i['run'] == 'success':
            if 'rescore' in i['feature'] and os.path.isfile( i['output_filename'] ):
                print i['output_filename'] + ' appears to have been successfully generated, do not run it again'
            # skip, either its not a rescore, or its a rescore that ran successfully
            continue
        
        # combine the individual relax runs
        silent_filenames = [j['output_filename'] for j in relax_commands if j['variant'] == i['variant'] and 'run' in j.keys() and j['run'] == 'success']
        # actually need to identify the combined_silent_filename, be sure the relax files have not already been merged
        # which variant
        target_variant = [j for j in task_summary['variants'].keys() if j.split( '_' )[-1] == i['variant'] and j.split( '_' )[0] in i['command']]
        if not target_variant:
            # its native
            combined_silent_filename = task_summary['other']['combined_native_silent_filename']
            combined_score_filename = task_summary['other']['combined_native_score_filename']
        elif len( target_variant ) > 1:
            raise Exception( '??? found more than on matching variant ???\n' + ', '.join( target_variant ) )
        else:
            # found it
            combined_silent_filename = task_summary['variants'][target_variant[0]]['combined_silent_filename']
            combined_score_filename = task_summary['variants'][target_variant[0]]['combined_score_filename']

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
        
        # else, its a rescore, check the corresponding relax
        # should already point to the proper combined_silent_filename
        if not 'rescore' in i['feature'] or not 'variant' in i.keys():
            raise Exception( '??? incomplete or erroneous command that is not a rescore ???\n' + str( i ) )
        else:
            # is a rescore AND its relax job it complete, run it
            run_local_commandline( i['command'] )    # not 'run_until_complete'?
            i['run'] = 'success'

    # return anything?
    # summary? updated with "run" status
    # should communicate by writing file

    # rewrite task summary...
    #if not isinstance( task_summary_filename , str ):
        #raise NotImplementedError( 'should input the task summary filename (not the summary itself)...' )
        # actually, is currently implemented...see above
    write_task_summary( task_summary , task_summary_filename )
    task_summary = load_task_summary( task_summary_filename )

    os.chdir( current_dir )
    
    return task_summary


def run_VIPUR_task_summaries_serially( task_summaries , single_relax = False , delete_intermediate_relax_files = True ):
    for i in xrange( len( task_summaries ) ):
        # tasks will check if the task summaries indicates they have already be run
        task_summaries[i] = run_task_commands_serially( task_summaries[i] , single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )


#####################
# MULTI-INPUT METHODS

# identify the targets and do pre-processing and post-processing together

def run_VIPUR_serially( pdb_filename = '' , variants_filename = '' ,
        out_path = '' , write_numbering_map = True ,
        single_relax = False , delete_intermediate_relax_files = True ,
        demo = False , rerun_preprocessing = False ):
    """
    Runs the VIPUR pipeline on the variants of the protein in  <pdb_filename>
    (a protein structure) that are described in  <variants_filename>
    
    Optionally specify the  <out_path>  (where to write output)*
    Optionally provide the  <sequence_filename>  you desire for the protein
        sequence extracted from  <pdb_filename>  (needed for input to PSIBLAST)
    Optionally specify the  <target_chain>  within the PDB
    Optionally  <write_numbering_map>  of how to convert 1-indexed positions
        (such as the default "pose numbering" of Rosetta) to the position
        numbers in  <pdb_filename>
        note: this is to avoid requiring all input PDBs be renumbered
    Optionally run  <sequence_only>  (requires an input sequence file in FASTA
        format for  <pdb_filename>, instead of a PDB file, does not generate
        variant structures, or run any structural analysis)
        many columns are indicated as "Not Applicable" in the output predictions
    
    
    VIPUR analyzes protein variants by considering conservation scores and
    structural scores to identify variants that are likely to disrupt protein
    function
    
    Conservation scores are derived from a PSSM of similar sequences found
    using PSIBLAST against the NCBI nr database
    
    Structural analysis is done using Rosetta to:
        consider variant structures by rapid structure optmization, allowing
        fast evaluation of approximate variant ddG values
        (using Rosetta ddg_monomer)
        and
        refining variant structures and considering the distribution of
        energies and structural scores (rms, gdtmm) across several low energy
        conformations (physically near the input conformation)
    
    Additional features are provided by an internal "aminochange" classification
    (crude similarity of amino acid properties) and the variant position
    surface area, evaluated using PROBE
    
    These analyses are combined with a learned Logistic Regression model to
    classify the input variants as "neutral" are "deleterious" and provide a
    structurally-informed hypothesis as to why variants are likely to
    disrupt the protein
    """
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


    # pre processing
    task_summaries = []
    for i in target_proteins:
        # check paths...er, done properly in preprocessing...

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

        task_summaries.append( task_summary_filename )

    # "sequence_only"
#    sequence_task_summaries = []
#    for i in fa_filenames:
#        task_summary_filename = run_preprocessing( i[0] , i[1] ,
#            out_path = out_path , write_numbering_map = write_numbering_map ,
#            single_relax = single_relax , sequence_only = True )

#        sequence_task_summaries.append( task_summary_filename )
#    raw_input( 'preprocessing' )

    # run them all
    run_VIPUR_task_summaries_serially( task_summaries , single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )

    # post processing
    for i in xrange( len( task_summaries ) ):
        # always okay to rerun post processing...should not make any difference
        sequence_only = target_proteins[i][2]
        task_summaries[i] = run_postprocessing( task_summaries[i] , sequence_only = sequence_only )

    return task_summaries


# sort the check_successful function
def determine_check_successful_function( command_dict , single_relax = False ):
    check_successful = lambda x : True    # assume successful
    if command_dict['feature'] == 'psiblast':
        # just the pssm
        # assumes blast output structure...
        check_successful = lambda x : check_psiblast_output( x['output_filename'] , PSIBLAST_OPTIONS['out']( x['output_filename'].replace( '.pssm' , '' ) ) )

    elif command_dict['feature'] == 'probe':
        check_successful = lambda x : check_probe_output( x['output_filename'] )

    elif command_dict['feature'] == 'ddg_monomer':
        check_successful = lambda x : check_ddg_monomer_output( x['output_filename'] )

    elif command_dict['feature'].replace( '_native' , '' ) == 'relax' and not 'rescore' in command_dict['feature']:
        check_successful = lambda x : check_relax_output( ROSETTA_RELAX_OPTIONS['out:file:scorefile']( x['output_filename'].replace( '.silent' , '' ) ) , single_relax = single_relax )

    return check_successful


