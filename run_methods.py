#!/usr/bin/env python
# :noTabs=true:

"""
"top level" run methods
combine smaller parts
arguably place in VIPUR.py but its getting cluttered...
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
# RUN METHODS

# different options/approaches for setting up job processing

def run_until_complete( command_dict , run_command , check_successful , max_tries ):
    tries = 0
    complete = False
    failure_summary = ''
    while tries < max_tries and not complete:
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


# simples task manager, just fire off sequentially + locally
def run_tasks_locally( task_summary_filename ,
        ddg_monomer_cleanup = True , max_tries = 2 ,
        single_relax = False , delete_intermediate_relax_files = False ):
    # alternate input type
    if isinstance( task_summary_filename , str ):
        task_summary = load_task_summary( task_summary_filename )
    else:
        # shrug? assume its okay
        task_summary = task_summary_filename
    
    # move into the target directory before continuing, for ddg_monomer local
    current_dir = os.getcwd()
    os.chdir( task_summary['out_path'] )
        
    # skip rescore, must wait until corresponding relax job is finished
    print 'launching jobs locally...\n'
    for i in task_summary['commands']:
        if 'rescore' in i['feature']:    # check for now...
            continue
    
        # run the command
        # replace with a while loop checking success...
        check_successful = lambda x : True    # assume successful
        if i['feature'] == 'psiblast':
            # just the pssm
            # assumes blast output structure...
            check_successful = lambda x : check_psiblast_output( x['output_filename'] , PSIBLAST_OPTIONS['out']( x['output_filename'].replace( '.pssm' , '' ) ) )

        elif i['feature'] == 'probe':
            check_successful = lambda x : check_probe_output( x['output_filename'] )

        elif i['feature'] == 'ddg_monomer':
            check_successful = lambda x : check_ddg_monomer_output( x['output_filename'] )

        elif i['feature'].replace( '_native' , '' ) == 'relax' and not 'rescore' in i['feature']:
            check_successful = lambda x : check_relax_output( ROSETTA_RELAX_OPTIONS['out:file:scorefile']( x['output_filename'].replace( '.silent' , '' ) ) , single_relax = single_relax )

        # alternate method, run until complete        
        completed , tries , failure_summary = run_until_complete( i , run_command = run_local_commandline , check_successful = check_successful , max_tries = max_tries )
    
        # optionally cleanup
        if ddg_monomer_cleanup and i['feature'] == 'ddg_monomer':#'ddg' in i['output_filename']:
            print 'ddg_monomer writes useless output files, deleting these now...'
            remove_intermediate_ddg_monomer_files()
        
        # check for complete? failed? how many tries?
        i['run'] = 'success'*completed + (str( tries ) +' tries;failure ' + failure_summary)*(not completed)

    # not do rescore, "but only if complete" (not currently checking)
    relax_commands = [i for i in task_summary['commands'] if i['feature'].replace( '_native' , '' ) == 'relax']
    for i in task_summary['commands']:
        if not 'rescore' in i['feature'] and 'run' in i.keys() and i['run'] == 'success':
            continue
        
        # combine the individual relax runs
        silent_filenames = [j['output_filename'] for j in relax_commands if j['variant'] == i['variant'] and 'run' in j.keys() and j['run'] == 'success']
        if not single_relax:
            if not len( silent_filenames ) == ROSETTA_RELAX_OPTIONS['nstruct']:
                raise Exception( '??? somehow the matching relax run(s) has failed ???\n' + str( i ) )
            score_filenames = [j.replace( '.silent' , '.sc' ) for j in silent_filenames]

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

            merge_rosetta_relax_output( silent_filenames , combined_silent_filename , score_filenames , combined_score_filename , delete_old_files = delete_intermediate_relax_files )
            # rescore already knows the proper filename
        else:
            # just a single match for each
            # output filename should be correct as is :)
            None
        
        # else, its a rescore, check the corresponding relax
        if not 'rescore' in i['feature'] or not 'variant' in i.keys():
            raise Exception( '??? incomplete or erroneous command that is not a rescore ???\n' + str( i ) )
        else:
            # is a rescore AND its relax job it complete, run it
            run_local_commandline( i['command'] )
            i['run'] = 'success'

    # return anything?
    # summary? updated with "run" status
    # should communicate by writing file

    # rewrite task summary...
    write_task_summary( task_summary , task_summary_filename )
    task_summary = load_task_summary( task_summary_filename )

    os.chdir( current_dir )
    
    return task_summary


# overall method, for now, only supports local
def run_VIPUR_in_stages( pdb_filename , variants_filename , prediction_filename = '' , out_path = '' ,
        target_chain = '' , sequence_filename = '' , write_numbering_map = True ,
        sequence_only = False , task_summary_filename = '' , single_relax = False , delete_intermediate_relax_files = True ):
    # preprocessing
    task_summary_filename = run_preprocessing( pdb_filename , variants_filename , prediction_filename , out_path = out_path ,
        target_chain = target_chain , sequence_filename = sequence_filename , write_numbering_map = write_numbering_map ,
        sequence_only = sequence_only , task_summary_filename = task_summary_filename , single_relax = single_relax )
    
    # run locally
    task_summary = run_tasks_locally( task_summary_filename ,
        single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )
    # "sequence_only" details are covered in the "task_summary_filename"
    
    # post processing
    run_postprocessing( task_summary , sequence_only = sequence_only )

    return task_summary

################################################################################
# MULTI-INPUT METHODS

# identify the targets and do pre-processing and post-processing together

def run_VIPUR_serially( pdb_filename = '' , variants_filename = '' ,
        out_path = '' , write_numbering_map = True ,
        single_relax = False , delete_intermediate_relax_files = True ,
        demo = False ):
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
        fa_filenames = [i for i in os.listdir( pdb_filename ) if get_file_extension( i ) == 'fa']
        fa_filenames = [(i , get_root_filename( i ) + variants_filename) for i in fa_filenames if os.path.isfile( pdb_filename +'/'+ get_root_filename( i ) + variants_filename )]

        print 'running VIPUR on all (.pdb,' + variants_filename + ') file pairs found in ' + pdb_filename
        # find .pdb files
        pdb_filenames = [i for i in os.listdir( pdb_filename ) if get_file_extension( i ) == 'pdb']

        # look for pairs
        pdb_filenames = [(i , get_root_filename( i ) + variants_filename) for i in pdb_filenames if os.path.isfile( pdb_filename +'/'+ get_root_filename( i ) + variants_filename )]
#        print [i for i in pdb_filenames if os.path.isfile( pdb_filename +'/'+ get_root_filename( i ) + variants_filename )]        
        print str( len( pdb_filenames ) ) + ' pairs found'

        print str( len( fa_filenames ) ) + ' pairs found'

        # go there...
        os.chdir( pdb_filename )

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
            pdb_filenames = [[pdb_filename , variants_filename]]
        else:
            fa_filenames = [[]]

    # combine all "filenames" to run into unified framework
    target_proteins = []#None]*(len( pdb_filenames ) + len( fa_filenames ))
    for i in pdb_filenames:
        target_proteins.append( i + [False , out_path] )
    for i in fa_filenames:
        target_proteins.append( i + [True , out_path] )

    # pre processing
    task_summaries = []
    for i in target_proteins:
        task_summary_filename = run_preprocessing( i[0] , i[1] ,
            sequence_only = i[2] , out_path = i[3] ,
            write_numbering_map = write_numbering_map , single_relax = single_relax )

        task_summaries.append( task_summary_filename )

    # "sequence_only"
#    sequence_task_summaries = []
#    for i in fa_filenames:
#        task_summary_filename = run_preprocessing( i[0] , i[1] ,
#            out_path = out_path , write_numbering_map = write_numbering_map ,
#            single_relax = single_relax , sequence_only = True )

#        sequence_task_summaries.append( task_summary_filename )

    # run them all
    for i in xrange( len( task_summaries ) ):
        task_summaries[i] = run_tasks_locally( task_summaries[i] , single_relax = single_relax , delete_intermediate_relax_files = delete_intermediate_relax_files )
        # works for both full and sequence only
#    for i in xrange( len( sequence_task_summaries ) ):
#        sequence_task_summaries[i] = run_tasks_locally( sequence_task_summaries[i] )    # those extra options are irrelevant
        
#    run_VIPUR_parallel( pdb_filename , variants_filename ,
#        prediction_filename = prediction_filename , out_path = out_path ,
#        target_chain = target_chain , write_numbering_map = write_numbering_map ,
#        sequence_only = sequence_only )
#    run_VIPUR_deprecated( pdb_filename , variants_filename ,
#        prediction_filename = prediction_filename , out_path = out_path ,
#        target_chain = target_chain , write_numbering_map = write_numbering_map ,
#        sequence_only = sequence_only )
#    run_preprocessing( pdb_filename , variants_filename ,
#        prediction_filename = prediction_filename , out_path = out_path ,
#        target_chain = target_chain , write_numbering_map = write_numbering_map ,
#        sequence_only = sequence_only )
#    run_VIPUR_in_stages( pdb_filename , variants_filename ,
#        prediction_filename = prediction_filename , out_path = out_path ,
#        target_chain = target_chain , write_numbering_map = write_numbering_map ,
#        sequence_only = sequence_only )

    # post processing
    for i in xrange( len( task_summaries ) ):
        sequence_only = target_proteins[i][2]
        task_summaries[i] = run_postprocessing( task_summaries[i] , sequence_only = sequence_only )
#    for i in xrange( len( sequence_task_summaries ) ):
#        sequence_task_summaries[i] = run_postprocessing( sequence_task_summaries[i] , sequence_only = True )

    task_summaries = task_summaries + sequence_task_summaries

    return task_summaries


#def run_VIPUR_serially( pdb_filenames ,
#        out_path = '' , write_numbering_map = True ,
#        single_realx = False ):
    # pre process all the things
#    for i in pdb_filenames:
    
    # run locally in serial
    
    # post process all the things


################################################################################
# COMPLETE METHODS

# integrate preprocessing, running, and postprocessing...
# single method

# single run method
# leave this here for testing
def run_VIPUR_deprecated( pdb_filename , variants_filename , prediction_filename = '' , out_path = '' ,
        target_chain = '' , sequence_filename = '' , relax_nstruct = '-1' , write_numbering_map = True ,
        sequence_only = False ):
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
    try:
        relax_nstruct = int( relax_nstruct )
    except:
        raise Exception( 'Illegal value for relax_nstruct: ' + str( relax_nstruct ) )

    # prepare output writing
    # support writing to  <out_path>
    debug_time = [('start' , time.time())]

    if not os.path.isfile( pdb_filename ):
        raise IOError( 'cannot find ' + pdb_filename + '!!!' )
    if not os.path.isfile( variants_filename ):
        raise IOError( 'cannot find ' + variants_filename + '!!!' )

    if out_path:
        if not os.path.isdir( out_path ):
            create_directory( out_path )

        # make copies of the input files
        # bad solution for now, will copy the input into output directory
        copy_file( os.path.abspath( pdb_filename ) , out_path )
        copy_file( os.path.abspath( variants_filename ) , out_path )
                
        pdb_filename = pdb_filename.split( '/' )[-1]
        variants_filename = variants_filename.split( '/' )[-1]

        # since ddg_monomer just writes into the directory its run in (facepalms)
        os.chdir( out_path )
    else:
        # more bad solutions...but this is just a temporary solution
        if '/' in pdb_filename:
            copy_file( os.path.abspath( pdb_filename ) , '.' )
            pdb_filename = pdb_filename.split( '/' )[-1]
        if '/' in variants_filename:
            copy_file( os.path.abspath( variants_filename ) , '.' )
            variants_filename = variants_filename.split( '/' )[-1]


    # extract the structure's sequence
    if sequence_only and not pdb_filename[-4:] == '.pdb':
        # alternatively accept FASTA input
        sequence = load_fasta( pdb_filename )
        
        if len( sequence ) > 1:
            print 'multiple sequences found in the input FASTA file,\ncontinuing with just the first entry:\n\t' + sequence[0][0]
            sequence = sequence[0][1]    # pairs
        
        # make a dummy residue map...if this doesn't work, the numbering is off
        residue_map = dict( [(str( i ) , i) for i in xrange( 1 , len( sequence ) + 1 )] )

    else:    # the normal protocol
        if not target_chain:    # also support int for index?
            target_chain = extract_chains_from_pdb( pdb_filename )
            target_chain = target_chain[0]    # always a list

        sequence , residue_map , sequence_filename = extract_protein_sequence_from_pdb( pdb_filename , out_filename = sequence_filename , target_chain = target_chain , write_numbering_map = write_numbering_map )
        # "residue_map" is a residue number map between PDB and 1-indexed positions
        # unfortunately many Rosetta programs will do this by default, such as
        # ddg_monomer, or otherwise output 1-indexed protein positions

    # load the variants
    variants = load_variants_file( variants_filename )
    
    # verify the input variants are proper
    variants , failed_variants = check_variants( variants , sequence , residue_map )
    if not variants:
        raise IOError( 'none of the variants in \"' + variants_filename + '\" are acceptable!!?!' )
    print '\ngenerating features for variants:\n\t' + '\n\t'.join( variants )
    sys.stdout.flush()

    # for storing the feature values
    variants = dict( [(i , {}) for i in variants] )
    
    # make the variant structures
    # this is currently done using pymol
#    create_variant_structures( pdb_filename , variants )
    if not sequence_only:
        # this case should never happen, should always know target_chain by here
        if not target_chain:    # also support int for index?
            target_chain = extract_chains_from_pdb( pdb_filename )
            target_chain = target_chain[0]    # always a list

        variant_structures = create_variant_protein_structures( pdb_filename , variants , target_chain )    
    
        # also store meta data, like filenames
        for i in variant_structures:
            for j in variants.keys():
                if ('_' + j + '.') in i:
                    variants[j]['variant structure filename'] = i
                    break
        # verify each structure has one
#        print variants
#        print variant_structures
        failed_structures = [i for i in variants.keys() if not 'variant structure filename' in variants[i].keys()]
        if failed_structures:
            raise IOError( 'something has gone wrong!!?! there are variants that do not have structures,\nsomething has likely gone wrong with variant structure generation' )
    else:
        print 'running in sequence only mode, no structural features calculated'
    
    ####################
    # feature generation
    
    # determine the "aminochange" values
    print '[[VIPURLOG]]evaluating the \"aminochange\" values'
    sys.stdout.flush()
    evaluate_aminochange = lambda nat , var : 2 - int( [i for i in AMINOCHANGE_GROUPS if nat in i][0] == [i for i in AMINOCHANGE_GROUPS if var in i][0] ) - int( nat == var )
    for i in variants.keys():
        variants[i]['aminochange'] = evaluate_aminochange( i[0] , i[-1] )

    debug_time.append( ('preprocessing' , time.time()) )
    
    # run PSIBLAST
    print '[[VIPURLOG]]running PSIBLAST...'
    sys.stdout.flush()
    psiblast_filename = run_psiblast( sequence_filename )
    psiblast_filename = sequence_filename.rstrip( '.fa' ) + '.pssm'
        
    # extract PSIBLAST feature values
    print '[[VIPURLOG]]extracting PSSM features from PSIBLAST output'
    sys.stdout.flush()
    # quality control? check to make sure the PSIBLAST output is proper?
    pssm = extract_pssm_from_psiblast_pssm( psiblast_filename )
    for i in variants.keys():
        nat = i[0]
        pos = i[1:-1]
        pos = residue_map[pos] + 1    # shifted from parsing
        var = i[-1]
        
        variants[i]['pssm_native'] = pssm[pos]['log-likelihood'][nat]
        variants[i]['pssm_variant'] = pssm[pos]['log-likelihood'][var]
        variants[i]['pssm_difference'] = pssm[pos]['log-likelihood'][nat] - pssm[pos]['log-likelihood'][var]
        variants[i]['pssm_information_content'] = pssm[pos]['information content']

    debug_time.append( ('psiblast' , time.time()) )


    # generate structure features
    if not sequence_only:
        # run PROBE
        print '[[VIPURLOG]]running PROBE and determining the ACCP'
        sys.stdout.flush()
        probe_output_filename , positions = run_probe( pdb_filename , variants )
    
        # extract the PROBE feature
        accp = extract_accp_from_probe( probe_output_filename )
        # positions is a parallel list of positions
        accp_dict = dict( [(positions[i] , accp[i]) for i in xrange( len( accp ) )] )
        for i in variants.keys():
            for j in accp_dict.keys():
                if j == i[1:-1]:
                    variants[i]['probe_accp'] = accp_dict[j]
                    continue
        debug_time.append( ('probe' , time.time()) )

    
        # run Rosetta ddg_monomer
        # only run once, NOT on each variant structure, ddg_monomer makes these mutations internally itself
        mut_filename = pdb_filename.rstrip( '.pdb' ) + '.mut'
        write_mut_file( variants , residue_map , mut_filename )

        print '[[VIPURLOG]]running Rosetta ddg_monomer...'
        sys.stdout.flush()
        ddg_monomer_out_filename = run_rosetta_ddg_monomer( pdb_filename , mut_filename )
    
        # extract the ddg_monomer score results
        # (remember the residue mapping)
        ddg_monomer_dict = extract_score_terms_from_ddg_monomer( out_filename = ddg_monomer_out_filename )
        # the residues from this will be "pose numbered", need to map back to the PDB numbering
        header = ddg_monomer_dict['description']
        for i in variants.keys():
            new_key = i[0] + str( residue_map[i[1:-1]] + 1 ) + i[-1]
            for j in xrange( len( header ) ):
                variants[i]['ddg_' + header[j]] = ddg_monomer_dict[new_key][j]
                # added "ddg_" for legacy compatability, is artibrary, make more informative
        debug_time.append( ('ddg_monomer' , time.time()) )

    
        # run Rosetta relax
        relax_new_options = {}
        if not relax_nstruct == -1:
            relax_new_options= {'nstruct' : relax_nstruct}

        print '[[VIPURLOG]]running Rosetta relax on the native and variant structures ...'
        sys.stdout.flush()

#        relax_cmds = {}
#        relax_cmds['__native__'] = run_rosetta_relax( pdb_filename, relax_new_options, run = True )    # reference for other structures
#        master_relax_task_file = open( 'MasterRelaxTaskFile' , 'w' )
#        master_relax_task_file.write( open( relax_cmds['__native__'][1] ).read() )
#        for v in variants:
            # note: this current version assumes 1:1 for variant structures and
            #    variants...could change in the future...
#            variant_structure = variants[v]['variant structure filename']
#            relax_cmds[v] = run_rosetta_relax( variant_structure , relax_new_options , run = True )
#            master_relax_task_file.write( open( relax_cmds[v][1] ).read() )
#        master_relax_task_file.close()

#        for relax_task , (command , tmp_filename , score_filename , silent_filename) in relax_cmds.iteritems():
#            print '>>>' , relax_task , command , tmp_filename , score_filename , silent_filename
        
        #run_local_commandline('parallel --line-buffer --wd %s -S 16/scda000,20/scda001,20/scda002,20/scda003 -a %s'%(out_path, mrtf.name))
        #TODO: in real use, we would doing something to allocated nodes, we fake it here
#        run_local_commandline( 'parallel --line-buffer --wd %s -S 10/scda000,10/linvm021,10/linvm022,10/linvm023,10/linvm024,10/linvm025 -a %s' % (out_path , master_relax_task_file.name) )
        # for each variant and the native structure, create a single silent and a single score file for the individual runs.
  
        # merge the partial output files.
#        for relax_task, (command , tmp_filename , score_filename, silentfn) in relax_cmds.iteritems():
#            score_file = open( score_filename , 'w' )
#            silent_file = open( silent_filename , 'w' )
#            first = True
#            for s in xrange( relax_nstruct ):
#                tag = '_%05d' % s # MUST BE CONSISTENT WITH rosetta_feature_generation.py
#                task_score_file = open( score_filename + tag ).readlines()
#                task_silent_file = open( silent_filename + tag ).readlines()
#                if first:
#                    first = False
#                else:
#                    task_score_file = task_score_file[1:]
#                    task_silent_file = task_silent_file[2:]
#                score_file.write( ''.join( tast_score_file ) )
#                silent_file.write( ''.join( task_silent_file ) )
#            score_file.close()
#            silent_file.close()

#        native_relax_filename = relax_cmds['__native__'][3]    # the silent file for this

##        native_relax_filename = run_rosetta_relax( pdb_filename )    # reference for other structures
#        native_score_filename = run_rosetta_rescore( native_relax_filename , native_filename = pdb_filename )
#        native_scorefile_dict = extract_scores_from_scorefile( native_score_filename )    # save time, only parse this once
#        debug_time.append( ('relax' , time.time()) )

        # old method, deprecated...but here for testing (for now)
        native_relax_filename = run_rosetta_relax_local( pdb_filename )    # reference for other structures
        native_score_filename = run_rosetta_rescore( native_relax_filename , native_filename = pdb_filename )
        native_scorefile_dict = extract_scores_from_scorefile( native_score_filename )    # save time, only parse this once
        debug_time.append( ('relax' , time.time()) )

    
        for i in variants.keys():
            # this file SHOULD exist for every variant
            # if it does not, there was a problem earlier
            variant_structure = variants[i]['variant structure filename']
            # note: this current version assumes 1:1 for variant structures and
            #    variants...could change in the future...
    
            print '[[VIPURLOG]]running Rosetta relax on ' + variant_structure + '...'
            sys.stdout.flush()
            relax_filename = run_rosetta_relax_local( variant_structure )    # old method
#            relax_filename = relax_cmds[i][3]    # the silent file
        
            debug_time.append( ('structure ' + variant_structure , time.time()) )

            # post process, extract scores and generate quartile scores
            print '[[VIPURLOG]]extracting the relax features...'
            sys.stdout.flush()

            # additional step, rescore with Rosetta
            # only required to add in 2 scores that SHOULD be output by default with relax
            # whatever, this is the easiest solution, just make an additional call
            score_filename = run_rosetta_rescore( relax_filename , native_filename = pdb_filename )

            # extract features, use the quartile method to extract comparisons
            # between the native and variant score distributions
            quartile_scores = extract_quartile_score_terms_from_scorefiles( score_filename , native_scorefile_dict )
            # make sure there is not overlap
            overlaps = [j for j in quartile_scores.keys() if j in variants[i].keys()]
            if overlaps:
                raise IOError( '!!?! somehow your scorefile has score terms that match the name(s) of features that already exist!!?! something has gone horribly wrong...' + str( overlaps ) )
            variants[i].update( quartile_scores )
        
    
    ############
    # prediction
    
    # combine and prune feature vector
    print 'classifying the variants as \"deleterious\" or \"neutral\"'
    
    # calculate scores using the logistic regression parameters
    predictions = {}
    for i in variants.keys():
        print 'Variant: ' + i
        if sequence_only:
            # just the sequence only classifier
            # standardize output, add in empty values
            interpretation = VIPUR_classifier.sequence_classifier.classify_instance( variants[i] )
            interpretation['structure label'] = 'NA'
            interpretation['structure P'] = 'NA'
            interpretation['sequence label'] = interpretation['label']
            interpretation['sequence P'] = interpretation['P']

            interpretation['exposure'] = 'NA'
            essential_score = 'NA'
            variants[i]['ddg_total'] = 'NA'
            interpretation['interpretation'] = 'NA'
            interpretation['explanation'] = ['NA']
        else:
            # full prediction
            interpretation = VIPUR_classifier.interpret_classification( variants[i] )
            provide_additional_interpretation( interpretation )    # adds directly into the dict
            essential_score = interpretation['P'] - interpretation['structure P']

        variants[i]['interpretation'] = interpretation
        
        # output append
        predictions[i] = '\t'.join( [str( j ) for j in [
            i ,
            interpretation['label'] ,
            interpretation['P'] ,
#            interpretation['P*'] ,    # not now
#            interpretation['score'] ,    # not now

            interpretation['structure label'] ,
            interpretation['structure P'] ,
            interpretation['sequence label'] ,
            interpretation['sequence P'] ,

            interpretation['exposure'] ,
#            interpretation['conservation'] ,
            essential_score ,
            variants[i]['ddg_total'] ,
            
            interpretation['interpretation'] ,
            ', '.join( interpretation['explanation'] )
#            + interpretation['terms']    # not now...
            ] ] )
    
    # write output predictions etc.
    if not prediction_filename:
        prediction_filename = pdb_filename.rstrip( '.pdb' ) + '.predictions'
    f = open( prediction_filename , 'w' )
    # very shameful, hardcoded indices again...
    temp = predictions.values()
#    print temp[0]
#    print [i.split( '\t' )[:5] for i in temp]
    temp = sorted( temp , key = lambda x : 1/float( x.split( '\t' )[2] ) )
    f.write( '\t'.join( PREDICTION_OUTPUT_HEADER ) +'\n'+ '\n'.join( temp ) )
    f.close()

    # debug putput
#    print variants
    
    debug_time.append( ('postprocessing' , time.time()) )
    # debug analyze time
    for i in xrange( len( debug_time ) - 1 ):
        dt = debug_time[i + 1][1] - debug_time[i][1]
        minutes = int( dt/60 )
        dt = dt - minutes*60
#        print debug_time[i + 1][0].ljust( 20 ) + str( round( debug_time[i + 1][1] - debug_time[i][1] , 4 ) ) +'s'
        print debug_time[i + 1][0].ljust( 40 ) + str( minutes ) + 'min ' + str( round( dt , 4 ) ) +'s'





# single run method
# currently unsupported, modifying to be "task" oriented
def run_VIPUR_parallel( pdb_filename , variants_filename , prediction_filename = '' , out_path = '' ,
        target_chain = '' , sequence_filename = '' , relax_nstruct = '-1' , write_numbering_map = True ,
        sequence_only = False ):
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
    try:
        relax_nstruct = int( relax_nstruct )
    except:
        raise Exception( 'Illegal value for relax_nstruct: ' + str( relax_nstruct ) )

    # prepare output writing
    # support writing to  <out_path>
    debug_time = [('start' , time.time())]

    if not os.path.isfile( pdb_filename ):
        raise IOError( 'cannot find ' + pdb_filename + '!!!' )
    if not os.path.isfile( variants_filename ):
        raise IOError( 'cannot find ' + variants_filename + '!!!' )

    if out_path:
        if not os.path.isdir( out_path ):
            create_directory( out_path )

        # make copies of the input files
        # bad solution for now, will copy the input into output directory
        copy_file( os.path.abspath( pdb_filename ) , out_path )
        copy_file( os.path.abspath( variants_filename ) , out_path )
                
        pdb_filename = pdb_filename.split( '/' )[-1]
        variants_filename = variants_filename.split( '/' )[-1]

        # since ddg_monomer just writes into the directory its run in (facepalms)
        os.chdir( out_path )
    else:
        # more bad solutions...but this is just a temporary solution
        if '/' in pdb_filename:
            copy_file( os.path.abspath( pdb_filename ) , '.' )
            pdb_filename = pdb_filename.split( '/' )[-1]
        if '/' in variants_filename:
            copy_file( os.path.abspath( variants_filename ) , '.' )
            variants_filename = variants_filename.split( '/' )[-1]


    # extract the structure's sequence
    if sequence_only and not pdb_filename[-4:] == '.pdb':
        # alternatively accept FASTA input
        sequence = load_fasta( pdb_filename )
        
        if len( sequence ) > 1:
            print 'multiple sequences found in the input FASTA file,\ncontinuing with just the first entry:\n\t' + sequence[0][0]
            sequence = sequence[0][1]    # pairs
        
        # make a dummy residue map...if this doesn't work, the numbering is off
        residue_map = dict( [(str( i ) , i) for i in xrange( 1 , len( sequence ) + 1 )] )

    else:    # the normal protocol
        if not target_chain:    # also support int for index?
            target_chain = extract_chains_from_pdb( pdb_filename )
            target_chain = target_chain[0]    # always a list

        sequence , residue_map , sequence_filename = extract_protein_sequence_from_pdb( pdb_filename , out_filename = sequence_filename , target_chain = target_chain , write_numbering_map = write_numbering_map )
        # "residue_map" is a residue number map between PDB and 1-indexed positions
        # unfortunately many Rosetta programs will do this by default, such as
        # ddg_monomer, or otherwise output 1-indexed protein positions

    # load the variants
    variants = load_variants_file( variants_filename )
    
    # verify the input variants are proper
    variants , failed_variants = check_variants( variants , sequence , residue_map )
    if not variants:
        raise IOError( 'none of the variants in \"' + variants_filename + '\" are acceptable!!?!' )
    print '\ngenerating features for variants:\n\t' + '\n\t'.join( variants )
    sys.stdout.flush()

    # for storing the feature values
    variants = dict( [(i , {}) for i in variants] )
    
    # make the variant structures
    # this is currently done using pymol
#    create_variant_structures( pdb_filename , variants )
    if not sequence_only:
        # this case should never happen, should always know target_chain by here
        if not target_chain:    # also support int for index?
            target_chain = extract_chains_from_pdb( pdb_filename )
            target_chain = target_chain[0]    # always a list

        variant_structures = create_variant_protein_structures( pdb_filename , variants , target_chain )    
    
        # also store meta data, like filenames
        for i in variant_structures:
            for j in variants.keys():
                if ('_' + j + '.') in i:
                    variants[j]['variant structure filename'] = i
                    break
        # verify each structure has one
#        print variants
#        print variant_structures
        failed_structures = [i for i in variants.keys() if not 'variant structure filename' in variants[i].keys()]
        if failed_structures:
            raise IOError( 'something has gone wrong!!?! there are variants that do not have structures,\nsomething has likely gone wrong with variant structure generation' )
    else:
        print 'running in sequence only mode, no structural features calculated'
    
    ####################
    # feature generation
    
    # determine the "aminochange" values
    print '[[VIPURLOG]]evaluating the \"aminochange\" values'
    sys.stdout.flush()
    evaluate_aminochange = lambda nat , var : 2 - int( [i for i in AMINOCHANGE_GROUPS if nat in i][0] == [i for i in AMINOCHANGE_GROUPS if var in i][0] ) - int( nat == var )
    for i in variants.keys():
        variants[i]['aminochange'] = evaluate_aminochange( i[0] , i[-1] )

    debug_time.append( ('preprocessing' , time.time()) )
    
    # run PSIBLAST
    print '[[VIPURLOG]]running PSIBLAST...'
    sys.stdout.flush()
    psiblast_filename = run_psiblast( sequence_filename )
    psiblast_filename = sequence_filename.rstrip( '.fa' ) + '.pssm'
        
    # extract PSIBLAST feature values
    print '[[VIPURLOG]]extracting PSSM features from PSIBLAST output'
    sys.stdout.flush()
    # quality control? check to make sure the PSIBLAST output is proper?
    pssm = extract_pssm_from_psiblast_pssm( psiblast_filename )
    for i in variants.keys():
        nat = i[0]
        pos = i[1:-1]
        pos = residue_map[pos] + 1    # shifted from parsing
        var = i[-1]
        
        variants[i]['pssm_native'] = pssm[pos]['log-likelihood'][nat]
        variants[i]['pssm_variant'] = pssm[pos]['log-likelihood'][var]
        variants[i]['pssm_difference'] = pssm[pos]['log-likelihood'][nat] - pssm[pos]['log-likelihood'][var]
        variants[i]['pssm_information_content'] = pssm[pos]['information content']

    debug_time.append( ('psiblast' , time.time()) )


    # run PROBE
    if not sequence_only:
        print '[[VIPURLOG]]running PROBE and determining the ACCP'
        sys.stdout.flush()
        probe_output_filename , positions = run_probe( pdb_filename , variants )
    
        # extract the PROBE feature
        accp = extract_accp_from_probe( probe_output_filename )
        # positions is a parallel list of positions
        accp_dict = dict( [(positions[i] , accp[i]) for i in xrange( len( accp ) )] )
        for i in variants.keys():
            for j in accp_dict.keys():
                if j == i[1:-1]:
                    variants[i]['probe_accp'] = accp_dict[j]
                    continue
        debug_time.append( ('probe' , time.time()) )

    
        # run Rosetta ddg_monomer
        # only run once, NOT on each variant structure, ddg_monomer makes these mutations internally itself
        mut_filename = pdb_filename.rstrip( '.pdb' ) + '.mut'
        write_mut_file( variants , residue_map , mut_filename )

        print '[[VIPURLOG]]running Rosetta ddg_monomer...'
        sys.stdout.flush()
        ddg_monomer_out_filename = run_rosetta_ddg_monomer( pdb_filename , mut_filename )
    
        # extract the ddg_monomer score results
        # (remember the residue mapping)
        ddg_monomer_dict = extract_score_terms_from_ddg_monomer( out_filename = ddg_monomer_out_filename )
        # the residues from this will be "pose numbered", need to map back to the PDB numbering
        header = ddg_monomer_dict['description']
        for i in variants.keys():
            new_key = i[0] + str( residue_map[i[1:-1]] + 1 ) + i[-1]
            for j in xrange( len( header ) ):
                variants[i]['ddg_' + header[j]] = ddg_monomer_dict[new_key][j]
                # added "ddg_" for legacy compatability, is artibrary, make more informative
        debug_time.append( ('ddg_monomer' , time.time()) )

    
        # run Rosetta relax
        relax_new_options = {}
        if not relax_nstruct == -1:
            relax_new_options= {'nstruct' : relax_nstruct}

        print '[[VIPURLOG]]running Rosetta relax on the native and variant structures ...'
        sys.stdout.flush()

        relax_cmds = {}
        relax_cmds['__native__'] = run_rosetta_relax( pdb_filename, relax_new_options, run = True )    # reference for other structures
        master_relax_task_file = open( 'MasterRelaxTaskFile' , 'w' )
        master_relax_task_file.write( open( relax_cmds['__native__'][1] ).read() )
        for v in variants:
            # note: this current version assumes 1:1 for variant structures and
            #    variants...could change in the future...
            variant_structure = variants[v]['variant structure filename']
            relax_cmds[v] = run_rosetta_relax( variant_structure , relax_new_options , run = True )
            master_relax_task_file.write( open( relax_cmds[v][1] ).read() )
        master_relax_task_file.close()

        for relax_task , (command , tmp_filename , score_filename , silent_filename) in relax_cmds.iteritems():
            print '>>>' , relax_task , command , tmp_filename , score_filename , silent_filename
        
        #run_local_commandline('parallel --line-buffer --wd %s -S 16/scda000,20/scda001,20/scda002,20/scda003 -a %s'%(out_path, mrtf.name))
        #TODO: in real use, we would doing something to allocated nodes, we fake it here
        run_local_commandline( 'parallel --line-buffer --wd %s -S 10/scda000,10/linvm021,10/linvm022,10/linvm023,10/linvm024,10/linvm025 -a %s' % (out_path , master_relax_task_file.name) )
        # for each variant and the native structure, create a single silent and a single score file for the individual runs.
  
        # merge the partial output files.
        for relax_task, (command , tmp_filename , score_filename, silentfn) in relax_cmds.iteritems():
            score_file = open( score_filename , 'w' )
            silent_file = open( silent_filename , 'w' )
            first = True
            for s in xrange( relax_nstruct ):
                tag = '_%05d' % s # MUST BE CONSISTENT WITH rosetta_feature_generation.py
                task_score_file = open( score_filename + tag ).readlines()
                task_silent_file = open( silent_filename + tag ).readlines()
                if first:
                    first = False
                else:
                    task_score_file = task_score_file[1:]
                    task_silent_file = task_silent_file[2:]
                score_file.write( ''.join( tast_score_file ) )
                silent_file.write( ''.join( task_silent_file ) )
            score_file.close()
            silent_file.close()

        native_relax_filename = relax_cmds['__native__'][3]    # the silent file for this

#        native_relax_filename = run_rosetta_relax( pdb_filename )    # reference for other structures
        native_score_filename = run_rosetta_rescore( native_relax_filename , native_filename = pdb_filename )
        native_scorefile_dict = extract_scores_from_scorefile( native_score_filename )    # save time, only parse this once
        debug_time.append( ('relax' , time.time()) )
    
        for i in variants.keys():
            # this file SHOULD exist for every variant
            # if it does not, there was a problem earlier
            variant_structure = variants[i]['variant structure filename']
            # note: this current version assumes 1:1 for variant structures and
            #    variants...could change in the future...
    
            print '[[VIPURLOG]]running Rosetta relax on ' + variant_structure + '...'
            sys.stdout.flush()
            relax_filename = relax_cmds[i][3]    # the silent file
        
            debug_time.append( ('structure ' + variant_structure , time.time()) )

            # post process, extract scores and generate quartile scores
            print '[[VIPURLOG]]extracting the relax features...'
            sys.stdout.flush()

            # additional step, rescore with Rosetta
            # only required to add in 2 scores that SHOULD be output by default with relax
            # whatever, this is the easiest solution, just make an additional call
            score_filename = run_rosetta_rescore( relax_filename , native_filename = pdb_filename )

            # extract features, use the quartile method to extract comparisons
            # between the native and variant score distributions
            quartile_scores = extract_quartile_score_terms_from_scorefiles( score_filename , native_scorefile_dict )
            # make sure there is not overlap
            overlaps = [j for j in quartile_scores.keys() if j in variants[i].keys()]
            if overlaps:
                raise IOError( '!!?! somehow your scorefile has score terms that match the name(s) of features that already exist!!?! something has gone horribly wrong...' + str( overlaps ) )
            variants[i].update( quartile_scores )
        
    
    ############
    # prediction
    
    # combine and prune feature vector
    print 'classifying the variants as \"deleterious\" or \"neutral\"'
    
    # calculate scores using the logistic regression parameters
    predictions = {}
    for i in variants.keys():
        print 'Variant: ' + i
        if sequence_only:
            # just the sequence only classifier
            # standardize output, add in empty values
            interpretation = VIPUR_classifier.sequence_classifier.classify_instance( variants[i] )
            interpretation['structure label'] = 'NA'
            interpretation['structure P'] = 'NA'
            interpretation['sequence label'] = interpretation['label']
            interpretation['sequence P'] = interpretation['P']

            interpretation['exposure'] = 'NA'
            essential_score = 'NA'
            variants[i]['ddg_total'] = 'NA'
            interpretation['interpretation'] = 'NA'
            interpretation['explanation'] = ['NA']
        else:
            # full prediction
            interpretation = VIPUR_classifier.interpret_classification( variants[i] )
            provide_additional_interpretation( interpretation )    # adds directly into the dict
            essential_score = interpretation['P'] - interpretation['structure P']

        variants[i]['interpretation'] = interpretation
        
        # output append
        predictions[i] = '\t'.join( [str( j ) for j in [
            i ,
            interpretation['label'] ,
            interpretation['P'] ,
#            interpretation['P*'] ,    # not now
#            interpretation['score'] ,    # not now

            interpretation['structure label'] ,
            interpretation['structure P'] ,
            interpretation['sequence label'] ,
            interpretation['sequence P'] ,

            interpretation['exposure'] ,
#            interpretation['conservation'] ,
            essential_score ,
            variants[i]['ddg_total'] ,
            
            interpretation['interpretation'] ,
            ', '.join( interpretation['explanation'] )
#            + interpretation['terms']    # not now...
            ] ] )
    
    # write output predictions etc.
    if not prediction_filename:
        prediction_filename = pdb_filename.rstrip( '.pdb' ) + '.predictions'
    f = open( prediction_filename , 'w' )
    # very shameful, hardcoded indices again...
    temp = predictions.values()
#    print temp[0]
#    print [i.split( '\t' )[:5] for i in temp]
    temp = sorted( temp , key = lambda x : 1/float( x.split( '\t' )[2] ) )
    f.write( '\t'.join( PREDICTION_OUTPUT_HEADER ) +'\n'+ '\n'.join( temp ) )
    f.close()

    # debug putput
#    print variants
    
    debug_time.append( ('postprocessing' , time.time()) )
    # debug analyze time
    for i in xrange( len( debug_time ) - 1 ):
        dt = debug_time[i + 1][1] - debug_time[i][1]
        minutes = int( dt/60 )
        dt = dt - minutes*60
#        print debug_time[i + 1][0].ljust( 20 ) + str( round( debug_time[i + 1][1] - debug_time[i][1] , 4 ) ) +'s'
        print debug_time[i + 1][0].ljust( 40 ) + str( minutes ) + 'min ' + str( round( dt , 4 ) ) +'s'


