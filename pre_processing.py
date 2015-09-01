#!/usr/bin/env python
# :noTabs=true:

"""
main preprocessing
"""

################################################################################
# IMPORT

# common modules
import os
import sys

# bigger modules

# custom modules
from helper_methods import load_variants_file , check_variants , extract_chains_from_pdb , get_file_extension , get_root_filename , create_directory , copy_file

from psiblast_feature_generation import load_fasta , extract_protein_sequence_from_pdb , run_psiblast
from probe_feature_generation import run_probe
from rosetta_feature_generation import create_variant_protein_structures , write_mut_file , run_rosetta_ddg_monomer , run_rosetta_relax_local , run_rosetta_rescore
from vipur_settings import ROSETTA_RELAX_OPTIONS

################################################################################
# MAIN PREPROCESSING

# prepare commands etc.
def run_preprocessing( pdb_filename , variants_filename , prediction_filename = '' , out_path = '' ,
        target_chain = '' , sequence_filename = '' , write_numbering_map = True ,
        sequence_only = False , task_summary_filename = '' ,
        single_relax = False , rosetta_relax_options = ROSETTA_RELAX_OPTIONS ,
        pymol_environment_setup = '' ):
    # prepare output writing
    # support writing to  <out_path>
    #debug_time = [('start' , time.time())]

    # check paths
    if not os.path.isfile( pdb_filename ):
        raise IOError( 'cannot find ' + pdb_filename + '!!!' )
    if not os.path.isfile( variants_filename ):
        raise IOError( 'cannot find ' + variants_filename + '!!!' )

    file_extension = get_file_extension( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    base_root_filename = root_filename.split( '/' )[-1]
    
#    print root_filename
#    raw_input( 'pause' )

    if not task_summary_filename:
        task_summary_filename = root_filename + '.task_summary'

    # where to store output?
    if out_path:
        if not os.path.isdir( out_path ):
#            print out_path
#            print os.getcwd()
            create_directory( out_path )

        # make copies of the input files
        # bad solution for now, will copy the input into output directory
        if not os.path.isfile( out_path +'/'+ base_root_filename +'.pdb' ):
            copy_file( os.path.abspath( pdb_filename ) , out_path )
        if not os.path.isfile( out_path +'/'+ base_root_filename + '.txt' ):
            copy_file( os.path.abspath( variants_filename ) , out_path )
        
        root_filename = out_path +'/'+ base_root_filename
        pdb_filename = root_filename +'.'+ file_extension
        variants_filename = root_filename +'.txt'
        
        # since ddg_monomer just writes into the directory its run in (facepalms)
#        current_dir = os.getcwd()
#        os.chdir( out_path )    # not happening here now err, for variant structures
    else:
        # otherwise, operate in the current directory
        if not os.path.isfile( pdb_filename ):
            copy_file( os.path.abspath( pdb_filename ) , '.' )
        if not os.path.isfile( variants_filename ):
            copy_file( os.path.abspath( variants_filename ) , '.' )

    # just the base filename
#    if '/' in pdb_filename:
#        pdb_filename = pdb_filename.split( '/' )[-1]
#    if '/' in variants_filename:
#        variants_filename = variants_filename.split( '/' )[-1]


    # extract the structure's sequence
    if sequence_only and not file_extension == 'pdb':
        # alternatively accept FASTA input
        sequence_filename = pdb_filename
        sequence = load_fasta( pdb_filename )
        
        if len( sequence ) > 1:
            print 'multiple sequences found in the input FASTA file,\ncontinuing with just the first entry:\n\t' + sequence[0][0]
            sys.stdout.flush()
        sequence = sequence[0][1]    # pairs (id , sequence), do this either way
        
        # make a dummy residue map...if this doesn't work, the numbering is off
        residue_map = dict( [(str( i ) , i - 1) for i in xrange( 1 , len( sequence ) + 1 )] )
                
        # write this map...however sad it is
        if isinstance( write_numbering_map , str ):
            numbering_map_filename = write_numbering_map
        else:
            numbering_map_filename = root_filename + '.numbering_map'

#        print numbering_map_filename
#        raw_input( 'moo' )
        f = open( numbering_map_filename , 'w' )
        keys = sorted( residue_map.keys() , key = lambda x : int( x ) )
        f.write( '\n'.join( ['\t'.join( [i , str( residue_map[i] + 1 ) , sequence[residue_map[i]]] ) for i in keys] ) )
        f.close()

        print 'wrote the numbering for ' + pdb_filename + ' to ' + numbering_map_filename
        sys.stdout.flush()
    else:    # the normal protocol
        if not target_chain:    # also support int for index?
#            print pdb_filename
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
    if not sequence_only:
        # this case should never happen, should always know target_chain by here
        if not target_chain:    # also support int for index?
            target_chain = extract_chains_from_pdb( pdb_filename )
            target_chain = target_chain[0]    # always a list

        variant_structures = create_variant_protein_structures( pdb_filename , variants , target_chain , pymol_environment_setup = pymol_environment_setup )    
    
        # also store meta data, like filenames
        for i in variant_structures:
            for j in variants.keys():
                if ('_' + j + '.') in i:
                    variants[j]['variant structure filename'] = i
                    break
        # verify each structure has one
        failed_structures = [i for i in variants.keys() if not 'variant structure filename' in variants[i].keys()]
        if failed_structures:
            raise IOError( 'something has gone wrong!!?! there are variants that do not have structures,\nsomething has likely gone wrong with variant structure generation' )
    else:
        print 'running in sequence only mode, no structural features calculated'
        sys.stdout.flush()


    # generate "commands" for feature generation

    # save for post processing
    # aminochange evaluation
    
    # run PSIBLAST
    print '[[VIPURLOG]]generating PSIBLAST run command'
    sys.stdout.flush()
    psiblast_command , psiblast_filename = run_psiblast( sequence_filename , run = False )
    # extract features from pssm


    # generate structure features
    if not sequence_only:
        # run PROBE
        print '[[VIPURLOG]]generating PROBE run command'
        sys.stdout.flush()
        probe_command , probe_output_filename , probe_positions = run_probe( pdb_filename , variants , run = False )
        # extract the PROBE feature

    
        # run Rosetta ddg_monomer
        # only run once, NOT on each variant structure, ddg_monomer makes these mutations internally itself
        mut_filename = pdb_filename.rstrip( file_extension ) + 'mut'    # "." still there
        write_mut_file( variants , residue_map , mut_filename )

        print '[[VIPURLOG]]generating Rosetta ddg_monomer run command'
        sys.stdout.flush()
        ddg_monomer_command , ddg_monomer_out_filename = run_rosetta_ddg_monomer( pdb_filename , mut_filename , out_path = out_path , run = False )
            # extract the ddg_monomer score results
        # (remember the residue mapping)

    
        # run Rosetta relax
        print '[[VIPURLOG]]generating Rosetta relax commands on the native and variant structures'
        sys.stdout.flush()

        # old method, deprecated...but here for testing (for now)
        nstruct = rosetta_relax_options['nstruct']
        jran = rosetta_relax_options['run:jran']
        native_relax_commands = []
        if single_relax:
            native_relax_command , native_relax_filename = run_rosetta_relax_local( pdb_filename , run = False )    # reference for other structures

            native_relax_commands.append( native_relax_command )
            native_relax_commands.append( native_relax_filename )
        else:
            # as separate runs, 1 per nstruct!
            native_relax_commands = []
            for j in xrange( nstruct ):
                extra_options = {
                    'nstruct' : 1 ,
                    'run:jran' : jran + j ,
                    'out:file:silent' : rosetta_relax_options['out:file:silent']( pdb_filename ).replace( '.' + file_extension + '.silent' , '_' + str( j + 1 ) + '.silent' ) ,
                    'out:file:scorefile' : rosetta_relax_options['out:file:scorefile']( pdb_filename ).replace( '.' + file_extension + '.sc' , '_' + str( j + 1 ) + '.sc' )
                    }
                native_relax_command , native_relax_filename = run_rosetta_relax_local( pdb_filename , run = False , extra_options = extra_options )
                
                    
                native_relax_commands.append( native_relax_command )
                native_relax_commands.append( native_relax_filename )
        
        # for records...
        combined_native_silent_filename = rosetta_relax_options['out:file:silent']( pdb_filename.replace( '.' + file_extension , '' ) )
        combined_native_score_filename = rosetta_relax_options['out:file:scorefile']( pdb_filename.replace( '.' + file_extension , '' ) )

        # rescore needs the combined filenames        
        native_score_command , native_score_filename = run_rosetta_rescore( combined_native_silent_filename , native_filename = pdb_filename , run = False )

        native_relax_commands.append( native_score_command )
        native_relax_commands.append( native_score_filename )


        # individual variant relax    
        variant_relax = {}
        for i in variants.keys():
            # this file SHOULD exist for every variant
            # if it does not, there was a problem earlier
            variant_structure = variants[i]['variant structure filename']
            # note: this current version assumes 1:1 for variant structures and
            #    variants...could change in the future...
    
            print '[[VIPURLOG]]generating Rosetta relax command on ' + variant_structure + '...'
            sys.stdout.flush()
            relax_commands = []
            if single_relax:
                relax_command , relax_filename = run_rosetta_relax_local( variant_structure , run = False )    # old method
#            relax_filename = relax_cmds[i][3]    # the silent file
        
                relax_commands.append( relax_command )
                relax_commands.append( relax_filename )
            else:
                # as separate runs, 1 per nstruct!
                relax_commands = []
                for j in xrange( nstruct ):
                    extra_options = {
                        'nstruct' : 1 ,
                        'run:jran' : jran + j ,
                        'out:file:silent' : rosetta_relax_options['out:file:silent']( variant_structure ).replace( '.' + file_extension + '.silent' , '_' + str( j + 1 ) + '.silent' ) ,
                        'out:file:scorefile' : rosetta_relax_options['out:file:scorefile']( variant_structure ).replace( '.' + file_extension + '.sc' , '_' + str( j + 1 ) + '.sc' )
                        }

                    relax_command , relax_filename = run_rosetta_relax_local( variant_structure , run = False , extra_options = extra_options )    # old method
                    
                    relax_commands.append( relax_command )
                    relax_commands.append( relax_filename )

            combined_silent_filename = rosetta_relax_options['out:file:silent']( variant_structure.replace( '.' + file_extension , '' ) )
            variants[i]['combined_silent_filename'] = combined_silent_filename
            combined_score_filename = rosetta_relax_options['out:file:scorefile']( variant_structure.replace( '.' + file_extension , '' ) )
            variants[i]['combined_score_filename'] = combined_score_filename

            # additional step, rescore with Rosetta
            # only required to add in 2 scores that SHOULD be output by default with relax
            # whatever, this is the easiest solution, just make an additional call
            #score_command , score_filename = run_rosetta_rescore( relax_filename , native_filename = pdb_filename , run = False )
            score_command , score_filename = run_rosetta_rescore( combined_silent_filename , native_filename = pdb_filename , run = False )

            # extract features, use the quartile method to extract comparisons
            variant_relax[i] = relax_commands + [score_command , score_filename]


    # write out a master summary file
    # root_filename| <root_filename>
    # out_path| # needed to create abs paths and cd for ddg_monomer (facepalms again)
    # options| chain: , sequence_only:
    # files| pdb_filename: , variants_filename: , prediction_filename: , sequence_filename:
    # (the rest of these can be multiple)
    # variant| name: , structure filename:
    # command| <command>

    summary_text = 'root_filename| ' + root_filename +'\n'

    if not out_path:
        # here we don't want it to be relative
        out_path = os.getcwd()
    summary_text += 'out_path| ' + out_path +'\n'

    summary_text += 'other| ' + 'target_chain:' + target_chain +','+ 'sequence_only:' + str( sequence_only ) +'\n'    # only these 2 for now...

    summary_text += 'files| ' + 'pdb_filename:' + pdb_filename +','+ 'variants_filename:' + variants_filename +','+ 'prediction_filename:' + prediction_filename +','+ 'sequence_filename:' + sequence_filename
    # lol...wtf
    summary_text += ',task_summary_filename:' + task_summary_filename
    if write_numbering_map:
        summary_text += ',numbering_map:' + root_filename +'.numbering_map'
    summary_text += '\n'

    # variants
    # well, lets have a record of these
    for i in failed_variants.keys():
        summary_text += 'variant| ' + 'name:' + root_filename +'_'+ i +','+ 'failed:' + failed_variants[i] +'\n'

    for i in variants.keys():
        summary_text += 'variant| ' + 'name:' + root_filename +'_'+ i
        # isn't really needed...for consistency? NO, should not make more complex for no reason...
        if not sequence_only:
            summary_text += ','+ 'structure_filename:' + variants[i]['variant structure filename']
            summary_text += ','+ 'combined_silent_filename:' + variants[i]['combined_silent_filename']
            summary_text += ','+ 'combined_score_filename:' + variants[i]['combined_score_filename']
        summary_text += '\n'

    # commands
    # psiblast, simplest, for the entire protein
    summary_text += 'command| ' + 'feature:psiblast' +','+ 'output_filename:' + psiblast_filename +','+ psiblast_command +'\n'
    if not sequence_only:
        # probe
        summary_text += 'command| ' + 'feature:probe' +','+ 'output_filename:' + probe_output_filename +','+ probe_command +'\n'
        summary_text += 'other| ' + 'probe_positions:' + ';'.join( [str( i ) for i in probe_positions] ) +'\n'
        # ddg_monomer
        summary_text += 'command| ' + 'feature:ddg_monomer' +','+ 'output_filename:' + ddg_monomer_out_filename +','+ ddg_monomer_command +'\n'
        # relax, ugh, these are variant specific

        if single_relax:
            summary_text += 'command| ' + 'feature:relax_native' +','+ 'output_filename:' + native_relax_filename +','+ 'variant:native,' + native_relax_command +'\n'
        else:
            # run as separate commands
            for j in xrange( 0 , len( native_relax_commands ) - 2 , 2 ):    # in pairs, skip the last 2, the rescore command
                summary_text += 'command| ' + 'feature:relax_native' +','+ 'output_filename:' + native_relax_commands[j + 1] +','+ 'variant:native,' + native_relax_commands[j] +'\n'
        summary_text += 'other| ' + 'combined_native_silent_filename:' + combined_native_silent_filename +'\n'
        summary_text += 'other| ' + 'combined_native_score_filename:' + combined_native_score_filename +'\n'
    
        summary_text += 'command| ' + 'feature:relax_native_rescore' +','+ 'output_filename:' + native_relax_commands[-1] +','+ 'variant:native,' + native_relax_commands[-2] +'\n'

        for i in variant_relax.keys():
            if len( variant_relax[i] ) == 4:
                # run as a batch
                summary_text += 'command| ' + 'feature:relax' +','+ 'output_filename:' + variant_relax[i][1] +','+ 'variant:' + i +','+ variant_relax[i][0] +'\n'
            else:
                # run as separate commands
                for j in xrange( 0 , len( variant_relax[i] ) - 2 , 2 ):    # in pairs, skip the last 2, the rescore command
                    summary_text += 'command| ' + 'feature:relax' +','+ 'output_filename:' + variant_relax[i][j + 1] +','+ 'variant:' + i +','+ variant_relax[i][j] +'\n'

            summary_text += 'command| ' + 'feature:relax_rescore' +','+ 'output_filename:' + variant_relax[i][-1] +','+ 'variant:' + i +','+ variant_relax[i][-2] +'\n'


    print '[[VIPURLOG]]writing task summary file to ' + task_summary_filename
    sys.stdout.flush()
    f = open( task_summary_filename , 'w' )
    f.write( summary_text.rstrip( '\n' ) )
    f.close()
    
    # testing   
#    load_task_summary( task_summary_filename )
    
#    os.chdir( current_dir )    # go back...
    
    return task_summary_filename


# need to communicate through simple text files
def load_task_summary( task_summary_filename ):
    # use "|" to determine greater variable grouping
    # remember rescore requires relax to be complete!!!
    f = open( task_summary_filename , 'r' )
    lines = [i.rstrip( '\n' ).split( '|' ) for i in f.xreadlines()]
    f.close()
    
    # sanity check
    what_lines = [i for i in lines if not len( i ) == 2]
    if what_lines:
        print what_lines
        raise Exception( 'something has gone wrong, \"|\" or \"\\n\" characters generated' )
    
    # specific variables
    root_filename = None
    out_path = None
    filenames = {}
    other = {}
    variants = {}
    commands = []
    for i in lines:
        if i[0] == 'root_filename':
            root_filename = i[1].strip()
        elif i[0] == 'out_path':
            out_path = i[1].strip()
        elif i[0] == 'other':
            other.update( dict( [j.split( ':' ) for j in i[1].strip().split( ',' )] ) )
        elif i[0] == 'files':
            filenames.update( dict( [j.split( ':' ) for j in i[1].strip().split( ',' )] ) )
        elif i[0] == 'variant':
            variant = dict( [j.split( ':' ) for j in i[1].strip().split( ',' )] )
            variants[variant['name']] = variant
        elif i[0] == 'command':
            command = i[1].strip().split( ',' )
            extra = dict( [j.split( ':' ) for j in command[:-1]] )
            # no, should be simple!
            # treat qsub command special, it uses ":" and ","
#            extra = [j.split( ':' ) for j in command[:-1]]
#            extra = dict( [(j[0] , ':'.join( j[1:] ))for j in extra] )
            extra['command'] = command[-1]
            commands.append( extra )

    # special cases
    if 'sequence_only' in other.keys():
        other['sequence_only'] = bool( other['sequence_only'] == 'True' )
    if 'probe_positions' in other.keys():
        other['probe_positions'] = other['probe_positions'].split( ';' )

    summary = {
        'root_filename' : root_filename , 
        'out_path' : out_path ,
        'other' : other ,
        'filenames' : filenames ,
        'variants' : variants ,
        'commands' : commands
        }
    return summary

# simple enough, wrapped for testing
def write_task_summary( task_summary , task_summary_filename ):
    summary_text = 'root_filename| ' + task_summary['root_filename'] +'\n'
    summary_text += 'out_path| ' + task_summary['out_path'] +'\n'

    summary_text += 'other| '
    for i in task_summary['other'].keys():
        if isinstance( task_summary['other'][i] , bool ):
            summary_text += i +':'+ str( task_summary['other'][i] ) +','
        elif isinstance( task_summary['other'][i] , str ):
            summary_text += i +':'+ task_summary['other'][i] +','
        elif isinstance( task_summary['other'][i] , list ):
            summary_text += i +':'+ ';'.join( task_summary['other'][i] ) +','
        else:
            raise Exception( 'not currently supported...' + i +' '+ str( task_summary['other'][i] ) )
    summary_text = summary_text.rstrip( ' ,' ) +'\n'

    # add "task_summary_filename" if it is absent?
    if not 'task_summary_filename' in task_summary['filenames'].keys():
        task_summary['filenames']['task_summary_filename'] = task_summary_filename    
    summary_text += 'files| ' + ','.join( [i +':'+ task_summary['filenames'][i] for i in task_summary['filenames'].keys()] ) +'\n'
    for i in task_summary['variants'].keys():
        # only add bridging "," if >1 keys
        summary_text += 'variant| ' + 'name:' + task_summary['variants'][i]['name'] +','*bool( len( task_summary['variants'][i] ) > 1 )+ ','.join( [j +':'+ task_summary['variants'][i][j] for j in task_summary['variants'][i].keys() if not j == 'name'] ) +'\n'
    for i in task_summary['commands']:
        summary_text += 'command| ' + ','.join( [j +':'+ i[j] for j in i.keys() if not j =='command'] ) +','+ i['command'] +'\n'

    # write it
    f = open( task_summary_filename , 'w' )
    f.write( summary_text.rstrip( '\n' ) )
    f.close()

#    return summary_text


