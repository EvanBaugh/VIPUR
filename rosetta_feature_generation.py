#!/usr/bin/env python
# :noTabs=true:

"""
ehb: 

#

Note: add documentation for release
"""

################################################################################
# IMPORT

# common modules
import os
import tempfile
from math import floor

# bigger modules

# custom modules
from vipur_settings import PATH_TO_ROSETTA_DDG_MONOMER , PATH_TO_ROSETTA_RELAX , PATH_TO_ROSETTA_SCORE , PATH_TO_PYMOL , USE_PYROSETTA , PATH_TO_VIPUR , ROSETTA_DDG_MONOMER_OPTIONS , ROSETTA_RELAX_OPTIONS , ROSETTA_SCORE_OPTIONS , ROSETTA_TERMS_TO_COMPARE , ROSETTA_RELAX_PARALLEL
from helper_methods import create_executable_str , run_local_commandline

################################################################################
# METHODS

# support PyRosetta OR PyMOL
def create_variant_protein_structures( pdb_filename , variants , chain , use_pyrosetta = USE_PYROSETTA , pymol_environment_setup = '' ):
    # optionally run environment setup
    if pymol_environment_setup:
        print 'setting up environment variables'
        run_local_commandline( pymol_environment_setup )

    # make sure the variants have been filtered
    if use_pyrosetta:
        # load the PDB as a pose
        pose = pose_from_pdb( pdb_filename )
        failed = {}
        root_filename = pdb_filename.rstrip( 'pdb' )
            
        # currently cannot handle multi-chain input
        # handle this before VIPUR
        if pose.chain( pose.total_residue() ) > 1:
            print 'CANNOT currently handle multi-chain PDBs (as pose), using PyMOL instead!'
            if PATH_TO_PYMOL:
                create_variant_protein_structures( pdb_filename , variants , chain , use_pyrosetta = False )
                return
            else:
                faulty = 'clean before VIPUR, cannot handle multi-chain PDBs'
                failed[faulty] = variants
        elif not pose.pdb_info().chain( 1 ) == chain:
            print '...not sure what it happening, you wanted chain ' + chain + ' but VIPUR found chain ' + pose.chain( 1 ) + ', skipping this entire sample!'
            faulty = 'clean before VIPUR, improper chain ID'
            failed[faulty] = variants
                
        # in case this condition is found:
        faulty = 'could not load position from PDB'
        variant_structures = []
        for variation in variants:
            # make a copy
            test_pose = Pose()
            test_pose.assign( pose )

            native = variation[0]
            position = variation[1:-1]
            mutant = variation[-1]
                
            # make sure the position was loaded
            icode = ' '    # default...this could cause problems...
            if not position[-1].isdigit():
                icode = position[-1]
                position = position[:-1]
                position = int( position )
            if not test_pose.pdb_info().pdb2pose( chain , position , icode ):
                if faulty in failed.keys():
                    failed[faulty].append( variation )
                else:
                    failed[faulty] = [variation]
                break    # stop the loop
            position = test_pose.pdb_info().pdb2pose( chain , position , icode )
                    
            # simple, use a mover to make the change
            # appears to have trouble with N terminal variants since it uses "replace_residue"
            # code is available that does not have this problem, however reloading into Rosetta with accurately determine the position of these atoms
            make_variant = MutateResidue( position , mutant )
            make_variant.apply( test_pose )

            # write out
            out_filename = self.root_filename +'.chain_'+ chain +'_'+ variation +'.pdb'
            variant_structures.append( out_filename )
            test_pose.dump_pdb( out_filename )
            print 'generated ' + variation + ' variant structure and wrote to ' + out_filename

        return variant_structures        
    else:
        # use the pymol script
        #for variants in self.variants['permissible']:
        # use default output naming
        # create command explicitly here, slightly different
        root_filename = pdb_filename.rstrip( '.pdb' )
        command = PATH_TO_PYMOL + ' -qcr ' + PATH_TO_VIPUR + '/pymol_make_variant_structure.py -- -p ' + pdb_filename + ' -m ' + ','.join( variants ) + ' -c ' + chain + ' -r ' + root_filename
#            print command
        if pymol_environment_setup:
            command = pymol_environment_setup +'\n\n'+ command
        run_local_commandline( command )
            
        # reconstruct the names
        variant_structures = [root_filename + '.chain_' + chain +'_'+ i +'.pdb' for i in variants]
        
        # verify they have been made
        if [None for i in variant_structures if not os.path.isfile( i )]:
            raise IOError( 'could not make variant protein structures,\ntry checking the input PDB file or the pymol script pymol_make_variant_structure.py' )
        
        return variant_structures

#############
# DDG_MONOMER

# simple format writing helper
def write_mut_file( variants , residue_map , mut_filename ):
    # ...any way to hangle multiple mutants?
    # need to consider multiple backgrounds, pairwise variants with different references
    text = 'total ' + str( len( variants ) ) +'\n'
    # must split by positions
        
    for i in variants:
        text += '1\n'    # hardcoded for now...
        pdb_position = i[1:-1]
        # 1-indexed, not 0-indexed
        # assumes that the sequence filtering is proper...may have problems on erroneously formatted PDBs
        pose_position = str( residue_map[pdb_position] + 1 )
                
        text += ' '.join( [i[0] , pose_position , i[-1]] ) +'\n'
            
        # ...unfamiliar with the "mut" file syntax...because documentation is poor

    # write it out
    f = open( mut_filename , 'w' )
    f.write( text.rstrip( '\n' ) )
    f.close()

# local
def run_rosetta_ddg_monomer( pdb_filename , mut_filename , out_filename = '' , out_path = '' , cleanup = True , run = True ):
    root_filename = os.path.abspath( pdb_filename ).rstrip( '.pdb' )
    # hardcoded...ddg_monomer is such a painful protocol...
    out_filename = ''
    if '/' in root_filename:
        out_filename += '/'.join( root_filename.split( '/' )[:-1] ) +'/'
    out_filename += 'ddg_predictions.out'
    # clear it out if it exists, otherwise it will be appended to...
    if os.path.exists( out_filename ):
        os.remove( out_filename )

    # collect the options, set the input, derive the output filenames
    ddg_monomer_options = {}
    ddg_monomer_options.update( ROSETTA_DDG_MONOMER_OPTIONS )
    ddg_monomer_options['in:file:s'] = pdb_filename
    ddg_monomer_options['ddg::mut_file'] = mut_filename
    for i in ddg_monomer_options.keys():
        if '__call__' in dir( ddg_monomer_options[i] ):
            ddg_monomer_options[i] = ddg_monomer_options[i]( root_filename )

    for i in ddg_monomer_options.keys():
        if isinstance( ddg_monomer_options[i] , str ) and os.path.isfile( ddg_monomer_options[i] ):
            ddg_monomer_options[i] = os.path.abspath( ddg_monomer_options[i] )
    
    command = ''
    # optionally move into the specific directory...
    if out_path:
        command += 'cd '+ out_path +'\n\n'
    
    command += create_executable_str( PATH_TO_ROSETTA_DDG_MONOMER , args = [] , options = ddg_monomer_options )

    if run:
        run_local_commandline( command )
    
        # optionally cleanup
        if cleanup:
            print 'ddg_monomer writes useless output files, deleting these now...'
            remove_intermediate_ddg_monomer_files()
        
        # the only output we need
        return out_filename
    else:
        return command , out_filename

# simple helper
def remove_intermediate_ddg_monomer_files():
    for i in os.listdir( '.' ):
        if i == 'wt_traj' or 'mutant_traj' == i[:11]:
            os.remove( i )

# simple, for now just check if empty or not
def check_ddg_monomer_output( ddg_monomer_output_filename ):
    # simple enough, for now just check if empty
    f = open( ddg_monomer_output_filename , 'r' )
    success = bool( f.read().strip() )    # load all of this!?
    f.close()
    
    # use the extract method? check if match desired positions?
    
    return success

# extract the score terms from the output and setup to match with residue numbers
def extract_score_terms_from_ddg_monomer( out_filename = 'ddg_predictions.out' , prefix = 'ddG:' ):
    f = open( out_filename , 'r' )
    lines = [i for i in f.xreadlines() if i.strip()]
    f.close()
    
    parse_ddg_monomer_line = lambda line : [i.strip() for i in line.lstrip( prefix ).split( ' ' ) if i.strip()]
    # added "ddg_" for legacy compatability, is artibrary, make more informative

    ddg_monomer_dict = [parse_ddg_monomer_line( i ) for i in lines]
    # will add a "description" entry for the header
    ddg_monomer_dict = dict( [(i[0] , i[1:]) for i in ddg_monomer_dict] )
    
    return ddg_monomer_dict


#######
# RELAX

# modified by njc
def run_rosetta_relax( pdb_filename , extra_options = {} , run = True , parallel = ROSETTA_RELAX_PARALLEL ):
    root_filename = pdb_filename.rstrip( '.pdb' )
    
    # collect the options, set the input, derive the output filenames
    relax_options = {}
    relax_options.update( ROSETTA_RELAX_OPTIONS )
    relax_options.update( extra_options )
    relax_options['s'] = pdb_filename
    relax_options['native'] = pdb_filename    # required to get gdtmm scores
    for i in relax_options.keys():
        if '__call__' in dir( relax_options[i] ):
            relax_options[i] = relax_options[i]( root_filename )

    # ...weird Rosetta append behavior...
#    if os.path.isfile( relax_options['out:file:silent'] ):
#        os.remove( relax_options['out:file:silent'] )
#    if os.path.isfile( relax_options['out:file:scorefile'] ):
#        os.remove( relax_options['out:file:scorefile'] )


    # for njc parallelization
    nstruct = int( relax_options.get( 'nstruct' , '0' ) )
    parallel = int( parallel )
    tmp_file = None
    if nstruct > 1 and parallel > 1:
        relax_options['nstruct'] = 1 #TODO: Add chunking option?
        score_filename = relax_options['out:file:scorefile']
        silent_filename = relax_options['out:file:silent']

        if 'run:jran' in relax_options:
            restoreJran = True
            jran = int( relax_options['run:jran'] )
        else:
            restoreJran = False
            jran = 123

        tmp_file = tempfile.NamedTemporaryFile( delete = False )
        print 'Parallel relax commands are in ' + tmp_file.name

        for s in xrange( nstruct ):
            tag = '_%05d' % s
            relax_options['run:jran'] = jran*nstruct + s
            relax_options['out:file:scorefile'] = score_filename + tag
            relax_options['out:file:silent'] = silent_filename + tag
            print >>tmp_file , create_executable_str( PATH_TO_ROSETTA_RELAX , args = [] , options = relax_options ) + " > %s 2>&1; echo '[[VIPURLOG]]' %s %d" % ((silent_filename + tag).replace( 'silent_' , 'log_' ) , pdb_filename , s + 1 )

        tmp_file.close()
        # the "find ... | xargs ..." idiom is used just in case nstruct is ever a *very* large number.
        command = '''\
parallel -j %d -a %s
find . -name '%s_[0-9]*[0-9]' | xargs cat | awk 'NR == 1 || $2 != "score" {print $0}' > %s
find . -name '%s_[0-9]*[0-9]' | xargs rm
find . -name '%s_[0-9]*[0-9]' | xargs cat | awk 'NR <= 2 || !($2 == "score" || $1 == "SEQUENCE:") {print $0}' > %s
find . -name '%s_[0-9]*[0-9]' | xargs rm
''' % (parallel , tmp_file.name , score_filename , score_filename , score_filename , silent_filename , silent_filename , silent_filename)
        print 'Parallel relax driver command:', command

        # restore option values
        relax_options['nstruct'] = str( nstruct )
        relax_options['out:file:scorefile'] = score_filename
        relax_options['out:file:silent'] = silent_filename
        if restoreJran:
            relax_options['run:jran'] = jran

        if run:
            return (command , tmp_file.name , score_filename , silent_filename)

        if tmp_file:
            os.unlink( tmp_file.name )
    else:
        command = create_executable_str( PATH_TO_ROSETTA_RELAX , args = [] , options = relax_options )

    if run:
        run_local_commandline( command )
    
#    command = create_executable_str( PATH_TO_ROSETTA_RELAX , args = [] , options = relax_options )

#    run_local_commandline( command )
    
    # the only output we need
#    return relax_options['out:file:scorefile']
    return relax_options['out:file:silent']

# local
def run_rosetta_relax_local( pdb_filename , extra_options = {} , run = True ):
    root_filename = os.path.abspath( pdb_filename ).replace( '.pdb' , '' )
    
    # collect the options, set the input, derive the output filenames
    relax_options = {}
    relax_options.update( ROSETTA_RELAX_OPTIONS )
    relax_options.update( extra_options )
    relax_options['s'] = pdb_filename
    relax_options['native'] = pdb_filename    # required to get gdtmm scores
    for i in relax_options.keys():
        if '__call__' in dir( relax_options[i] ):
            relax_options[i] = relax_options[i]( root_filename )

    for i in relax_options.keys():
        if isinstance( relax_options[i] , str ) and os.path.isfile( relax_options[i] ):
            relax_options[i] = os.path.abspath( relax_options[i] )

    # ...weird Rosetta append behavior...
    if os.path.isfile( relax_options['out:file:silent'] ):
        os.remove( relax_options['out:file:silent'] )
    if os.path.isfile( relax_options['out:file:scorefile'] ):
        os.remove( relax_options['out:file:scorefile'] )
    
    command = create_executable_str( PATH_TO_ROSETTA_RELAX , args = [] , options = relax_options )

    if run:
        run_local_commandline( command )
    
        # the only output we need
        return relax_options['out:file:silent']
    else:
        return command , relax_options['out:file:silent']

# simple, for now just check if empty or not
def check_relax_output( relax_score_filename , target_number_of_trajectories = ROSETTA_RELAX_OPTIONS['nstruct'] , header_lines = 1 , single_relax = False ):
    # simple enough, for now just check if empty
    f = open( relax_score_filename , 'r' )
    trajectories = len( f.readlines() ) - header_lines
    f.close()
    
    # optionally split into individual jobs
    if not single_relax:
        target_number_of_trajectories = 1
    
    # use the silent file instead? so much bulkier...
    success = (trajectories == int( target_number_of_trajectories ))
    
    return success , trajectories

# crude crude method
def merge_rosetta_relax_output( silent_filenames , combined_silent_filename , score_filenames , combined_score_filename , delete_old_files = False ):
    # ??? just combine all the text
    text = ''
    for i in silent_filenames:
        f = open( i , 'r' )
        new_text = f.read()
        f.close()
        
        if text and not text[-1] == '\n':
            text += '\n'

        text += new_text

    f = open( combined_silent_filename , 'w' )
    f.write( text )
    f.close()
    
    # optionally delete the old files
    if delete_old_files:
        for i in silent_filenames:
            os.remove( i )    # should all be abspath files...

    text = ''
    for i in score_filenames:
        f = open( i , 'r' )
        new_text = f.readlines()
        f.close()
        
        # remove headers, unless its the first one
        if text:
            new_text = new_text[1:]
        new_text = ''.join( new_text )
        
        # 'glue' with newlines
        if text and not text[-1] == '\n':
            text += '\n'

        text += new_text

    f = open( combined_score_filename , 'w' )
    f.write( text )
    f.close()

    # optionally delete the old files
    if delete_old_files:
        for i in score_filenames:
            os.remove( i )    # should all be abspath files...


#######
# SCORE

# added to make sure additional scores are in the output file
def run_rosetta_rescore( silent_filename , native_filename , score_filename = '' , run = True ):
    """
    Performs extraction of individual PDB structures from  <silent_filename>
    to  <out_dir>  (default to current location) using the "score" protocol
    of Rosetta (built against 3.5)
    
    Optionally specify  <extra_options>
    """
    root_filename = os.path.abspath( silent_filename ).rstrip( '.silent' )
    
    score_options = {}
    score_options.update( ROSETTA_SCORE_OPTIONS )
    score_options['in:file:silent'] = silent_filename
    score_options['in:file:native'] = native_filename    # required to get gdtmm scores
    for i in score_options.keys():
        if '__call__' in dir( score_options[i] ):
            score_options[i] = score_options[i]( root_filename )

    # necessary...
    if 'out:file:scorefile' in score_options.keys() and not 'rescore.sc' in score_options['out:file:scorefile']:
        score_options['out:file:scorefile'] = score_options['out:file:scorefile'].replace( '.sc' , '_rescore.sc' )

    for i in score_options.keys():
        if isinstance( score_options[i] , str ) and os.path.isfile( score_options[i] ):
            score_options[i] = os.path.abspath( score_options[i] )

    # ...weird Rosetta append behavior...
    if os.path.isfile( score_options['out:file:scorefile'] ):
        os.remove( score_options['out:file:scorefile'] )
        
    # default options
    command = create_executable_str( PATH_TO_ROSETTA_SCORE , args = [] , options = score_options )

    if run:
        run_local_commandline( command )
    
        return score_options['out:file:scorefile']
    else:
        return command , score_options['out:file:scorefile']


####################
# FEATURE EXTRACTION

# makes a dict summarizing the scores in the scorefile, divided by score term (column in the scorefile)
def extract_scores_from_scorefile( scorefilename , header = 0 , hit = 'SCORE: ' , as_float = ROSETTA_TERMS_TO_COMPARE ):
    # load it
    f = open( scorefilename , 'r' )
    lines = f.readlines()
    f.close()

    # find the score terms, as dict for easier parsing
    score_terms = [i.strip() for i in lines[header].replace( hit , '' ).split( ' ' ) if i.strip()]
    scores = dict( [(i , []) for i in score_terms] )

    # load the terms on each line
    for i in lines[header + 1:]:
        # skip lines that are not "hits"
        if not i[:len( hit )] == hit:
            continue

        # split into the columns
        line = [j.strip() for j in i.replace( hit , '' ).split( ' ' ) if j.strip()]

        if not len( line ) == len( scores.keys() ):
            raise IOError( '??!? wrong number of columns (' + str( len( line ) ) + ', should be ' + str( len( score_terms ) ) + ') found !!?!\n\n' + i )

        # add values to each distribution
        for j in xrange( len( line ) ):
            # convert the specified terms to float
            if score_terms[j] in as_float:
                line[j] = float( line[j] )
            scores[score_terms[j]].append( line[j] )

    return scores

# find or calculate the value of the xth quartile e.g. Q2=.5 on a distribution (the value at Q2)
def determine_quartile_value( quartile , distribution , tolerance = 1e-7 ):
    # sort, just in case
    distribution.sort()
    
    total = len( distribution ) - 1    # largest index possible
    
    if quartile <= 0:
        # return the smallest value observed
        quartile_value = distribution[0]
    elif quartile > 1:
        # return the largest value observed
        quartile_value = distribution[-1]
    else:
        # interpolate:
        # determine the index of the closest value (under) the target
        # scale the remaining desired quartile value by the difference bewteen the observed values (interpolation)
        base_index = quartile*total    # actually the "perceived index" for now
        remaining = base_index - floor( base_index )    # on "quartile scale"
        if remaining < tolerance:
            # well, if no remainder, just take what we found!
            quartile_value = distribution[int( base_index )]
        else:
            base_index = int( floor( base_index ) )    # now its the base index, the closest value no exceeding the target
            quartile_value = distribution[base_index] + remaining*( distribution[base_index + 1] - distribution[base_index] )

    return quartile_value

# find or calculate the quartile corresponding to a particular value on distribution e.g. what quartile is 290?
def determine_quartile( quartile_value , distribution , tolerance = 1e-7 ):
    # sort just in case
    distribution.sort()

    total = len( distribution ) - 1
    
    if quartile_value <= distribution[0]:
        quartile = 0    # too low, not even in the range
    elif quartile_value >= distribution[-1]:
        quartile = total    # too big, not even in the range
    else:
        # find the exact value in the distribution
        base_index = [i for i in xrange( len( distribution ) ) if distribution[i] >= quartile_value][0]
        # or, the closest value greater than the target
        if abs( distribution[base_index] - quartile_value ) <= tolerance:
            quartile = base_index    # note, adapted from code that steps down from here, however this is ALREADY done by the search above
        # linearly interpolate using the two closest points
        # calculate what fraction of the difference to the target (from the closest lower point) scaled by the difference between adjacent points
        else:
            quartile = base_index - 1 + float( quartile_value - distribution[base_index - 1] )/( distribution[base_index] - distribution[base_index - 1] )

    # scale by the maximum value
    quartile = float( quartile )/total
    
    return quartile

# helper for mapping between distributions
def determine_quartile_from_quartile_value_of_another_distribution( quartile , query_distribution , reference_distribution ):
    # get the quartile value on the query distribution
    quartile_value = determine_quartile_value( quartile , query_distribution )
    
    # place it on the reference distribution
    quartile2 = determine_quartile( quartile_value , reference_distribution )
    
    return quartile2

# loads both score files and compares the distributions
def extract_quartile_score_terms_from_scorefiles( variant_distribution , native_distribution , quartiles = {'Q1' : .25 , 'Q2' : .5 , 'Q3' : .75} , terms = ROSETTA_TERMS_TO_COMPARE ):
    # if lazy and input score filename
    if isinstance( variant_distribution , str ):
        variant_distribution = extract_scores_from_scorefile( variant_distribution )
    if isinstance( native_distribution , str ):
        native_distribution = extract_scores_from_scorefile( native_distribution )
    
    # add as unique_terms
    quartile_comparisons = {}
    
    for term in terms:
        for quartile in quartiles.keys():
            # using the legacy names of these features
            # rename these to be less cumbersome
            feature_name = 'quartile_' + term + quartile
            quartile_comparisons[feature_name] = determine_quartile_from_quartile_value_of_another_distribution( quartiles[quartile] , variant_distribution[term] , native_distribution[term] )

    return quartile_comparisons


