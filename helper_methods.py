#!/usr/bin/env python
# :noTabs=true:

"""
Simple "helper methods" used by VIPUR
wrappers for file and directory manipulation and running programs with subprocess

includes methods for loading and quality-checking variants

Currently other pre-processing methods that are feature-specific, such as variant structure generation
are contained within the script associated with those features
"""

################################################################################
# IMPORT

# common modules
import os
import shutil
import subprocess

# bigger modules

# custom modules
from vipur_settings import PROTEIN_LETTERS

################################################################################
# METHODS

get_file_extension = lambda in_filename: in_filename.split( '.' )[-1]
get_file_extension.__doc__ = 'Returns the file extension of  <in_filename>\n\nin_filename.split( \'.\' )[-1]'

# hacky version
get_root_filename = lambda in_filename: in_filename[:-len( get_file_extension( in_filename ) ) - 1]
get_root_filename.__doc__ = 'Returns the \"root filename\" of  <in_filename>  (pre file extension)\n\nin_filename[:len( in_filename.split( \'.\' )[-1] ) - 1]\na little hacky...'
# better version
#get_root_filename = lambda in_filename: ''.join( [i for i in in_filename.split( '.' )[:-1]] )


# helper for creating a directory, checks and delets existing name
def create_directory( dir_name , tagline = ' to sort the data' ):
    """
    Creates the directory  <dir_name>
    
    WARNING: this will delete the directory and its contents if it already
    exists!
    
    Optionally output something special in  <tagline>
    """
    # check if it exists
    print 'Creating a new directory ' + os.path.relpath( dir_name ) + tagline
    if os.path.isdir( dir_name ):
        print 'a directory named ' + os.path.relpath( dir_name ) + ' already exists, deleting it now...'
        shutil.rmtree( dir_name )
    os.mkdir( dir_name )

# copy helper
def copy_file( filename , destination , display = False ):
    """
    Copy  <filename>  to/into  <destination>
    
    just a cp wrapper...what?
    """
    if display:    # optional
        if os.path.isdir( destination ):
            print 'placing a copy of ' + os.path.relpath( filename ) + ' into the ' + os.path.relpath( destination ) + ' directory'
        elif os.path.isfile( destination ):
            print 'copying ' + os.path.relpath( filename ) + ' to ' + os.path.relpath( destination )
    shutil.copy( filename , destination )


#####################
# subprocess wrappers

# helper for making bash scripts
def create_executable_str( executable , args = [] , options = {} , out_filename = '' , args_first = False , extra_options = '' , append = False ):
    """
    Returns a str for a commandline call to  <executable>  (input as a string)
    using a list og  <args>  and a dict of  <options>  (keys as options, values
    as option values, empty str for flags)
    
    Optionally write stout to  <out_filename>  (only if specified)
    Optionally provide the arguments to the method first  (<args_first>)
    
    Add additional text  <extra_options>  to the call
    """
    if not isinstance( executable , str ):
        raise IOError( 'variable  <executable>  must be a string!' )
    if not isinstance( args , list ) and not isinstance( args , tuple ):
        raise IOError( 'variable  <args>  must be a list! (or otherwise iterable)' )
    if not isinstance( options , dict ):
        raise IOError( 'variable  <executable>  must be a string!' )
        
    # setup the options and args
    options = ''.join( [' ' + ('-')*bool( isinstance( i , str ) and len( i ) > 0 ) + str( i ) + ( ' ' + str( options[i] ) )*bool( options[i] ) for i in options.keys()] )
    args = ''.join( [' ' + str( i ) for i in args] )
    
    # choose general format
    if args_first:
        perform = executable + args + options + extra_options
    else:
        perform = executable + options + args + extra_options
    
    # optionally write stdout
    if out_filename:
        if append:
            perform += ' >> ' + out_filename
        else:
            perform += ' > ' + out_filename

    return perform

# runs a commandline, usually combined with create_executable_str above
def run_local_commandline( command , collect_stdout = False ):
    """
    Runs the explicit  <command>  by performing a system call with subprocess
    """
    # get the output
    print '\n'+ '='*80 + '\nPerforming system call:\n' + command + '\n' + '='*80 +'\n'
    if not collect_stdout:
        subprocess.call( command , shell = True )    # just the command, no output piping

    # older call that pipes the output into Python for manipulation
    else:
        #stdout = subprocess.Popen( command , shell = True , stdout = subprocess.PIPE , stdin = subprocess.PIPE , stderr = subprocess.STDOUT ).communicate()[0].strip()
        # shell = True, do NOT .split the command, = False, DO .split the command
        stdout = subprocess.Popen( command , shell = True , stdout = subprocess.PIPE , stdin = subprocess.PIPE , stderr = subprocess.STDOUT ).communicate()[0]#.strip()    # for consistency
        return stdout

#######################
# preprocessing methods

# simple file loading
def load_variants_file( variants_filename ):
    """
    Returns a list of the variants found in  <variants_filename>
    the variants are currently supported only as a list (successive lines) of
    variations/mutations in the form:
        native amino acid + protein position + variant amino acid
    
    ex.
        E17H
        E17G
        F42Y
    """
    f = open( variants_filename , 'r' )
    variants = f.read().strip( ',\n' )
    f.close()

    # "\n" delimited variants, each variant is on a new line
    variants = variants.split( '\n' )

    return variants

# scan for faulty mutations, input can be messy
def check_variants( variants , sequence , residue_map , protein_letters = PROTEIN_LETTERS ):
    """
    Returns a list of variants contained in  <variants>  which are properly
    formatted
    Requires the  <sequence>  of the native protein and a  <residue_map>
    converting the PDB residue numbers into the sequence positions*
    
    currently checks for:
        improper amino acids (not in  <protein_letters>)
        positions not in the structure (not in  <residue_map>)
        identical amino acids (e.g. not an actual variant)
        improper native amino acids (do not match the position in  <sequence>)

    *Note: since many protein structures (in PDB format) map to only a region
        of the protein of interest, VIPUR will often be run only on this region
        the  <residue_map>  ensures the actual protein positions are maintained
        during analysis (since many of these programs would consider the start
        of the modeled region as the first residue in the protein)
    """
    permissible_variants = []
    failed_variants = {}
    for variation in variants:
        native = variation[0]
        position = variation[1:-1]
        mutant = variation[-1]
            
        # check for nonsense variation
        if not mutant in protein_letters or not native in protein_letters:
            exit_message = '\"' + variation + '\" is improperly formatted; unaccepted amino acid single letter code'
            print exit_message
            failed_variants[variation] = exit_message
            continue
        
        # check the position
        if not position in residue_map.keys():
            exit_message = 'residue \"' + position + '\" is not in the PDB; skipping \"' + variation + '\"'
            print exit_message
            failed_variants[variation] = exit_message
            continue

        # check for redundancy
        if native == mutant:
            exit_message = 'native and variant are identical (?); skipping \"' + variation + '\"'
            print exit_message
            failed_variants[variation] = exit_message
            continue
        
        # check for native errors
        pdb_aa = sequence[residue_map[position]]
        if not native == pdb_aa:
            exit_message = 'wrong native amino acid! found \"' + pdb_aa + '\" in the pdb; not \"' + native + '\"'
            print exit_message
            failed_variants[variation] = exit_message
            continue
        
        # it passed!
        permissible_variants.append( variation )

    return permissible_variants , failed_variants

# extract the chains from the PDB
# only considers ATOM lines, mainly for use with clean_nucleic_acid_lines_from_pdb
def extract_chains_from_pdb( pdb_filename , only = ['ATOM'] ):
    """
    Returns the chains found in  <pdb_filename>
    
    Only consider lines starting with  <only>
    """
    pdb_filename = os.path.abspath( pdb_filename )
    if os.path.exists( pdb_filename ):
        # load the data
        f = open( pdb_filename , 'r' )
        data = [i for i in f.xreadlines() if i[:6].strip() in only]
        f.close()
    
        # find unique chains
        chains = []
        for i in data:
            if not i[21] in chains:
                chains.append( i[21] )

        return chains
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False


