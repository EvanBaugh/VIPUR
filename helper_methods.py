#!/usr/bin/env python
# :noTabs=true:

"""
Simple "helper methods" used by VIPUR
wrappers for file and directory manipulation and running programs with subprocess
"""

################################################################################
# IMPORT

# common modules
import os
import shutil
import subprocess

# bigger modules

# custom modules

################################################################################
# METHODS

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
    options = ''.join( [' -' + str( i ) + ( ' ' + str( options[i] ) )*bool( options[i] ) for i in options.keys()] )
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
def run_local_commandline( command ):
    """
    Runs the explicit  <command>  by performing a system call with subprocess
    """
    # get the output
    print '\n'+ '='*80 + '\nPerforming system call:\n' + command + '\n' + '='*80 +'\n'
    subprocess.call( command , shell = True )    # just the command, no output piping

    # older call that pipes the output into Python for manipulation
    #stdout = subprocess.Popen( command , shell = True , stdout = subprocess.PIPE , stdin = subprocess.PIPE , stderr = subprocess.STDOUT ).communicate()[0].strip()


