#!/usr/bin/env python
# :noTabs=true:

"""
Methods for preparing input for running PROBE and parsing the output
for VIPUR, PROBE is only run on the relevant protein positions (variant
positions) so additional information on the variants is required to map these
ACCP values to protein positions

Note: the run options for PROBE are setup in run_probe and set in
settings.py in PROBE_OPTIONS
"""

################################################################################
# IMPORT

# common modules
import os

# bigger modules

# custom modules
from vipur_settings import AMINO_ACID_CODES , PATH_TO_PROBE , PROBE_OPTIONS
from helper_methods import create_executable_str , run_local_commandline

################################################################################
# METHODS

# local
def run_probe( pdb_filename , variants , probe_output_filename = '' , run = True ):
    """
    Runs PROBE on  <pdb_filename>  on the positions found among  <variants>
    using the default options in PROBE_OPTIONS and writes the output to
    <probe_output_filename>  (also returns this output filename)
    """
    if not probe_output_filename:
        probe_output_filename = os.path.abspath( pdb_filename ).rstrip( '.pdb' ) + '.probe_out'

    # get the unique variant positions
    positions = list( set( [i[1:-1] for i in variants] ) )
    positions.sort()
    
    # generate the commands to run
#    command = '#!/bin/sh\nrm ' + probe_output_filename + '\ntouch ' + probe_output_filename + '\n'
    command = 'rm ' + probe_output_filename + ';touch ' + probe_output_filename + ';'
    # delete any prior copy since we will append to it
    
    for i in positions:
        probe_options = {}
        probe_options.update( PROBE_OPTIONS )
            
        probe_options['out'] = pdb_filename
        probe_options['Q'] = str( i )

        command += create_executable_str( PATH_TO_PROBE , [] , probe_options , probe_output_filename , append = True ) +';'#'\n'

    # run PROBE, store the output
    if run:
        run_local_commandline( command )

        return probe_output_filename , positions
    else:
        # the command, well, get positions etc. too
        return command , probe_output_filename , positions

# simple, for now just check if empty or not
def check_probe_output( probe_output_filename ):
    # simple enough, for now just check if empty
    f = open( probe_output_filename , 'r' )
    success = bool( f.read().strip() )    # load all of this!?
    f.close()
    
    # use the extract method? check is match desired positions?
    
    return success

# simple parsing, only run PROBE on variant positions, rank ordered
def extract_accp_from_probe( probe_output_filename ):
    """
    Returns a list of the  ACCP (Protein ACCessible surface area) values
    found in  <probe_output_filename>
    
    Note: VIPUR only runs PROBE on the relevant protein positions (variant
    positions) so additional information on the variants is required to map
    these ACCP values to protein positions
    """
    f = open( probe_output_filename , 'r' )
    lines = f.readlines()
    f.close()
        
    # ...or keep it simple...
#    sasa = [float( i.split( ':' )[-1].strip().split( ' ' )[0] ) for i in lines if 'accessible surface area' in i]
    contact_area = [float( i.split( ':' )[-1].strip().split( ' ' )[0] ) for i in lines if 'contact surface area' in i]
    potential_area = [float( i.split( ':' )[-1].strip().split( ' ' )[0] ) for i in lines if 'potential area' in i]
        
    # inhereted, unsure if this is the "best" way to do this
    accp = [round( contact_area[i]/potential_area[i]*100 , 2 ) for i in xrange( len( contact_area ) )]

    return accp


