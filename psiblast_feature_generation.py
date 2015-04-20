#!/usr/bin/env python
# :noTabs=true:

"""
Methods for preparing input for PSIBLAST (sequence file), running PSIBLAST,
and parsing the output PSSM

Note: the run options for PSIBLAST are setup in run_psiblast and set in
settings.py in PSIBLAST_OPTIONS
"""

################################################################################
# IMPORT

# common modules

# bigger modules

# custom modules
from settings import AMINO_ACID_CODES , PATH_TO_PSIBLAST , PSIBLAST_OPTIONS , PROTEIN_LETTERS
from helper_methods import create_executable_str , run_local_commandline

################################################################################
# METHODS

# for sequence searches
def extract_protein_sequence_from_pdb( pdb_filename , out_filename = '' , target_chain = 'A' , write_numbering_map = True ):
    """
    Returns the protein sequence found in  <pdb_filename>  matching
    <target_chain>  and writes the sequence (in FASTA format) to  <out_filename>
    
    Optionally  <write_numbering_map>  that maps the PDB residue numbering to
        the sequence (0-indexed) numbering*
    
    Note: although the sequence file is only required for input to PSIBLAST,
        the sequence and residue numbering map are useful for identifying
        improper/malformed variant inputs

    *Note: since many protein structures (in PDB format) map to only a region
        of the protein of interest, VIPUR will often be run only on this region
        the  <residue_map>  ensures the actual protein positions are maintained
        during analysis (since many of these programs would consider the start
        of the modeled region as the first residue in the protein)
    """
    # default output
    root_filename = pdb_filename.rstrip( '.pdb' )
    if not out_filename:
        out_filename = root_filename + '.fa'

    # load the ATOM lines
    f = open( pdb_filename , 'r' )
    lines = [i for i in f.xreadlines() if i[:4] == 'ATOM']
    f.close()
        
    sequence = ''
    residues = []
    last_resi = None
    for i in lines:
        if not i[21:22] == target_chain:
            continue
        
        resi = i[22:27].strip()
        if not resi == last_resi:
            # a new residue, add the amino acid
            sequence += AMINO_ACID_CODES[i[17:20]]
            residues.append( resi )
#            print resi
            last_resi = resi
    
    print 'extracted protein sequence from ' + pdb_filename + ' chain ' + target_chain
        
    # optionally write the numbering map
    if write_numbering_map:
        # allow it to be a str
        if isinstance( write_numbering_map , str ):
            numbering_map_filename = write_numbering_map
        else:
            numbering_map_filename = pdb_filename.rstrip( '.pdb' ) + '.numbering_map'

        f = open( numbering_map_filename , 'w' )
        f.write( '\n'.join( ['\t'.join( [residues[i] , str( i ) , sequence[i]] ) for i in xrange( len( residues ) )] ) )
        f.close()

        print 'wrote the numbering for ' + pdb_filename + ' chain ' + target_chain + ' to ' + numbering_map_filename
    
    # convert into a dict
    residues = dict( [(residues[i] , i) for i in xrange( len( residues ) )] )

    # "optionally" write out the sequence
    # currently a default, not an option - will determine an appropriate  <out_filename>  if not provided as an argument
    if out_filename:
        # write it to output directory
        f = open( out_filename , 'w' )
        f.write( '>' + root_filename +'_chain_'+ target_chain +'\n'+ sequence )
        f.close()
        print 'wrote the sequence for ' + pdb_filename + ' chain ' + target_chain + ' to ' + out_filename

    return sequence , residues , out_filename

# contingency method used if run in "sequence only" mode
def load_fasta( fasta_filename ):
    f = open( fasta_filename , 'r' )
    lines = f.readlines()
    f.close()
    
    sequences = []
    for i in lines:
        if i[0] == '>':
            sequences.append( [i.lstrip( '>' ) , ''] )
        else:
            sequences[-1][1] += i.strip()    # remove other characters
    
    # summarize mutliple sequences?
    return sequences

# local
def run_psiblast( sequence_filename ):
    """
    Runs PSIBLAST on  <sequence_filename>  using the default options in
    PSIBLAST_OPTIONS and returns the relevant output file: "out_ascii_pssm"
    """
    root_filename = sequence_filename.rstrip( '.fa' )
    
    # collect the options, set the input, derive the output filenames
    psiblast_options = {}
    psiblast_options.update( PSIBLAST_OPTIONS )
    psiblast_options['query'] = sequence_filename
    for i in psiblast_options.keys():
        if '__call__' in dir( psiblast_options[i] ):
            psiblast_options[i] = psiblast_options[i]( root_filename )
    
    command = create_executable_str( PATH_TO_PSIBLAST , args = [] , options = psiblast_options )

    run_local_commandline( command )
    
    # the only output we need
    return psiblast_options['out_ascii_pssm']

# copied parsing method
def extract_pssm_from_psiblast_pssm( pssm_filename , headers = 3 , trailers = 7 , columns = len( PROTEIN_LETTERS ) , first_columns_width = 3 , second_columns_width = 4 ):
    """
    Returns a dict summarizing the contents of PSIBLAST output  <pssm_filename>
    
    several variables help determines the file format, indicating the number of
    <headers>  and  <trailers>  (lines), the number of  <columns>  and the
    <first_columns_width>  and  <second_columns_width>  (PSIBLAST PSSM output
    is effectively two columns for each amino acid, however they are written
    "in bulk" for human readability)
    """
    # load the raw text
    f = open( pssm_filename , 'r' )
    lines = [i.strip( '\n' ) for i in f.xreadlines()]
    f.close()

    # parse the header
    header = lines[2][9:]    # empty spacing
    # cheating, assume empty space rather than counting characters
    first_columns_keys = [i.strip() for i in header[:first_columns_width*columns].split( ' ' ) if i.strip()]
    second_columns_keys = [i.strip() for i in header[first_columns_width*columns:].split( ' ' ) if i.strip()]    # nothing else

    # remove useless
    lines = lines[headers:-1*trailers]

    pssm_dict = {}
    for i in lines:
        # look for basic info
        front = [j for j in i[:9].split( ' ' ) if j.strip()]
        position = int( front[0] )
        query_identity = front[1]
        back = i[9:]

        # first set
        first_columns_values = [int( back[j*first_columns_width:(j + 1)*first_columns_width].strip() ) for j in xrange( columns )]
        second_columns_values = [float( back[first_columns_width*columns + 1 + j*second_columns_width:first_columns_width*columns + 1 + (j + 1)*second_columns_width].strip() )/100 for j in xrange( columns )]    # as "percentages"...arg!

        back = [float( j ) for j in back[columns*(first_columns_width + second_columns_width) + 1:].split( ' ' ) if j.strip()]

        # turn into dicts
        first_columns_dict = dict( [(first_columns_keys[j] , first_columns_values[j]) for j in xrange( len( first_columns_keys ) )] )
        second_columns_dict = dict( [(second_columns_keys[j] , second_columns_values[j]) for j in xrange( len( second_columns_keys ) )] )

        # all we need
        pssm_dict[position] = {
            'position' : position ,
            'query identity' : query_identity ,
            # are these proper?
            'log-likelihood' : first_columns_dict ,
            'approximate frequencies' : second_columns_dict ,
            'information content' : back[0] ,
            '?' : back[1]
            }

    return pssm_dict


