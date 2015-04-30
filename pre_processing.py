#!/usr/bin/env python
# :noTabs=true:

"""
VIPUR methods for loading and quality-checking variants

Currently other pre-processing methods that are feature-specific, such as variant structure generation
are contained within the script associated with those features
"""

################################################################################
# IMPORT

# common modules
import os

# bigger modules

# custom modules
from settings import PROTEIN_LETTERS

################################################################################
# METHODS

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
    for variation in variants:
        native = variation[0]
        position = variation[1:-1]
        mutant = variation[-1]
            
        # check for nonsense variation
        if not mutant in protein_letters or not native in protein_letters:
            print '\"' + variation + '\" is improperly formatted, unaccepted amino acid single letter code'
            continue
        
        # check the position
        if not position in residue_map.keys():
            print 'residue \"' + position + '\" is not in the PDB, skipping \"' + variation + '\"'
            continue

        # check for redundancy
        if native == mutant:
            print 'native and variant are identical (?), skipping \"' + variation + '\"'
            continue
        
        # check for native errors
        pdb_aa = sequence[residue_map[position]]
        if not native == pdb_aa:
            print 'wrong native amino acid! found \"' + pdb_aa + '\" in the pdb, not \"' + native + '\"'
            continue
        
        # it passed!
        permissible_variants.append( variation )

    return permissible_variants

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


