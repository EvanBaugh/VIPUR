#!/usr/bin/env python
# :noTabs=true:

"""
ehb: 12/3

PyMOL script that creates protein structures from the list of variations

INPUT:


Note: finish documentation for release, currently the preferred method of variant structure generation
...although this introduces dependence on PyMOL (currently MUCH more practical than relying on Rosetta or PyRosetta)
"""

################################################################################
# IMPORT

# common modules

# bigger modules
import pymol
import optparse    # for commandline options

# custom modules

################################################################################
# METHODS

AMINO_ACID_CODES = {
    'A' : 'ALA' ,
    'C' : 'CYS' ,
    'D' : 'ASP' ,
    'E' : 'GLU' ,
    'F' : 'PHE' ,
    'G' : 'GLY' ,
    'H' : 'HIS' ,
    'I' : 'ILE' ,
    'K' : 'LYS' ,
    'L' : 'LEU' ,
    'M' : 'MET' ,
    'N' : 'ASN' ,
    'P' : 'PRO' ,
    'Q' : 'GLN' ,
    'R' : 'ARG' ,
    'S' : 'SER' ,
    'T' : 'THR' ,
    'V' : 'VAL' ,
    'W' : 'TRP' ,
    'Y' : 'TYR' ,
    }

# use pymol to create a mutation
def mutate_residue( selection , mutation , out_filename = '' , mutant_selection_name = 'mutant' , reference = 'reference' ):
    if not selection in pymol.cmd.get_names() and os.path.isfile( selection ):
        # assume its an input filename
        pymol.cmd.load( selection , reference )
        selection = reference
    # handle the input mutation
    if isinstance( mutation , str ):
        if '.' in mutation:
            mutation = mutation.split( '.' )
        else:
            # just make it a list
            mutation = [mutation]
    elif not isinstance( mutation , list ):
        raise IOError( 'not currently supported' )
    
    # make a copy
    pymol.cmd.copy( mutant_selection_name , selection )
    
    pymol.cmd.wizard( 'mutagenesis' )
    for i in mutation:
        # split it
        native = i[0]
        mutant = AMINO_ACID_CODES[i[-1]]
        if not native in AMINO_ACID_CODES.keys():
            # assume the native was not specified
            position = i[:-1]
        else:
            # verify the positions
            native = AMINO_ACID_CODES[native]
            position = i[1:-1]
            # use iterate to get the resn and chain
            # need the chain for selecting the proper position
            pymol.cmd.iterate( mutant_selection_name + ' and resi ' + position + ' and name CA' , '(pymol.stored.native_resn , pymol.stored.native_chain) = (resn , chain)' )
            if not pymol.stored.native_resn == native:
                raise IOError( 'input native as ' + i[:-1] + ' but it is actually ' + pymol.stored.native_resn + position + ' (chain ' + pymol.stored.native_chain +')' )
#            elif pymol.stored.native_resn == mutant:
#                raise IOError( 'input mutation as ' + i + ' and position ' + position + ' is already ' + pymol.stored.native_resn + ' (chain ' + pymol.stored.native_chain +')' )

        # make the changes
        pymol.cmd.get_wizard().set_mode( mutant )
        pymol.cmd.edit( '/' + mutant_selection_name +'//'+ pymol.stored.native_chain +'/'+ position + '/CA' )
#        pymol.cmd.edit( mutant_selection_name + ' and chain ' + pymol.stored.native_chain + ' and resi ' + position + ' and name CA' )
        pymol.cmd.get_wizard().do_pick( 0 )
        pymol.cmd.get_wizard().apply()
    
    pymol.cmd.set_wizard()
    print 'created ' + mutant_selection_name + ' from ' + selection + ' by changing: ' + ', '.join( mutation )

    # write out?
    if out_filename:
        pymol.cmd.save( out_filename , mutant_selection_name )
pymol.cmd.extend( 'mutate_residue' , mutate_residue )

# alternate interface, load a file and apply multiple changes etc., optionally write out
def mutate_pdb( filename , mutation , root_filename = True , chain = 'A' , mutant_selection_name = 'mutant' , reference = 'reference' ):
    # load it first
    pymol.cmd.load( filename , reference )
    if chain:
        # ...just delete everything else
        pymol.cmd.remove( reference + ' and not chain ' + chain )
    # interpret mutation input
    if isinstance( mutation , str ):
        if ',' in mutation:
            mutation = [i.split( '.' ) for i in mutation.split( ',' )]
        elif '.' in mutation:
            mutation = mutation.split( '.' )
        else:
            # just make it a list
            mutation = [mutation]
    elif not isinstance( mutation , list ):
        raise IOError( 'not currently supported' )
    # output
    if root_filename and not isinstance( root_filename , str ):
        root_filename = filename.rstrip( '.pdb' )    # must be ".pdb"

    for i in mutation:
        if isinstance( i , list ):
            out_filename = root_filename + '.chain_' + chain +'_'+ '_'.join( i ) +'.pdb'
        else:
            out_filename = root_filename + '.chain_' + chain +'_'+ i +'.pdb'
        # write out by default
#        print out_filename
        mutate_residue( reference , i , out_filename , mutant_selection_name = mutant_selection_name )
pymol.cmd.extend( 'mutate_pdb' , mutate_pdb )

################################################################################
# MAIN

if __name__ == '__main__':
    # parser object for managing input options
    parser = optparse.OptionParser()
    # essential data
    parser.add_option( '-p' , dest = 'pdb_filename' ,
        default = '' ,
        help = 'the PDB structure of the native protein' )
    parser.add_option( '-m' , dest = 'mutations' ,
        default = '' ,
        help = 'the mutations to make, separated by \",\" between variants and by \".\" for multiple changes to the same variant' )

    parser.add_option( '-r' , dest = 'root_filename' ,
        default = True ,
        help = 'prefix for output variant PDB structures, defaults to \"True\" which uses the input PDB filename' )
    parser.add_option( '-c' , dest = 'chain' ,
        default = 'A' ,
        help = 'which chain to mutate' )

    parser.add_option( '-s' , dest = 'mutant_selection_name' ,
        default = 'mutant' ,
        help = 'internal PyMOL name for the variant structure' )
    parser.add_option( '-f' , dest = 'reference' ,
        default = 'reference' ,
        help = 'internal PyMOL name for the native structure' )

    (options,args) = parser.parse_args()

    # check inputs
    pdb_filename = options.pdb_filename
    # don't bother formatting input here, handled by functions
    mutations = options.mutations
    
    root_filename = options.root_filename
    chain = options.chain
    mutant_selection_name = options.mutant_selection_name
    reference = options.reference
    
    mutate_pdb( pdb_filename , mutations , root_filename = root_filename , chain = chain , mutant_selection_name = mutant_selection_name , reference = reference )


