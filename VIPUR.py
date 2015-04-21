#!/usr/bin/env python
# :noTabs=true:

"""
What is VIPUR?
The VIPUR pipeline analyzes protein variants by considering conservation scores
and structural scores to identify variants that are likely to disrupt protein
function
Using Rosetta, these structural features allow you to interpret what is causing
the variant to be identified as deleterious

Conservation scores are derived from a PSSM of similar sequences found
using PSIBLAST against the NCBI nr database

Structural analysis is done using Rosetta to:
    consider variant structures by rapid structure optmization, allowing
    fast evaluation of approximate variant ddG values (Rosetta ddg_monomer)
        and
    refine variant structures and consider the distribution of energies and
    structural scores (rms, gdtmm) across several low energy conformations
    (physically near the input conformation, Rosetta relax)

Additional features are provided by an internal "aminochange" classification
(crude similarity of amino acid properties) and the variant position
surface area, evaluated using PROBE

These analyses are combined with a learned Logistic Regression model to
classify input variants as "neutral" or "deleterious" and provide a
structurally-informed hypothesis as to why variants are likely to
disrupt the protein



##################
(rewrite all this)
VIPUR INSTALLATION
- download the VIPUR module
- verify that you have the necessary software
 running VIPUR requires
    PSIBLAST
    PROBE
    Rosetta
    PyMOL (or PyRosetta)
-add paths to settings.py


DOWNLOAD PSIBLAST
PSIBLAST is a sequence search tool based on iterative alignment profile scans
PSIBLAST is part of the free NCBI BLAST+ distribution that can be found at:

ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

There are many tutorials on downloading, installing, and running PSIBLAST,
we suggest the NCBI usage book:

http://www.ncbi.nlm.nih.gov/books/NBK52640/

VIPUR has been benchmarked using BLAST 2.2.25+
with the NCBI nr database (downloaded Feb 2012)

Camacho, C. et al. BLAST+: architecture and applications.
BMC Bioinformatics 10(421) (2009).
	

DOWNLOAD PROBE
PROBE is a tool for analyzing contacts within protein structures (PDB format)
Here, we use PROBE as an accurate method for calculating the
ACCessible surface area of Protein amino acids (ACCP), the fraction of
potential surface area in contact with other residues
PROBE is freely available for download:

http://kinemage.biochem.duke.edu/software/probe.php

VIPUR has been benchmarked using probe.2.12.071128

Word, et. al. Visualizing and Quantifying Molecular Goodness-of-Fit:
Small-probe Contact Dots with Explicit Hydrogens.
J. Mol. Biol. 285, 1711-1733 (1999).


DOWNLOAD ROSETTA
Rosetta is a software suite capable of modeling and designing proteins
and other biomacromolecules
Here, we use the well established Rosetta protocols "relax" and "ddg_monomer"
for refining variant protein structures and evaluating protein energetics
Rosetta is free for academic use with the license and download tutorial at:

https://www.rosettacommons.org/software/license-and-download

Compiling Rosetta may introduce additional dependencies depending on your system
Here, we will use the default "release" compilation settings
(e.g. no "mode=MPI" etc.)
Please consult the Rosetta user documenation and forums if you encounter
complications during setup:

https://www.rosettacommons.org/docs/latest/Build-Documentation.html

Note: both relax and ddg_monomer perform stochastic searches and may output
slightly different energy values on evaluation with different random seeds
VIPUR has been benchmarked using the Rosetta 3.4 release version (54167)

Leaver-Fay, A. et al. ROSETTA3: An object-oriented software suite for
the simulation and design of macromolecules.
Methods in Enzymology 487, 548-574 (2011).


DOWNLOAD PYMOL
PyMOL is a tool for molecular visualization and analysis
PyMOL has several versions with advanced features, however we only require
simple PyMOL functionality
Please consult the PyMOL website to determine which license works best for you:

http://pymol.org/educational/

PyMOL can be downloaded at:

http://pymol.org/dsc/

Note: PyMOL will be used to create variant structures, PyRosetta can
alternatively be used for this task, however you only need one of these programs

The PyMOL Molecular Graphics System, Version 1.5.0.4 Schrodinger, LLC.


Note: you only need a copy of PyMOL or PyRosetta to run VIPUR
(for creating variant structures)
DOWNLOAD PYROSETTA
PyRosetta is an interactive Python interface to Rosetta
PyRosetta is free for academic use but on a separate license from Rosetta:

http://c4c.uwc4c.com/express_license_technologies/pyrosetta

The PyRosetta software and a download tutorial can be found at:

http://www.pyrosetta.org/dow

Chaudhury, S. et al. PyRosetta: a script-based interface for implementing
molecular modeling algorithms using Rosetta.
Bioinformatics 26(5), 689-691 (2010).



###########
USING VIPUR
VIPUR is written as a simple Python module
Once you have VIPUR and the required software setup, you can use VIPUR.py as
a Python script or a library (from VIPUR import run_VIPUR)
other scripts that make up VIPUR outline the specific feature generation steps,
including input file parsing and output file analysis

Currently, VIPUR supports a commandline interface to the options necessary to
run a single protein variant, specifying the input structure file (-p)
<pdb_filename>  and variants (-v)  <variant_filename>

ex.
python VIPUR.py -p example_input/2C35.pdb -v example_input/2C35.txt

VIPUR can also run on a directory of (PDB, .txt) file pairs by inputting a path
to the directory containing the files (-p), defaults to the current directory

ex.
python VIPUR.py -p example_input/

Please consult the help (-h) for more details on run options, including:
    -d      filename of the output predictions file
    -o      path for output to be written (several intermediate files)
    -c      chain of the input PDB that contains the native sequence (if not A)
    -s      filename of the intermediate (native) sequence file (FASTA format)
    -w      option to write out the numbering map (PDB to 1-indexed)
    -q      option to run in "sequence only" mode, no structural analysis

An example of VIPUR input files is provided in "example_input"
and their expected output is provided in "example_output_reference"
You can automatically run VIPUR on this demo with the --demo option

ex.
python VIPUR.py --demo


VARIANT INPUT FORMAT
VIPUR currently takes variant input files as plain text files with one variant
per line, expecting a reference amino acid for the native protein

ex.
E14R
R84P
A101W

the native reference amino acid is necessary to ensure the input numbering
matches the input protein structure (so you don't have to worry about constantly
renumbering your indices or PDB files)
currently, any input variants that do not have correct positions or native
amino acids will be skipped


INTERPRETING RESULTS
VIPUR outputs a ".predictions" file containing a summary of the predictions and
analysis of input variants. The default output file is tab delimited ('\t')
with one line per variant indicating:
    variant                         the variant (ex. P335R)
    predicted label                 "deleterious" or "neutral"*
    prediction confidence           P(deleterious|analysis), deleterious score*
    structure-only label            prediction of the structure-only classifier
    structure-only confidence       structure-only confidence score
    sequence-only label             prediction of the sequence-only classifier
    sequence-only confidence        sequence-only confidence score
    exposure                        "surface" or "interior" for the position**
    "essential" score               can indicate conserved positions***
    ddG prediction                  prediction of ddG in kcal/mol (approx)
    interpretation                  simple interpretation of the effect****
    explanation                     top contributions to the interpretation****

*VIPUR was trained on a dataset of natural variants, pseudomutations, and
protein variants from mutagenesis experiments. Variants are curated as
"deleterious" if they have literature or UniProt annotations indicating loss of
an essential protein function. Pseudomutations are from HumDiv and assumed to be
"neutral" (though not all neutral examples come from pseudomutations).
Please see the  <main text>  for a full description of the VIPUR training set.
Note: this binary label refers to PROTEIN function which may, or may not,
indicate disease association (a more complicated phenotype)
variants that are predicted "deleterious" are expected to lack an essential
protein function (in many cases, the protein variant is misfolded or
too unstable to maintain its native fold, as indicated by the structure-based
features), however they may display a severely reduced function or
maintain other functions
variants that are predicted "neutral" are expected to function effectively as
the native protein, however they may have slightly reduced function
(nearly-neutral, not predicted as not well-curated data exists)

Note: our "deleterious" predictions are NOT synonymous with the impact
on stability, deleterious variants may minimally change or even stabilize
a protein (ex. can make an enzyme too rigid) and neutral variants can
reasonably destabilize a protein as long as it still folds,
please consult the ddG prediction provided by the ddg_monomer protocol to
interpret changes in stability
Note: our "deleterious" predictions are NOT synonymous with "disease assocition"
while many disease associated/causal variants are loss-of-function changes,
the disruption of protein function is itself insufficient to indicate the cause
of a disease (this comes from knowledge of the protein's function or prior
variant association)
the goal of VIPUR is to provide a clear prediction of variants that disrupt
protein function to help INTERPRET variants that already have some correlated
label (e.g. disease or other phenotype)
not all deleterious predictions of VIPUR necessarily influence disease, however
confident deleterious predictions of variants known to be disease associated
suggests a strong effect worth investigating

Since VIPUR contains a binary classifier, the output confidence metric (the
learned conditional probability) indicates confidence of the binary prediction.
This confidence score is effectively a "deleterious" score, with higher values
indicating increased probability of being deleterious.
Lower values similarly indicate an increased probability of being neutral
(since the labels are binary).
You can filter your results to contain the most confident deleterious
predictions (ex. identify the 5 protein variants with the highest
deleterious confidence etc.).

**"exposure" indicates the local environment of the variant amino acid position
Here, we restrict "exposure" to "surface" and "interior", indicating if the
position is (approximately) on the protein surface or inside the protein.
The VIPUR training set does not suggest there is a significant difference
in performance for positions on the surface or interior of the protein,
although available data is imbalanced (many surface variants are annotated
"neutral", many interior variants are annotated "deleterious").

***Disagreement between the overall VIPUR classifier (sequence+structure) and
the structure-only classifier indicates a strong sequence signal without
apparent destabilization of the monomer structure. In some cases, this is due to
the inadequacy of the monomeric structure to capture the protein energetics,
suggesting this amino acid is important of interactions (can be with other
proteins, metals, ligands, nucleic acids etc.).
When analyzing VIPUR output, consider surface variants with
high (>=.8) deleterious scores and low structure-only scores (difference >=.2)
as "potential interaction sites".
Please see the  <main text>  for more information on VIPUR analysis and scores.

****VIPUR automatically identifies which features contribute to
deleterious predictions. The structure-based features are directly interpretable
since they represent specific physical interactions (hydrogen bonding,
solvent favorability, dilsufide bond strain, etc.) and expected distributions
based on other proteins (backbone conformational strain,
side chain interactions, etc.).
The top structure features are included as an "explanation" (3 by default,
can be set to include more).
We also include a reduced "interpretation" of the variant deleterious effect as:

structural conservation
The variant destabilizes the protein structure, likely due to improper size or
surface area. This includes unfit deviations in backbone conformation
(e.g. G and P positions) and packing configurations
(e.g. V, L, I, M differences). In some cases, it is seen that binding sites
have more restricted conformations, detectable by destabilizing variations even
in the absence of the binding partner (e.g. structurally conserved and the
variant cannot attain the necessary conformation).

disrupted chemical interactions
This variant eliminates an important and/or stabilizing chemical interaction of
the protein. This includes hydrogen bonding partners, conserved hydroxyl groups,
salt bridges, and disulfide bonds.

potential interaction site
As noted above*** high deleterious scores with low confidence structure-only
scores can indicate conservation where there is no structural evidence and
can indicate potential interaction sites on the protein surface.

other
The remaining features are difficult to directly interpret. They indicate
sequence conservation and structural conservation, but cannot clearly suggest
why.


PREPARING STRUCTURES FOR ROSETTA
Rosetta requires structural models that are cleanly readable. For
many applications, this includeds removal of waters, ligands, nucleic acids,
and metals. Rosetta has methods for handling all of these inputs, however
our benchmark uses protein structures stripped of all these additional
coordinates (note: while Rosetta can handle these cases, many proteins lack
models of docked conformations, ligands, etc. and we wanted to ensure the same
information was available for all samples). There are several methods available
for cleaning structures for Rosetta, including simply removing all HETATM and
nucleic acid ATOM coordinates (note: Rosetta can handle input structures with
missing densities or regions ex. removed a non-canonical amino acid).
The script available at:

https://github.com/Olvikon/miscellaneous_scripts/blob/master/process_pdb.py

outputs a directory containing the monomeric structure cleaned for use with
Rosetta (either ".clean.pdb" or ".protein.pdb" if nucleic acid lines are
removed).
More suggestions on how to clean structures for Rosetta is available at:

https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/dd/da1/preparing_structures.html

note: additional refinement steps of initial models are unnecessary since
Rosetta relax is run as part of feature generation, and the benchmark
performance is using structures that have not been pre-relaxed
"""

################################################################################
# IMPORT

# common modules
import os
import shutil
import subprocess
import optparse
import time    # for debugging only

# bigger modules

# custom modules
from settings import *
from helper_methods import create_directory , copy_file
from pre_processing import *

from psiblast_feature_generation import *
from probe_feature_generation import *
from rosetta_feature_generation import *

from classification import *

################################################################################
# MAIN

# use disclaimer
# still not sure how we are going to gate this
# currently making code available...
disclaimer = '='*80 + '\n\nRosetta is freely available to academic and government laboratories,\nA license must first be obtained through the University of Washington through the Express Licensing Program.\nThese agreements have standard terms and conditions that allow for rapid licensing by companies, organizations, or individuals.\nCommercial licesenses are also available and Rosetta Commercial users have a higher priority for support assistance from developers.\n\nYou can obtain a license for Rosetta at:\n\thttps://www.rosettacommons.org/software/license-and-download\n\nDO NOT use Rosetta code if you have not obtained a license, it is AGAINST THE LAW\nif you are a commercial user, Rosetta licenses are available\n\n'
print disclaimer
#if not LICENSE:
#    r = raw_input( 'have you obtained a Rosetta license? (yes/no): ' )
#    if not r == 'yes':
#        raise IOError( 'you MUST obtain a Rosetta license before using the Rosetta software' )

# credits to PROBE, BLAST+, Rosetta)
credits = 'VIPUR also relies on:\nPROBE\nWord, et. al. Visualizing and Quantifying Molecular Goodness-of-Fit:\nSmall-probe Contact Dots with Explicit Hydrogens.\nJ. Mol. Biol. 285, 1711-1733 (1999).\n\nBLAST+\nCamacho, C. et al. BLAST+: architecture and applications.\nBMC Bioinformatics 10(421) (2009).'
print credits


# single run method
def run_VIPUR( pdb_filename , variants_filename , prediction_filename = '' , out_path = '' ,
        target_chain = 'A' , sequence_filename = '' , write_numbering_map = True ,
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
        sequence , residue_map , sequence_filename = extract_protein_sequence_from_pdb( pdb_filename , out_filename = sequence_filename , target_chain = target_chain , write_numbering_map = write_numbering_map )
        # "residue_map" is a residue number map between PDB and 1-indexed positions
        # unfortunately many Rosetta programs will do this by default, such as
        # ddg_monomer, or otherwise output 1-indexed protein positions

    # load the variants
    variants = load_variants_file( variants_filename )
    
    # verify the input variants are proper
    variants = check_variants( variants , sequence , residue_map )
    if not variants:
        raise IOError( 'none of the variants in \"' + variants_filename + '\" are acceptable!!?!' )
    print '\ngenerating features for variants:\n\t' + '\n\t'.join( variants )

    # for storing the feature values
    variants = dict( [(i , {}) for i in variants] )
    
    # make the variant structures
    # this is currently done using pymol
#    create_variant_structures( pdb_filename , variants )
    if not sequence_only:
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
    print 'evaluating the \"aminochange\" values'
    evaluate_aminochange = lambda nat , var : 2 - int( [i for i in AMINOCHANGE_GROUPS if nat in i][0] == [i for i in AMINOCHANGE_GROUPS if var in i][0] ) - int( nat == var )
    for i in variants.keys():
        variants[i]['aminochange'] = evaluate_aminochange( i[0] , i[-1] )

    debug_time.append( ('preprocessing' , time.time()) )
    
    # run PSIBLAST
    print 'running PSIBLAST...'
    psiblast_filename = run_psiblast( sequence_filename )
    psiblast_filename = sequence_filename.rstrip( '.fa' ) + '.pssm'
        
    # extract PSIBLAST feature values
    print 'extracting PSSM features from PSIBLAST output'
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
        print 'running PROBE and determining the ACCP'
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

        print 'running Rosetta ddg_monomer...'
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
        print 'running Rosetta relax on the native structure...'
        native_relax_filename = run_rosetta_relax( pdb_filename )    # reference for other structures
        native_score_filename = run_rosetta_rescore( native_relax_filename , native_filename = pdb_filename )
        native_scorefile_dict = extract_scores_from_scorefile( native_score_filename )    # save time, only parse this once
        debug_time.append( ('relax' , time.time()) )
    
        for i in variants.keys():
            # this file SHOULD exist for every variant
            # if it does not, there was a problem earlier
            variant_structure = variants[i]['variant structure filename']
            # note: this current version assumes 1:1 for variant structures and
            #    variants...could change in the future...
    
            print 'running Rosetta relax on ' + variant_structure + '...'
            relax_filename = run_rosetta_relax( variant_structure )
        
            debug_time.append( ('structure ' + variant_structure , time.time()) )

            # post process, extract scores and generate quartile scores
            print 'extracting the relax features...'

            # additional step, rescore with Rosetta
            # only required to add in 2 scores that SHOULD be output by default with relax
            # whatever, this is the easiest solution, just make an additional call
            score_filename = run_rosetta_rescore( relax_filename , native_filename = pdb_filename )

            # extract features, use the quartile method to extract comparisons
            # between the native and variant score distributions
            quartile_scores = extract_quartile_score_terms_from_scorefiles( score_filename , native_scorefile_dict )
            # make sure there is not overlap
            if [None for j in quartile_scores.keys() if j in variants[i].keys()]:
                raise IOError( '!!?! somehow your scorefile has score terms that match the name(s) of features that already exist!!?! something has gone horribly wrong...' )
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


if __name__ == '__main__':
    # parser object for managing input options
    parser = optparse.OptionParser()
    # essential data
    parser.add_option( '-p' , dest = 'pdb_filename' ,
        default = '' ,
        help = 'the PDB file containing the native structure (and inherently the native sequence)' )
    parser.add_option( '-v' , dest = 'variants_filename' ,
        default = '' ,
        help = 'the file containing line-separated variants to predict' )

    parser.add_option( '-d' , dest = 'prediction_filename' ,
        default = '' ,
        help = 'desired filename for the output predictions' )
    parser.add_option( '-o' , dest = 'out_path' ,
        default = '' ,
        help = 'desired path (or directory) for output (derived from input PDB by default)' )

    parser.add_option( '-c' , dest = 'target_chain' ,
        default = 'A' ,
        help = 'the chain ID of the native protein within the PDB, assumes \"A\" by default (!!!)' )
    parser.add_option( '-s' , dest = 'sequence_filename' ,
        default = '' ,
        help = 'desired filename for the protein sequence file (FASTA format) of the input PDB file (written out as default)' )

    # options you wont use...but options none the less
    parser.add_option( '-w' , dest = 'write_numbering_map' ,
        default = True , action = 'store_false' ,
        help = 'boolean (default True), write the numbering map among the output files (provide -w for False)' )
    parser.add_option( '-q' , dest = 'sequence_only' ,
        default = False , action = 'store_true' ,
        help = 'boolean (default False), run in \"sequence only\" mode, no structural features' )
    parser.add_option( '--demo' , dest = 'demo' ,
        default = False , action = 'store_true' ,
        help = 'optional flag that overrides normal behaviour and runs the demo' )

    # optionally allow specification of input paths
    # well...do this another day...

    (options,args) = parser.parse_args()
    
    pdb_filename = options.pdb_filename
    variants_filename = options.variants_filename
    prediction_filename = options.prediction_filename
    out_path = options.out_path
    target_chain = options.target_chain
    sequence_filename = options.sequence_filename
    write_numbering_map = bool( options.write_numbering_map )
    sequence_only = bool( options.sequence_only )
    demo = bool( options.demo )

    # for the example input
    if demo:
        pdb_filename = PATH_TO_VIPUR + '/example_input/2C35.pdb'
        variants_filename = PATH_TO_VIPUR + '/example_input/2C35.txt'

        prediction_filename = ''
        out_path = PATH_TO_VIPUR + '/example_output'
        target_chain = 'A'

    # alternatively, run on an entire directory
    if not pdb_filename and not variants_filename:
        print 'no input provided, assuming you want to run on every (.pdb,.txt) file pair found in the current directory'
        pdb_filename = os.getcwd()
    if os.path.isdir( pdb_filename ) and not variants_filename:
        variants_filename = '.txt'
    if os.path.isdir( pdb_filename ) and variants_filename[0] == '.':
        # instead, run on the directory
        print 'running VIPUR on all (.pdb,' + variants_filename + ') file pairs found in ' + pdb_filename
        # find .pdb files
        pdb_filenames = [i for i in os.listdir( pdb_filename ) if i[-4:] == '.pdb']
        # look for pairs
        os.chdir( pdb_filename )
#        print pdb_filename , pdb_filenames , os.getcwd() , variants_filename , pdb_filenames[0][:-4] + variants_filename
        pdb_filenames = [(i , i[:-4] + variants_filename) for i in pdb_filenames if os.path.isfile( i[:-4] + variants_filename )]
        print str( len( pdb_filenames ) ) + ' pairs found'
        if not pdb_filenames:
            raise IOError( '!!! no (.pdb,' + variants_filename + ') file pairs found in ' + pdb_filename + '!!?!' )
        
        # run them all
        for i in pdb_filenames:
            # assume all other values are defauls
            run_VIPUR( i[0] , i[1] )
        
        exit()    # so it does not error-out after the loop
    
        
    run_VIPUR( pdb_filename , variants_filename ,
        prediction_filename = prediction_filename , out_path = out_path ,
        target_chain = target_chain , write_numbering_map = write_numbering_map ,
        sequence_only = sequence_only )


