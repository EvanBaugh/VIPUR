#!/usr/bin/env python
# :noTabs=true:

"""
Contains the paths and settings for VIPUR
Modify the paths or method settings as needed

When setting up VIPUR, make sure you provide all the paths needed
"""

################################################################################
# PATHS

PATH_TO_VIPUR = '/home/evan/VIPUR_pipeline/VIPUR'    # remember, this should be blank each push
#PATH_TO_VIPUR = '/scratch/ehb250/VIPUR_pipeline/VIPUR'    # remember, this should be blank each push
PATH_TO_VIPUR_EXECUTABLES = '/home/evan/VIPUR_pipeline/VIPUR_feature_executables' #/VIPUR_feature_executables'

PATH_TO_PSIBLAST = PATH_TO_VIPUR_EXECUTABLES + '/psiblast_blast+2.2.25'
#PATH_TO_PSIBLAST = '/home/ehb250/MEC_files/ncbi-blast-2.2.25+/bin/psiblast'#PATH_TO_VIPUR_EXECUTABLES + '/psiblast_blast+2.2.25'
PATH_TO_BLAST_DATABASE = '/home/evan/bio/databases/blast/nr/nr'
#PATH_TO_BLAST_DATABASE = '/scratch/ehb250/nr_db/nr'#/home/evan/bio/databases/blast/nr/nr'
PATH_TO_PROBE = PATH_TO_VIPUR_EXECUTABLES + '/probe'
#PATH_TO_PROBE = '/home/ehb250/MEC_files/probe/probe'#PATH_TO_VIPUR_EXECUTABLES + '/probe'

LICENSE = ''

PATH_TO_ROSETTA_DATABASE = PATH_TO_VIPUR_EXECUTABLES + '/rosetta_database'
#PATH_TO_ROSETTA_DATABASE = '/scratch/ehb250/rosetta-3.4/rosetta_database'#PATH_TO_VIPUR_EXECUTABLES + '/rosetta_database'
PATH_TO_ROSETTA_DDG_MONOMER = PATH_TO_VIPUR_EXECUTABLES + '/ddg_monomer_r54167.64bit.linuxgccrelease'
#PATH_TO_ROSETTA_DDG_MONOMER = '/scratch/ehb250/rosetta-3.4/rosetta_source/bin/ddg_monomer.linuxgccrelease'#PATH_TO_VIPUR_EXECUTABLES + '/ddg_monomer_r54167.64bit.linuxgccrelease'
PATH_TO_ROSETTA_RELAX = PATH_TO_VIPUR_EXECUTABLES + '/relax_r54167.64bit.linuxgccrelease'
#PATH_TO_ROSETTA_RELAX = '/scratch/ehb250/rosetta-3.4/rosetta_source/bin/relax.mpi.linuxgccrelease'#PATH_TO_VIPUR_EXECUTABLES + '/relax_r54167.64bit.linuxgccrelease'
PATH_TO_ROSETTA_SCORE = PATH_TO_VIPUR_EXECUTABLES + '/score_r54167.64bit.linuxgccrelease'
#PATH_TO_ROSETTA_SCORE = '/scratch/ehb250/rosetta-3.4/rosetta_source/bin/score.linuxgccrelease'#PATH_TO_VIPUR_EXECUTABLES + '/score_r54167.64bit.linuxgccrelease'

# alternate method for making variant structures, needs both paths
PATH_TO_PYMOL = 'pymol'
#PATH_TO_PYMOL = '/share/apps/pymol/1.5.0.1/bin/pymol'

# check for PyRosetta, used for making mutant structures
USE_PYROSETTA = False
if USE_PYROSETTA:
    try:
        from rosetta import init
        init()
        USE_PYROSETTA = True
    except:
        None
if not USE_PYROSETTA:
    # anything in response?
    if not PATH_TO_PYMOL:
        raise IOError( 'please edit settings.py to include your PATH_TO_PYMOL' )

# classifier files
# note: parameter values could be hardcoded to reduce file IO
FINAL_CLASSIFIER_WEIGHTS = PATH_TO_VIPUR + '/VIPUR_trained_model_parameters/final_classifier.weights'
STRUCTURE_ONLY_MODEL_WEIGHTS = PATH_TO_VIPUR + '/VIPUR_trained_model_parameters/structure_only.weights'
SEQUENCE_ONLY_MODEL_WEIGHTS = PATH_TO_VIPUR + '/VIPUR_trained_model_parameters/sequence_only.weights'

TRAINING_SET_FEATURE_MEANS = PATH_TO_VIPUR + '/VIPUR_trained_model_parameters/training_set_feature.means'
TRAINING_SET_FEATURE_STDS = PATH_TO_VIPUR + '/VIPUR_trained_model_parameters/training_set_feature.stds'

################################################################################
# AMINO ACID DETAILS

# for verifying input variants match the input structure
# the PDB format lists the amino acids differently, need to convert to sequence
PROTEIN_LETTERS = 'ACDEFGHIKLMNPQRSTVWY'
AMINO_ACID_CODES = {
    'ALA' : 'A' ,
    'CYS' : 'C' ,
    'ASP' : 'D' ,
    'GLU' : 'E' ,
    'PHE' : 'F' ,
    'GLY' : 'G' ,
    'HIS' : 'H' ,
    'ILE' : 'I' ,
    'LYS' : 'K' ,
    'LEU' : 'L' ,
    'MET' : 'M' ,
    'ASN' : 'N' ,
    'PRO' : 'P' ,
    'GLN' : 'Q' ,
    'ARG' : 'R' ,
    'SER' : 'S' ,
    'THR' : 'T' ,
    'VAL' : 'V' ,
    'TRP' : 'W' ,
    'TYR' : 'Y' ,
    # pseudo-standard 3 letter codes for the standard aa
    'CYD' : 'C' ,
    'CYZ' : 'C' ,
    'HID' : 'H' ,
    'HIE' : 'H' ,
    'HIP' : 'H'
    }


# feature generation settings

# amino acid groups for the "aminochange" feature
AMINOCHANGE_GROUPS = [
    'GP' ,
    'AVIL' ,
    'MFWY' ,
    'HKR' ,
    'NQ' ,
    'DE' ,
    'CST'
    ]

################################################################################
# FEATURE GENERATION SETTINGS

PSIBLAST_OPTIONS = {
    'db' : PATH_TO_BLAST_DATABASE ,

#    'num_threads' : 40 ,    # messes up my PBS submission
    'num_iterations' : 2 ,
    'pseudocount' : 2 ,

    'evalue' : 1 ,
    'inclusion_ethresh' : .001 ,
    'comp_based_stats' : 1 ,

    'outfmt' : 7 ,
    'num_descriptions' : 3000 ,
    'num_alignments' : 300 ,

    # still keeping output for potential features later
    'out' : lambda x : x + '.pb' ,
    'out_ascii_pssm' : lambda x : x + '.pssm' ,
    'out_pssm' : lambda x : x + '.cp' ,
    'export_search_strategy' : lambda x : x + '.ss'
    }

PROBE_OPTIONS = {
    'rad1.4' : '' ,    # sets the radius for "sphere rolling"
    'C' : ''    # ?
    }

ROSETTA_DDG_MONOMER_OPTIONS = {
    'database' : PATH_TO_ROSETTA_DATABASE ,

    'ddg::dump_pdbs' : 'false' ,
    'ddg::output_silent' : 'false' ,
    'ddg::suppress_checkpointing' : 'true' ,
    'ddg::iterations' : 5 ,
    'ddg::weight_file' : 'soft_rep_design' ,
    'ddg::local_opt_only' : 'true' ,
    'ddg::min_cst' : 'false' ,
    'ddg::mean' : 'true' ,
    'ddg::min' : 'false' ,
    'ddg::sc_min_only' : 'false' ,
    'ddg::ramp_repulsive' : 'false' ,
    'ddg::opt_radius' : 8.0 ,
    'linmem_ig' : 10 ,
    'run:ignore_zero_occupancy' : 'false' ,
#    'mute' : 'all'

    'run:constant_seed' : '' ,    # for reproducibility
    'run:jran' : 17 ,
    }

ROSETTA_RELAX_OPTIONS = {
    'database' : PATH_TO_ROSETTA_DATABASE ,

    'nstruct' : 50 ,    # 50!
#    'nstruct' : 50 ,
    'relax:fast' : '' ,
    'evaluation:gdtmm' : 'true' ,
#    'in:file:native' : lambda x : x + '.pdb' ,
#    'bGDT' : 'true' ,
    'run:ignore_zero_occupancy' : 'false' ,
#    'jd2:mpi_file_buf_job_distributor' : 'false' ,
#    'run:multiple_processes_writing_to_one_directory' : '' ,

    'run:constant_seed' : '' ,    # for reproducibility
    'run:jran' : 17 ,
    
    'out:file:silent' : lambda x : x + '.silent' ,
    'out:file:scorefile' : lambda x : x + '.sc' ,
    
#    'parallel' : 40 ,    # a Rosetta option?
    }
ROSETTA_RELAX_PARALLEL = 40#False    # OPTIONS should be reserved for explicit options to Rosetta

ROSETTA_SCORE_OPTIONS = {
    'database' : PATH_TO_ROSETTA_DATABASE ,
    'in:file:fullatom' : '' ,
    'rescore:verbose' : '' ,
    'out:file:scorefile' : lambda x : x + '.sc'
#    'score:weights' : 'test.wts'
    }

################################################################################
# PBS QUEUE SYSTEM INTERACTION

# specific to the author's queuing system at NYU

PBS_USER = 'ehb250'

# any commands that are needed to setup the environment for preprocessing
PBS_ENVIRONMENT_SETUP = 'module load pymol'

PBS_QUEUE_QUOTA = 20    # how many jobs can be in the queue simultaneously (excludes "R"unning jobs, that quota is set elsewhere now...)
PBS_QUEUE_MONITOR_DELAY = 60    # seconds, how long to wait between checking the queue

# simply add the flanking text necessary
PBS_BASH_SCRIPT = lambda x : '#!/bin/bash\n\n' + x.replace( ';' , '\n\n' ) +'\n\n'


PBS_SERIAL_JOB_OPTIONS = {
    'q' : 's48' ,
    'l' : 'nodes=1:ppn=1,mem=6gb,walltime=12:00:00' ,
    'o' : lambda x : x.replace( '.sh' , '.log.out' ) ,
    'e' : lambda x : x.replace( '.sh' , '.log.err' ) ,
#    'V' : '' ,    # necessary? seems to be case specific for me...
    }

PBS_PARALLEL_JOB_OPTIONS = {
    'q' : 'p24' ,
    'l' : 'nodes=3:ppn=12,mem=46gb,walltime=4:00:00' ,
    'o' : lambda x : x.replace( '.sh' , '.log.out' ) ,
    'e' : lambda x : x.replace( '.sh' , '.log.err' ) ,
    }

################################################################################
# SLURM QUEUE SYSTEM INTERACTION

SLURM_USER = 'ebaugh'

SLURM_QUEUE_QUOTA = 10
SLURM_QUEUE_MONITOR_DELAY = 60    # less then this?

SLURM_BASH_SCRIPT = lambda x : '#!/bin/bash\n\n' + x.replace( ';' , '\n\n' ) +'\n\n'

SLURM_SERIAL_JOB_OPTIONS = {
#    'n' : ,    # how many to request? is it actually threaded? just do 1 for now, later we can combine multiple for arrays or array like implementation
    # lol, these two are the same
    'o' : lambda x : x.replace( '.sh' , '.log.out' ) ,
    'e' : lambda x : x.replace( '.sh' , '.log.err' ) ,
    }

# no parallel options for now, Rosetta does not optimize this anyway...so SLURM's interface should make Rosetta MPI obsolete other than resource sharing...

################################################################################
# POST PROCESSING

# for feature extraction
ROSETTA_TERMS_TO_COMPARE = [
    'score' ,
    'fa_atr' ,
    'fa_rep' ,
    'fa_sol' ,
    'fa_intra_rep' ,
    'pro_close' ,
    'fa_pair' ,
    'hbond_sr_bb' ,
    'hbond_lr_bb' ,
    'hbond_bb_sc' ,
    'hbond_sc' ,
    'dslf_ss_dst' ,
    'dslf_cs_ang' ,
    'dslf_ss_dih' ,
    'dslf_ca_dih' ,
    'rama' ,
    'omega' ,
    'fa_dun' ,
    'p_aa_pp' ,
    'ref' ,
    'rms' ,
    'gdtmm' ,
    'maxsub' ,
#    'time' ,
#    'description'
    'maxsub2.0' ,
    'gdtmm7_4' ,
    'gdtmm4_3' ,
    'gdtmm3_3' ,
    'gdtmm2_2' ,
    'irms' ,
    'gdtmm1_1' ,
    'allatom_rms'
    ]

################################################################################
# CLASSIFICATION SETTINGS

# order of features read in, must match the order in the weights file
LABEL_DESCRIPTION = [
    'aminochange' ,
    'pssm_variant' ,
    'pssm_native' ,
    'pssm_difference' ,
    'pssm_information_content' ,
    'probe_accp' ,
    'quartile_scoreQ1' ,
    'quartile_scoreQ2' ,
    'quartile_scoreQ3' ,
    'quartile_fa_atrQ1' ,
    'quartile_fa_atrQ2' ,
    'quartile_fa_atrQ3' ,
    'quartile_fa_repQ1' ,
    'quartile_fa_repQ2' ,
    'quartile_fa_repQ3' ,
    'quartile_fa_solQ1' ,
    'quartile_fa_solQ2' ,
    'quartile_fa_solQ3' ,
    'quartile_fa_intra_repQ1' ,
    'quartile_fa_intra_repQ2' ,
    'quartile_fa_intra_repQ3' ,
    'quartile_pro_closeQ1' ,
    'quartile_pro_closeQ2' ,
    'quartile_pro_closeQ3' ,
    'quartile_fa_pairQ1' ,
    'quartile_fa_pairQ2' ,
    'quartile_fa_pairQ3' ,
    'quartile_hbond_sr_bbQ1' ,
    'quartile_hbond_sr_bbQ2' ,
    'quartile_hbond_sr_bbQ3' ,
    'quartile_hbond_lr_bbQ1' ,
    'quartile_hbond_lr_bbQ2' ,
    'quartile_hbond_lr_bbQ3' ,
    'quartile_hbond_bb_scQ1' ,
    'quartile_hbond_bb_scQ2' ,
    'quartile_hbond_bb_scQ3' ,
    'quartile_hbond_scQ1' ,
    'quartile_hbond_scQ2' ,
    'quartile_hbond_scQ3' ,
    'quartile_dslf_ss_dstQ1' ,
    'quartile_dslf_ss_dstQ2' ,
    'quartile_dslf_ss_dstQ3' ,
    'quartile_dslf_cs_angQ1' ,
    'quartile_dslf_cs_angQ2' ,
    'quartile_dslf_cs_angQ3' ,
    'quartile_dslf_ss_dihQ1' ,
    'quartile_dslf_ss_dihQ2' ,
    'quartile_dslf_ss_dihQ3' ,
    'quartile_dslf_ca_dihQ1' ,
    'quartile_dslf_ca_dihQ2' ,
    'quartile_dslf_ca_dihQ3' ,
    'quartile_ramaQ1' ,
    'quartile_ramaQ2' ,
    'quartile_ramaQ3' ,
    'quartile_omegaQ1' ,
    'quartile_omegaQ2' ,
    'quartile_omegaQ3' ,
    'quartile_fa_dunQ1' ,
    'quartile_fa_dunQ2' ,
    'quartile_fa_dunQ3' ,
    'quartile_p_aa_ppQ1' ,
    'quartile_p_aa_ppQ2' ,
    'quartile_p_aa_ppQ3' ,
    'quartile_refQ1' ,
    'quartile_refQ2' ,
    'quartile_refQ3' ,
    'quartile_allatom_rmsQ1' ,
    'quartile_allatom_rmsQ2' ,
    'quartile_allatom_rmsQ3' ,
    'quartile_gdtmmQ1' ,
    'quartile_gdtmmQ2' ,
    'quartile_gdtmmQ3' ,
    'quartile_gdtmm1_1Q1' ,
    'quartile_gdtmm1_1Q2' ,
    'quartile_gdtmm1_1Q3' ,
    'quartile_gdtmm2_2Q1' ,
    'quartile_gdtmm2_2Q2' ,
    'quartile_gdtmm2_2Q3' ,
    'quartile_gdtmm3_3Q1' ,
    'quartile_gdtmm3_3Q2' ,
    'quartile_gdtmm3_3Q3' ,
    'quartile_gdtmm4_3Q1' ,
    'quartile_gdtmm4_3Q2' ,
    'quartile_gdtmm4_3Q3' ,
    'quartile_gdtmm7_4Q1' ,
    'quartile_gdtmm7_4Q2' ,
    'quartile_gdtmm7_4Q3' ,
    'maxsub' ,    # no longer necessary/used, leave for compatability
    'maxsub2.0' ,
    'ddg_total' ,
    'ddg_fa_atr' ,
    'ddg_fa_rep' ,
    'ddg_fa_sol' ,
    'ddg_pro_close' ,
    'ddg_fa_pair' ,
    'ddg_hbond_sr_bb' ,
    'ddg_hbond_lr_bb' ,
    'ddg_hbond_bb_sc' ,
    'ddg_hbond_sc' ,
    'ddg_dslf_ss_dst' ,
    'ddg_dslf_cs_ang' ,
    'ddg_dslf_ss_dih' ,
    'ddg_dslf_ca_dih' ,
    'ddg_fa_dun' ,
    'ddg_p_aa_pp' ,
    'ddg_ref'
    ]

# how many of the top ranked features go into the output
TOP_FEATURES_TO_INCLUDE = 3

# used to qualitatively classify "interior" vs "surface" from the "probe_accp" term
PROBE_ACCP_INTERIOR_CUTOFF = 12.5
# from secondary structure/exposure rules, use the lower bound, shrinks "buried" to include only the lowest values

# skip these terms during summary, not very useful
SKIP_DURING_EXPLANATION = [
    'pssm_information_content' ,
    'pssm_variant' ,
    'pssm_difference' ,
    'pssm_native' ,
    'ddg_ref' ,
    'quartile_refQ1' ,
    'quartile_refQ2' ,
    'quartile_refQ3'
    ]

# this map enables
# for qualitative description
DELETERIOUS_PREDICTION_DESCRIPTION = {
    'aminochange' : 'unfavorable change in amino acid chemical properties' ,
    'pssm_variant' : 'conserved?' ,
    'pssm_native' : 'conserved?' ,
    'pssm_difference' : 'conserved?' ,
    'pssm_information_content' : 'conserved?' ,
    'probe_accp' : 'improper side chain surface area or volume' ,
    'quartile_scoreQ1' : 'highly destabilizing' ,
    'quartile_scoreQ2' : 'highly destabilizing' ,
    'quartile_scoreQ3' : 'highly destabilizing' ,
    'quartile_fa_atrQ1' : 'side chains cannot \"pack\" against each other (clashes)' ,    #'poor packing (clashes)' ,
    'quartile_fa_atrQ2' : 'side chains cannot \"pack\" against each other (clashes)' ,
    'quartile_fa_atrQ3' : 'side chains cannot \"pack\" against each other (clashes)' ,
    'quartile_fa_repQ1' : 'side chains cannot \"pack\" against each other (clashes)' ,
    'quartile_fa_repQ2' : 'side chains cannot \"pack\" against each other (clashes)' ,
    'quartile_fa_repQ3' : 'side chains cannot \"pack\" against each other (clashes)' ,
    'quartile_fa_solQ1' : 'unfavorable solvation (exposed hydrophobic or buried hydrophilic)' ,
    'quartile_fa_solQ2' : 'unfavorable solvation (exposed hydrophobic or buried hydrophilic)' ,
    'quartile_fa_solQ3' : 'unfavorable solvation (exposed hydrophobic or buried hydrophilic)' ,
    'quartile_fa_intra_repQ1' : 'unfavorable amino acid conformation' ,
    'quartile_fa_intra_repQ2' : 'unfavorable amino acid conformation' ,
    'quartile_fa_intra_repQ3' : 'unfavorable amino acid conformation' ,
    'quartile_pro_closeQ1' : 'change to/from proline not tolerated, likely conserved' ,
    'quartile_pro_closeQ2' : 'change to/from proline not tolerated, likely conserved' ,    #'disrupts a delicate proline, likely conserved' ,
    'quartile_pro_closeQ3' : 'change to/from proline not tolerated, likely conserved' ,
    'quartile_fa_pairQ1' : 'disrupted side chain contact (ex. salt bridge)' ,
    'quartile_fa_pairQ2' : 'disrupted side chain contact (ex. salt bridge)' ,
    'quartile_fa_pairQ3' : 'disrupted side chain contact (ex. salt bridge)' ,
    'quartile_hbond_sr_bbQ1' : 'disrupted backbone hydrogen bonding' ,
    'quartile_hbond_sr_bbQ2' : 'disrupted backbone hydrogen bonding' ,
    'quartile_hbond_sr_bbQ3' : 'disrupted backbone hydrogen bonding' ,
    'quartile_hbond_lr_bbQ1' : 'disrupted backbone hydrogen bonding' ,
    'quartile_hbond_lr_bbQ2' : 'disrupted backbone hydrogen bonding' ,
    'quartile_hbond_lr_bbQ3' : 'disrupted backbone hydrogen bonding' ,
    'quartile_hbond_bb_scQ1' : 'disrupted backbone - side chain hydrogen bonding' ,
    'quartile_hbond_bb_scQ2' : 'disrupted backbone - side chain hydrogen bonding' ,
    'quartile_hbond_bb_scQ3' : 'disrupted backbone - side chain hydrogen bonding' ,
    'quartile_hbond_scQ1' : 'disrupted hydrogen bond between side chains (likely conserved)' ,
    'quartile_hbond_scQ2' : 'disrupted hydrogen bond between side chains (likely conserved)' ,
    'quartile_hbond_scQ3' : 'disrupted hydrogen bond between side chains (likely conserved)' ,
    'quartile_dslf_ss_dstQ1' : 'destabilized disulfide' ,
    'quartile_dslf_ss_dstQ2' : 'destabilized disulfide' ,
    'quartile_dslf_ss_dstQ3' : 'destabilized disulfide' ,
    'quartile_dslf_cs_angQ1' : 'destabilized disulfide' ,
    'quartile_dslf_cs_angQ2' : 'destabilized disulfide' ,
    'quartile_dslf_cs_angQ3' : 'destabilized disulfide' ,
    'quartile_dslf_ss_dihQ1' : 'destabilized disulfide' ,
    'quartile_dslf_ss_dihQ2' : 'destabilized disulfide' ,
    'quartile_dslf_ss_dihQ3' : 'destabilized disulfide' ,
    'quartile_dslf_ca_dihQ1' : 'destabilized disulfide' ,
    'quartile_dslf_ca_dihQ2' : 'destabilized disulfide' ,
    'quartile_dslf_ca_dihQ3' : 'destabilized disulfide' ,
    'quartile_ramaQ1' : 'unfavorable backbone conformation' ,
    'quartile_ramaQ2' : 'unfavorable backbone conformation' ,
    'quartile_ramaQ3' : 'unfavorable backbone conformation' ,
    'quartile_omegaQ1' : 'unfavorable backbone conformation' ,
    'quartile_omegaQ2' : 'unfavorable backbone conformation' ,
    'quartile_omegaQ3' : 'unfavorable backbone conformation' ,
    'quartile_fa_dunQ1' : 'unfavorable side chain conformation' ,
    'quartile_fa_dunQ2' : 'unfavorable side chain conformation' ,
    'quartile_fa_dunQ3' : 'unfavorable side chain conformation' ,
    'quartile_p_aa_ppQ1' : 'unfavorable amino acid?' ,
    'quartile_p_aa_ppQ2' : 'unfavorable amino acid?' ,
    'quartile_p_aa_ppQ3' : 'unfavorable amino acid?' ,
    'quartile_refQ1' : 'reference?' ,
    'quartile_refQ2' : 'reference?' ,
    'quartile_refQ3' : 'reference?' ,
    'quartile_allatom_rmsQ1' : 'improper fold, cannot accomodate the variant amino acid' ,
    'quartile_allatom_rmsQ2' : 'improper fold, cannot accomodate the variant amino acid' ,
    'quartile_allatom_rmsQ3' : 'improper fold, cannot accomodate the variant amino acid' ,
    'quartile_gdtmmQ1' : 'improper fold, cannot accomodate the variant amino acid' ,
    'quartile_gdtmmQ2' : 'improper fold, cannot accomodate the variant amino acid' ,
    'quartile_gdtmmQ3' : 'improper fold, cannot accomodate the variant amino acid' ,
    'quartile_gdtmm1_1Q1' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm1_1Q2' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm1_1Q3' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm2_2Q1' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm2_2Q2' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm2_2Q3' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm3_3Q1' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm3_3Q2' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm3_3Q3' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm4_3Q1' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm4_3Q2' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm4_3Q3' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm7_4Q1' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm7_4Q2' : 'improper fold, cannot accomodate the variant' ,
    'quartile_gdtmm7_4Q3' : 'improper fold, cannot accomodate the variant' ,
    'maxsub' : 'improper fold, cannot accomodate the variant amino acid' ,
    'maxsub2.0' : 'improper fold, cannot accomodate the variant amino acid' ,
    'ddg_total' : 'highly destabilizing' ,
    'ddg_fa_atr' : 'side chains cannot \"pack\" against each other (clashes)' ,
    'ddg_fa_rep' : 'side chains cannot \"pack\" against each other (clashes)' ,
    'ddg_fa_sol' : 'unfavorable solvation (hydrophobic exposed or buried hydrophilic)' ,
    'ddg_pro_close' : 'change to/from proline not tolerated, likely conserved' ,
    'ddg_fa_pair' : 'disrupted side chain contact (ex. salt bridge)' ,
    'ddg_hbond_sr_bb' : 'disrupted backbone hydrogen bonding' ,
    'ddg_hbond_lr_bb' : 'disrupted backbone hydrogen bonding' ,
    'ddg_hbond_bb_sc' : 'disrupted backbone - side chain hydrogen bonding' ,
    'ddg_hbond_sc' : 'disrupted hydrogen bond between side chains (likely conserved)' ,
    'ddg_dslf_ss_dst' : 'destabilized disulfide' ,
    'ddg_dslf_cs_ang' : 'destabilized disulfide' ,
    'ddg_dslf_ss_dih' : 'destabilized disulfide' ,
    'ddg_dslf_ca_dih' : 'destabilized disulfide' ,
    'ddg_fa_dun' : 'unfavorable side chain conformation' ,
    'ddg_p_aa_pp' : 'unfavorable amino acid?' ,
    'ddg_ref' : 'reference?'
    }

# how different must the confidence scores be to be considered an "essential" position
ESSENTIAL_POSITION_SCORE_DIFFERENCE = .2

# simpler terms...
DELETERIOUS_PREDICTION_CRUDE_DESCRIPTION = {
    'structural conservation' : [
        'quartile_scoreQ1' ,
        'quartile_scoreQ2' ,
        'quartile_scoreQ3' ,
        'quartile_fa_atrQ1' ,    #'poor packing (clashes)' ,
        'quartile_fa_atrQ2' ,
        'quartile_fa_atrQ3' ,
        'quartile_fa_repQ1' ,
        'quartile_fa_repQ2' ,
        'quartile_fa_repQ3' ,
        'quartile_fa_solQ1' ,
        'quartile_fa_solQ2' ,
        'quartile_fa_solQ3' ,
        'quartile_fa_intra_repQ1' ,
        'quartile_fa_intra_repQ2' ,
        'quartile_fa_intra_repQ3' ,
        'quartile_pro_closeQ1' ,
        'quartile_pro_closeQ2' ,    #'disrupts a delicate proline, likely conserved' ,
        'quartile_pro_closeQ3' ,
        'quartile_hbond_sr_bbQ1' ,
        'quartile_hbond_sr_bbQ2' ,
        'quartile_hbond_sr_bbQ3' ,
        'quartile_hbond_lr_bbQ1' ,
        'quartile_hbond_lr_bbQ2' ,
        'quartile_hbond_lr_bbQ3' ,
        'quartile_hbond_bb_scQ1' ,
        'quartile_hbond_bb_scQ2' ,
        'quartile_hbond_bb_scQ3' ,
        'quartile_dslf_ss_dstQ1' ,
        'quartile_dslf_ss_dstQ2' ,
        'quartile_dslf_ss_dstQ3' ,
        'quartile_dslf_cs_angQ1' ,
        'quartile_dslf_cs_angQ2' ,
        'quartile_dslf_cs_angQ3' ,
        'quartile_dslf_ss_dihQ1' ,
        'quartile_dslf_ss_dihQ2' ,
        'quartile_dslf_ss_dihQ3' ,
        'quartile_dslf_ca_dihQ1' ,
        'quartile_dslf_ca_dihQ2' ,
        'quartile_dslf_ca_dihQ3' ,
        'quartile_ramaQ1' ,
        'quartile_ramaQ2' ,
        'quartile_ramaQ3' ,
        'quartile_omegaQ1' ,
        'quartile_omegaQ2' ,
        'quartile_omegaQ3' ,
        'quartile_fa_dunQ1' ,
        'quartile_fa_dunQ2' ,
        'quartile_fa_dunQ3' ,
        'quartile_allatom_rmsQ1' ,
        'quartile_allatom_rmsQ2' ,
        'quartile_allatom_rmsQ3' ,
        'quartile_gdtmmQ1' ,
        'quartile_gdtmmQ2' ,
        'quartile_gdtmmQ3' ,
        'quartile_gdtmm1_1Q1' ,
        'quartile_gdtmm1_1Q2' ,
        'quartile_gdtmm1_1Q3' ,
        'quartile_gdtmm2_2Q1' ,
        'quartile_gdtmm2_2Q2' ,
        'quartile_gdtmm2_2Q3' ,
        'quartile_gdtmm3_3Q1' ,
        'quartile_gdtmm3_3Q2' ,
        'quartile_gdtmm3_3Q3' ,
        'quartile_gdtmm4_3Q1' ,
        'quartile_gdtmm4_3Q2' ,
        'quartile_gdtmm4_3Q3' ,
        'quartile_gdtmm7_4Q1' ,
        'quartile_gdtmm7_4Q2' ,
        'quartile_gdtmm7_4Q3' ,
        'maxsub' ,
        'maxsub2.0' ,
        'ddg_total' ,
        'ddg_fa_atr' ,
        'ddg_fa_rep' ,
        'ddg_fa_sol' ,
        'ddg_pro_close' ,
        'ddg_hbond_sr_bb' ,
        'ddg_hbond_lr_bb' ,
        'ddg_hbond_bb_sc' ,
        'ddg_dslf_ss_dst' ,
        'ddg_dslf_cs_ang' ,
        'ddg_dslf_ss_dih' ,
        'ddg_dslf_ca_dih' ,
        'ddg_fa_dun'
        ] ,
    'disrupted chemical interactions' : [
        'aminochange' ,
        'probe_accp' ,
        'quartile_fa_pairQ1' ,
        'quartile_fa_pairQ2' ,
        'quartile_fa_pairQ3' ,
        'quartile_hbond_scQ1' ,
        'quartile_hbond_scQ2' ,
        'quartile_hbond_scQ3' ,
        'ddg_fa_pair' ,
        'ddg_hbond_sc'
        ] ,
    # ignore these for now
    'other' : [
        'pssm_variant' ,
        'pssm_native' ,
        'pssm_difference' ,
        'pssm_information_content' ,
        'quartile_p_aa_ppQ1' ,
        'quartile_p_aa_ppQ2' ,
        'quartile_p_aa_ppQ3' ,
        'quartile_refQ1' ,
        'quartile_refQ2' ,
        'quartile_refQ3' ,
        'ddg_p_aa_pp' ,
        'ddg_ref'
        ]
    }

# header for the output predictions
PREDICTION_OUTPUT_HEADER = [
    'variant' ,
    'label prediction' ,
    'confidence score' ,
    
    'structure-only label prediction' ,
    'structure-only confidence score' ,
    'sequence-only label prediction' ,
    'sequence-only confidence score' ,
    
    'exposure' ,
    'essential score' ,    # score difference
    'ddG prediction ~kcal/mol' ,
    
    'interpretation' ,
    'top ' + str( TOP_FEATURES_TO_INCLUDE ) + ' ranked features'
    ]

