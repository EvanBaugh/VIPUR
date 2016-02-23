#!/usr/bin/env python
# :noTabs=true:

"""
amusingly, when I built this as 3 sections, this was also in the pre_processing file
"""

################################################################################
# IMPORT

# common modules
import os

# bigger modules

# custom modules
from vipur_settings import AMINOCHANGE_GROUPS , PREDICTION_OUTPUT_HEADER

from pre_processing import load_task_summary

from psiblast_feature_generation import load_numbering_map , extract_pssm_from_psiblast_pssm
from probe_feature_generation import extract_accp_from_probe
from rosetta_feature_generation import extract_score_terms_from_ddg_monomer , extract_scores_from_scorefile , extract_quartile_score_terms_from_scorefiles

from classification import VIPUR_classifier , provide_additional_interpretation

################################################################################
# MAIN POSTPROCESSING

# ...not much it can do if the features are incomplete...
# though I suppose that is the point, feature generation should take care of this

def run_postprocessing( task_summary_filename , sequence_only = False ):
    # load the tasks
    if isinstance( task_summary_filename , str ) and os.path.isfile( task_summary_filename ):
        task_summary = load_task_summary( task_summary_filename )
    else:
        # input is what we want
        task_summary = task_summary_filename

    # residue map
    residue_map = {}
    if 'numbering_map' in task_summary['filenames'].keys():
#        filename = task_summary['out_path'] +'/'+ task_summary['filenames']['numbering_map']
        filename = task_summary['filenames']['numbering_map']
        residue_map = load_numbering_map( filename )
    
    
    # extract the features
    
    # determine the "aminochange" values
    evaluate_aminochange = lambda nat , var : 2 - int( [i for i in AMINOCHANGE_GROUPS if nat in i][0] == [i for i in AMINOCHANGE_GROUPS if var in i][0] ) - int( nat == var )

    # do this mapping once, rather than locate each task when needed
    important_tasks = {}
    for i in task_summary['commands']:
        if i['feature'] == 'psiblast':
            if 'psiblast' in important_tasks.keys():
                raise Exception( '??? duplicate psiblast task ???' )
            important_tasks['psiblast'] = i
        elif i['feature'] == 'probe':
            if 'probe' in important_tasks.keys():
                raise Exception( '??? duplicate probe task ???' )
            important_tasks['probe'] = i
        elif i['feature'] == 'ddg_monomer':
            if 'ddg_monomer' in important_tasks.keys():
                raise Exception( '??? duplicate ddg_monomer task ???' )
            important_tasks['ddg_monomer'] = i
        elif i['feature'] == 'relax_native_rescore':
            if 'relax_rescore_task' in important_tasks.keys():
                raise Exception( '??? duplicate relax_native_rescore task ???' )
            important_tasks['relax_native_rescore'] = i
        elif i['feature'] == 'relax_rescore':
            if 'relax_rescore_' + i['variant'] in important_tasks.keys():
                raise Exception( '??? duplicate relax_rescore_' + i['variant'] +' task ???' )
            important_tasks['relax_rescore_' + i['variant']] = i
        #else:
            # misc/relax runs...
    
    # psiblast
#    psiblast_task = [i for i in task_summary['commands'] if i['feature'] == 'psiblast']
#    if not psiblast_task or not 'run' in psiblast_task[0].keys() or not 'success' in psiblast_task[0]['run']:
    if not 'psiblast' in important_tasks.keys() or not 'run' in important_tasks['psiblast'].keys() or not 'success' in important_tasks['psiblast']['run']:
        raise Exception( 'psiblast did not complete successfully!!!' )

    pssm = extract_pssm_from_psiblast_pssm( important_tasks['psiblast']['output_filename'] )
    if not residue_map:
        residue_map = dict( [(str( i ) , str( i )) for i in pssm.keys()] )

    for i in task_summary['variants'].keys():
        if 'failed' in task_summary['variants'][i].keys():
            continue
    
        nat = i.split( '_' )[-1]
        pos = nat[1:-1]
        pos = int( residue_map[pos] ) + 1    # shifted from parsing
        var = nat[-1]
        nat = nat[0]
        
        task_summary['variants'][i]['features'] = {
            'aminochange' : evaluate_aminochange( nat , var ) ,
            'pssm_native' : pssm[pos]['log-likelihood'][nat] ,
            'pssm_variant' : pssm[pos]['log-likelihood'][var] ,
            'pssm_difference' : pssm[pos]['log-likelihood'][nat] - pssm[pos]['log-likelihood'][var] ,
            'pssm_information_content' : pssm[pos]['information content'] ,
            }

    if not sequence_only:
        # probe
#        probe_task = [i for i in task_summary['commands'] if i['feature'] == 'probe']
#        if not probe_task or not 'run' in probe_task[0].keys() or not 'success' in probe_task[0]['run']:
        if not 'probe' in important_tasks.keys() or not 'run' in important_tasks['probe'].keys() or not 'success' in important_tasks['probe']['run']:
            raise Exception( 'probe did not complete successfully!!!' )

        if not 'probe_positions' in task_summary['other']:
            raise Exception( 'task summary was not written properly' )
        positions = task_summary['other']['probe_positions']
        
        # extract the PROBE feature
        accp = extract_accp_from_probe( important_tasks['probe']['output_filename'] )
        # positions is a parallel list of positions
        accp_dict = dict( [(positions[i] , accp[i]) for i in xrange( len( accp ) )] )
        for i in task_summary['variants'].keys():
            if 'failed' in task_summary['variants'][i].keys():
                continue

            pos = i.split( '_' )[-1][1:-1]
            if pos in accp_dict.keys():
                task_summary['variants'][i]['features']['probe_accp'] = accp_dict[pos]
            else:
                continue

        # ddg monomer
#        ddg_monomer_task = [i for i in task_summary['commands'] if i['feature'] == 'ddg_monomer']
#        if not ddg_monomer_task or not 'run' in ddg_monomer_task[0].keys() or not 'success' in ddg_monomer_task[0]['run']:
        if not 'ddg_monomer' in important_tasks.keys() or not 'run' in important_tasks['ddg_monomer'].keys() or not 'success' in important_tasks['ddg_monomer']['run']:
            raise Exception( 'ddg_monomer did not complete successfully!!!' )
        ddg_monomer_dict = extract_score_terms_from_ddg_monomer( out_filename = important_tasks['ddg_monomer']['output_filename'] )
        # the residues from this will be "pose numbered", need to map back to the PDB numbering
        header = ddg_monomer_dict['description']
        for i in task_summary['variants'].keys():
            if 'failed' in task_summary['variants'][i].keys():
                continue

            mutation = i.split( '_' )[-1]
            new_key = mutation[0] + str( int( residue_map[mutation[1:-1]] ) + 1 ) + mutation[-1]
            for j in xrange( len( header ) ):
#                print mutation
#                print new_key
#                print ddg_monomer_dict.keys()
                task_summary['variants'][i]['features']['ddg_' + header[j]] = float( ddg_monomer_dict[new_key][j] )
                # added "ddg_" for legacy compatability, is artibrary, make more informative
        
        # relax
        # get native relax reference scores
#        native_task = [i for i in task_summary['commands'] if i['feature'] == 'relax_native_rescore']# and i['variant'] == 'native']
#        if not native_task or not 'run' in native_task[0].keys() or not 'success' in native_task[0]['run']:
        if not 'relax_native_rescore' in important_tasks.keys() or not 'run' in important_tasks['relax_native_rescore'].keys() or not 'success' in important_tasks['relax_native_rescore']['run']:
            raise Exception( 'relax for the native structure reference did not complete successfully!!!' )
        native_scorefile_dict = extract_scores_from_scorefile( important_tasks['relax_native_rescore']['output_filename'] )    # save time, only parse this once
        
        for i in task_summary['variants'].keys():
            if 'failed' in task_summary['variants'][i].keys():
                continue

            # this file SHOULD exist for every variant
            # if it does not, there was a problem earlier
            variant_structure = task_summary['variants'][i]['structure_filename']
            # note: this current version assumes 1:1 for variant structures and
            #    variants...could change in the future...

            mutation = i.split( '_' )[-1]
#            rescore_task = [j for j in task_summary['commands'] if j['feature'] == 'relax_rescore' and j['variant'] == mutation]
#            if not rescore_task or not 'run' in rescore_task[0].keys() or not 'success' in rescore_task[0]['run']
            if not 'relax_rescore_' + mutation in important_tasks.keys() or not 'run' in important_tasks['relax_rescore_' + mutation].keys() or not 'success' in important_tasks['relax_rescore_' + mutation]['run']:
                raise Exception( 'rescore (or relax?) did not complete successfully!!!' )
    
            # extract features, use the quartile method to extract comparisons
            # between the native and variant score distributions
            quartile_scores = extract_quartile_score_terms_from_scorefiles( important_tasks['relax_rescore_' + mutation]['output_filename'] , native_scorefile_dict )
            # make sure there is not overlap
            overlaps = [j for j in quartile_scores.keys() if j in task_summary['variants'][i]['features'].keys()]
            if overlaps:
                raise IOError( '!!?! somehow your scorefile has score terms that match the name(s) of features that already exist!!?! something has gone horribly wrong...' + str( overlaps ) )
            task_summary['variants'][i]['features'].update( quartile_scores )
    
    
    ############
    # run the classifier
    
    # combine and prune feature vector
    print 'classifying the variants as \"deleterious\" or \"neutral\"'
    
    # calculate scores using the logistic regression parameters
    predictions = {}
    for i in task_summary['variants'].keys():
        if 'failed' in task_summary['variants'][i].keys():
            # add as a failure
            predictions[i] = '\t'.join( [i.split( '/' )[-1] , i , 'failed: ' + task_summary['variants'][i]['failed'] , str( 0 )] )    # added feb 2016 to match below
            continue

        print 'Variant: ' + i.split( '/' )[-1]
        if sequence_only:
            # just the sequence only classifier
            # standardize output, add in empty values
            interpretation = VIPUR_classifier.sequence_classifier.summarize_classification( task_summary['variants'][i]['features'] )
            interpretation['structure label'] = 'NA'
            interpretation['structure P'] = 'NA'
            interpretation['sequence label'] = interpretation['label']
            interpretation['sequence P'] = interpretation['P']

            interpretation['exposure'] = 'NA'
            essential_score = 'NA'
            task_summary['variants'][i]['features']['ddg_total'] = 'NA'
            interpretation['interpretation'] = 'NA'
            interpretation['explanation'] = ['NA']
        else:
            # full prediction
            interpretation = VIPUR_classifier.interpret_classification( task_summary['variants'][i]['features'] )
            provide_additional_interpretation( interpretation )    # adds directly into the dict
            essential_score = interpretation['P'] - interpretation['structure P']

        task_summary['variants'][i]['interpretation'] = interpretation
        
        # output append
        predictions[i] = '\t'.join( [str( j ) for j in [
            i.split( '/' )[-1] ,    # added jan 2016
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
            task_summary['variants'][i]['features']['ddg_total'] ,
            
            interpretation['interpretation'] ,
            ', '.join( interpretation['explanation'] )
#            + interpretation['terms']    # not now...
            ] ] )
        
    # write output predictions etc.
    prediction_filename = task_summary['filenames']['prediction_filename']
    if not prediction_filename:
        prediction_filename = task_summary['filenames']['pdb_filename']
        prediction_filename = prediction_filename.replace( '.pdb' , '' ).replace( '.fa' , '' ).replace( '.fasta' , '' ) + '.predictions'
        task_summary['filenames']['prediction_filename'] = prediction_filename

    f = open( prediction_filename , 'w' )
    # very shameful, hardcoded indices again...
    temp = predictions.values()
    temp = sorted( temp , key = lambda x : 1/(float( x.split( '\t' )[2] ) + 1e-7) )    # make sure no /0
    f.write( '\t'.join( PREDICTION_OUTPUT_HEADER ) +'\n'+ '\n'.join( temp ) )
    f.close()

    # debug analyze time
#    for i in xrange( len( debug_time ) - 1 ):
#        dt = debug_time[i + 1][1] - debug_time[i][1]
#        minutes = int( dt/60 )
#        dt = dt - minutes*60
##        print debug_time[i + 1][0].ljust( 20 ) + str( round( debug_time[i + 1][1] - debug_time[i][1] , 4 ) ) +'s'
#        print debug_time[i + 1][0].ljust( 40 ) + str( minutes ) + 'min ' + str( round( dt , 4 ) ) +'s'

    # write out task summary? any new details added?
    # meh, should be derivative analysis using post_processing anyway...
    
    return task_summary


