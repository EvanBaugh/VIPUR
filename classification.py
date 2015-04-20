#!/usr/bin/env python
# :noTabs=true:

"""
The VIPUR Classifier
built an a classifier object to allow easy testing of alternative values
and more complication analyses of output terms etc.

Currently constructed as a VIPURClassifier with three LogisticRegressionClassifiers
that represent different parameter sets for:
    final classifier        the combined method, gives the best performance
    structure only          using only structural analysis, precise but not highly generalized
    sequence only           using only sequence analysis, useful when no sufficient model is available

NOTE: the parameter values are NOT hardcoded and instead live in the VIPUR code release,
as several text files
these text files currently have a hardcoded, unspecified order, to properly use them,
the LABEL_DESCRIPTION list must match the feature order


ehb: this code was grafted from elsewhere, it should be cleaned up for release
"""

################################################################################
# IMPORT

# common modules
import math

# bigger modules

# custom modules
from settings import LABEL_DESCRIPTION , DELETERIOUS_PREDICTION_DESCRIPTION, FINAL_CLASSIFIER_WEIGHTS, STRUCTURE_ONLY_MODEL_WEIGHTS, SEQUENCE_ONLY_MODEL_WEIGHTS, TRAINING_SET_FEATURE_MEANS, TRAINING_SET_FEATURE_STDS , PROBE_ACCP_INTERIOR_CUTOFF , SKIP_DURING_EXPLANATION , TOP_FEATURES_TO_INCLUDE , DELETERIOUS_PREDICTION_CRUDE_DESCRIPTION , ESSENTIAL_POSITION_SCORE_DIFFERENCE

################################################################################
# METHODS

# single classifier
# evaluation of logistic regression is fast and it is informative to compare our combined results to the individual results from sequence-only and structure-only methods

class LogisticRegressionClassifier():
    def __init__( self , weights_filename , means_filename , stds_filename , description = LABEL_DESCRIPTION , scale = 1 , rescale = True , from_liblinear = False ):
        # load the weights
        f = open( weights_filename , 'r' )
        if from_liblinear:
            # specific to liblinear
            self.weights = [float( i.strip() ) for i in f.readlines()[6:]]
            self.weights = [self.weights[-1]] + self.weights[:-1]
        else:
            self.weights = [float( i.strip() ) for i in f.xreadlines()]
        f.close()

        # load the means
        f = open( means_filename , 'r' )
        self.means = [float( i.strip() ) for i in f.xreadlines()]
        f.close()
    
        # load the stds
        f = open( stds_filename , 'r' )
        self.stds = [[float( k ) for k in [j for j in i.strip().split( ' ' ) if j.strip()] if float( k )][0] for i in f.xreadlines()]
        f.close()

        # save this for later
        if description == None:
            description = ['']*len( self.weights )
        # ...pretty sure this is a Python trap, having a default be a list in a constructor...
        self.description = description
        # should be indexed THE SAME as the weights
        self.scale = scale    # only use this for display :)
        self.rescale = rescale

        # calculate the bias
        self.bias = self.weights[0]
        self.weights = self.weights[1:]
        self.relevant_indices = [i for i in xrange( len( self.weights ) ) if self.weights[i]]

        # adjust the bias for the normalized features
#        temp = [self.weights[i]*self.means[i]*self.stds[i] for i in xrange( len( self.means ) )]
#        self.bias -= sum( [temp[i] for i in self.relevant_indices] )
        # shouldn't need relevant_indices, others are 0 already!
        # NO! leave the weights "as is", transform the input vector

        # and the weights
#        self.weights = [self.weights[i]*self.stds[i] for i in xrange( len( self.weights ) )]
        
        # reset the scale
        self.scale = min( [abs( self.weights[i] ) for i in self.relevant_indices] )


        # internal display method :)
        self.summarize_classifier( self.rescale )


    # simple enough, I like doing things this way okay <:]
    def __str__( self , rescale = None ):
        if rescale == None:
            rescale = self.rescale
        # for consistency
        scale = 1
        if rescale:
            scale = self.scale

        text = 'bias'.ljust( 30 ) + str( self.bias/scale )
#        print self.description
#        print self.weights
#        print len( self.description ) , len( self.weights )
        for i in self.relevant_indices:
            text += '\n'+ self.description[i].ljust( 30 ) + str( self.weights[i]/scale )
        
        return text

    # display the classifier properties
    def summarize_classifier( self , rescale = None ):
        print self.__str__( rescale )
    
    def reorder_input_vector( self , input_vector ):
        if isinstance( input_vector , dict ):
            # for now, assume it will work, however fragile this is
#            input_vector = [float( input_vector[i] ) for i in self.description]
            input_vector = [float( input_vector[i] )  if i in input_vector.keys() else  0 for i in self.description]
        return input_vector
    
    # wrapped, in case you need it for other purposes
    def normalize_input( self , input_vector ):
        if isinstance( input_vector , dict ):
            input_vector = self.reorder_input_vector( input_vector )
        return [(input_vector[i] - self.means[i])*self.stds[i]  if i in self.relevant_indices else  0 for i in xrange( len( input_vector ) )]
    
    # actual classification, input_vector should be indexed the same as weights
    def classify_instance( self , input_vector , base = math.e , normalize = True ):
        if isinstance( input_vector , dict ):
            input_vector = self.reorder_input_vector( input_vector )
        # normalize the input vector
        # ...can save a bunch of useless computation
        # no need to even do this outside "relevant_indices"
#        input_vector = [(input_vector[i] - self.means[i])*self.stds[i]  if i in self.relevant_indices else  0 for i in xrange( len( input_vector ) )]
        if normalize:
            input_vector = self.normalize_input( input_vector )
    
        # calculate all terms
        summary = {'terms' : [self.bias] + [self.weights[i]*input_vector[i] for i in self.relevant_indices]}
        
        # sum it, convert to prob
        summary['score'] = sum( summary['terms'] )
        summary['P'] = (1 + base**(-summary['score']))**-1
        summary['P*'] = 2*abs(summary['P'] - .5)

        # oh, and the label too!
        summary['label'] = 'deleterious'*(summary['P'] >= .5) + 'neutral'*(summary['P'] < .5)
        
        return summary

    # classify it AND summarize the output plz
    def summarize_classification( self , input_vector , base = math.e , rescale = None ):
        if rescale == None:
            rescale = self.rescale
        scale = 1
        if rescale:
            scale = self.scale
        
        # actually classify
        summary = self.classify_instance( input_vector , base )

        # don't skip the bias
        terms = summary['terms']
        descriptions = ['*bias'] + [self.description[i] for i in self.relevant_indices]
        ranked = sorted( range( len( terms ) ) , key = lambda x : -terms[x] )
        # so...these are currently in "terms" indices
        # make them into the overal indices -err- no need :)
        positive = [(terms[i]/scale , descriptions[i]) for i in ranked if abs( terms[i] ) > 1e-7 and terms[i] >= 0]
        negative = [(terms[i]/scale , descriptions[i]) for i in ranked if abs( terms[i] ) > 1e-7 and terms[i] < 0]

        # hmmm...display
        print '='*3 + ' deleterious ' + '='*3
        text = ''
        for i in positive:
            print i[1].ljust( 30 ) + str( i[0] )
#        print '='*10 +' '+ str( self.bias/scale ) +' '+ '='*10 +' '+ str( summary['score'] )
        print '='*10 +' '+ str( summary['score'] ) +' '+ '='*10 +' '+ str( summary['label'] ) +', '+ str( summary['P'] )
        for i in negative:
            print i[1].ljust( 30 ) + str( i[0] )
        print '='*3 + ' neutral ' + '='*3
    
        return summary


# merge all the classifiers
# create a single object to hold all three

class VIPURClassifier():
    # give default paths
    def __init__( self , combined_weights_filename = FINAL_CLASSIFIER_WEIGHTS ,
            structure_only_weights_filename = STRUCTURE_ONLY_MODEL_WEIGHTS , sequence_only_weights_filename = SEQUENCE_ONLY_MODEL_WEIGHTS ,
            means_filename = TRAINING_SET_FEATURE_MEANS , stds_filename = TRAINING_SET_FEATURE_STDS ,
            description = LABEL_DESCRIPTION , scale = 1 , rescale = False , from_liblinear = False , top_features = TOP_FEATURES_TO_INCLUDE ):

        # load the individual classifiers
        print 'loading VIPUR weights'
        self.combined_classifier = LogisticRegressionClassifier( combined_weights_filename , means_filename = means_filename , stds_filename = stds_filename , description = description , scale = scale , rescale = rescale , from_liblinear = from_liblinear )
        print 'loading structure-only weights'
        self.structure_classifier = LogisticRegressionClassifier( structure_only_weights_filename , means_filename = means_filename , stds_filename = stds_filename , description = description , scale = scale , rescale = rescale , from_liblinear = from_liblinear )
        print 'loading sequence-only weights'
        self.sequence_classifier = LogisticRegressionClassifier( sequence_only_weights_filename , means_filename = means_filename , stds_filename = stds_filename , description = description , scale = scale , rescale = rescale , from_liblinear = from_liblinear )

        # keep doing it this way?
        self.scale = scale
        self.rescale = rescale
        
        self.top_features = top_features

    # not currently supported
#    def reorder_input_vector( self , input_vector ):
#        if isinstance( input_vector , dict ):
            # for now, assume it will work, however fragile this is
#            input_vector = [float( input_vector[i] ) for i in self.description]
#            input_vector = [float( input_vector[i] )  if i in input_vector.keys() else  0 for i in self.description]
#        return input_vector
    
    # classify it AND summarize the output plz
    def combined_classification( self , input_vector , base = math.e , rescale = None ):
        if rescale == None:
            rescale = self.rescale
        scale = 1
        if rescale:
            scale = self.scale

        # not currently supported
#        if isinstance( input_vector , dict ):
#            input_vector = self.reorder_input_vector( input_vector )
        
        # support sequence only classification?
        
        # actually classify
        summary = self.combined_classifier.classify_instance( input_vector , base )
        # ...saving this output?
        structure = self.structure_classifier.classify_instance( input_vector , base )
        sequence = self.sequence_classifier.classify_instance( input_vector , base )
        
        return summary

    # add properties etc.
    def interpret_classification( self , input_vector , top_features = None , deleterious_description = DELETERIOUS_PREDICTION_DESCRIPTION , skip_during_explanation = SKIP_DURING_EXPLANATION , display = True ):
        scale = self.scale    # for now
        if not top_features:
            top_features = self.top_features

        # not currently supported
#        if isinstance( input_vector , dict ):
#            input_vector = self.reorder_input_vector( input_vector )

        # make the predictions
        # ...redundant with the above method...
        # support sequence only classification?
        combined_prediction = self.combined_classifier.classify_instance( input_vector )
        structure_prediction = self.structure_classifier.classify_instance( input_vector )
        sequence_prediction = self.sequence_classifier.classify_instance( input_vector )
        
        combined_prediction['structure label'] = structure_prediction['label']
        combined_prediction['structure P'] = structure_prediction['P']
        combined_prediction['sequence label'] = sequence_prediction['label']
        combined_prediction['sequence P'] = sequence_prediction['P']
    
        # add in surface etc.? separate analysis
        # just use ACCP 12.5 as the cutoff
        exposure = 'surface'
        if isinstance( input_vector , list ) and len( input_vector ) > 5 and input_vector[5] > PROBE_ACCP_INTERIOR_CUTOFF:    # I LOATHE THIS! magic number AND hardcoded index!!!
            exposure = 'interior'
        elif isinstance( input_vector , dict ) and float( input_vector['probe_accp'] ) > PROBE_ACCP_INTERIOR_CUTOFF:
            exposure = 'interior'
            
        combined_prediction['exposure'] = exposure

        # save this analysis for the protocol itself...
        # consider pairs
        # for now, just combined but not structure
#        conservation = ''
#        if combined_prediction['label'] == 'deleterious' and not structure_prediction['label'] == 'deleterious':
            # sequence terms make it deleterious, not detected by structure
#            if exposure == 'surface':
#                conservation = 'potential interaction site'
#            else:
#                conservation = 'restricted chemical identity but does not destabilize the monomer (conserved, predicted deleterious, and not on the surface)'
#        combined_prediction['conservation'] = conservation
    
        # take the top few ranked terms
        terms = combined_prediction['terms']
        descriptions = ['*bias'] + [self.combined_classifier.description[i] for i in self.combined_classifier.relevant_indices]
        ranked = sorted( range( len( terms ) ) , key = lambda x : -terms[x] )
        # so...these are currently in "terms" indices
        # make them into the overal indices -err- no need :)
        positive = [(terms[i]/scale , descriptions[i]) for i in ranked if abs( terms[i] ) > 1e-7 and terms[i] >= 0]
        negative = [(terms[i]/scale , descriptions[i]) for i in ranked if abs( terms[i] ) > 1e-7 and terms[i] < 0]
        
        combined_prediction['explanation'] = ['appears neutral']
        if combined_prediction['label'] == 'deleterious':
#            combined_prediction['explanation'] = [deleterious_description[positive[i][1]] for i in  if not 'bias' in positive[i][1]]
            combined_prediction['explanation'] = [deleterious_description[i[1]] + ' [' + i[1] +']' for i in positive if not 'bias' in i[1] and not i[1] in skip_during_explanation][:top_features]

        # present the details
        if display:
            print combined_prediction['label'] +', '+ str( combined_prediction['P'] ) + ' confidence'
            print combined_prediction['exposure']# +', '*bool( conservation ) + conservation
            if 'explanation' in combined_prediction.keys():
                for i in combined_prediction['explanation']:
                    print '\t'+ i
#            print '='*10
            print    # for interactive display

        return combined_prediction

################################################################################
# ADDITIONAL INTERPRETATION

# adds more details into the prediction dict
def provide_additional_interpretation( prediction ):
    # special conservation case
    interpretation = '?'    # default, for debugging
    if not prediction['label'] == 'deleterious':
        interpretation = 'neutral'
    elif prediction['structure label'] == 'deleterious' and prediction['exposure'] == 'surface' and (prediction['P'] - prediction['structure P']) >= ESSENTIAL_POSITION_SCORE_DIFFERENCE:
        interpretation = 'potential interaction site'
    elif prediction['explanation']:
        interpretation = [i for i in DELETERIOUS_PREDICTION_CRUDE_DESCRIPTION.keys() if prediction['explanation'][0].split( ' ' )[-1].strip( '[]' ) in DELETERIOUS_PREDICTION_CRUDE_DESCRIPTION[i]][0]    # shoule only be 1 match
    
    prediction['interpretation'] = interpretation
    
    # any other details?

################################################################################

# create the actual instance
VIPUR_classifier = VIPURClassifier()


