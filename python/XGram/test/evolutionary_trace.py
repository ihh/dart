import string 

import XGram

from XGram.Model import *
from XGram.Generator.Prebuilt import ProteinSites

if __name__ ==  "__main__":

    filename_parameters = XGram.PATH_DATA + "/evolutionary_trace.eg"

    input_model= XGram.Parser.parseGrammar( open(filename_parameters).readlines() )
    
    num_bins = 3
    ## build a map of annotations:
    annotate_terminals = {}
    labels = string.letters.upper()
    for x in range(num_bins):
        t = "T%i" % x
        annotate_terminals[ t ] = (Annotation( row = "RATES",
                                               column = t,
                                               label = labels[x % len(labels)]),)
    
    ## generate the model
    generated_model = ProteinSites.buildModel( num_bins = num_bins,
                                               annotate_terminals = annotate_terminals )

    ## build a map of names in evolutionary_trace.eg to the new values
    map_this2other = {}
    letters = map(string.upper, generated_model.mAlphabet.getTokens())    
    # map frequencies: piA -> pA
    for l in letters:
        map_this2other[ 'p%s' % l ] = 'pi%s' % l
    
    # map transitions: rAB -> rAB: nothing to be done
    
    # map rates:
    # s1 -> K0 and 
    # w1 -> wT0
    for x in range(num_bins):
        map_this2other[ 'K%i' % x ] = 's%i' % (x+1)
        map_this2other[ 'wT%i' % x ] = 'w%i' % (x+1)
    
    ## set initial values according to the input model
    generated_model.mGrammar.copyParameters( input_model.mGrammar, 
                                             map_this2other = map_this2other,
                                             ignore_missing = True )

    print generated_model.getGrammar()
    

    
    