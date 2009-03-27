"""Test system using Kimura's two parameter model.
"""

import sys, time, string, re

from XGram.Model import Annotation
import XGram

if __name__ == "__main__":
    
    input_grammar_file = XGram.PATH_DATA + "/kimura2.eg"
    input_data_file1 = XGram.PATH_DATA + "/dpse_dmel.stk"
    input_data_file2 = XGram.PATH_DATA + "/cluster_5750.stk"

    input_model = XGram.Parser.parseGrammar( open(input_grammar_file, "r").readlines() )
    xgram = XGram.XGram()

    for data in input_data_file1, input_data_file2:        
        
        print "data=", data
        for alpha,beta in ( (4.0, 1.0), (3.0, 1.0), (1.0, 1.0) ):
        
            input_model.mGrammar.setParameter( "alpha", alpha )
            input_model.mGrammar.setParameter( "beta", beta )        

            result = xgram.train( input_model, data = data )
        
            # print result.getModel().getGrammar()
    
            print "iterations=", result.getNumIterations(), \
                "runtime=", result.getRunTime(), \
                "likelihood=", result.getLogLikelihood(), \
                "params=", result.getModel().mGrammar.getParameters(), \
                "input=", input_model.mGrammar.getParameters()
            
    
