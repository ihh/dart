"""test the PFold grammar."""

import sys, time, string, re

import XGram
from XGram.Generator.Prebuilt import RNA
import XGram.Run
import XGram.Parser

if __name__ ==  "__main__":

    t0 = time.time()

    ## file with initial values for pfold
    parameters = XGram.PATH_DATA + "/pfold.eg"
    data = XGram.PATH_DATA + "/ncrna.FBtr0081064.stk"
    symmetric = False
    
    xgram = XGram.XGram()
    
    generated_model = RNA.buildPFold()
    
    data_model = XGram.Parser.parseGrammar( open( parameters, "r").readlines() )

    ## load PFold parameterization
    generated_model.mGrammar.copyParameters( data_model.mGrammar,
                                             remove_missing = True)
    
    print data_model.getGrammar()    
    print generated_model.getGrammar()
    
    result = xgram.annotate( generated_model, data,
                             tree_model = generated_model )
    
    print "".join(result.mLog)
    print "".join(result.mData)
    