import sys, time, string, re

import XGram
from XGram.Generator.Prebuilt import RNA
from XGram.Run import Iteration

if __name__ ==  "__main__":

    t0 = time.time()

    data =  XGram.PATH_DATA + "/rf00013_withtree.stk"    

    data =  XGram.PATH_DATA + "/rf00013.stk"
    data =  XGram.PATH_DATA + "/rf00503.stk"
    data =  XGram.PATH_DATA + "/rf00001.stk"
        
    ref_data = XGram.PATH_DATA + "/pfold.eg"
    tree_model = XGram.PATH_DATA + "/pfold.eg"
    
    input_model = RNA.buildPFold( strand_symmetric=True, 
                                  weak_strong = True, 
                                  copy_parameters=ref_data)
    
    xgram = XGram.XGram()
    xgram.setDebug( True )
    
    iterator = Iteration( xgram = xgram )
    iterator.setTreeModel( tree_model )
    iterator.setLogLevel( 5 )
    iterator.setDumpPatternModel( "/home/aheger/iterations/rf00001/%s.eg" )
    iterator.setDumpPatternLog( "/home/aheger/iterations/rf00001/%s.log" )
    iterator.setDumpPatternData( "/home/aheger/iterations/rf00001/%s.data" )    
    
    results = iterator.run( input_model, data )
    
    print "# retrieved %i results" % len(results)
    
    print "# final grammar:"
    
    print results[-1].getModel().getGrammar()
            