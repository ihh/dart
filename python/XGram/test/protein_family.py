"""test script for building a protein family model.
"""

import sys, time, string, re

import XGram
from XGram.Generator.Prebuilt import ProteinFamily

if __name__ ==  "__main__":

    t0 = time.time()

    data = XGram.PATH_DATA + "/pf00017_withtree.stk"
    data = XGram.PATH_DATA + "/pf00017_truncated.stk"    

    model = ProteinFamily.buildProteinFamily( data )
    
    print model.getGrammar()
    
    sys.stdout.flush()
    
    xgram = XGram.XGram()
    xgram.setMinIncrement( 0.01 )
    xgram.setDebug( True )
    
    result = xgram.train( model, data )
    
    print result.getModel().getGrammar()
    
    print "\n".join(result.getData())
    
    print "\n".join(result.getLog())
    
    print "execution lasted %i seconds." % result.getRunTime()
    