import XGram.Parser
import XGram

import sys

if __name__ ==  "__main__":
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = XGram.PATH_DATA + "/sn.eg"
        
    lines= open( filename, "r" ).readlines()
    
    ## test the parser
    result = XGram.Parser.parseGrammar(lines)

    print result.getGrammar()
    
    print "# testing if I can evaluate frequencies:"
    print result.evaluateTerminalFrequencies()
    
    print result.evaluateRateMatrix()