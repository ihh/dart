import sys, time, string, re

import XGram
from  XGram.Generator.Prebuilt import RNA

if __name__ ==  "__main__":

    generated_model = XGram.Generator.RNA.buildNullRNA()
    
    print generated_model.getGrammar()
    
    generated_model = XGram.Generator.RNA.buildPFold()
    
    print generated_model.getGrammar()