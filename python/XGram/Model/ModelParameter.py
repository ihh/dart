import string, re, types
from XGram.Model import BuildExpression

class ParameterGroup:
    """a parameter group.

    This class remembers certain flags for each parameter.
    1. Is it a probability group or a rate?
    2. Is it fixed or is it variable?
    3. Is it part of pgroup/rate or param/const block.
    4. Is it a dummy variable? Dummy variables are place holder
    for parameters used in a chain, but defined in the grammar section.
    """
    def __init__(self, 
                 param,
                 is_const = False,
                 is_explicit = False ):
        
        ## set to true if 
        self.mIsConst = is_const
        self.mIsExplicit = is_explicit
        self.mIsDummy = False
        self.mIsRate = False
        
        ## key,value pairs of actual parameters
        self.mParameters = {}

        if type(param) == types.StringType:
            self.mIsDummy = True
            self.mParameters[param] = None
        elif type(param[0]) == types.StringType:
            self.mIsRate = True
            self.mParameters[param[0]] = param[1]
        else:
            for x in param:
                self.mParameters[x[0]] = x[1]
                
    def isConst(self):
        """return True if parameter is const."""
        return self.mIsConst
    
    def setConst(self, const=True):
        self.mIsConst = const

    def isRate(self):
        """return True if paramter is variable."""
        return self.mIsRate
        
    def isExplicit(self):
        """return True if parameter is part of pgroup/rate block."""
        return self.mIsExplicit
    
    def isDummy(self):
        return self.mIsDummy
    
    def items(self):
        return self.mParameters.items()
    
    def __getitem__(self, key):
        return self.mParameters[key]
         
    def __setitem__(self, key, value):
        self.mParameters[key] = value
        
    def __delitem__(self, key):
        del self.mParameters[key]
        
    def __contains__(self, key):
        return key in self.mParameters
    
    def __len__(self):
        return len(self.mParameters)
    
    def __str__(self):
        return "\n".join(self.getGrammarLines())
        
    def getGrammarLines(self):
        keys = self.mParameters.keys()
        keys.sort()
        lines = []
        if not self.mIsRate: lines.append("(")
        for x in keys:
            lines.append( "(%s %s)" % (x, self.mParameters[x]) )
        if not self.mIsRate: lines.append(")")
        return lines
