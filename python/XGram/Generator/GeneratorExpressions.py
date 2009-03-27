#######################################################################
## Generator for values of expressions
class Expressions:
    """return an expression for use in transitions.
    
    default generator. Returns 1 for any rate."""
    def __init__(self):
        pass
    
    def getRate(self, k, b=None):
        """returns an intial value for rate k."""
        return (1.0,)

class ExpressionsDictionary(Expressions):
    """return values for rates from a dictionay."""
    def __init__(self, dictionary):
        Expression.__init__(self)
        self.mExpressionMap = dictionary
    
    def getRate(self, a, b = None ):
        """returns an intial value for rate k."""        
        if b:
            return self.mExpressionMap[a][b]
        else:
            return self.mExpressionMap[a]            

    
