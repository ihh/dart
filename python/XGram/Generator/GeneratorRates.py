#######################################################################
## Generator for initial values for rates
class Rates:
    """Generator for rates in the param/const section.
    
    default generator. Returns 1 for any rate.
    """
    def __init__(self):
        pass
    
    def getRate(self, k, b=None):
        """returns an intial value for rate k."""
        return 1.0

class RatesDictionary(Rates):
    """return values for rates from a dictionay."""
    def __init__(self, dictionary):
        Rates.__init__(self)
        self.mRateMap = dictionary
    
    def getRate(self, a, b = None ):
        """returns an intial value for rate k."""        
        if b:
            return self.mRateMap[a][b]
        else:
            return self.mRateMap[a]            

class RatesSmall(Rates):
    """Generator for rates in the param/const section.
    
    Returns 0.01 for any rate.
    """
    def __init__(self):
        Rates.__init__(self)
    
    def getRate(self, k, b=None):
        """returns an intial value for rate k."""
        return 0.01
