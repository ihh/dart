class ChainComponent:
    """base class for generators used to build a chain.
    They all share methods to talk to the chain generator."""
    def __init__(self, 
                 chain_generator = None, 
                 is_const = False, 
                 prefix = "", 
                 shared_frequencies = False, 
                 shared_rates = False):
        self.mChainGenerator = chain_generator
        self.mIsConst = is_const
        self.mPrefix = prefix
        self.mSharedRates = shared_rates
        self.mSharedFrequencies = shared_frequencies
        
    def registerChainGenerator(self, generator): 
       """register master generator.
       """
       self.mChainGenerator = generator
        
    def registerParameters(self, chain):
        """register parameters with chain."""
        pass

    def registerTerminals(self, chain):
        """register any new terminals with chain."""
        pass
    
    def setIsConst(self, flag):
        self.mIsConst = flag
    
    def buildGrammar(self, chain):
        """build grammar."""
        self.registerParameters( chain )
        self.registerTerminals(chain)
        
    def quoteFrequency(self, name):
        """quote a frequency.
        
        return name + a prefix if shared_frequencies is set 
        """
        if not self.mSharedFrequencies and self.mPrefix:
            return "%s_%s" % (self.mPrefix, name)
        else:
            return name
        
    def quoteRate(self, name):
        """quote a frequency.
        
        return name + a prefix if shared_rates is set 
        """
        if not self.mSharedRates and self.mPrefix:
            return "%s_%s" % (self.mPrefix, name)
        else:
            return name
