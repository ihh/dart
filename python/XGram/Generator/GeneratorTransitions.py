from XGram.Generator import ChainComponent, Rates, Expressions

#######################################################################        
## Generators of transitions 
class Transitions(ChainComponent):
    def __init__(self, 
                 rate_generator = None,
                 expression_generator = None,
                 *args, **kwargs ):
        """prefix: prefix all rate variables with this prefix."""
        ChainComponent.__init__( self, *args, **kwargs )

        if not rate_generator:
            self.mRateGenerator = Rates()
        else:
            self.mRateGenerator = rate_generator

        if not expression_generator:
            self.mExpressionGenerator = Expressions()
        else:
            self.mExpressionGenerator = expression_generator

        
    def buildGrammar(self, chain):
        """build grammar.
        
        This stub generically iterates over the alphabet and
        calls self.addTransitions for each pair of letters.
        """
        ChainComponent.buildGrammar( self, chain )
        alphabet = self.mChainGenerator.getAlphabet()
        l = len(alphabet )
        for l1 in range(0,l-1):
            a = alphabet[l1]
            for l2 in range(l1 + 1, l):                
                self.addTransitions( chain, a, alphabet[l2] )

class TransitionsConst( Transitions ):
    """default generator for transitions. 
    
    Writes uniform rates for all possible transitions within the alphabet. 
    """ 
    def __init__(self, 
                 *args, **kwargs ):
        """prefix: prefix all rate variables with this prefix."""
        Transitions.__init__( self, *args, **kwargs )
        
    def addTransitions(self, chain, l1, l2 ):
        """add transitions for l1 -> l2 and l2 -> l1."""
        chain.addTransition( l1, l2, ("rate", self.mExpressionGenerator.getRate(l1,l2)) )
        chain.addTransition( l2, l1, ("rate", self.mExpressionGenerator.getRate(l2,l1)) )        
        
    def getAlphabet(self):
        return self.mChainGenerator.getAlphabet()                    

    def registerTerminals(self, chain ):
        """register terminals.
        """
        t = []
        alphabet = self.mChainGenerator.getAlphabet()
        if len(alphabet[0]) > 1 or not self.mPrefix:            
            for x in range(len(alphabet[0])):
                t.append( "%sA%x" % (self.mPrefix, x) )
            chain.setTerminals( tuple(t) )
        else:
            chain.setTerminals( (self.mPrefix,) )
            
    def registerParameters(self, chain):
        """register parameters."""
        pass
    
    def getInitialProbability(self, emissions ):
        """return the initial probability for a set of emissions.
        """
        return 1.0 / len( self.mChainGenerator.getAlphabet() )


class TransitionsScaledConst(Transitions):
    """writes a parametric transition matrix of
    the form K * PAB, where A and B are given by
    the alphabet and PAB is a constant rate.
    
    The variable K can be specified with the optional
    "rate" option.
    """
    def __init__(self, rate = 'K' , **kwargs ):
        Transitions.__init__(self, **kwargs)
        self.mRate = rate
        
    def addTransitions(self, chain, a1, a2):
        r = "%s * r%s%s" 
        
        chain.addTransition( a1, a2, 
                             ("rate", 
                              r % (self.mRate, a1.upper(), a2.upper()) ))
        chain.addTransition( a2, a1, 
                             ("rate", 
                              r % (self.mRate, a2.upper(), a1.upper()) ))                             
        
    def registerParameters(self, chain):
        chain.addVariable( (self.mRate, self.mRateGenerator.getRate( self.mRate ) ) )
        alphabet = self.mChainGenerator.getAlphabet()            
        l = len(alphabet)
        for a in alphabet:
            for b in alphabet:
                if a == b: continue
                rate = "r%s%s" % (a.upper(), b.upper())
                if self.mIsConst:
                    chain.addConst( (rate, self.mRateGenerator.getRate(rate) ))
                else:
                    chain.addVariable( (rate, self.mRateGenerator.getRate(rate) ))
                
