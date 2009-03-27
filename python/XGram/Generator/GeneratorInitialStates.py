from XGram.Generator import ChainComponent

#######################################################################
## Generators for the initial states
class InitialStates(ChainComponent):
    def __init__(self, *args, **kwargs):
        ChainComponent.__init__(self, *args, **kwargs )
                
    def buildGrammar(self, chain):
        """build grammar for chain by adding initial states."""        
        ChainComponent.buildGrammar(self,chain)


class InitialStatesConst(InitialStates):
    """generates initial states and probability distributions
    over them. The default generator writes uniform probabilities
    for all states.
    """
    def __init__(self, *args, **kwargs):
        InitialStates.__init__( self, *args, **kwargs )
        
    def buildGrammar(self, chain):
        InitialStates.buildGrammar(self, chain)
        alphabet = self.mChainGenerator.getAlphabet()
        prob = 1.0 / len( alphabet )
        for l1 in alphabet:
            chain.addInitialState( l1, ("prob", prob) )   
       
class InitialStatesParametric(InitialStates):
    """Generate parametrized initial distribution over emissions.
    
    Initial frequencies are given as P%s where %s is the alphabet character.
    """
    
    def __init__(self, *args, **kwargs):
        InitialStates.__init__( self, *args, **kwargs )
        
    def buildGrammar(self, chain):
        """build grammar for chain by adding initial states."""        
        InitialStates.buildGrammar(self, chain)
        alphabet = self.mChainGenerator.getAlphabet()
        for l1 in alphabet:
            chain.addInitialState( l1,
                                   ("prob", self.quoteFrequency( "p%s" % l1.upper()) ))

    def registerParameters(self, chain):
        """register parameters with chain."""
        alphabet = self.mChainGenerator.getAlphabet()
        parameters = []
        p = 1.0/ len(alphabet)
        for l1 in alphabet:
            parameters.append( (self.quoteFrequency("p%s" % l1.upper()), p))
    
        chain.addParameter( parameters, 
                            is_const = self.mIsConst,
                            is_explicit = True)
        
class InitialStatesFromTransitions(InitialStates):
    """Generate parametrized initial distribution over emissions.
    
    Objects of this class talk to a GeneratorTransitions object to
    ask for probabilities.
    """
    
    def __init__(self, *args, **kwargs):
        InitialStates.__init__( self, *args, **kwargs )
        
    def buildGrammar(self, chain):
        """build grammar for chain by adding initial states."""        
        InitialStates.buildGrammar(self, chain)

        gt = self.mChainGenerator.getGeneratorTransitions()
        alphabet = self.mChainGenerator.getAlphabet()
        for l1 in alphabet:
            chain.addInitialState( l1,
                                   ("prob", gt.getInitialProbability( l1 ) ))

