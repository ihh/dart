## Generator for chains
import XGram.Generator

import XGram.Exceptions
import XGram.Model

## Generators of chains             
class Chain:
    """write a chain with uniform rates and probabilities."""
    
    def __init__(self, 
                 alphabet = None,
                 generator_initial_states = None,
                 generator_transitions = None):
        
        if not alphabet: 
            self.mAlphabet = Bio.Alphabet.IUPAC.unambiguous_dna.letters.lower()
        else:
            self.mAlphabet = alphabet
        
        if generator_transitions:
            self.mGeneratorTransitions = generator_transitions
        else:
            self.mGeneratorTransitions = XGram.Generator.TransitionsConst()
        if generator_initial_states:
            self.mGeneratorInitialStates = generator_initial_states
        else:
            self.mGeneratorInitialStates = XGram.Generator.InitialStatesConst()
        
        # register the chain generator
        self.mGeneratorTransitions.registerChainGenerator( self )
        self.mGeneratorInitialStates.registerChainGenerator( self )        
        self.mName = "full"
            
    def getAlphabet(self):
        return self.mAlphabet
    
    def setAlphabet(self, alphabet):
        self.mAlphabet = alphabet
    
    def getGeneratorTransitions(self):
        return self.mGeneratorTransitions
    
    def getGeneratorInitialStates(self):
        return self.mGeneratorInitialStates
    
    def setGeneratorTransitions(self, generator):
        self.mGeneratorTransitions = generator
        self.mGeneratorTransitions.registerGenerator( self )    
    
    def setGeneratorInitialStates(self, generator):
        self.mGeneratorInitialStates = generator
        self.mGeneratorInitialStates.registerGenerator( self )            

    
    def buildGrammar(self, chain):
        """output the grammar."""
        self.mGeneratorInitialStates.buildGrammar( chain )
        self.mGeneratorTransitions.buildGrammar( chain )
        chain.setName( self.mName )
                   
if __name__ == "__main__":
    
    chain = XGram.Model.Chain()
    
    generator = ChainCodonML( prefix = "C1_" )
    
    generator.buildGrammar( chain )
    
    print chain.getGrammar()
    