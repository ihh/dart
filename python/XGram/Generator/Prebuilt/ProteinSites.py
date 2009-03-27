## Build evolutionary trace like models.

import sys, re, string
from types import *
import Bio.Alphabet.IUPAC

import XGram.Model
from XGram.Generator import Chain, \
    TransitionsScaledConst, InitialStatesParametric, \
    GrammarLinearSequenceMultipleChains, \
    AlphabetProtein
    
def buildModel( num_bins = 3, annotate_terminals = None ):
    """build evolutionary trace like model.
    
    This model is parameterized by the number of bins. Each
    bin corresponds to a column evolving at a differnent rate.

    X chains {X: 0..num_bins-1} are created, in which transitions are parameterized as
    KX * PAB, where PAB is the rate for going from A to B (amino acid alphabet),
    and KX is rate constant for chain X.
    
    If initial_model is given, initial frequencies and rate estimates
    are taken from this model. Otherwise, uniform frequencies are chosen.
    Note, that the naming scheme in the initial model has to follow the
    one chosen here.
    """
    
    model = XGram.Model.Model()

    grammar = XGram.Model.Grammar()
    alphabet = XGram.Model.Alphabet()
    
    for x in range(num_bins):
        chain = XGram.Model.Chain()            
        generator_chain = Chain( \
            alphabet = Bio.Alphabet.IUPAC.IUPACProtein.letters.lower(),
            generator_initial_states = InitialStatesParametric( is_const = True), 
            generator_transitions = TransitionsScaledConst( rate = "K%i" % x, is_const = True)
            )
        generator_chain.buildGrammar( chain )
        chain.setTerminals( ("T%i" % x,) )
        grammar.addChain( chain )
    
    # add annotaions
    generator_grammar = GrammarLinearSequenceMultipleChains( annotate_terminals = annotate_terminals )  
    generator_grammar.buildGrammar( grammar )     

    grammar.setName( "XGramProteinSites" )
    
    generator_alphabet = AlphabetProtein()
    generator_alphabet.buildGrammar(alphabet)
     
    model.setAlphabet( alphabet )
    model.setGrammar( grammar )
    
    return model


    