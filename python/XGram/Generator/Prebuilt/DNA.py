"""This module provides standard nucleotide based models.
"""

import sys, re, string
from types import *
import Bio.Alphabet.IUPAC

from XGram.Generator import Transitions, TransitionsConst, \
    InitialStates, InitialStatesConst, InitialStatesParametric, \
    InitialStatesFromTransitions, Chain, \
    GrammarLinearSequence, \
    GrammarLinearSequenceBlocks, \
    GrammarLinearSequenceMultipleChains, \
    AlphabetDNA
    
import XGram.Model
from XGram import Exceptions

import Bio
import Bio.Alphabet
import Bio.Alphabet.IUPAC

class ChainDNA( Chain ):
    """generator for codon models."""
    def __init__( self,
                  generator_initial_states = None,
                  generator_transitions = None,
                  *args, **kwargs ):

        if not generator_transitions:
            generator_transitions = Transitions( *args, **kwargs )
        
        # the alphabet is given by the codon table, but the states
        # are separated.
        alphabet = ('a', 'c', 'g', 't' )
        Chain.__init__(self, 
                       alphabet = alphabet,
                       generator_initial_states = generator_initial_states,
                       generator_transitions = generator_transitions, 
                       )
    
        self.mName = "dna"
        
    def buildGrammar(self, chain):
        """output the grammar."""
        Chain.buildGrammar( self, chain )
        chain.setUpdatePolicy( "parametric" )
        
class TransitionsK80( Transitions ):
    """transitions for Kimura80 2-parameter model.
    
    Transitions happen with rate alpha.
    Transversions happen with rate beta.
    """
    
    def __init__(self, *args, **kwargs):
        Transitions.__init__( self, *args, **kwargs )

    def addTransitions(self, chain, aa1, aa2):
        """decides, if a transition is to be added and adds it.
        Both directions will be added.
        """
        is_transition = (aa1,aa2) in (("a","g"), ("g","a"), ("t","c"), ("c","t"))        

        if is_transition:
            rate = "alpha"
        else:
            rate = "beta"

        chain.addTransition( aa1, aa2, ("rate", self.quoteRate(rate)))  
        chain.addTransition( aa2, aa1, ("rate", self.quoteRate(rate)))          

    def registerParameters(self, chain ):
        """register parameters for the chain to built.
        
        Use small values, because xrate has better convergence if it starts
        from parameter values that are too small.
        """

        chain.addParameter( (self.quoteRate("alpha"), 0.1), is_explicit = True )
        chain.addParameter( (self.quoteRate("beta"), 0.1), is_explicit = True )        

    def registerTerminals(self, chain):
        chain.setTerminals( (self.mPrefix,) )

class TransitionsGTR( Transitions ):
    """transitions for GTR 9-parameter model.

    see Felsenstein 2004, pp205
    """
    
    mRates = { ("a","c") : "beta",
               ("a","g") : "alpha",
               ("a","t") : "gamma",
               ("g","c") : "delta",
               ("g","t") : "epsilon",
               ("c","t") : "theta" }
    
    def __init__(self, *args, **kwargs):
        Transitions.__init__( self, *args, **kwargs )

    def addTransitions(self, chain, aa1, aa2):
        """decides, if a transition is to be added and adds it.
        Both directions will be added.
        """

        if (aa1,aa2) in self.mRates:
            rate = self.quoteRate( self.mRates[(aa1,aa2)])
        else:
            rate = self.quoteRate( self.mRates[(aa2,aa1)])
            
        chain.addTransition( aa1, aa2, ("rate", rate + " * %s" % (self.quoteFrequency( "p%s" % aa2.upper()))))  
        chain.addTransition( aa2, aa1, ("rate", rate + " * %s" % (self.quoteFrequency( "p%s" % aa1.upper()))))          

    def registerParameters(self, chain ):
        """register parameters for the chain to built.
        
        Use small values, because xrate has better convergence if it starts
        from parameter values that are too small.
        """

        chain.addParameter( (self.quoteRate( "alpha"), 0.1), is_explicit = True )
        chain.addParameter( (self.quoteRate( "beta"), 0.1), is_explicit = True )        
        chain.addParameter( (self.quoteRate( "gamma"), 0.1), is_explicit = True )        
        chain.addParameter( (self.quoteRate( "delta"), 0.1), is_explicit = True )        
        chain.addParameter( (self.quoteRate( "epsilon"), 0.1), is_explicit = True )        
        chain.addParameter( (self.quoteRate( "theta"), 0.1), is_explicit = True )        

    def registerTerminals(self, chain):
        chain.setTerminals( (self.mPrefix,) )

def buildModel( substitution_model = "jc69",
                grammar_type = "linear",
                num_blocks = 1,
                shared_frequencies = False,
                shared_rates = False,
                annotate_terminals = False,
                explicit_extension = False,
                ):
    """build a standard nucleotide models.
    
    model: can be either one of
        jc69: Jukes-Cantor model
        k80: Kimura 2-parameter model
        gtr: general time reversible model (9 parameters)
        rev: general 12 parameter model
        
    exclicit_extension: add an explicit extension probability
    (called ext and not_ext) to the grammar.

    grammar_type can be one of "linear", "linear-blocks".
            linear: simple linear sequence with gaps
            linear-blocks: blocks of linear sequence

    num_blocks: the number of blocks for the "linear-blocks" grammar.
            
    shared_frequencies: the nucleotide frequencies are shared between blocks.
        The default (false) has separate nucleotide frequencies for each block.
            
    shared_rates: the rates between blocks are shared. The default (false)
        has separate rates for each block.
            
    annotate_terminals:
        The linear-blocks grammar allows to annotate the input
        alignments. Provide a hash mapping the emitted states to
        Annotations. 
                    
    """
    
    model = XGram.Model.Model()
    grammar = XGram.Model.Grammar()
    alphabet = XGram.Model.Alphabet()
            
    if substitution_model == "jc69":
        generator_initial_states = InitialStatesConst
        generator_transitions = TransitionsConst
    elif substitution_model == "k80":
        generator_initial_states = InitialStatesConst
        generator_transitions = TransitionsK80
    elif substitution_model == "rev":
        generator_initial_states = InitialStatesConst
        generator_transitions = TransitionsConst                   
    elif substitution_model == "gtr":
        generator_initial_states = InitialStatesParametric
        generator_transitions = TransitionsGTR              
    else:
        raise Exceptions.UsageError( "model %s unknown" % substitution_model )
        

    generator_grammar = GrammarLinearSequence( explicit_extension = explicit_extension )      

    if grammar_type == "linear":
        chain = XGram.Model.Chain()     
        generator_chain = ChainDNA( generator_initial_states(),
                                    generator_transitions() )

        
        generator_chain.buildGrammar( chain )

        if substitution_model == "rev":
            chain.setUpdatePolicy( "rev" )

        grammar.addChain( chain )
        generator_grammar = GrammarLinearSequence( explicit_extension = explicit_extension )  
        
    elif grammar_type in ("linear-blocks", "multiple-blocks" ) :
        for x in range( num_blocks ):
            chain = XGram.Model.Chain()    
            prefix = "B%x" % x
            generator_chain = ChainDNA( generator_initial_states( shared_rates = shared_rates,
                                                                  shared_frequencies = shared_frequencies,
                                                                  prefix = prefix ),
                                        generator_transitions( shared_rates = shared_rates,
                                                               shared_frequencies = shared_frequencies,
                                                               prefix=prefix ),
                                        prefix = prefix )
            if substitution_model == "rev":
                chain.setUpdatePolicy( "rev" )

            generator_chain.buildGrammar( chain )
            grammar.addChain( chain )
        if grammar_type == "multiple-blocks":
            generator_grammar = GrammarLinearSequenceMultipleChains( annotate_terminals = annotate_terminals )            
        elif grammar_type == "linear-blocks":
            generator_grammar = GrammarLinearSequenceBlocks( annotate_terminals = annotate_terminals )
    else:
        raise XGram.Exceptions.UsageError( 'unknown grammar %s' % grammar_type )
        
    generator_grammar.buildGrammar( grammar )
    grammar.setName( substitution_model )
    grammar.addComment( "grammar built by DNA.py " )
    grammar.addComment( "  shared rates = %s" % shared_rates)
    grammar.addComment( "  shared frequencies = %s" % shared_frequencies)
    grammar.addComment( "  type = %s" % grammar_type)
    grammar.addComment( "  nblocks = %i" % num_blocks)
    grammar.addComment( "  substitution model = %s" % substitution_model)
            
    generator_alphabet = AlphabetDNA()
    generator_alphabet.buildGrammar(alphabet)
     
    model.setAlphabet( alphabet )
    model.setGrammar( grammar )
    
    return model
