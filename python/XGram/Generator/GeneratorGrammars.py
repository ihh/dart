## Generator for grammars

import Bio
import Bio.Alphabet
import Bio.Alphabet.IUPAC
import Bio.Data.CodonTable

from types import * 
import XGram.Model

class Grammar:
    
    def __init__(self, 
                 pseudo_terminals = None, 
                 annotate_terminals = None,
                 explicit_extension = False,
                 *args, **kwargs):
        """if pseudo_terminals is not set, the pseudo_terminals are obtained 
        from the associated chains.
        
        annotate_terminals is a map of terminals to annotations.
        
        explicit_extension: use explicit parameter "ext" for extension probability
        """
        
        self.mPseudoTerminals = pseudo_terminals
        self.mAnnotateTerminals = annotate_terminals
        self.mExplicitExtension = explicit_extension
        
    def buildGrammar(self, grammar):
        """stub for building grammars.
        
        The base class takes care of update-rules/update-rates.
        """
        grammar.setName( "linear_sequence" )
        grammar.setUpdateRules( True )
        
        if self.mExplicitExtension:
            grammar.setIsParametric( True )        
        
        for chain in grammar.getChains().values():
            if chain.isMutable():
                grammar.setUpdateRates( True )
        
class GrammarLinearSequence(Grammar):
    """generator for grammars of linear sequences with gaps."""
    
    def __init__(self, no_gaps = False, *args, **kwargs):
        Grammar.__init__(self, *args, **kwargs)
        self.mNoGaps = no_gaps
    
    def buildGrammar(self, grammar ):
        """build simple grammar for pairwise or multiple alignment."""
        Grammar.buildGrammar(self, grammar)
        
        if self.mPseudoTerminals:    
            pseudo_terminals = " ".join(self.mPseudoTerminals)
        else:
            chains = grammar.getChains()
            if len(chains) != 1:
                raise "more than one chain given."
            else:
                ## retrieve pseudo terminals of first (and only) chain
                pseudo_terminals = " ".join(chains.keys()[0])

        if self.mExplicitExtension:
            rate_ext = ("prob", "ext")
            rate_stop = ("prob", "not_ext" )
            grammar.addProbability( (("ext", 0.5),
                                     ("not_ext", 0.5)), 
                                     is_const = False )
        else:
            rate_ext = ("prob", 0.5)
            rate_stop = ("prob", 0.5)
            
        grammar.addRule( XGram.Model.Rule ( ("START",), 
                                ("EMIT", ),
                                rate = rate_ext ) )
        grammar.addRule( XGram.Model.Rule ( ("START",), 
                                ("",), 
                                rate = rate_stop ) )
        if self.mNoGaps:
            grammar.addRule( XGram.Model.Rule ( ("EMIT",), 
                                    ("%s EMIT*" % pseudo_terminals,) , 
                                    ))
        else:
            grammar.addRule( XGram.Model.Rule ( ("EMIT",), 
                                    ("%s EMIT*" % pseudo_terminals,) , 
                                    minlen = 1,
                                    gapsok = True ))
            
        grammar.addRule( XGram.Model.Rule ( ("EMIT*",), 
                                ("",), 
                                rate = rate_stop ) )
        grammar.addRule( XGram.Model.Rule ( ("EMIT*",), 
                                ("EMIT",), 
                                rate = rate_ext ) )

    
class GrammarLinearSequenceMultipleChains(Grammar):
    """build a grammar with multiple chains which commute freely
    between each other.
    
    There is some problem in xREI when the emission is pre-ceded
    by a transformation rule.
    
    """    
    def __init__(self, *args, **kwargs):
        Grammar.__init__(self, *args, **kwargs)
        
    def buildGrammar(self, grammar):
        
        Grammar.buildGrammar( self, grammar )
        
        if self.mPseudoTerminals:    
            pseudo_terminals = " ".join(self.mPseudoTerminals)
        else:
            chains = grammar.getChains()
            ## retrieve pseudo terminals of first (and only) chain
            pseudo_terminals = chains.keys()
            pseudo_terminals.sort()

        # uniform probabilities
        prob = 1.0 / len(pseudo_terminals)
        for x in pseudo_terminals:
            T = " ".join(x)
            grammar.addRule( XGram.Model.Rule( ("Start", ),
                                   ( "NT_%s" % T, ),
                                   ( "prob", prob )) )
            

        prob = 1.0 / (len(pseudo_terminals) + 1)
        for x in pseudo_terminals:
            T = " ".join(x)
            ## produces a terminal
            if self.mAnnotateTerminals and \
                T in self.mAnnotateTerminals:
                annotations = self.mAnnotateTerminals[T]
            else:
                annotations = []

            grammar.addRule( XGram.Model.Rule( 
                                              start = ("NT_%s" % T, ),
                                              end = ("%s NT_%s*" % (T, T), ),
                                              gapsok = True,
                                              minlen = 1.0,
                                              annotations = annotations ))

            for y in pseudo_terminals:
                TY = " ".join(y)
                grammar.addRule( XGram.Model.Rule( 
                                                  start = ("NT_%s*" % T, ),
                                                  end = ("NT_%s" % TY, ),
                                                  rate = ("prob", prob) ) )

            grammar.addRule( XGram.Model.Rule( 
                                              start = ("NT_%s*" % T, ),
                                              end = (),
                                              rate = ("prob", prob) ) )
                
        grammar.setIsParametric( False )
            
class GrammarLinearSequenceBlocks(Grammar):
    """build a grammar with multiple chains in blocks. 
    
    In this grammar, you can only go forward in the block
    sequence. The grammar and start or terminate from any
    block.
    
    If not pseudo-terminals are given when this generator
    is called, the pseudo-terminals of the chains are taken
    and sorted alphabetically.
    """
    def __init__(self, *args, **kwargs):
        Grammar.__init__(self, *args, **kwargs)

    def addBlock(self, grammar, term1, term2 = None):
        """add emission and transition to end state for block term1.
        
        If term1 is not None, a transition to the next block is added
        (given by term2)
        """

        T1 = "".join( term1 )
        E1 = " ".join( term1 )
        
        ## produces a terminal
        if self.mAnnotateTerminals and \
            term1 in self.mAnnotateTerminals:
            annotations = self.mAnnotateTerminals[term1]
        else:
            annotations = []
            
        grammar.addRule( XGram.Model.Rule( 
                                          start = ("NT_%s" % T1, ),
                                          end = ("%s NT_%s*" % (E1, T1), ),
                                          gapsok = True,
                                          minlen = 1.0,
                                          annotations = annotations ))

        if term2:
            p = 1.0 / 3.0
        else:
            p = 1.0 / 2.0
        ## go to end 
        grammar.addRule( XGram.Model.Rule( start = ("NT_%s*" % T1, ),
                                           end = ("",),
                                           rate = ("prob", p) ) )

        ## or cycle in block
        grammar.addRule( XGram.Model.Rule( start = ("NT_%s*" % T1, ),
                                           end = ("NT_%s" % T1,),
                                           rate = ("prob", p) ) )

        ## or go to next block
        if term2:
            T2 = "".join( term2 )
            grammar.addRule( XGram.Model.Rule( start = ("NT_%s*" % T1, ),
                                               end = ("NT_%s" % T2,),
                                               rate = ("prob", p) ) )

        
    def buildGrammar(self, grammar):
        
        Grammar.buildGrammar( self, grammar )
        
        if self.mPseudoTerminals:    
            pseudo_terminals = " ".join(self.mPseudoTerminals)
        else:
            chains = grammar.getChains()
            ## retrieve pseudo terminals of first (and only) chain
            pseudo_terminals = chains.keys()
            ## make sure that they are block order
            pseudo_terminals.sort()
            
        ## start can go into any block state:
        p = 1.0 / len(pseudo_terminals)
        for x in pseudo_terminals:
            T = "".join(x)
            
            grammar.addRule( XGram.Model.Rule( 
                                              ("Start", ),
                                              ( "NT_%s" % T, ),
                                              ( "prob", p )))
            
            
        ## it is only possible to move to the next block
        ## (or ()) from any block
        for x in range(len(pseudo_terminals) - 1):
            
            self.addBlock( grammar, 
                           pseudo_terminals[x],
                           pseudo_terminals[x+1])

        self.addBlock( grammar, pseudo_terminals[-1] )

        grammar.setIsParametric( False )
