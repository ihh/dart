"""This module provides methods to build models for
protein families. The method here is very simple.
"""

import XGram.Utils 

from XGram import Model, Generator

from XGram.Model import Rule

import Bio.Alphabet.IUPAC

import sys, string, re

def addEmissionMatch( grammar, column, generator_chain,
                      t, nt, pnt):

    ## emission from match state
    grammar.addRule( Rule( start = ( nt, ),
                           end = (t, pnt),
                           rate = ("prob", 1.0 ) ))
    
    ## insert chain for this state
    chain = Model.Chain()
    generator_chain.buildGrammar(chain)
    chain.setTerminals( (t,) )
    chain.setUpdatePolicy( "rev" )        
    grammar.addChain( chain )

def buildProteinFamily( data, 
                        tree_model = None,
                        name = "proteinfamily",
                        local = False):
    """build protein family for alignment in data.
    
    tree_model: model used to estimate tree. If it is not given,
        the tree is assumed to be part of the alignment file.
    
    1. A tree is estimated data using a default model.
    2. Each state is assigned to be match or indel state
        (50% aligned?) -> could use a grammar for this?
        What about the parametric gap model?
    3. equilibrium frequencies are determined via EM
    """
    
    mali = XGram.Utils.Mali()
    
    mali.readFromFile( open(data), format = "stockholm" )
    
    model = Model.Model()
    alphabet = Model.Alphabet()
    grammar = Model.Grammar()
    
    start_state = "start"
    end_state = ""

    previous_state = "col0"
    previous_t = "t_%s" % previous_state    
    previous_nt = "nt_%s" % previous_state
    previous_pnt = "nt_%s*" % previous_state

    ## the chain generator to use
    generator_chain = Generator.Chain( \
         alphabet = Bio.Alphabet.IUPAC.IUPACProtein.letters.lower(),
         )
    
    columns = mali.getColumns()

    ## add transition form start to first state
    grammar.addRule( Rule( start = ( start_state, ),
                           end = ( previous_nt, ),
                           rate = ("prob", 1.0 ) ))

    addEmissionMatch( grammar, columns[0], generator_chain,
                      previous_t, previous_nt, previous_pnt )

    
    for c in range( 1, len(columns) ):
        
        ## add transition from start state for local alignments
        if local:
            grammar.addRule( Rule( start = ( start_state, ),
                                   end = ( current_nt, ),
                                   rate = ("prob", 1.0 ) ))
        
        column = columns[c]
        current_state = "col%i" % c
        current_t = "t_%s" % current_state        
        current_nt = "nt_%s" % current_state
        current_pnt = "nt_%s*" % current_state

        ## transition from previous to this state
        grammar.addRule( Rule( start = ( previous_pnt, ),
                               end = ( current_nt, ),
                               rate = ("prob", 1.0 ) ))


        addEmissionMatch( grammar, column, generator_chain,
                          current_t, current_nt, current_pnt)
        
        ## add transition to end state for local alignments
        if local:
            grammar.addRule( Rule( start = ( current_pnt, ),
                                   end = ( end_state, ),
                                   rate = ("prob", 1.0 ) ))
                         
        previous_t = current_t                        
        previous_nt = current_nt                         
        previous_pnt = current_pnt
        
    ## add final transition , if local alignment
    if not local:
        grammar.addRule( Rule( start = ( current_pnt, ),
                               end = ( end_state, ),
                               rate = ("prob", 1.0 ) ))
        
    grammar.setUpdateRules(True)
    grammar.setUpdateRates(True)
    grammar.setName( name )

    generator_alphabet = Generator.AlphabetProtein()
    generator_alphabet.buildGrammar(alphabet)
     
    model.setAlphabet( alphabet )
    model.setGrammar( grammar )
    
    return model
    
    
    