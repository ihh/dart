################################################################################
#   XGram parser
#
#   $Id: Parser.py,v 1.14 2008/01/30 14:07:24 grepall Exp $
#
#   Copyright (C) 2006 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
import os, sys

from types import *

import Model
import Exceptions

# Import parser engine
import ply.yacc
# Get the token map from the lexer.  This is required.
from Lexer import tokens

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
##
## The parser works by filling a dictionary with values according to the
## different components for grammars, alphabets and chains. The dictionary
## is translated into an object manually in p_grammar, p_chain and p_alphabet.
##
def JoinComponents( p ):
    """This function is called by all xxx_components functions to update
    the dictionary.

    If a field does not exist, it is created.
    If a field does exist, it is converted to a List (if necessary) and
    the new values are appended.
    """
    
    p[0] = {}    
    if len(p) > 2:
        if p[1]:
            p[0] = p[1]

        ## aggregate components
        ## if more than one component of the same type is present,
        ## the results will be appended.
        if p[3]:
            for k, v in p[3].items():
                if k not in p[0]:
                    p[0][k] = v
                else:
                    if type(p[0][k]) == DictType:
                        for kk, vv in v.items():
                            p[0][k][kk] = vv
                    else:
                        if type(p[0][k]) != ListType:
                            p[0][k] = [p[0][k],]
                        if v == ListType:
                            p[0][k] += v
                        else:
                            p[0][k].append( v )
    

###########################################################################
## Start of rule definitions
###########################################################################
def p_xgram(p):
    """xgram :  xgram_alphabet_first
    | xgram_grammar_first
    """
    # A grammar consists of both alphabet and grammar. 
    p[0] = p[1]
    
def p_xgram_alphabet_first(p):
    """xgram_alphabet_first : alphabet grammar
    """
    p[0] = Model.Model()
    p[0].mAlphabet = p[1]
    p[0].mGrammar = p[2]

def p_xgram_grammar_first(p):
    """xgram_grammar_first : grammar alphabet
    """
    p[0] = Model.Model()
    p[0].mGrammar = p[1]
    p[0].mAlphabet = p[2]

###########################################################################
###########################################################################
###########################################################################
## Start of grammar section
###########################################################################    
def p_grammar(p):
    """grammar : BO GRAMMAR grammar_components BC
    """
    g = Model.Grammar()

    # build the grammar from the properties defined
    # down the tree

    for k,v in p[3].items():
        if k == 'name':
            g.setName(v)
        elif k == 'chain':
            if type(v) == ListType:
                for gg in v:
                    g.addChain( gg )
            else:
                g.addChain( v )
        elif k == 'params':
            for vv in v:
                if len(vv) == 1:
                    g.addVariable( vv[0] )
                else:
                    g.addVariable( vv )
        elif k == 'rgroup':
            for vv in v:
                g.addRate( vv, is_const = False )
        elif k == 'const_rgroup':
            for vv in v:
                g.addRate( vv, is_const = True )
        elif k == 'pgroup':
            for vv in v:
                g.addProbability( vv, is_const = False )
        elif k == 'const_pgroup':
            for vv in v:
                g.addProbability( vv, is_const = True )
        elif k == 'const':      
            for vv in v:
                if len(vv) == 1:
                    g.addConst( vv[0] )
                else:
                    g.addConst( vv )
        elif k == 'rule':
            for vv in v:
                g.addRule( vv )     
        elif k == 'nonterminal':
            for vv in v:
                g.addNonTerminal( vv )
        elif k == 'update-rules':
            g.setUpdateRules( v )
        elif k == 'update-rates':
            g.setUpdateRates( v )
        elif k == 'parametric':
            g.setIsParametric( v )
        elif k == 'pseudo-counts':
            g.setPseudoCounts( v )
        elif k == 'observed-counts':
            g.setObservedCounts( v )
        elif k == 'observed-chain':
            # TODO: save observed counts here
            pass
        else: 
            raise "unknown component %s" % k

    # add Parameters to chain as declared in Grammar
    for t, c in g.getChains().items():
        c.addParameters(g)
        
    p[0] = g

def p_grammar_components(p):
    """
    grammar_components : grammar_components BO grammar_component BC
    | empty
    """
    JoinComponents( p )
    
def p_grammar_component(p):
    """grammar_component : name
    | update_rules
    | update_rates
    | nonterminal
    | parametric
    | rule
    | parameters
    | consts
    | pgroup
    | const_pgroup
    | rgroup
    | const_rgroup
    | chain
    | hybridchain
    | pseudocounts
    | observedcounts
    | observedchaincounts
    """
    p[0] = p[1]

def p_update_rules(p):
    """update_rules : UPDATERULES FLOAT"""
    p[0] = { 'update-rules' : p[2] != "0" }
        
def p_update_rates(p):
    """update_rates : UPDATERATES FLOAT"""
    p[0] = { 'update-rates' : p[2] != "0" }

def p_nonterminal(p):
    """nonterminal : NONTERMINAL BO NAME ATOM BC"""
    p[0] = { 'nonterminal' : p[4] }
    
def p_parametric(p):
    """parametric : PARAMETRIC"""
    p[0] = { p[1] : True }

def p_rule(p):
    """rule : TRANSFORM rule_components"""

    rule = Model.Rule( p[2]['start'], p[2]['end'] )

    for k, v in p[2].items():
        if k == 'gapsok':
            rule.setGapsOk( True )
        elif k == 'minlen':
            rule.setMinLen( v )
        elif k == 'annotation':
            if type(v) == ListType:
                for vv in v:
                    rule.addAnnotation( vv )
            else:
                rule.addAnnotation( v )
        elif k == 'rate':
            rule.setRate( v )
    
    p[0] = { 'rule' : rule}    
    
def p_rule_components(p):
    """rule_components : rule_components BO rule_component BC 
    | empty
    """
    JoinComponents(p)
    
def p_rule_component(p):
    """rule_component : rule_from
    | rule_to
    | rule_rate
    | rule_annotation
    | rule_gapsok
    | rule_minlen
    """
    p[0] = p[1]
    
def p_rule_rate(p):
    """rule_rate : probability
    | rate
    """
    p[0] = {'rate' : p[1] }

def p_rule_from(p):
    """rule_from : FROM BO states BC"""
    p[0] = {'start' : p[3] }

def p_rule_to(p):
    """rule_to : TO BO states BC"""
    p[0] = {'end' : p[3] }

def p_rule_minlen(p):
    """rule_minlen : MINLEN FLOAT"""
    p[0] = {'minlen': int(p[2]) }
    
def p_rule_gapsok(p):
    """rule_gapsok : GAPSOK"""
    p[0] = {'gapsok' : True }
    
def p_rule_annotation(p):
    """rule_annotation : ANNOTATE annotation_components"""
    
    row, column, label = None, None, None
    for k,v in p[2].items():
        if k == 'row':
            row = v
        elif k == 'column':
            column = v
        elif k == 'label':
            label = v

    annotation = Model.Annotation( row = row, 
                                         column = column,
                                         label = label )
    
    p[0] = {'annotation' : annotation }

def p_annotation_components(p):
    """annotation_components : annotation_components BO annotation_component BC 
    | empty
    """
    JoinComponents(p)
    
def p_annotation_component(p):
    """annotation_component : ATOM ATOM
    | ROW ATOM
    | LABEL ATOM
    """
    p[0] = { p[1] : p[2] }

def p_consts(p):
    """consts : CONST parameter_list"""
    p[0] = { 'const' : p[2] }

def p_parameters(p):
    """parameters : PARAMS parameter_list
    """
    p[0] = { 'params': p[2] }

def p_pgroup(p):
    """pgroup : PGROUP parameter_list
    """
    p[0] = { 'pgroup': p[2]  }
    
def p_const_pgroup(p):
    """const_pgroup : CONSTPGROUP parameter_list
    """
    p[0] = {'const_pgroup' : p[2] }
    
def p_rgroup(p):
    """rgroup : RATE parameter_list
    """
    p[0] = { 'rgroup': p[2] }

def p_const_rgroup(p):
    """const_rgroup : CONSTRGROUP parameter_list
    """
    p[0] = { 'const_rgroup': p[2] }
    
def p_parameter_list(p):
    """parameter_list : BO parameter_group BC parameter_list
    | BO parameter BC parameter_list
    | empty
    """
    if len(p) == 5:
        p[0] = [p[2],] + p[4]        
    else:
        p[0] = []
        
def p_parameter_group(p):
    """parameter_group : BO parameter BC parameter_group
    | empty
    """
    if len(p) > 2:
        p[0] = [p[2],] + p[4]
    else:
        p[0] = []

def p_parameter(p):
    """parameter : ATOM FLOAT
    """
    p[0] = (p[1], float(p[2]))
        
###########################################################################
###########################################################################
###########################################################################
## Rules for getting counts from results file
## Rate parameters have waiting times associated with them, thus these blocks
## are not 100% equivalent to the paramter blocks.
###########################################################################
def p_pseudocounts(p):
    """pseudocounts : PSEUDOCOUNTS counts_list
    """
    p[0] = { 'pseudo-counts': p[2] }

def p_observedcounts(p):
    """observedcounts : OBSERVEDCOUNTS counts_list
    """
    p[0] = { 'observed-counts' : p[2] }
    
def p_counts_list(p):
    """counts_list : counts_list BO counts_group BC 
    | empty
    """
    JoinComponents( p )
    
def p_counts_group(p):
    """counts_group : ATOM FLOAT
    | ATOM FLOAT FLOAT
    | empty
    """
    if len(p) == 3:
        p[0] = { p[1] : { 'counts': float(p[2]) } }
    elif len(p) == 4:
        p[0] = { p[1] : { 'counts': float(p[2]), 'time': float(p[3]) } }        
    else:
        p[0] = {}

###########################################################################
###########################################################################
###########################################################################
## Rules for observed chain counts.
###########################################################################
def p_observedchaincounts(p):
    """observedchaincounts : OBSERVEDCHAINCOUNTS observed_chains
    """
    p[0] = {'observed-chain' : p[2] }
    

def p_observed_chains(p):
    """observed_chains : observed_chains BO observed_chain BC
    | empty """
    JoinComponents( p )
    
def p_observed_chain( p ):
    """observed_chain : observed_chain_components"""
    
    c = Model.Chain()

    for k,v in p[1].items():
        if k == 'terminal':
            c.setTerminals( v )
        elif k == 'transition':
            for start, end, rate in v:
                c.addTransition( start, end, rate )
        elif k == 'initial_state':
            for state, rate in v:
                c.addInitialState( state, rate )
        elif k == 'wait_state':
            pass
#            for state, rate in v:
#                c.addInitialState( state, rate )
            
    p[0] = { c.getTerminals(): c}


def p_observed_chain_components(p):
    """
    observed_chain_components : observed_chain_components BO observed_chain_component BC
    | empty
    """
    JoinComponents( p )
    
def p_observed_chain_component(p):
    """observed_chain_component : terminal
    | initial_state
    | transition
    | wait_state
    """
    p[0] = p[1]

    
###########################################################################
###########################################################################
###########################################################################
## Rules for the chain
###########################################################################
def p_chain(p):
    """chain : CHAIN chain_components"""
    c = Model.Chain()

    for k,v in p[2].items():
        if k == 'terminal':
            c.setTerminals( v )
        elif k == 'transition':
            for start, end, rate in v:
                c.addTransition( start, end, rate )
        elif k == 'initial_state':
            for state, rate in v:
                c.addInitialState( state, rate )
        elif k == 'update-policy':
            c.setUpdatePolicy( v )
            
    p[0] = {'chain': c }

def p_chain_components(p):
    """
    chain_components : chain_components BO chain_component BC
    | empty
    """
    JoinComponents( p )
    
def p_chain_component(p):
    """chain_component : terminal
    | update_policy
    | initial_state
    | transition
    | hidden_class
    """
    p[0] = p[1]
    
def p_terminal(p):
    """terminal : TERMINAL BO states BC"""
    p[0] = {'terminal' : p[3] }

def p_update_policy(p):
    """update_policy : UPDATEPOLICY ATOM
     | UPDATEPOLICY PARAMETRIC"""
    p[0] = {'update-policy' : p[2] }

def p_initial_state(p):
    """initial_state : INITIAL BO ATOM BO states BC BC BO probability BC
     | INITIAL BO ATOM BO states BC BC BO counts BC"""
    p[0] = { 'initial_state' : (p[5], p[9] ) }

def p_wait_state(p):
    """wait_state : WAIT BO ATOM BO states BC BC BO time BC"""
    p[0] = { 'wait_state' : (p[5], p[9] ) }
    
def p_transition(p):
    """transition : MUTATE BO FROM BO states BC BC BO TO BO states BC BC BO rate BC
     | MUTATE BO FROM BO states BC BC BO TO BO states BC BC BO counts BC
    """
    p[0] = { 'transition': (p[5], p[11], p[15]) }

def p_hidden_class(p):
    """hidden_class : HIDDENCLASS hidden_class_options_list"""
    pass

def p_hidden_class_options_list(p):
    """hidden_class_options_list : hidden_class_options_list BO hidden_class_options BC
    | empty
    """
    pass

def p_hidden_class_options( p ):
    """hidden_class_options : ROW ATOM
    | LABEL BO states BC
    """
    pass

###########################################################################
###########################################################################
###########################################################################
## Rules for hybrid chains
###########################################################################
def p_hybridchain(p):
    """hybridchain : HYBRIDCHAIN hybridchain_parts
    """
    return { 'hybrid-chain' : None }

def p_hybridchain_parts(p):
    """
    hybridchain_parts : hybridchain_parts BO hybridchain_part BC
    | empty
    """
    JoinComponents( p )
    
def p_hybridchain_part(p):
    """hybridchain_part : terminal
    | row
    | components
    """
    p[0] = p[1]

def p_row( p ):
    """row : ROW ATOM
    """
    return { 'row' : p[2] }

def p_components( p ):
    """components : COMPONENTS components_list"""
    pass

def p_components_list( p ):
    """components_list : components_list BO component BC
    | empty""" 
    pass

def p_component( p ):
    """component : BO LABEL ATOM BC BO TERMINAL BO ATOM BC BC"""
    pass

###########################################################################
###########################################################################
###########################################################################
## Rules for the alphabet
###########################################################################
def p_alphabet( p ):
    """alphabet : BO ALPHABET alphabet_components BC"""

    g = Model.Alphabet()

    # build the grammar from the properties defined
    # down the tree

    for k,v in p[3].items():
        if k == 'name':
            g.setName( v )
        elif k == 'tokens':
            g.setTokens( v )
        elif k == 'complements':
            g.setComplements( v )
        elif k == 'extend':     
            for kk, vv in v.items():
                g.addExtension( kk, vv)       
        elif k == 'wildcard':
            g.setWildcard( v )
            
    p[0] = g

def p_alphabet_components(p):
    """
    alphabet_components : alphabet_components BO alphabet_component BC
    | empty
    """
    JoinComponents( p )

def p_alphabet_component(p):
    """alphabet_component : name
    | tokens
    | complements
    | extension
    | wildcard
    """
    p[0] = p[1]

def p_tokens(p):
    """tokens : TOKEN BO token_list BC"""
    p[0] = {'tokens' : p[3] }

def p_complements(p):
    """complements : COMPLEMENT BO token_list BC"""
    p[0] = { 'complements' : p[3] }

def p_token_list(p):
    """token_list : ATOM token_list
    | empty
    """
    ## Token list are tuples so that they are hashable
    if len(p) == 3:
        p[0] = (p[1],) + p[2]
    else:
        p[0] = ()
        
def p_extension(p):
    """extension : EXTEND BO TO ATOM BC from_list"""
    p[0] = {'extend' : { p[4] :  p[6] } }
    
def p_from_list(p):
    """from_list : BO FROM ATOM BC from_list
    | empty
    """
    if len(p) == 6:
        p[0] = (p[3],) + p[5]
    else:
        p[0] = ()
    
def p_wildcard(p):
    """wildcard : WILDCARD ATOM
                | WILDCARD TIMES"""
    p[0] = { 'wildcard': p[2] }

###########################################################################
###########################################################################
###########################################################################
## General rules for shared productions
###########################################################################
def p_name(p):
    """name : NAME ATOM
    """
    p[0] = { 'name' : p[2] }

def p_probability(p):
    """probability : PROB expression
    """
    p[0] = (p[1], p[2])
    
def p_rate(p):
    """rate : RATE expression
    """
    p[0] = (p[1], p[2])

def p_counts(p):
    """counts : COUNT FLOAT
    """
    p[0] = (p[1], p[2])

def p_time(p):
    """time : TIME expression
    """
    p[0] = (p[1], p[2])

## do not parse expressions, just save the tokens
def p_expression(p):
    """
    expression : expression PLUS expression
           | expression MINUS expression   
           | expression TIMES expression   
           | expression DIVIDE expression  
           | BO expression BC               
           | FLOAT                          
           | ATOM
           | POUND ATOM
    """
    if len(p) == 4:
        if p[1] == "(":
            p[0] = [p[1]] + p[2] + [p[3]]
        else:
            p[0] = p[1] + [p[2]] + p[3]
    elif len(p) == 2:
        p[0] = [p[1]]
    elif len(p) == 3:
        p[0] = [p[1] + " " + p[2]]

def p_states(p):
    """states : STATE states
    | ATOM states
    | NULLNONTERMINAL states
    | empty
    """
    # States are tuples of characters. This way they can be hashed.
    if len(p) > 2:
        p[0] = (p[1],) + p[2]        
    else:
        p[0] = ()
    
###########################################################################
## Empty rule
def p_empty(p):
    'empty :'
    pass

###########################################################################
# Error rule for syntax errors
def p_error(p):

    raise Exceptions.ParsingError( "parsing error at %s; state=%s, next token=%s" %\
                                   (p.value, str(p), ply.yacc.token()), str(p.lineno))
    
# Build the parser
ply.yacc.yacc()

###########################################################################
###########################################################################
## End of parser definition
###########################################################################

###########################################################################
###########################################################################
# Use this if you want to build the parser using LALR(1) instead of SLR
# yacc.yacc(method="LALR")

def parseGrammar( lines ):
    """parse a grammar given in a list of lines."""
    
    data = []
    # remove comments
    for line in lines:
        try:
            x = line.index( ";" )
            line = line[:x]
        except ValueError:
            pass
        data.append( line )
    data = "".join(data)        
    return ply.yacc.parse(data)
