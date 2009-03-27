################################################################################
#   XGram parser
#
#   $Id: Lexer.py,v 1.9 2007/10/15 09:20:04 grepall Exp $
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

import Exceptions
import ply.lex

# List of token names.   This is always required
tokens = (
   'BO',
   'BC',
   'GRAMMAR',
   'NAME',
   'FROM',
   'TO',
   'TRANSFORM',
   'UPDATERULES',
   'UPDATERATES',
   'PARAMETRIC',
   'GAPSOK',
   'MINLEN',
   'RATE',
   'PROB',
   'PARAMS',
   'CONST',
   'FLOAT',
   'ATOM',
   'MUTATE',
   'INITIAL',
   'CHAIN',
   'ROW',
   'TERMINAL',
   'UPDATEPOLICY',
   'HIDDENCLASS',
   'ALPHABET',
   'EXTEND',
   'TOKEN',
   'COMPLEMENT',
   'WILDCARD',
   'ANNOTATE',
   'LABEL',
   'PLUS',
   'MINUS',
   'TIMES',
   'DIVIDE',
   'STATE',
   'NULLNONTERMINAL',
   'PGROUP',
   'CONSTPGROUP',
   'CONSTRGROUP',
   'PSEUDOCOUNTS', 
   'OBSERVEDCOUNTS',
   'OBSERVEDCHAINCOUNTS', 
   'HYBRIDCHAIN',
   'COMPONENTS',
   'NONTERMINAL',
   'WAIT',
   'TIME',
   'COUNT',
   'POUND',
)

## List of reserved words. Those not matching are atoms
reserved = {
    'grammar' : 'GRAMMAR',
    'name' : 'NAME',
    'from' : 'FROM',
    'to': 'TO',
    'transform': 'TRANSFORM',
    'update-rules': 'UPDATERULES',
    'update-rates': 'UPDATERATES',
    'parametric' : 'PARAMETRIC',
    'gaps-ok': 'GAPSOK',
    'minlen': 'MINLEN',
    'rate' : 'RATE',
    'prob' : 'PROB',
    'params' : 'PARAMS',
    'const': 'CONST',
    'pgroup' : 'PGROUP',
    'const-pgroup': 'CONSTPGROUP',
    'const-rate' : 'CONSTRGROUP',  
    'mutate': 'MUTATE',
    'initial': 'INITIAL',
    'chain': 'CHAIN',
    'terminal': 'TERMINAL',
    'update-policy': 'UPDATEPOLICY',
    'hidden-class': 'HIDDENCLASS',
    'label' : 'LABEL',
    'row' : 'ROW',
    'alphabet' : 'ALPHABET',
    'extend' : 'EXTEND',
    'token' : 'TOKEN',
    'complement' : 'COMPLEMENT',
    'wildcard' : 'WILDCARD',
    'annotate' : 'ANNOTATE',
    'pseudocounts' : 'PSEUDOCOUNTS', 
    'observed-counts' : 'OBSERVEDCOUNTS', 
    'observed-chain-counts' : 'OBSERVEDCHAINCOUNTS', 
    'hybrid-chain' : 'HYBRIDCHAIN', 
    'components' : 'COMPONENTS',
    'nonterminal' : 'NONTERMINAL',
    'count' : 'COUNT',
    'wait' : 'WAIT',
    'time' : 'TIME', 
    '+' : 'PLUS',
    '-' : 'MINUS',
    '*' : 'TIMES',
    '/' : 'DIVIDE',
    '#' : 'POUND',
    }

# Regular expression rules for simple tokens
t_PLUS    = r'\+'
t_MINUS   = r'-'
t_TIMES   = r'\*'
t_DIVIDE  = r'/'

t_POUND = r'\#'
t_STATE=r'[a-zA-Z][a-zA-Z0-9]*'
t_NULLNONTERMINAL=r'[a-zA-Z][a-zA-Z0-9]*[*]'
t_FLOAT = r'[0-9.]+[0-9.+-Ee]*'
t_BO = r'[(]'
t_BC = r'[)]'
t_ignore  = ' \t\n'
t_NAME = 'name'
t_PARAMETRIC = 'parametric'
t_UPDATERATES = 'update-rates'
t_UPDATERULES = 'update-rules'
t_FROM = 'from'
t_TO = 'to'
t_TRANSFORM = 'transform'
t_GAPSOK = 'gaps-ok'
t_MINLEN = 'minlen'
t_RATE = 'rate'
t_PROB = 'prob'
t_PARAMS = 'params'
t_CONST = 'const'
t_GRAMMAR = 'grammar'
t_MUTATE = 'mutate'
t_INITIAL = 'initial'
t_CHAIN = 'chain'
t_TERMINAL = 'terminal'
t_UPDATEPOLICY = 'update-policy'
t_HIDDENCLASS = 'hidden-class'
t_LABEL = 'label'
t_ROW = 'row'
t_ALPHABET = 'alphabet'
t_EXTEND = 'extend'
t_TOKEN = 'token'
t_COMPLEMENT = 'complement'
t_WILDCARD = 'wildcard'
t_ANNOTATE = 'annotate'
t_HYBRIDCHAIN = 'hybrid-chain'
t_PSEUDOCOUNTS = 'pseudocounts'
t_OBSERVEDCOUNTS = 'observed-counts'
t_OBSERVEDCHAINCOUNTS = 'observed-chain-counts'
t_TIME = 'time'
t_WAIT = 'wait'
t_COMPONENTS = 'components'
t_NONTERMINAL = 'nonterminal'
t_COUNT = 'count'

## Anything that is not a reserved word and not a number
def t_ATOM(t):
    r'[a-zA-Z_/*<>0-9?~][a-zA-Z_0-9-/*<>\-.]*'
    try:
        float(t.value)
        t.type = ('FLOAT')
    except ValueError:
        t.type = reserved.get(t.value,'ATOM')    # Check for reserved words
        
    ## translate ~ to not_
    if t.value[0] == "~":
        t.value = "not_%s" % t.value[1:]
    
    return t

## Error handling rule
def t_error(t):
    raise Exceptions.ParsingError("Illegal character '%s'" % (t.value[0]), str(t.lineno))

def t_newline(t):
    r'\n+'
    t.lineno += len(t.value)

## Build the lexer
ply.lex.lex( )


