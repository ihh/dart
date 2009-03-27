"""Maps an XGram model definition into an XGram model object.

The XGram model definition includes various components like
an alphabets, chains, grammars, etc. The classes in this
module map these definitions onto objects.

An XGram model definition is mapped onto the Model class.
The Model class contains an alphabet (Alphabet) and a
grammar (Grammar). The grammar contains one or more
of the following:

    * Rule: a transformation rule in the stochastic context
        free grammar. Rules might contain annotations (Annotation)
    
    * Chain: a markov chain describing the mutational process    
    
Both Chain and Grammar are derived from ProbabilisticModel, which
manages variables and constants defined in a model.

Each component has a getGrammar() method that returns the model
definition in XGram's LISP-S format.

Note: not everything is modelled. For example rates, expressions
and other concepts are moslty past around as values or tuples.
"""

import sys, string, re
from types import *

import XGram.Exceptions

def PrintStructure( ss ):
    """print a list/tuple of lists/tuples as a LISP type hierarchy of brackets.
    The last level of list/tuple is simply concatenated."""
    entries = []
    for s in ss:
        if type(s) in (ListType, TupleType):
            entries.append( PrintStructure( s ) )
        else:           
            entries.append( str(s) )
    return "(" + ' '.join(entries) + ")"

def PrintMembers( for_class ):
    """print members of a class."""
    members = for_class.__dict__
    member_keys = list(members.keys())
    member_keys.sort()
    rows = []
    for member in member_keys:
        if member[0] == 'm':
            v = members[member]            
            if type(v) == DictType:
                rows.append( "%-40s: Dictionary" % member )                
                for k, vv in v.items():
                    s = str(vv)
                    if s.split("\n") > 1:
                        rows.append( "%s:" % str(k))
                        rows.append( "%s" % s )
                    else:
                        rows.append( "%s -> %s" % (str(k), s))
            else:
                rows.append( "%-40s: %s" % (member, str(v)))
    return "\n".join(rows)

def BuildExpression( expression ):
    """returns a string for expression.
    
    The input is a tuple of (kind,expression),
    where kind is the type of expression (rate|prob)
    and expresssion is either a string or a tuple/list.
    """
    kind, expr = expression
    if type(expr) in (TupleType, ListType):
        expr = " ".join(map(str, expr))
    return "%s %s" % (kind, expr)
    
from ModelAlphabet import *
from ModelAnnotation import *
from ModelProbabilisticModel import *
from ModelChain import *
from ModelGrammar import *
from ModelModel import *
from ModelRule import *
from ModelParameter import *


        