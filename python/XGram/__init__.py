"""Python interface to XGram.

The module provides the following packages/modules:

Model: Maps an XGram grammar onto a XGram model containing
        chains, alphabets, etc.

Parser: A module to parse XGram grammars and
        return a model. Uses the Lexer class.
 
Generator: A package to build XGram models from scratch using
        a set of predefined components.
   
Run: An interface to the XGram executable with methods to
        train models and annotate alignments.
        
Exceptions: A set of exceptions thrown in this module.
"""
## path towards to the data files

import os
PATH_DATA=os.path.abspath( __file__ + "/../data" )

import Generator
import Model
import Parser
import Lexer
from Exceptions import *
from XGram import XGram     
   


