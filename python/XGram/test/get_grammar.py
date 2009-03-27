USAGE="""generic module to retrieve pre-built grammars.
"""

import XGram.Generator

from XGram.Generator.Prebuilt import RNA, ProteinSites, Codons

import sys, string, os, re, optparse

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: get_grammar.py,v 1.2 2006/10/04 16:39:39 grepall Exp $", usage = USAGE)

    parser.add_option("-m", "--model", dest="model", type="string",
                      help="model to build."  )

    parser.set_defaults(
                        model = None,
                        )

    (options, args) = parser.parse_args()
    
    if not options.model:
        raise "please supply a model."

    options.model ="RNA.buildPFold( strand_symmetric=True, weak_strong = True, copy_parameters='/home/aheger/workspace_dart/XGram/XGram/data/pfold.eg')"
    # options.model ="Codons.buildCodonML( codon_model='f3x4-four', grammar_type = 'linear-blocks' )"
    grammar = eval(options.model)
    
    print grammar.getGrammar()
    