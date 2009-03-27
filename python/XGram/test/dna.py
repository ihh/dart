import sys, time, string, re, optparse

USAGE="""test module for calculating nucleotide substitution rates.
"""

from XGram.Generator.Prebuilt import DNA
from XGram.Model import Annotation
import XGram.Run

if __name__ ==  "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: dna.py,v 1.1 2007/09/19 12:49:45 grepall Exp $", usage = USAGE)

    parser.add_option("-v", "--verbose", dest="loglevel", type="int", 
                      help="loglevel to use."  )

    parser.add_option("-m", "--model", dest="models", type="choice", action="append",
                      choices=("f3x4-two", "f3x4-four"),
                      help="models to use."  )

    parser.add_option("-w", "--write", dest="write", type="choice", action="append",
                      choices=("input_fixed", "trained_fixed", "input_variable", "trained_variable", "all" ),
                      help="output sections to write."  )

    parser.add_option("-o", "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern for output files." )
    
    parser.add_option("-t", "--test", dest="test", type="choice",
                      choices=("models", "blocks" ),
                      help="perform test on models or block/non-block models." )

    parser.add_option("-i", "--insert-frequencies", dest="insert_frequencies", action="store_true",
                      help="insert frequencies from input file." )


    parser.set_defaults(
                        loglevel = 1,
                        model = "jc69",
                        test = True,
                        write = [],
                        output_pattern = "%s.eg",
                        stdout = sys.stdout,
                        stdlog = sys.stdout,
                        value_format = "%6.4f",
                        )

    (options, args) = parser.parse_args()
    
    xgram = XGram.XGram()

    model = DNA.buildModel( substitution_model = options.model )

    # print model.getGrammar()

    if len(args) > 0:
        data = args[0]
    else:
        data = XGram.PATH_DATA + "/dpse_dmel.stk"
    
    if options.test:
        
        xgram.setDebug()
        data = XGram.PATH_DATA + "/dpse_dmel.stk"    
        
        # print trained_model.getGrammar()
        print "result according to %s" % options.model
        
        if options.model in ( "k80", ):

            result = xgram.train( model, data )
            trained_model = result.getModel()
            
        elif options.model in ('jc69', ):
            result = xgram.buildTree( model, data )
        
        if options.model == "k80":
            alpha, beta = \
                ( trained_model.mGrammar.getParameter( 'alpha' ),
                  trained_model.mGrammar.getParameter( 'beta' ) )
            # this assumes that the branch length in the input is normalized to 1
            # which is true for dpse_dmel.stk
            distance = 2 * beta + alpha
            kappa = alpha / beta
            
            print "\t".join( ("alpha", "beta", "kappa", "distance", "logL" ) )
            print "\t".join( map(lambda x: options.value_format % x, (alpha, beta, kappa, distance, result.getLogLikelihood()) ) )
            
        elif options.model == "jc69":
            
            tree = result.getTree()
            ## multiply distance by tree, as rates are set to 1 and 
            ## thus the matrix is scaled by a factor of 3
            distance = 3.0 * float(re.search("\(\S+:([0-9.]+)\)", tree).groups()[0])
            print "distance"
            print options.value_format % distance

    