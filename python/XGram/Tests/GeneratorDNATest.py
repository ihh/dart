import unittest, string

import XGram, XGram.Parser, XGram.Exceptions
from XGram.Generator.Prebuilt import DNA
from XGram.Model import Annotation

class GeneratorCodonsTest(unittest.TestCase):
    """
    A test class for testing a grammar
    """

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.mInputFile= XGram.PATH_DATA+"/dpse_dmel.stk"

        self.mXgram = XGram.XGram()
        self.mXgram.setDebug()

    def tearDown(self):
        """
        tear down any data used in tests
        tearDown is called after each test function execution.
        """
        pass

    def testModelK80(self):
        """test the k80 model
        """
        
        self.buildAndCheckModel( "k80" )

    def testModelREV(self):
        """test the rev model
        """
        self.buildAndCheckModel( "rev" )

    def testModelGTR(self):
        """test the gtr model
        """
        
        self.buildAndCheckModel( "gtr" )


    def buildAndCheckModel(self, substitution_model):
        """build various models checking parameter settings."""

        print "##### %s : default ##########" % (substitution_model )    
        model = DNA.buildModel( substitution_model = substitution_model )
        self.checkModel( model )
        
        print "##### %s : explicit ##########" % (substitution_model )    
        model = DNA.buildModel( substitution_model = substitution_model,
                                explicit_extension = True )
        self.checkModel( model )       

        num_blocks = 2
        for grammar in ("linear-blocks", "multiple-blocks"):

            print "##### %s : %s : shared rates ##########" % (substitution_model, grammar )    

            model = DNA.buildModel( substitution_model = substitution_model,
                                    grammar_type = grammar,
                                    shared_rates = True,
                                    shared_frequencies = False,
                                    num_blocks = num_blocks )
            self.checkModel( model )       

            print "##### %s : %s : shared freqs ##########" % (substitution_model, grammar )    
            
            model = DNA.buildModel( substitution_model = substitution_model,
                                    grammar_type = grammar,
                                    shared_rates = False,
                                    shared_frequencies = True,
                                    num_blocks = num_blocks )
            self.checkModel( model )       

            print "##### %s : %s : shared all ##########" % (substitution_model, grammar )    
            model = DNA.buildModel( substitution_model = substitution_model,
                                    grammar_type = grammar,
                                    shared_rates = True,
                                    shared_frequencies = True,
                                    num_blocks = num_blocks )
            self.checkModel( model )       

            print "##### %s : %s : shared all with annotations ##########" % (substitution_model, grammar )                
            ## test model with annotations
            ## build annotation
            labels = string.letters.upper()
            annotate_terminals = {}
            for x in range(num_blocks):
                annotations = []
                key = "B%i" % x
                annotations.append( Annotation( row = "STATE",
                                                column = key,
                                                label = labels[x % len(labels)] ))

                annotate_terminals[ key ] = annotations
                
            model = DNA.buildModel( substitution_model = substitution_model,
                                    grammar_type = grammar,
                                    shared_rates = True,
                                    shared_frequencies = True,
                                    num_blocks = num_blocks,
                                    annotate_terminals = annotate_terminals )
            self.checkModel( model )       
                
    def checkModel(self, model ):
        """check a model."""
        
        print "###################################################################"
        print model.getGrammar()
        print "###################################################################"
        # frequencies = model.evaluateTerminalFrequencies()
        # matrix = model.evaluateRateMatrix()
        
if __name__ == '__main__':
    unittest.main()
        

