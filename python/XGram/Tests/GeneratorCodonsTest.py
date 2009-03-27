import unittest
import string
from pprint import pprint

import XGram, XGram.Parser, XGram.Exceptions
from XGram.Generator.Prebuilt import Codons
from XGram.Model import Annotation
import Bio.Data.CodonTable

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

    def testModelF3X4Two(self):
        """test f3x4-two model.
        """
        self.buildAndCheckModel( "f3x4-two" )

    def testModelF3X4Four(self):
        """test f3x4-four model.
        """
        self.buildAndCheckModel( "f3x4-four" )

    def testModelF3X4Four(self):
        """test f3x4-four model.
        """
        self.buildAndCheckModel( "f3x4-fourproducts" )

    def testModelCodonsTwo(self):
        """test codons-two model
        """        
        self.buildAndCheckModel( "codons-two" )
        
        codons = Bio.Data.CodonTable.standard_dna_table.forward_table
        codon_frequencies = {}
        n = 1
        f = 61 * 62 / 2
        for codon in Bio.Data.CodonTable.standard_dna_table.forward_table:
            codon_frequencies[codon] = float(n)/f
            n += 1
        
        self.buildAndCheckModel( "codons-four", codon_frequencies = codon_frequencies )                

    def testModelCodonsFour(self):
        """test codons-four model
        """        
        self.buildAndCheckModel( "codons-four" )
        
        codons = Bio.Data.CodonTable.standard_dna_table.forward_table
        codon_frequencies = {}
        n = 1
        f = 61 * 62 / 2
        for codon in Bio.Data.CodonTable.standard_dna_table.forward_table:
            codon_frequencies[codon] = float(n)/f
            n += 1
        
        self.buildAndCheckModel( "codons-four", codon_frequencies = codon_frequencies )                
        
    def buildAndCheckModel(self, codon_model, **kwargs):
        """build various models checking parameter settings."""
        
        model = Codons.buildCodonML(codon_model = codon_model,
                                    **kwargs )
        self.checkModel( model )
        
        model = Codons.buildCodonML(codon_model = codon_model,
                                    fix_kappa = True,
                                    **kwargs )
        self.checkModel( model )
        
        model = Codons.buildCodonML(codon_model = codon_model,
                                    fix_omega = True,
                                    **kwargs )
        self.checkModel( model )

        model = Codons.buildCodonML(codon_model = codon_model,
                                    fix_omega = True,
                                    fix_kappa = True,
                                    **kwargs )
        self.checkModel( model )

        model = Codons.buildCodonML( codon_model,
                                     num_blocks=2,
                                     grammar_type="linear-blocks",
                                     shared_frequencies = False,
                                     shared_rates = False,                                     
                                     **kwargs )
        
        self.checkModel(model)
        
        num_blocks = 2
        model = Codons.buildCodonML( codon_model,
                                     num_blocks=num_blocks,
                                     grammar_type="linear-blocks",
                                     shared_frequencies = True,
                                     shared_rates = False,
                                     **kwargs)
        self.checkModel(model)

        num_blocks = 2
        model = Codons.buildCodonML( codon_model,
                                     num_blocks=num_blocks,
                                     grammar_type="linear-blocks",
                                     shared_frequencies = False,
                                     shared_rates = True,
                                     **kwargs)
        self.checkModel(model)

        num_blocks = 2
        model = Codons.buildCodonML( codon_model,
                                     num_blocks=num_blocks,
                                     grammar_type="linear-blocks",
                                     shared_frequencies = True,
                                     shared_rates = True,
                                     **kwargs)
        self.checkModel(model)

        ## test model with annotations
        ## build annotation
        labels = string.letters.upper()
        annotate_terminals = {}
        for x in range(num_blocks):
            annotations = []
            key = []

            for c in range( 0,3 ):
                t = "B%i_COD%i" % (x, c)
                key.append(t)
                annotations.append( Annotation( row = "STATE",
                                                column = t,
                                                label = labels[x % len(labels)] ))

            annotate_terminals[ tuple(key) ] = annotations

        model = Codons.buildCodonML( codon_model,
                                     num_blocks=2,
                                     grammar_type="linear-blocks",
                                     shared_frequencies = True,
                                     annotate_terminals = annotate_terminals,
                                     **kwargs )
                                      
        # print model.getGrammar()
        self.checkModel(model)
        
    def checkModel(self, model ):
        """check a model."""
        model.getGrammar()
        frequencies = model.evaluateTerminalFrequencies()
        matrix = model.evaluateRateMatrix()
        
if __name__ == '__main__':
    unittest.main()
        

