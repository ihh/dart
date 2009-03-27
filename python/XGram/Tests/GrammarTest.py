import unittest
import glob
from pprint import pprint

import XGram, XGram.Parser, XGram.Exceptions
from XGram.Generator.Prebuilt import Codons

class SimpleGrammarTest(unittest.TestCase):
    """
    A test class for testing a grammar
    """

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.mInputFile= XGram.PATH_DATA+"/kimura2.eg"
        lines= open( self.mInputFile, "r" ).readlines()        
        self.mModel = XGram.Parser.parseGrammar(lines)
    
    def tearDown(self):
        """
        tear down any data used in tests
        tearDown is called after each test function execution.
        """
        pass

    def testSetParameter(self):
        """test setting the value of an existing parameter
        """
        new_value = 2.0
        old_value = 4.0
        grammar = self.mModel.mGrammar
        
        return_value = grammar.setParameter( "alpha", new_value )
        self.assertEqual( grammar.getParameter( "alpha" ), new_value)
        self.assertEqual( old_value, return_value )
        self.mModel.evaluateTerminalFrequencies()
        self.mModel.evaluateRateMatrix()

    def testRenameParameter(self):
        """test renaming a parameter.
        """
        grammar = self.mModel.mGrammar
        old_value = grammar.getParameter( "beta" ) 
        grammar.renameParameter( "beta", "delta" )
        self.assertEqual( old_value, grammar.getParameter("delta") )
        self.mModel.evaluateTerminalFrequencies()
        self.mModel.evaluateRateMatrix()
        
    def testRemoveParameter(self):
        """test removing a parameter.
        """
        old_value = 4.0
        grammar = self.mModel.mGrammar
        
        return_value = grammar.removeParameter( "alpha" )
        self.assertEqual( old_value, return_value )
        # the model is now incomplete, so there will be an error.
        # TODO: remove transitions/states involving a deleted parameter 

class ComplexGrammarTest(unittest.TestCase):
    """
    A test class for testing a grammar
    """

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.mModel = Codons.buildCodonML( "codons-four",
                                           num_blocks=2,
                                           grammar_type="linear-blocks",
                                           shared_frequencies = False,
                                           shared_rates = False )                                     
    
    def tearDown(self):
        """
        tear down any data used in tests
        tearDown is called after each test function execution.
        """
        pass

    def testRenameParameter(self):
        """test renaming a parameter.
        """
        grammar = self.mModel.mGrammar
        old_value = grammar.getParameter( "B0_Rsi" ) 
        grammar.renameParameter( "B0_Rsi", "Rsi" )
        self.assertEqual( old_value, grammar.getParameter("Rsi") )

        old_value = grammar.getParameter( "B1_Rsi" )         
        grammar.renameParameter( "B1_Rsi", "Rsi" )
        self.assertEqual( old_value, grammar.getParameter("Rsi") )
        
        lines = grammar.getGrammar()
        self.mModel.evaluateTerminalFrequencies()
        self.mModel.evaluateRateMatrix()
        
if __name__ == '__main__':
    unittest.main()