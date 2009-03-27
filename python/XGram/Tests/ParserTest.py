import unittest
import glob
from pprint import pprint

import XGram, XGram.Parser, XGram.Exceptions

class ParserTest(unittest.TestCase):
    """
    A test class for the feedparser module.
    """

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.mInputFiles= glob.glob( XGram.PATH_DATA+"/*.eg" )
        if len(self.mInputFiles) == 0:
            raise "could not find any test grammars in %s" % XGram.PATH_DATA
    
    def tearDown(self):
        """
        tear down any data used in tests
        tearDown is called after each test function execution.
        """
        pass

    def testParser(self):
        """
        test parsing
        """
        print "ParserTest.testParser()"
        nerrors = 0
        for filename in self.mInputFiles:
            print filename, 
            lines= open( filename, "r" ).readlines()
    
            ## test the parser
            try:
                result = XGram.Parser.parseGrammar(lines)
            except XGram.Exceptions.ParsingError, msg:
                nerrors += 1
                print msg
                continue

            freqs = result.evaluateTerminalFrequencies()
            try:
                matrix = result.evaluateRateMatrix()
            except XGram.Exceptions.UsageError, msg:
                print msg
                continue
            
            print "ok"
            # pprint(freqs)
            
        if nerrors > 0:
            self.fail( "Error: %i grammar(s) could be parsed." % nerrors)
                
if __name__ == '__main__':
    unittest.main()