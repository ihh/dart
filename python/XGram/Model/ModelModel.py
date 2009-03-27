from types import *
from XGram.Model import Grammar, Alphabet
from XGram import Exceptions

def IndentString( src, start = 0, increase = 4 ):
    """indent a mult-line string.
    
    start at indent start and increase by increase.
    """
    
    ss = src.split("\n")
    new_ss = []
    indent = start
    for s in ss:
        no = s.count( "(" )
        nc = s.count( ")" )
        if no > nc:
            new_ss.append( " " * indent + s)        
            indent += increase
        elif nc > no:
            indent -= increase
            new_ss.append( " " * indent + s)                                        
        else:
            new_ss.append( " " * indent + s)                                        
                        
    return "\n".join(new_ss)

class Model:
    def __init__(self):
        self.mGrammar = Grammar()
        self.mAlphabet = Alphabet()
        
    def __str__(self):
        return 'grammar:\n%s\nalphabet:\n%s' % (self.mGrammar, self.mAlphabet)
    
    def setAlphabet(self, alphabet):
        self.mAlphabet = alphabet
    
    def setGrammar(self, grammar):
        self.mGrammar = grammar
    
    def getGrammar(self):
        s = "\n".join( self.getGrammarLines())
        return IndentString( s )
    
    def getGrammarLines(self):
        """write the grammar for this chain."""
        
        lines = []
        lines += self.mGrammar.getGrammarLines()
        lines += self.mAlphabet.getGrammarLines()
        
        return lines
 
    def renameParameter(self, old, new):
        """rename a parameter in the model."""
        self.mGrammar.renameParameter( old, new )
 
    def evaluateRateMatrix(self, terminals = None):
        """returns the evaluated rate Matrix for the given terminals.
        
        This method uses the parameter section of the grammar to evaluate
        the expressions in the rate matrix.
        """
        
        if terminals == None:
            tt = self.mGrammar.getTerminals()
        else:
            tt = [terminals,]
                
        results = {}
        
        for terminals in tt:
            chain = self.mGrammar.getChain( terminals )
            lines = self.mGrammar.getParameterCode()
                    
            transitions = chain.getTransitions()
            
            for start, ends in transitions.items():
                for end, rate in ends.items():
                    if type(rate[1]) == StringType:
                        expr = rate[1]
                    else:
                        n = []
                        for x in rate[1]:
                            if x[:2] == "# ":
                                n.append( x[2:])
                            else:
                                n.append( x )
                        expr = " ".join( tuple(n) )
                
                    lines.append( "matrix[%s][%s] = %s" % ( (str(start), str(end), expr) ))
    
            # build the full rate matrix
            matrix = {}
            states = chain.getInitialStates().keys()
            for x in states:
                s = {}
                for y in states:
                    s[y] = 0
                matrix[x] = s
    
            ## fill the matrix with values
            try:
                exec("\n".join(lines))
            except SyntaxError, msg:
                raise Exceptions.ProgrammingError( 'error "%s" while computing rate matrix:\n%s\n' % \
                                                   (msg,"\n".join(lines)))
            
            ## set the diagonal elements
            for x in states:
                matrix[x][x] = 0
                s = sum( matrix[x].values() )
                matrix[x][x] = -s
            
            results[terminals] = matrix
            
        return results
            
    def evaluateTerminalFrequencies(self, terminals = None):
        """returns frequencies for states given the model.
        
        Initial frequencies might be given by parametric expressions and can
        thus not be read directly from the grammar but need to be computed
        via the parameters.
        """
        
        return self.mGrammar.evaluateTerminalFrequencies( terminals )    
    
    def swopConstantsVariables(self):
        """swop constants and variables in the model."""
        self.mGrammar.swopConstantsVariables()