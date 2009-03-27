import string, re
from XGram.Model import BuildExpression

class Rule:
    """a transformational rule in the grammar."""
    def __init__(self, 
                 start, end, 
                 rate = None, 
                 minlen = None,
                 gapsok = None,
                 annotations = None):
        self.mStart, self.mEnd = start, end
        self.mRate = rate
        self.mMinLen = minlen
        self.mGapsOk = gapsok
        if annotations:
            self.mAnnotations = annotations
        else:
            self.mAnnotations = []
        
    def getStart(self):
        return self.mStart
    def getEnd(self):
        return self.mEnd
    
    def setRate(self, rate):
        self.mRate = rate
    def setMinLen(self, minlen):
        self.mMinLen = minlen
    def setGapsOk(self, flag ):
        self.mGapsOk = flag
    def setAnnotations(self, annotations):
        self.mAnnotations = annotations
    def addAnnotation(self, annotation):
        self.mAnnotations.append( annotation )
        
    def getGrammar(self):
        """build expression for a rule."""
        options = []
        if self.mRate:
            options.append( "( %s )" % BuildExpression( self.mRate )  )

        if self.mGapsOk != None and self.mGapsOk:
            options.append( "( gaps-ok )" )
        
        if self.mMinLen != None:
            options.append( "( minlen %i )" % self.mMinLen)
            
        for annotation in self.mAnnotations:
            options.append(annotation.getGrammar())
            
        if options:
            return "( transform ( from ( %s ) ) ( to ( %s ) ) \n%s\n )" %\
                          (" ".join(self.mStart), " ".join(self.mEnd), "\n".join( options ))
        else:
            return "( transform ( from ( %s ) ) ( to ( %s ) ) )" %\
                      (" ".join(self.mStart), " ".join(self.mEnd) )

