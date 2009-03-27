from XGram.Model import PrintMembers

class Alphabet:

    def __init__(self):
        self.mTokens = None
        self.mComplements = None
        self.mExtensions = {}
        self.mName = ""
        self.mWildcard = None
        
    def setName(self, name):
        self.mName = name
        
    def __str__(self):
        return PrintMembers( self )
    
    def addExtension(self, start, ends):
        self.mExtensions[start] = ends

    def setTokens(self, tokens):
        self.mTokens = tokens
        
    def getTokens(self):
        return self.mTokens

    def setWildcard(self, wildcard):
        self.mWildcard = wildcard
        
    def setComplements(self, complements):
        self.mComplements = complements
        
    def getGrammar(self):
        return "\n".join( self.getGrammarLines())
    
    def getGrammarLines(self):
        """write the grammar for this chain."""
        
        lines = []
        lines.append( "( alphabet" )
        lines.append( "( name %s )" % self.mName )       
        
        lines.append( "( token ( %s ) )" % (" ".join(self.mTokens ) ) )
        if self.mComplements:
            lines.append( "( complement ( %s ) )" % (" ".join(self.mComplements ) ) )    
        
        for start, ends in self.mExtensions.items():
            lines.append( "( extend (to %s) %s)" % ( start, " ".join(map( lambda x: "(from %s) " % x, ends))))
        if self.mWildcard:
            lines.append( "( wildcard %s)" % self.mWildcard )
        
        lines.append( ") ;; finished alphabet %s" % self.mName )
        
        return lines
    
