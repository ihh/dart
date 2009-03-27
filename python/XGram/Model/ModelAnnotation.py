class Annotation:
    """an annotation for a rule."""
    def __init__(self, row, column, label):
        self.mRow = row
        self.mColumn = column
        self.mLabel = label
    
    def getGrammar(self):
        return "( annotate ( row %s ) ( column %s ) ( label %s ) )" %\
            ( self.mRow, self.mColumn, self.mLabel)
            
    def setRow(self, row):
        self.mRow = row
        
    def setColumn(self, column):
        self.mColumn = column
        
    def setLabel(self, label):
        self.mLabel = label
        
