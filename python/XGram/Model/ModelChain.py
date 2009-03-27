import string, re
from types import *
from XGram import Exceptions
from XGram.Model import ProbabilisticModel, BuildExpression, PrintMembers

class Chain( ProbabilisticModel ):
    
    def __init__(self):
        ProbabilisticModel.__init__(self)
        self.mTerminals = ("A0",)
        self.mTransitions = {}
        self.mInitialStates = {}
        self.mUpdatePolicy = None
        self.mIsMutable = True
    
    def setTerminals(self, terminals):
        """convert terminals to a tuple, so that it is hashable"""
        self.mTerminals = tuple(terminals)
    
    def isMutable(self):
        """returns whether chain is mutable. 
        
        A chain is mutable, if it has parameters."""
        for x in self.mParameters:
            if not x.isConst():
                return True
        else:
            return False
    
    def getTransitions(self):
        """returns the transition matrix."""
        return self.mTransitions
    
    def getInitialStates(self):
        return self.mInitialStates
            
    def setUpdatePolicy(self, policy):
        self.mUpdatePolicy = policy
        
    def getTerminals(self):
        return self.mTerminals
        
    def addTransition( self, start, end, rate ):
        """adds a transition to the chain. 
        
        Rate can be either string or expression. Strings are automatically
        converted to expressions."""
        if start not in self.mTransitions:
            self.mTransitions[start] = {}
        x = self.mTransitions[start] 
        if end not in x: x[end] = {}
        ## this is patch until I write an object handling expressions
        ## and rates.
        t, r = rate
        if type(r) == StringType:
            r = tuple(r.split(" "))
        x[end] = (t,r)

    def addInitialState(self, state, probability):
        self.mInitialStates[state] = probability
    
    def getGrammar(self):
        return "\n".join( self.getGrammarLines())
    
    def getGrammarLines(self):
        """write the grammar for this chain."""
        
        lines = []
        
        terminal_token = " ".join(self.mTerminals)
            
        lines.append( "( chain ")
        lines.append( "( terminal ( %s ) )" % terminal_token )
        
        if self.mUpdatePolicy:
            lines.append( "( update-policy %s )" % self.mUpdatePolicy)
            
        lines.append( ";; initial states" ) 
        for state, prob in self.mInitialStates.items():
            expr = BuildExpression( prob )
            lines.append( "( initial ( state ( %s ) ) ( %s ) )" % (" ".join(state), expr)) 
            
        lines.append( ";; substitution rates" )
        for start, tos in self.mTransitions.items():
            for to, rate in tos.items():
                expr = BuildExpression( rate )
                lines.append( "( mutate ( from ( %s ) ) ( to ( %s ) ) ( %s ) )" %\
                              ( " ".join(start), " ".join(to), expr ) ) 
        lines.append( ") ;; end of chain for terminals %s" % " ".join(self.mTerminals) )
        return lines
    
    def guessParameters(self):    
        """retrieve a list of parameters defined in the transitions.
        """
        n = set()
        rx = re.compile( "^[0-9+\-*/.Ee()]+$")
        for start, ends in self.mTransitions.items():
            for end, rate in ends.items():
                for x in rate[1]:
                    if type(x) != StringType: continue
                    ## ignore fix prefix
                    if x[:2] == "# ": x = x[2:]
                    if not rx.match( x ):
                        n.add(x)
        return n
    
    def checkContent(self, grammar ):
        """check content of chain and make consistent.
        
        If fix is set to true, names for used but undeclared
        parameters are added to the variables section. This
        is needed for parsing, because the chain is parsed
        apart from the grammar, it is not aware of any parameters
        declared elsewhere.
        
        The rate is set to 0.0.
        """
        # check for parameters
        defined_parameters = set(grammar.getParameters().keys())
        # add parameters defined in transitions as dummy parameters
        params = self.guessParameters() 
        for p in params: self.addParameter( p )         
        used_parameters = set(self.getParameters().keys())
        extra_used = used_parameters.difference(defined_parameters)
        ## removed some code here to remove extra_defined parameters
        ## don't think it necessary any more.
        
        if extra_used:
            raise Exceptions.UsageError( \
                 "transitions used in chain, but not declared in grammar: %s" %\
                   ",".join(tuple(extra_used)))
        
    def __str__(self):
        return PrintMembers( self )        

    def renameParameter(self, old, new ):
        """rename a parameter from old to new."""

        n = set()
        for start, ends in self.mTransitions.items():
            for end, rate in ends.items():
                if type(rate[1]) == StringType:
                    ends[end] = (rate[0], re.sub( old, new, rate[1] ))
                else:
                    new_rate = []
                    for v in rate[1]:
                        if type(v) != StringType: 
                            new_vals.append(v)
                        else:
                            new_rate.append( re.sub( old, new, v) )
                    ends[end] = (rate[0], new_rate) 

        for state, vals in self.mInitialStates.items():
            if type(vals[1]) == StringType:
                self.mInitialStates[state] = (vals[0], re.sub( old, new, vals[1] ))
            else:
                new_val = []
                if type(vals[1]) in (TupleType, ListType):
                    for v in vals[1]:
                        if type(v) != StringType: 
                            new_val.append( v )
                        else:
                            new_val.append( re.sub( old, new, v) )
                    self.mInitialStates[state] = ( vals[0], new_val )
                elif type(vals[1]) == StringType:
                    self.mInitialStates[state] = ( vals[0], re.sub(old,new,vals[1]) )
                    

        ProbabilisticModel.renameParameter( self, old, new )
        
    def copyParameters(self, other, 
                       copy = True,
                       ignore_missing = False,
                       remove_missing = False):
        """copy transition rates and initial states from other 
        to this chains.
        """
        xx = self.mTransitions.keys()
        for x in xx:
            if x in other.mTransitions:
                yy = self.mTransitions[x].keys()
                for y in yy:
                    if y in other.mTransitions[x]:
                        if copy:
                            self.mTransitions[x][y] = other.mTransitions[x][y]
                    elif remove_missing:
                        del self.mTransitions[x][y]
                    elif not ignore_missing:
                        raise \
                        XGram.Exceptions.UsageError( "transition from %s to %s in %s not in %s" % \
                                                     (x, y, 
                                                      " ".join(self.getTerminals()),
                                                      " ".join(other.getTerminals())))

            elif remove_missing:
                del self.mTransitions[x]
            elif not ignore_missing:
                raise XGram.Exceptions.UsageError( "transition from %s in %s not in %s" % \
                                                   (x, 
                                                    " ".join(self.getTerminals()),
                                                    " ".join(other.getTerminals())))
        
        xx = self.mInitialStates.keys()
        for x in xx:
            if x in other.mInitialStates:
                if copy:
                    self.mInitialStates[x] = other.mInitialStates[x]
            elif remove_missing:
                del self.mInitialStats[x]
            elif not ignore_missing:
                raise XGram.Exceptions.UsageError( "initial state %s in %s not in %s" % \
                                                   (x, 
                                                    " ".join(self.getTerminals()),
                                                    " ".join(other.getTerminals())))                
        
