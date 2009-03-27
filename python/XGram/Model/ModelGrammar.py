from types import *
import string, re
from XGram.Model import ProbabilisticModel, PrintStructure, PrintMembers
from XGram import Exceptions

class Grammar( ProbabilisticModel ):
    
    def __init__(self):
        ProbabilisticModel.__init__(self)
        self.mChains = {}
        self.mName = "unnamed"
        self.mRules = {}
        self.mOrderedRules = []
        self.mUpdateRules = None
        self.mUpdateRates = None
        self.mIsParametric = None
        self.mStartState = False
        self.mPseudoCounts = None
        self.mObservedCounts = None
        self.mNonTerminals = []
        self.mComments = []
        
    def setUpdateRules(self, flag ):
        self.mUpdateRules = flag

    def setUpdateRates(self, flag ):
        self.mUpdateRates = flag
        
    def setIsParametric(self, flag):
        self.mIsParametric = flag
        
    def getName(self):
        return self.mName
    
    def setName(self, name):
        self.mName = name
    
    def setPseudoCounts(self, counts):
        """set pseudocounts to counts."""
        self.mPseudoCounts = counts
        
    def setObservedCounts(self, counts):
        """set observed counts."""
        self.mObservedCounts = counts
        
    def addNonTerminal(self, name ):
        """add a non-terminal symbol."""
        if name not in self.mNonTerminals:
            self.mNonTerminals.append(name)
            
    def addComment(self, text):
        """add a comment to the grammar."""
        self.mComments.append( text )
        
    def addRule(self, rule):
        """add a rule from start to end with optional rate/probability and options.
        
        The order in which rules are added is preserved as the order of rules
        defines precedence.
        """
        start, end =rule.getStart(), rule.getEnd()
        if not self.mStartState:
            self.mStartState = start
        
        if start not in self.mRules:
            self.mRules[start] = {}
            
        self.mRules[ start ][end] = rule
        self.mOrderedRules.append( rule )
        
    def addChain(self, chain):
        self.mChains[chain.getTerminals()] = chain
        self.mParameters += chain.mParameters 

    def renameParameter(self, old, new):
        """rename a parameter from old to new."""

        ProbabilisticModel.renameParameter( self, old, new )

        for terminals, chain in self.getChains().items():
            chain.renameParameter( old, new )

    def getChain(self, terminals = None):
        """returns a chain specified by the terminals. If not terminals
        are given, the first is returned, if it is the only one."""
        if terminals:
            chain = self.mChains[terminals]
        else:
            if len(self.mChains) > 1:
                raise Exceptions.UsageError( "more than one chain - please supply terminals.")
            chain = self.mChains.values()[0]
        return chain

    def getChains(self):
        """return chains organized via their terminals."""
        return self.mChains
        
    def getTerminals(self):
        """return a list of all terminals."""
        return self.mChains.keys()
        
    def __str__(self):
        return PrintMembers( self )
    
    def getGrammar(self):
        return "\n".join( self.getGrammarLines())
        
    def checkContent(self):
        """check contents of chain and make consistent."""
        used_parameters = set()
        for terminal, chain in self.mChains.items():
            chain.checkContent( self )
            used_parameters.update( set(chain.getParameters().keys()))
#        make sure used parameters also includes the grammar!!!
#
#        defined_parameters = set(self.getParameters().keys())
#
#        extra_defined = defined_parameters.difference(used_parameters)
#        extra_used = used_parameters.difference(defined_parameters)
#
#        if extra_defined:
#            print extra_defined
#            for x in extra_defined:
#                self.removeParameter( x )
#
#        if extra_used:
#            raise XGram.Exceptions.UsageError( \
#                  "transitions used in chain, but not declared in grammar: %s" %\
#                   ",".join(tuple(extra_used)))
     
    def getGrammarLines(self):
        """write the grammar for this chain."""
        
        self.removeRedundantParameters()
        self.checkContent()
        
        lines = []
        lines.append( "( grammar " )
        lines.append( "( name %s )" % self.mName )       

        lines += [ ";;; %s" % x for x in self.mComments ]

        lines.append(";; flags")
        
        if self.mUpdateRules != None:
            if self.mUpdateRules:
                lines.append("( update-rules 1 )"  )
            else:
                lines.append("( update-rules 0 )"  )                
        
        if self.mUpdateRates != None:
            if self.mUpdateRates:
                lines.append("( update-rates 1 )"  )
            else:
                lines.append("( update-rates 0 )"  )                
                
        if self.mIsParametric == True:
            lines.append("( parametric )" )

        if self.mNonTerminals:
            lines.append(";; Nonterminal declarations and modifiers" )
            for t in self.mNonTerminals:
                lines.append( "( nonterminal (name %s) )" % t)
        
        lines.append(";; rules")
        ## start with start state
        if not self.mStartState:
            raise XGram.Exceptions.UsageError( "no start state defined." )
        
        for rule in self.mOrderedRules:
            lines.append( rule.getGrammar() )       
        
        lines.append(";; parameters")

        def printParameters( lines, section, is_const, is_explicit, is_rate ):
            v = filter( lambda x: x.isConst() == is_const and \
                                x.isExplicit() == is_explicit and \
                                x.isRate() == is_rate and \
                                x.isDummy() == False, \
                                self.mParameters )
            if len(v) == 0: return
            lines.append("( %s" % section)

            for p in v: lines += p.getGrammarLines()
                
            lines.append( ") ;; end of %s" % section)
                

        printParameters( lines, "param", is_const=False, is_explicit=False, is_rate=False) 
        printParameters( lines, "param", is_const=False, is_explicit=False, is_rate=True)         
        printParameters( lines, "const", is_const=True, is_explicit=False, is_rate=False) 
        printParameters( lines, "const", is_const=True, is_explicit=False, is_rate=True)                          
        printParameters( lines, "pgroup", is_const=False, is_explicit=True, is_rate=False) 
        printParameters( lines, "const-pgroup", is_const=True, is_explicit=True, is_rate=False)                                                   
        printParameters( lines, "rate", is_const=False, is_explicit=True, is_rate=True) 
        printParameters( lines, "const-rate", is_const=True, is_explicit=True, is_rate=True)                                                   

        lines.append( ";; end of params" )
        
        if self.mPseudoCounts:
            lines.append( "( pseudocounts ")
            for key, val in self.mPseudoCounts.items():
                if 'time' in val:
                    lines.append( "( %s %s %s )" % (key, str(val['counts']), str(val['time']) ))
                else:
                    lines.append( "( %s %s )" % (key, str(val['counts'])))
            lines.append( ") ;; end pseudocounts" )

        if self.mObservedCounts:
            lines.append( "( observed-counts ")
            for key, val in self.mObservedCounts.items():
                if 'time' in val:
                    lines.append( "( %s %s %s )" % (key, str(val['counts']), str(val['time']) ))
                else:
                    lines.append( "( %s %s )" % (key, str(val['counts'])))
            lines.append( ") ;; end observed counts" )
            
        lines.append(";; chains")         
         ## write chains
        for terminal, chain in self.mChains.items():
            lines += chain.getGrammarLines()
        
        lines.append( ") ;; end of grammar %s" % (self.mName) )
            
        return lines
        
    def copyParameters( self, 
                        other, 
                        copy = True,
                        map_this2other = {} ,
                        ignore_missing = False,
                        remove_missing = False,
                        skip_variables = False,
                        skip_constants = False,
                        skip_chains = False,
                        skip_rules = False):
        """copy rates from another model.
        
        map_this2other is a dictionary, that maps states in this grammar
        to states in the other grammar.
        
        If ignore_missing is True, rates will be left unchanged for those
        with entry in other. 
        
        If remove_missing is True, the missing transitions
        will be removed.
        
        If copy is False: values are not copied. Instead, if remove_missing is
        set to True, the transitions/params/const/... are filtered.
        
        This routine transfers paramers. It considers
        * the params and const section
        * each chain, if there is a matching chain for the same terminals in other.
        * the grammar rules

        Note: nothing is added to the grammar nor any chains.
        """
    
        def __updateParameter(p):
            """update a parameter. 
            
            Either copy values or delete keys.
            """
            
            for key, value in p.items():
                if key in map_this2other:
                    other_key = map_this2other[key]
                else:
                    other_key = key
                    
                if other_key in other.getParameters():
                    if copy:
                        p[key] = other.getParameter( other_key )
                    else:
                        del p[key]
                else:
                    if not ignore_missing:
                        raise IndexError, "parameter %s->%s not in model %s" % (key, other_key,
                                                                                other.getName())                    
                        
        ## transfer parameters
        for x in self.mParameters:
            if skip_constants and x.isConst(): continue
            if skip_variables and x.isVariable(): continue
            __updateParameter( x )
                
        ## tansfer rates for each chain
        if not skip_chains:
            for terminals, chain in self.mChains.items():
                if terminals in other.mChains:
                    chain.copyParameters( other.mChains[terminals],
                                          copy = copy, 
                                          ignore_missing = ignore_missing,
                                          remove_missing = remove_missing )
            
        ## transfer rates for the production rules
        if not skip_rules:
            xx = self.mRules.keys()
            for x in xx:
                if x in other.mRules:
                    yy = self.mRules[x].keys()
                    for y in yy:
                        if y in other.mRules[x]:
                            self.mRules[x][y] = other.mRules[x][y]
                        elif remove_missing:
                            del self.mRules[x][y]
                        elif not ignore_missing:
                            raise \
                            XGram.Exceptions.UsageError( "rule from %s to %s in %s not in %s" % \
                                                         (x, y, 
                                                          " ".join(self.getTerminals()),
                                                          " ".join(other.getTerminals())))
    
                elif remove_missing:
                    del self.mRules[x]
                elif not ignore_missing:
                    raise Exceptions.UsageError( "rule from %s in src, but not in other" % \
                                                 (x ))
            
    def evaluateTerminalFrequencies(self, terminals = None):
        """returns frequencies for states given the model.
        
        Initial frequencies might be given by parametric expressions and can
        thus not be read directly from the grammar but need to be computed
        via the parameters.
        
        The terminal frequencies are normalized to be 1. They not necessarily sum up
        to 1, for example, if codon frequencies are determined by nucleotide frequencies,
        but stop-codons are disallowed. 
        """
        
        if terminals == None:
            tt = self.getTerminals()
        else:
            tt = [terminals,]
                
        results = {}
        
        for terminals in tt:
            chain = self.getChain( terminals )
            lines = self.getParameterCode()
                    
            # this works by building a little code block and executing it.
            # All parameter variables are added into the name space and then
            # the expressions are added explicitely from the model.
            for states, expr in chain.getInitialStates().items():
                if type(expr[1]) == StringType:
                    expr = expr[1]
                elif type(expr[1]) == FloatType:
                    expr = str(expr[1])
                else:            
                    expr = " ".join(expr[1])
                
                lines.append( "result[%s] = %s" % (str(states), expr))
                
            result = {}
            exec("\n".join(lines))
            t = sum(result.values())
            for k,v in result.items():
                result[k] /= t
                
            results[terminals] = result
            
        return results
    
