from types import *

import XGram
from XGram import Exceptions
from XGram.Model import ModelParameter
    
class ProbabilisticModel:
    """base class for probabilistic models.
    
    This class manages parameters and constants.
    """
    def __init__(self):
        
        self.mParameters = []
        
        self.mName = "unnamed"
        self.mParameterMap = {}
        
    def setName(self, name):
        self.mName = name
    
    def addRate( self, param, is_const = False ):
        """add a rate parameter for this chain."""
        self.mParameters.append( ModelParameter.ParameterGroup( param,
                                                 is_const = is_const, 
                                                 is_explicit = True) )
    
    def addProbability(self, param, is_const = False):
        """add a probability parameter for this chain.
        
        Parameter needs to be a tuple of values.
        """
        self.mParameters.append( ModelParameter.ParameterGroup( param,
                                                 is_const = is_const, 
                                                 is_explicit = True) )

    def addParameter(self, param, is_rate = False, is_const = False,
                     is_explicit=False, normalize = False):
        """master function for adding parameters."""

        # normalize parameters, if they are probabilities                
        # normalize parameters, if they are probabilities        
        if normalize and type(param[0]) != StringType:
            total = sum( map( lambda x: x[1], param) )
            param = map( lambda x: (x[0], x[1] / total), param )
            
        self.mParameters.append(  ModelParameter.ParameterGroup(param,
                                                                is_const = is_const, 
                                                                is_explicit = is_explicit ) )
        
    def addVariable(self, param, normalize = False):
        """register variable parameters for this chain.
        
       Parametes are tuples of (identifier, value).
        A list of tuples can also be submitted for probabilities.
        [ (identifier, value), (identifier, value), ... )
        
        If they are probabilities, they are normalized
        so that they sum to one (if normalize is set to True).
        """
        self.addParameter( param, 
                           is_const = False,
                           is_explicit = False,
                           normalize = normalize )
        
    def addConst(self, param, normalize = False):
        """register constant parameters for this chain.
        """        
        self.addParameter( param, 
                           is_const = False,
                           is_explicit = False,
                           normalize = normalize )

    def removeRedundantParameters( self ):
        """remove all redundant parameters.
        """
        self.mParameterMap = {}
        new_parameters = []
        found = set()
        
        for x in self.mParameters:
            nredundant = 0
            for key, value in x.items():
                if key in found: nredundant += 1
                found.add(key)
            if nredundant == 0:
                new_parameters.append( x )
            elif nredundant != len(x):
                raise XGram.Exceptions.UsageError("incomplete redundant group: %s" %\
                                                   str(x.keys()))

            
        self.mParameters = new_parameters
        
    def getParameterMap(self, only_consts = False, only_variables = False):
        """return a map of parameters to values.
        
        This function extracts this information from all the parameter groups
        that have been defined.
        """
        ## note: better might be a flag to denote consts/var
        parameter_map = {}
            
        for x in self.mParameters:
            if only_consts and not x.isConst(): continue
            if only_variables and x.isConst(): continue
            
            for key, value in x.items():
                parameter_map[key] = value
        
        return parameter_map
    
    def buildParameterMap(self, only_consts = False, only_variables = False):
        """build a flat dictionary of parameters to values."""
        self.mParameterMap = self.getParameterMap( only_consts, only_variables )
                
    def getParameter(self, parameter):
        """returns the value of a given parameter."""
        if not self.mParameterMap:
            self.buildParameterMap()
        return self.mParameterMap[parameter]
    
    def getParameters(self):
        """returns a dictionary of all parameters.
        
        both variable and const parameters are returned.
        """
        if not self.mParameterMap:
            self.buildParameterMap()
        return self.mParameterMap

    def removeParameter(self, parameter):
        """remove parameter from this model."""
        
        if parameter not in self.getParameters():
            raise XGram.Exceptions.UsageError("parameter %s not declared." % parameter)        
        
        for x in range( len(self.mParameters) ):
            param_group = self.mParameters[x]
            if parameter in param_group:
                old_value = param_group[parameter]
                del param_group[parameter]
                if self.mParameterMap:
                    del self.mParameterMap[parameter]
                if len(param_group) == 0:
                    del self.mParameters[x]
                return old_value
    
    def renameParameter(self, old, new):
        """rename a parameter from old to new."""

        ## do not throw an error, as this function is called
        ## by derived classes. 
        if old not in self.getParameters():
            return 
                
        for x in range( len(self.mParameters) ):
            param_group = self.mParameters[x]
            if old in param_group:
                if new in param_group:
                    raise XGram.Exceptions.UsageError("parameter new already exists in group %s." %\
                                                      (new, str(param_group)))        
                if self.mParameterMap:
                    self.mParameterMap[new] = self.mParameterMap[old] 
                    del self.mParameterMap[old]
                param_group[new] = param_group[old]
                del param_group[old]
                
        ## somewhere I have a bug with parameters related to
        ## update of a grammar map and chain maps
        self.buildParameterMap()
        
    def setParameter(self, parameter, new_value):
        """set a parameter to a new value.
        
        The parameter has to exist. The old value is returned.
        """

        if parameter not in self.getParameters():
            raise XGram.Exceptions.UsageError("parameter %s not declared." % parameter)        
        
        for param_group in self.mParameters:
            if parameter in param_group:
                old_value = param_group[parameter]
                param_group[parameter] = new_value
                if self.mParameterMap:
                    self.mParameterMap[parameter] = new_value
                return old_value
        return None
        
    def getVariables(self):
        """return (possibly nested) list with variables."""
        return map( lambda x: x.getParameter(), filter( lambda x: not x.isConst(), self.mParameters ))

    def getConsts(self):
        """return (possibly nested) list with constants."""
        return map( lambda x: x.getParameter(), filter( lambda x: not x.isConst(), self.mParameters ))
    
    def getParameterCode(self):
        """returns a code block (lines) with the parameters of this model.
        
        Note: this method is maybe not the most efficient way to do this and
        not the safest either ;=)
        """
        lines = []

        for parameter, value in self.getParameters().items():
            lines.append( "%s = %s" % (parameter, str(value)))
        return lines
 
    def guessParameters(self):
        """returns a list of identifiers, that this object 
        regards as parameters."""
        return []
    
    def addParameters(self, other):
        """adds parameters from another object."""
        missing_parameters = self.guessParameters()

        for p in missing_parameters:
            self.addParameter( p )
  
    def swopConstantsVariables(self):
        """swop constants and variables. What has been a constant is now
        variable and vice versa.
        """
        for x in self.mParameters:
            if x.isConst():
                x.setConst( const = False )
            else:
                x.setConst( const = True )
        self.mParameterMap = {}
        
    def moveVariableToConst(self, parameter ):
        """move a parameter from variables to consts.
        
        This will apply to the whole group.
        """
        for x in self.mParameters:
            if parameter in x:
                x.setConst( const = True )
