"""Wrapper around XGram executable.

This module defines the XGram class, which is an interface
for the XGram executable. XGram has two methods - train
and annotate, which return an XGramResult type object.

Note: maybe later this class can make use of the XGram C++ 
libraries directly.
"""## Xgram class
## Wrapper around the xgram executable

import subprocess, tempfile
import os.path, sys, shutil, re, string, time

from types import *

import Exceptions
import Parser
import Utils.Stk

class XGramResult:
    """container for XGram result."""
    def __init__( self, 
                 model = None, 
                 log = None, 
                 data = None, 
                 ):
        self.mModel = model
        self.mLog = log
        self.mData = data
        self.mDataParsed = None
        
        self.__parseLog()
        self.mIterations = []
        
        self.mRunTime = -1
        
    def getModel( self ):
        """return the model."""
        return self.mModel
    
    def getLog( self ):
        """return logfile."""
        return self.mLog
    
    def getRunTime(self):
        """return runtime."""
        return self.mRunTime

    def setRunTime(self, runtime):
        """return runtime."""
        self.mRunTime = runtime
    
    def getAnnot(self, pLoc, pKey):
        """return annotation."""
        annot = None
        if pLoc == 0:
            annot = self.mDataParsed.gf[pKey]
        else:
            annot = self.mDataParsed.gc[pKey]
        return annot

    def getData( self ):
        """return data."""
        return self.mData
    
    def getTree(self):
        """return tree."""
        if not self.mDataParsed: self.parseData()
        if 'NH' in self.mDataParsed.gf:
            return self.mDataParsed.gf['NH']
        else:
            return None
    
    def getDataParsed( self ):
        """return parsed data."""
        return self.mDataParsed
    
    def getIterations( self ):
        """return a list of iterations."""
        if not self.mIterations: self.__parseLog()
        return self.mIterations

    def getNumIterations( self ):
        """return number of iterations."""
        if not self.mIterations: self.__parseLog()        
        return len( self.mIterations )
    
    def getLogLikelihood( self ):
        """get best log likelihood."""
        if not self.mIterations: self.__parseLog()
        return self.mLogLikelihood

    def getPostIterationLogLikelihood( self ):
        """get post-iteration log likelihood."""        
        if not self.mIterations: self.__parseLog()
        return self.mPostIterationLogLikelihood

    def parseData( self ):
        """parse the datafile."""
        self.mDataParsed = Utils.Stk.Stk(pList=self.mData)

    def __parseLog( self ):
        """parse the logfile."""
        self.mIterations = []
        rx_iter = re.compile( "EM iteration #\d+: log-likelihood = ([-0-9.]+) bits" )
        rx_pi_logl = re.compile( "Post-iteration log-likelihood = ([-0-9.]+) bits" )
        rx_best_logl = re.compile( "Best log-likelihood: ([-0-9.]+)" )        
        
        for line in self.mLog:
            
            if rx_iter.match( line ):
                self.mIterations.append( float( rx_iter.search( line ).groups()[0] ) )
            elif rx_best_logl.match( line ):
                self.mLogLikelihood = float( rx_best_logl.search( line ).groups()[0] )
            elif rx_pi_logl.match( line ):
                self.mPostIterationLogLikelihood = float( rx_pi_logl.search( line ).groups()[0] )

class XGram:
    """Wrapper around the XGram executable.
    
    Options provided are
    --train
    --tree
    --log
    
    Other options can be added via the addOption method.
    """
    
    mExecutable = "xgram"
    mLogLevel = 6
    
    mMinIncrement = 0.0
    mMaxIterations = 0
    
    def __init__( self, pDartDir=None ):
        self.mOptions = []
        self.mDebug = False
        self.mExtraOptions = []
        
        ## runtime of last job
        self.mRunTime = -1
        if pDartDir: 
            self.dartDir = pDartDir
            self.dartBinDir = self.dartDir + '/bin/'
            self.mExecutable = self.dartBinDir + self.mExecutable
        
    def addOption( self, option ):
        """adds an options to the XGram call."""
        self.mExtraOptions.append( option )
        
    def setLogLevel( self, loglevel ):
        """sets the XGram loglevel."""
        self.mLogLevel = loglevel
    
    def getRunTime(self):
        """returns execution duration."""
        return self.mRunTime

    def setMinIncrement(self, mininc):
        """minimum (relative) increment between iterations.
        
        A value of 0 is uses the xgram default.
        """
        self.mMinIncrement = mininc
        
    def setMaxIterations(self, max_iterations):
        """maximum number of iterations.
        
        A value of 0 means no limit.
        """
        self.mMaxIterations = max_iterations
    
    def setDebug(self, debug = True):
        """sets/unsets debug flag. 
        
        If debug is true, temporary files won't be deleted and
        a message will be printed about where they can be found and
        the statement that has been executed.
        """
        self.mDebug = debug
    
    def __runStatement( self, statement ):
        """run a statement."""
        
        if self.mDebug:
            print "# executing statement %s" % statement
            print "# working directory: %s" % self.mTempdir
            sys.stdout.flush()        
            
        self.mFilenameOutput = self.mTempdir + "/stdout"
        self.mFilenameError = self.mTempdir + "/stderr"        
            
        t0 = time.time()
        s = subprocess.Popen( statement, 
                              shell = True, 
                              stdout = open(self.mFilenameOutput, "w"), 
                              stderr = open(self.mFilenameError, "w"),
                              cwd = self.mTempdir, 
                              close_fds = True )

        s.communicate()

        out = open(self.mFilenameOutput).readlines()
        err = open(self.mFilenameError).readlines()

        if s.returncode != 0:
            raise Exceptions.UsageError, "Error in running XGram\nstatement=%s\n%s\nTemporary directory in %s" % \
            ( statement, err, self.mTempdir )

        t1 = time.time()
        
        self.mRunTime = t1 - t0

        out = open(self.mFilenameOutput).readlines()
        err = open(self.mFilenameError).readlines()

        return out, err

    def setupModel( self, model, template = "input" ):
        """setups the file with a model.
        
        template is the filename for the model (usually
        input or tree).
        """
        
        if type( model ) == StringType:
            filename = model
        else:
            filename = "%s/%s.eg" % ( self.mTempdir, template )
            outfile = open( filename, "w" )
            outfile.write( model.getGrammar() )
            outfile.close()

        return filename 
    
    def setupData( self, data ):
        """setup data file.
        
        If len(data) is less than 200, it is assumed to be a filename.
        Otherwise, a file is created in the temporary directory.
        """
        if len(data) < 200:
            self.mFilenameData = os.path.abspath( data )            
            if not os.path.exists( self.mFilenameData ):
                raise Exceptions.UsageError("data file %s not found." % data )
        else:
            self.mFilenameData = "%s/data.stk" % (self.mTempdir)
            outfile = open(self.mFilenameData,"w")
            outfile.write( data + "\n" )
            outfile.close()

    def setupInputModel( self, model ):
        """create temporary directory and write file with input model."""
        
        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameInputModel = self.setupModel( model, "input" )
    
    def startUp(self, model, data, tree_model = None):
        """common startup code."""
        self.mOptions = []
        self.setupInputModel(model)
        self.setupData(data)

        if tree_model:
            if tree_model == "auto":
                self.mFilenameTreeModel = self.mFilenameInputModel
            else:
                self.mFilenameTreeModel = self.setupModel( tree_model, "tree" )
            self.mOptions.append( "--tree %s" % self.mFilenameTreeModel )
    
        if self.mMaxIterations:
            self.mOptions.append( "--maxrounds %i" % self.mMaxIterations)
        
        if self.mMinIncrement:
            self.mOptions.append( "--mininc %f" % self.mMinIncrement)            
    
    def cleanUp( self ):
        """remove temporary directory."""
        shutil.rmtree( self.mTempdir )
    
    def annotate( self, model, data, tree_model = None):
        """run XGram and return data annotated by the model."""
        
        self.startUp(model, data, tree_model)

        statement = "%s --grammar %s -log %i %s %s %s" % \
            ( self.mExecutable, 
              self.mFilenameInputModel, 
              self.mLogLevel, 
               " ".join( self.mExtraOptions ),           
               " ".join( self.mOptions ), 
               self.mFilenameData )
        
        out, err = self.__runStatement( statement )
        
        if not self.mDebug:
            self.cleanUp()
            
        result = XGramResult( model = None, 
                              data = out, 
                              log = err )
        
        result.setRunTime( self.getRunTime() )
        
        return result                    
            
    def train( self, model, data, 
               tree_model = None,
               no_annotate = False ):
        """run XGram and return the input model trained on data.
        
        if no_annotate is True, no re-annotation will be attempted.
        """
        
        self.startUp(model, data, tree_model)        
        
        filename_trained_model = "%s/trained.eg" % self.mTempdir
        
        if no_annotate:
            self.mOptions.append( "--noannotate" )
            
        statement = "%s --train %s --grammar %s -log %i %s %s %s" % \
        ( self.mExecutable, 
          filename_trained_model, 
          self.mFilenameInputModel, 
          self.mLogLevel, 
           " ".join( self.mExtraOptions ),           
           " ".join( self.mOptions ), 
           self.mFilenameData )
    
        out, err = self.__runStatement( statement )
        
        lines = open( filename_trained_model, "r" ).readlines()
        try:
            output_model = Parser.parseGrammar( lines )
        except Exceptions.ParsingError, msg:
            print "parsing failed with msg:", msg
            print "input file are in directory %s." % self.mTempdir
            sys.exit( 1 )

        if not self.mDebug:
            self.cleanUp()
        
        result = XGramResult( model = output_model, 
                              data = out, 
                              log = err )
        
        result.setRunTime( self.getRunTime() )
                
        return result
    
    def buildTree( self, model, data ):
        """run XGram to build a tree.
        
        The tree is stored in the data-section of the
        result.
        """
        self.startUp(model, data)        

        statement = "%s --tree %s -log %i %s %s %s" % \
        ( self.mExecutable, 
          self.mFilenameInputModel, 
          self.mLogLevel, 
           " ".join( self.mExtraOptions ),           
           " ".join( self.mOptions ), 
           self.mFilenameData )
    
        out, err = self.__runStatement( statement )
        
        if not self.mDebug:
            self.cleanUp()
        
        result = XGramResult( model = None, 
                              data = out, 
                              log = err )
        
        result.setRunTime( self.getRunTime() )
                
        return result
    
        
