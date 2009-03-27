"""Module for running XGram iteratively.

"""
import string, os, re, time, sys

import XGram

class Iteration:
    
    mLogLevel = 0
    
    def __init__(self,
                 xgram = None,
                 max_iterations = 0,
                 min_log_likelihood_difference = 1e-3 ,
                 tree_model = None ):
        """setup run. 
        
        Use xgram as your xgram executable. If None is given,
        use xgram with default settings.
        
        niterations/min_log_likelihood_difference determine the
        number of iterations performed. Niterations is the maximum
        bound of iterations. If set to 0, an unlimited number of iterations
        are run or until the likelihood difference between the current
        and previous run is less than min_log_likelihood_difference.
        
        If tree_model is None, a tree is assumed to be part of the data file.
        if tree_model is "auto", then the same grammar for training is used
        to build the tree. Otherwise, the model in tree_model is used to built
        a tree.
        """
        
        if xgram:
            self.mXGram = xgram
        else:
            self.mXGram = XGram.XGram()
            
        self.mMaxIterations = max_iterations
        self.mMinLogLikelihoodDifference = min_log_likelihood_difference
    
        self.mResults = []

        self.mDumpLevel = 0
    
        ## pattern for dumped files
        self.mDumpPatternModel = None
        self.mDumpPatternLog = None
        self.mDumpPatternData = None
        
        self.mTreeModel = tree_model
        
    def setLogLevel(self, loglevel ):
        self.mLogLevel = loglevel

    def setTreeModel(self, tree_model = "auto" ):
        """sets the tree model."""
        self.mTreeModel = tree_model

    def setDumpPatternModel(self, pattern = None):
        """pattern to dump models. Should contain one %s.
        
        Do not dump, if set to None.
        """
        self.mDumpPatternModel = pattern

    def setDumpPatternLog(self, pattern = None):
        """pattern to dump log files. Should contain one %s.
        
        Do not dump, if set to None.        
        """        
        self.mDumpPatternLog = pattern

    def setDumpPatternData(self, pattern = None):
        """pattern to dump data. Should contain one %s.
        
        Do not dump, if set to None.        
        """        
        self.mDumpPatternData = pattern

    def run(self, model, data, update_model = None):
        """run xgram iteratively. 
        
        After every iteration, the callback function is called with
        the trained grammar as a parameter.
        
        update_model is the callback function. If callback is not given, 
        the var and const section of the model are swopped.
        """
        xgram = XGram.XGram()
        
        iteration = 1
        
        results = []
        
        training_model = model

        ## check if we need to build a tree first:
        if self.mTreeModel:
            if self.mTreeModel != "auto":
                t1 = time.time()
                result = xgram.buildTree( self.mTreeModel, data )
                ## build tree and add it to data
                
                t2 = time.time()
                if self.mLogLevel >= 2:
                    print "# tree estimation: seconds=%i" % (t2 - t1)
                    sys.stdout.flush()
                # tree has been added to data, no recomputation of tree
                work_data = result.getData()
                tree_model = None
            else:
                ## re-estimation of tree at every step
                work_data = data
                tree_model = "auto"
        else:
            work_data = data
            tree_model = None
            
        last_ll = 0
        t0 = time.time()
        while 1:
            
            t1 = time.time()
            result = self.mXGram.train( training_model, 
                                        data = work_data,
                                        tree_model = tree_model)
            
            ll = result.getLogLikelihood()
            n = result.getNumIterations()
            results.append(result)
            
            t2 = time.time()

            ## various dump options
            if self.mLogLevel >= 1:
                print "# iteration %i: seconds=%i, iterations=%i, log-likelihood=%f" % \
                    (iteration, t2 - t1, n, ll)
                sys.stdout.flush()
                
            if self.mDumpPatternModel:
                outfile = open(self.mDumpPatternModel % str(iteration), "w" )
                outfile.write( result.getModel().getGrammar() + "\n")
                outfile.close()
                
            if self.mDumpPatternData:
                outfile = open(self.mDumpPatternData % str(iteration), "w" )
                outfile.write( result.getData() + "\n")
                outfile.close()                               
                
            if self.mDumpPatternLog:
                outfile = open(self.mDumpPatternLog % str(iteration), "w" )
                outfile.write( result.getLog() + "\n")
                outfile.close()
                
            ## check if iteration shall continue
            if last_ll:
                d = abs( (ll - last_ll) / last_ll)
                if d < self.mMinLogLikelihoodDifference:
                    if self.mLogLevel >= 2:
                        print "# leaving loop because log-likelihood difference %f smaller than minimum" % d
                    break

            last_ll = ll

            if self.mMaxIterations and iteration == self.mMaxIterations:
                if self.mLogLevel >= 2:
                    print "# leaving loop because iteration %i at maximum" % iteration
                break

            ## prepare for next iteration
            training_model = result.getModel()        
            if update_model:
                update_model( training_model )
            else:
                training_model.swopConstantsVariables()

            iteration += 1

        return results
        