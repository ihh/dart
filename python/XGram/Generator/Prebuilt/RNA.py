## Build models for RNA

from XGram.Generator import \
    AlphabetRNA, Transitions, InitialStatesFromTransitions, Chain, Grammar
from XGram.Model import Rule, Annotation
import XGram.Parser


class TransitionsSymmetricStem(Transitions):
    """class for generating a transition matrix with symmetric
    rates such that AC -> AT has the same rate as CA -> TA.
    """
    
    def __init__(self, fix_frequencies = False, *args,**kwargs):
        Transitions.__init__(self,*args,**kwargs)
        self.mFixFrequencies = fix_frequencies       
         
    def getRate(self, l1, l2):
        """returns the rate variable for emissions l1 and l2.
        
        returns ABCD for A:B C:D and B:A D:C
        returns ABCC for A:B C:C and B:A C:C         
        returns AACD for A:A C:D and A:A D:C
        """
        a,b = map(lambda x: x.upper(), l1)
        c,d = map(lambda x: x.upper(), l2)
        
        ## sort by lexicographical order of a and b
        if a > b: 
            a, b = b, a
            c, d = d, c
        ## but if they are the same, sort by c and d
        elif a == b and c > d:
            c, d = d, c
                
        return "r%s%s%s%s" % (a,b,c,d)

    def getProbability(self, l1):
        return "p%s%s" % tuple( map(lambda x:x.upper(), l1))
    
    def getInitialProbability(self, emissions):
        return self.getProbability( emissions )
    
    def registerParameters(self, chain):
        """register parameters for this chain."""
        alphabet = self.mChainGenerator.getAlphabet()            
        l = len(alphabet)
        written = {}
        for x in range(0, l):
            l1 = ("".join(alphabet[x])).upper()
            for y in range(0, l):
                if x == y: continue
                l2 = ("".join(alphabet[y])).upper()
                rate = self.getRate( l1, l2)
                if rate in written: continue
                written[rate] = True
                if self.mIsConst:
                    chain.addConst( (rate, self.mRateGenerator.getRate(rate) ))
                else:
                    chain.addVariable( (rate, self.mRateGenerator.getRate(rate) ))
        
        probs = []
        for x in alphabet:
            probs.append( (self.getProbability(x), self.mRateGenerator.getRate(x) ) )
            
        if self.mFixFrequencies:
            chain.addConst( probs )
        else:
            chain.addVariable( probs )
            
    def addTransitions(self, chain, l1, l2 ):
        chain.addTransition( l1, l2, 
                             ('rate', 
                              (self.getRate(l1, l2), 
                               "*",
                               self.getProbability(l2))))
        chain.addTransition( l2, l1, 
                             ('rate', 
                              (self.getRate(l2, l1),
                               "*",
                               self.getProbability(l1))))                               

class TransitionsSymmetricStemWeakStrong(TransitionsSymmetricStem):
    """class for generating a transition matrix with symmetric
    rates such that AC -> AT has the same rate as CA -> TA.
    
    Transition rates are further parameterized by a rate
    factor kXX, where XX is WW, WS, SW, SS with W and S being
    weak and strong base pairs, respectively.
    
    Strong base pairs are Watson Crick basepairs and
    the Wobble base pair G:U.
    
    The rates are initially in the const section.
    """
    
    def __init__(self, *args,**kwargs):
        TransitionsSymmetricStem.__init__(self,*args,**kwargs)

    def registerParameters(self, chain):
        """register parameters for this chain.
        
        adds kWW, kWS, etc. as consts.
        """
        TransitionsSymmetricStem.registerParameters(self, chain)

        for x in ('WW','WS', 'SW', 'SS'):
            rate = self.mPrefix + "k%s" % x
            chain.addConst( (rate, self.mRateGenerator.getRate( rate)) )

    def getStrength(self, l1):
        """returns the strength (S or W) of base pair l1."""
        bp = "".join( sorted(map(lambda x: x.upper(), l1)))
        if bp in ('AU', 'CG', 'GU' ):
            return "S"
        else:
            return "W"

    def addTransitions(self, chain, l1, l2 ):
        """add transitions to chain between states l1 and l2.
        
        A transition is of the type:
        rABCD * pCD * kSS
        rABCD: base line transition rate between pair AB -> CD
        pCD: equilibrium frequency of base pair CD
        kSS: interaction strength of base pair AB and CD
        S: strong, W: weak
        """

        strength = self.mPrefix + "k%s%s" % (self.getStrength(l1),
                                             self.getStrength(l2))
        
        chain.addTransition( l1, l2, 
                             ('rate', 
                              (self.getRate(l1, l2), 
                               "*",
                               self.getProbability(l2),
                               "*",
                               strength)))
        chain.addTransition( l2, l1, 
                             ('rate', 
                              (self.getRate(l2, l1),
                               "*",
                               self.getProbability(l1),
                               "*",
                               strength)))
            

            
class GeneratorGrammarPFold( Grammar ):
    def __init__(self, *args, **kwargs):
        Grammar.__init__(self, *args, **kwargs)
    
    def buildGrammar(self, grammar):
        """build PFold like grammar."""
        
        Grammar.buildGrammar( self, grammar)
        
        # state pfoldS
        grammar.addRule( Rule( ("pfoldS",),
                               ("pfoldL",),
                               ("prob", (1.0,)) ) )

        grammar.addRule( Rule( ("pfoldS",),
                               ("pfoldB",),
                               ("prob", (1.0,)) )  )

        # state pfoldF
        annotation1 = Annotation( row = "PFOLD",
                                  column = "LNUC",
                                  label = "<" )
        annotation2 = Annotation( row = "PFOLD",
                                  column = "RNUC",
                                  label = ">" )

        
        grammar.addRule( Rule( ("pfoldF",),
                               ("LNUC", "pfoldF*", "RNUC"),
                               gapsok = True,
                               annotations = [annotation1,annotation2] ) )
        
        grammar.addRule( Rule( ("pfoldF*",),
                               ("pfoldF",),
                               ("prob", (1.0,) ) ) )
        
        grammar.addRule( Rule( ("pfoldF*",),
                               ("pfoldB",),
                               ("prob", (1.0,) ) ))
        
        ## state pfoldB
        grammar.addRule( Rule( ("pfoldL",),
                               ("pfoldF",),
                               ("prob", (1.0,)) ) )

        grammar.addRule( Rule( ("pfoldL",),
                               ("pfoldU",),
                               ("prob", (1.0,)) ) )
        
        ## state pfoldB
        grammar.addRule( Rule( ("pfoldB",),
                               ("pfoldL", "pfoldS") ) )

        annotation3 = Annotation( row = "PFOLD",
                                  column = "NUC",
                                  label = "_" )
        
        ## state pfoldB
        grammar.addRule( Rule( ("pfoldU",),
                               ("NUC", "pfoldU*"), 
                               gapsok = True,
                               annotations = [annotation3,]))
        
        grammar.addRule( Rule( ("pfoldU*",),
                               (),
                               ("prob", (1.0,)) ) )


def buildNullRNA( ):
    """simple RNA Model."""
    
    alphabet = XGram.Model.Alphabet()
    Alphabets.AlphabetRNA().buildGrammar( alphabet )
    
    chain = XGram.Model.Chain()
    Chains.Chain( alphabet = alphabet.getTokens() ).buildGrammar( chain )
    chain.setUpdatePolicy( "rev" )
    
    grammar = XGram.Model.Grammar()
    grammar.addChain( chain )
    grammar.setName( "nullrna")
    
    Grammars.GrammarLinearSequence( no_gaps = True).buildGrammar( grammar )
    
    model = XGram.Model.Model()
    model.setGrammar( grammar )
    model.setAlphabet( alphabet )
    
    return model

def buildPFold( strand_symmetric = False,
                fix_frequencies = False,
                weak_strong = False,
                copy_parameters = None ):
    """build the PFold grammar.
    
    if copy_parameters is set to a filename of a grammar, then parameters from
    the grammar are copied to the generated chain.
    
    if weak_strong is set to true, than the strand-symmetric/weak-strong
    model is built. Note: for this option to be in effect strand_symmetric 
    has to be set to True .
    """
    alphabet = XGram.Model.Alphabet()
    XGram.Generator.AlphabetRNA().buildGrammar( alphabet )
    
    # build di-nucleotide alphabet
    aa_mono = []
    for x in alphabet.getTokens(): 
        aa_mono.append( (x,) )

    chain1 = XGram.Model.Chain()
    XGram.Generator.Chain( alphabet = aa_mono ).buildGrammar( chain1 )
    chain1.setUpdatePolicy( "rev" )
    chain1.setTerminals( ('NUC',) )

    aa_di = []
    for x in alphabet.getTokens(): 
        for y in alphabet.getTokens():
            aa_di.append( (x,y) )

    chain2 = XGram.Model.Chain()
        
    if strand_symmetric:
        if weak_strong:
            generator_transitions = \
                TransitionsSymmetricStemWeakStrong( fix_frequencies = fix_frequencies )            
        else:
            generator_transitions = \
                TransitionsSymmetricStem( fix_frequencies = fix_frequencies )
        generator_initial_states = InitialStatesFromTransitions()
        Chain( generator_transitions = generator_transitions,
               generator_initial_states = generator_initial_states,
               alphabet = aa_di ).buildGrammar( chain2 )
        chain2.setUpdatePolicy( "parametric" )
    else:
        Chain( alphabet = aa_di ).buildGrammar( chain2 )        
        chain2.setUpdatePolicy( "rev" )        
        
    chain2.setTerminals( ('LNUC', 'RNUC') )

    grammar = XGram.Model.Grammar()    
    grammar.addChain( chain1 )
    grammar.addChain( chain2 )    
    
    ## build Grammar for PFold

    GeneratorGrammarPFold().buildGrammar( grammar )

    if strand_symmetric:
        if weak_strong:
            grammar.setName( "pfold_strand_symmetric_weak_strong" )
        else:
            grammar.setName( "pfold_strand_symmetric")
    else:
        grammar.setName( "pfold")        

    ## if desired, copy data from another chain
    if copy_parameters:
        other = XGram.Parser.parseGrammar( open(copy_parameters, "r").readlines() )
        if strand_symmetric:
            ## first copy rules over
            grammar.copyParameters( other.mGrammar, 
                                    remove_missing = True,
                                    skip_rules = False,
                                    skip_variables = True,
                                    skip_constants = True,
                                    skip_chains = True )
            ## filter states from the chains
            grammar.copyParameters( other.mGrammar, 
                                    remove_missing = True,
                                    copy = False,
                                    skip_rules = True,
                                    skip_variables = True,
                                    skip_constants = True,
                                    skip_chains = False )
            grammar.checkContent()
        else:
            grammar.copyParameters( other.mGrammar,
                                    remove_missing = True)


    model = XGram.Model.Model()
    model.setGrammar( grammar )
    model.setAlphabet( alphabet )
    return model 

   