"""This module provides codon based models.
"""

import sys, re, string
from types import *
import Bio.Alphabet.IUPAC

from XGram.Generator import Transitions, InitialStates, \
    InitialStatesFromTransitions, Chain, Rates, RatesDictionary, \
    GrammarLinearSequence, GrammarLinearSequenceBlocks, \
    AlphabetDNA
    
import XGram.Model

import Bio
import Bio.Alphabet
import Bio.Alphabet.IUPAC
import Bio.Data.CodonTable

class TransitionsCodons( Transitions ):
    """transitions for codon based model."""
    
    def __init__(self, codon_table = None, *args, **kwargs):
        self.mCodonTable = codon_table
        
        self.mDefaultRates = {}
        Transitions.__init__( self, 
                              rate_generator = RatesDictionary( self.mDefaultRates ),
                              *args, 
                              **kwargs )

    def registerTerminals(self, chain):
        terminals = []
        for x in range(0,3):
            terminals.append( "%sCOD%i" % (self.mPrefix, x))
        chain.setTerminals( terminals )

class ChainCodons( Chain ):
    """generator for codon models."""
    def __init__( self,
                  codon_table = Bio.Data.CodonTable.standard_dna_table, 
                  generator_initial_states = None,
                  generator_transitions = None ):

        if not generator_transitions:
            generator_transitions = TransitionsCodons( codon_table = codon_table )
        
        self.mCodonTable = codon_table
        # the alphabet is given by the codon table, but the states
        # are separated.
        alphabet = map( lambda x: tuple(x.lower()), codon_table.forward_table.keys())
        Chain.__init__(self, 
                       alphabet = alphabet,
                       generator_initial_states = generator_initial_states,
                       generator_transitions = generator_transitions, 
                       )
    
        self.mName = "codons"

class CodonFrequencyGenerator:
    """generate nucleotide frequencies.
    
    The default implementation generates uniform constant codon frequencies.
    """
    
    def __init__(self, 
                 codon_table, 
                 prefix = "",
                 fix_frequencies = False ):
        self.mCodonTable = codon_table
        self.mPrefix = prefix
        self.mFixFrequencies = fix_frequencies
        self.mDefaultRates = {}
        self.mRateGenerator = RatesDictionary( self.mDefaultRates )

    def getInitialProbability( self, emissions):
        """return initial probability uniform frequencies
        """        
        return 1.0 / 61.0
        
    def getPj(self, codon ):
        """return pj for uniform frequencies
        """
        return 1.0 / 61.0  
        
    def registerParameters(self, chain ):
        """register parameters for the chain to built."""
        pass

class CodonFrequencyGeneratorF3X4( CodonFrequencyGenerator ):
    """generate nucleotide frequencies according to f3x4 model
    
        The normalization factor (nfactor) is necessary to convert from the full
        64x64 matrix that xrate uses to the 61x61 matrix that Goldman&Yang employ.
    
        The nfactor is 64.0/61.0 for uniform codon frequencies and in general is
        the factor necessary to make codon frequencies scale to one. The default
        is 64.0/61.0.
    """
    def __init__(self, **kwargs):
        CodonFrequencyGenerator.__init__(self, **kwargs)
        self.mNormalizationFactor = "%snfactor * "        

        ## set default uniform rates for f3x4 model
        na = self.mCodonTable.nucleotide_alphabet.letters
                
        for x in range(0,3):
            for n in na:
                self.mDefaultRates["%sp%s%i" % (self.mPrefix, n.lower(), x)] = 0.25

    def getInitialProbability( self, emissions):
        """return initial probability for f3x4 model.
        """
        
        f3x4 = self.mNormalizationFactor % self.mPrefix +\
            " * ".join ([ "%sp%s%i" % (self.mPrefix, emissions[x].lower(),x) for x in range(0,3)] )                    
        return f3x4
        
    def getPj(self, codon ):
        """return pj for f3x4 model.
        """
        
        f3x4 = self.mNormalizationFactor % self.mPrefix +\
            " * ".join([ "%sp%s%i" % (self.mPrefix, codon[x].lower(), x) for x in range(0,3)] )                    
        return f3x4        
        
    def registerParameters(self, chain ):
        """register parameters for the chain to built."""

        na = self.mCodonTable.nucleotide_alphabet.letters.lower()
                            
        for x in range(0,3):
            probs = []
            for n in na:
                rate = "%sp%s%i" % (self.mPrefix, n, x)
                probs.append( (rate, float(self.mRateGenerator.getRate(rate)) ) )
            if self.mFixFrequencies:
                chain.addParameter( probs, 
                                    is_const = True,
                                    is_explicit = True,
                                    normalize = True )
            else:
                chain.addParameter( probs,
                                    is_const = False,
                                    is_explicit = True, 
                                    normalize = True )
        
        nfactor = "%snfactor" % (self.mPrefix)
        chain.addParameter( (nfactor, 64.0/61.0 ),
                            is_explicit = True,
                            is_const = True,
                            is_rate = True )


class CodonFrequencyGeneratorCodons( CodonFrequencyGenerator ):
    """Nucleotide frequencies for codon model a la Goldman & Yang (1994).

    The nucleotide frequencies have to be given by the parameter codons. If
    not given, codon frequencies are set to uniform (1/61).
    """
    
    def __init__(self, codon_frequencies = None, **kwargs):
        CodonFrequencyGenerator.__init__(self, **kwargs)
        self.mCodonFrequencies = codon_frequencies
        
        if self.mCodonFrequencies:
            for codon in self.mCodonTable.forward_table:
                if codon.upper() not in self.mCodonFrequencies:
                    raise "missing codon %s in codon frequencies" % (codon)
        
    def getInitialProbability( self, emissions):
        """return initial probability for f3x4 model.
        """                
        if self.mCodonFrequencies:
            return self.mCodonFrequencies["".join(emissions).upper()]
        else:
            return 1/61.0

    def getPj(self, codon ):
        """return pj."""        
        if self.mCodonFrequencies:
            return self.mCodonFrequencies["".join(codon).upper()]
        else:
            return 1/61.0
         
    def registerParameters(self, chain ):
        """register parameters for the chain to built."""
        pass
        
class TransitionsCodonML( TransitionsCodons ):
    """transitions for Goldman & Yang (1994) like models."""
    def __init__(self,
                 codon_frequency_generator,
                 fix_omega = False, 
                 fix_kappa = False,
                 fix_frequencies = False,
                 shared_frequencies = False,
                 shared_rates = False,
                 **kwargs):
        TransitionsCodons.__init__( self, **kwargs )
        self.mFixOmega = fix_omega
        self.mFixKappa = fix_kappa
        self.mFixFrequencies = fix_frequencies
        self.mSharedRates = shared_rates
        self.mCodonFrequencyGenerator = codon_frequency_generator 
        
        if self.mSharedRates:
            self.mRatePrefix =""
        else:
            self.mRatePrefix = self.mPrefix
                
    def getInitialProbability(self, emissions):
        return self.mCodonFrequencyGenerator.getInitialProbability( emissions )
    
    def buildGrammar(self, chain):
        """five different rates:
            0 if more than two nucleotide changes between codons. 
            transitions/transversions x synonymous/non-synonymous.
            Each rate is multiplied by pi_j, the equilibrium frequency.
            
            Note that the alphabet is given as tuples of emissions (lower case),
            while the codon table uses upper case characters.
            """
        TransitionsCodons.buildGrammar( self, chain)
        
        chain.setUpdatePolicy( "parametric" )

    def registerParameters(self, chain):
        """register parameters."""
        self.mCodonFrequencyGenerator.registerParameters( chain )

class TransitionsCodonMLConst(TransitionsCodonML):
    """transitions for Goldman & Yang (1994) like models."""
    def __init__(self, 
                 **kwargs):
        TransitionsCodonML.__init__( self, **kwargs )

    def getInitialProbability( self, emissions):
        return 1.0/61.0
        
    def registerParameters(self, chain ):
        """register parameters for the chain to built."""

        if self.mFixOmega:
            num_rates = 2
        else:
            num_rates = 4
            
        for x in range( 0, num_rates):
            rate = "%sK%i" % (self.mRatePrefix, x)
            chain.addVariable( (rate, self.mRateGenerator.getRate(rate) ))

    def addTransitions(self, chain, codon1, codon2):
        """decides, if a transition is to be added and adds it.
        Both directions will be added.
        """
        changes = []
        for x in range(0,3):
            if codon1[x] != codon2[x]:
                changes.append( (x, codon1[x], codon2[x]) )
        
        aa1 = self.mCodonTable.forward_table[("".join(codon1)).upper()]
        aa2 = self.mCodonTable.forward_table[("".join(codon2)).upper()]
         
        # uniform probability:
        p = 0.004
        
        if len(changes) == 1:
            change = ( changes[0][1], changes[0][2] )
            is_transition = change in (("a","g"), ("g","a"), ("t","c"), ("c","t"))
            if self.mFixOmega:
                is_synonymous = True
            else:
                is_synonymous = aa1 == aa2
                
            if is_synonymous:            
                if is_transition: 
                    rate = "%sK0" % self.mRatePrefix
                else:
                    rate = "%sK1" % self.mRatePrefix
            else:
                if is_transition:
                    rate = "%sK2" % self.mRatePrefix
                else:
                    rate = "%sK3" % self.mRatePrefix
            
            chain.addTransition( codon1, codon2, 
                                 ("rate", "%s * %f" %\
                                  (rate, p ) ) )
            chain.addTransition( codon2, codon1, 
                                 ("rate", "%s * %f" %\
                                  (rate, p ) ))

class TransitionsCodonMLTwoParameters( TransitionsCodonML ):
    """Transitions for F3X4 model a la Goldman & Yang (1994).
    
    This model is parameterized by kappa (Transition/Transversion ratio)
    and omega (ka/ks-ratio)

    | rate    |  meaning                         | in PAML |
    | Rsi     |  synonymous & transition         | kappa   |
    | Rsv     |  synonymous & transversion       | 1       |
    | Rni     |  non-synonymous & transition     | kappa * omega   |
    | Rnv     |  non-synonymous & transversion   | omega       |
    """
    
    def __init__(self, **kwargs):
        TransitionsCodonML.__init__(self, **kwargs)

    def addTransitions(self, chain, codon1, codon2):
        """decides, if a transition is to be added and adds it.
        Both directions will be added.
        """
        changes = []
        for x in range(0,3):
            if codon1[x] != codon2[x]:
                changes.append( (x, codon1[x], codon2[x]) )
        
        aa1 = self.mCodonTable.forward_table[("".join(codon1)).upper()]
        aa2 = self.mCodonTable.forward_table[("".join(codon2)).upper()]

        if len(changes) == 1:
            change = ( changes[0][1], changes[0][2] )
            is_transition = change in (("a","g"), ("g","a"), ("t","c"), ("c","t"))
            if self.mFixOmega:
                is_synonymous = True
            else:
                is_synonymous = aa1 == aa2
            
            if is_synonymous:            
                if is_transition: 
                    rate = "%skappa" % self.mRatePrefix
                else:
                    rate = "1" 
            else:
                if is_transition:
                    rate = "%skappa * %somega" % (self.mRatePrefix, self.mRatePrefix)
                else:
                    rate = "%somega" % self.mRatePrefix
            
            chain.addTransition( codon1, codon2, 
                                 ("rate", "%s * %s" %\
                                  (rate, self.mCodonFrequencyGenerator.getPj( codon2 ) )))
            chain.addTransition( codon2, codon1, 
                                 ("rate", "%s * %s" %\
                                  (rate, self.mCodonFrequencyGenerator.getPj( codon1 ) )))
            
    def registerParameters(self, chain ):
        """register parameters for the chain to built."""

        TransitionsCodonML.registerParameters( self, chain )

        if not self.mFixKappa:
            rate = "%skappa" % self.mRatePrefix
            self.mDefaultRates[rate]  = 2.0
            chain.addParameter( (rate, self.mRateGenerator.getRate(rate) ),
                                is_const = False,
                                is_explicit = True)            
        
        if not self.mFixOmega:
            rate = "%somega" % self.mRatePrefix            
            self.mDefaultRates[rate]  = 0.3
            chain.addParameter( (rate, self.mRateGenerator.getRate(rate) ),
                                is_const = False,
                                is_explicit = True)
        
class TransitionsCodonMLFourParameters( TransitionsCodonML ):
    """Transitions for model a la Goldman & Yang (1994).

    Parameterize a GY grammar using up to four parameters for the
    four kinds of transitions
    
    | rate    |  meaning                         | in PAML |
    | Rsi     |  synonymous & transition         | k       |
    | Rsv     |  synonymous & transversion       | 1       |
    | Rni     |  non-synonymous & transition     | k * w   |
    | Rnv     |  non-synonymous & transversion   | w       |

    This grammar is equivalent to akaksgc.eg.
    
    If fix_kappa or fix_omega are set to true, then the number of rates
    is reduced to two (kappa, not_kappa) and (omega, not_omega).
    
    Note that the following holds:
    
    1 = (Rsi * Rnv) / (Rsv * Rni)
    
    """
    
    def __init__(self, **kwargs):
        TransitionsCodonML.__init__(self, **kwargs)

    def addTransitions(self, chain, codon1, codon2):
        """decides, if a transition is to be added and adds it.
        Both directions will be added.
        """
        changes = []
        for x in range(0,3):
            if codon1[x] != codon2[x]:
                changes.append( (x, codon1[x], codon2[x]) )
        
        aa1 = self.mCodonTable.forward_table[("".join(codon1)).upper()]
        aa2 = self.mCodonTable.forward_table[("".join(codon2)).upper()]
        
        if self.mFixOmega and self.mFixKappa:
            rsi, rsv, rni, rnv = "r", "r", "r", "r"
        elif self.mFixOmega:
            rsi, rsv, rni, rnv = "kappa", "not_kappa", "kappa", "not_kappa"
        elif self.mFixKappa:
            rsi, rsv, rni, rnv = "omega", "omega", "not_omega", "not_omega"
        else:
            rsi, rsv, rni, rnv = "Rsi", "Rsv", "Rni", "Rnv"
            
        if codon1 == codon2:
            chain.addTransition( codon1, codon2 )
        
        elif len(changes) == 1:
            change = ( changes[0][1], changes[0][2] )
            is_transition = change in (("a","g"), ("g","a"), ("t","c"), ("c","t"))
            is_synonymous = aa1 == aa2
            
            if is_synonymous:            
                if is_transition: 
                    rate = "%s%s" % (self.mRatePrefix, rsi)
                else:
                    rate = "%s%s" % (self.mRatePrefix, rsv)
            else:
                if is_transition:
                    rate = "%s%s" % (self.mRatePrefix, rni)
                else:
                    rate = "%s%s" % (self.mRatePrefix, rnv)

            chain.addTransition( codon1, codon2, 
                                 ("rate", "%s * %s" %\
                                  (rate, self.mCodonFrequencyGenerator.getPj( codon2 ) )))
            chain.addTransition( codon2, codon1, 
                                 ("rate", "%s * %s" %\
                                  (rate, self.mCodonFrequencyGenerator.getPj( codon1 ) )))

    def registerParameters(self, chain):
        """register rates of the model."""

        TransitionsCodonML.registerParameters( self, chain )
        
        ## set default rates by adding them into the dictionary
        ## mDefaultRates that is used by the RateGenerator to
        ## lookup the rates
        if self.mFixOmega and self.mFixKappa:
            rates = ("r" )
            self.mDefaultRates['%sr' % self.mRatePrefix]  = 1.0
        elif self.mFixOmega:
            rates = ("kappa" , "not_kappa")
            ## default rates are set for a kappa of 1.0 and a ds of 0.3
            self.mDefaultRates.update( { "%skappa" % self.mRatePrefix: 0.1, 
                                        "%snot_kappa" % self.mRatePrefix : 0.1 } )
        elif self.mFixKappa:
            rates = ("omega", "not_omega")
            self.mDefaultRates.update( { "%somega" % self.mRatePrefix : 1.0, 
                                        "%snot_omega" % self.mRatePrefix : 1.0 } )
        else:
            rates = ( "Rsi", "Rsv", "Rni", "Rnv" )
            ## default rates are set for omega = 0.1 kappa = 2.0 ds = 0.1
            ## rsi=4.7951, rsv=2.3975, rni=0.4795, rnv=0.2398
            self.mDefaultRates.update( { "%sRsi" % self.mRatePrefix: 4.80, 
                                         "%sRsv" % self.mRatePrefix: 2.4, 
                                         "%sRni" % self.mRatePrefix: 0.48, 
                                         "%sRnv" % self.mRatePrefix: 0.24 } )
            
        for r in rates:
            rate = "%s%s" % (self.mRatePrefix, r)
            chain.addParameter( (rate, self.mRateGenerator.getRate(rate) ),
                                is_explicit = True,
                                is_rate = True )

class TransitionsCodonMLFourProducts( TransitionsCodonML ):
    """Transitions for model a la Goldman & Yang (1994).

    Parameterize a GY grammar using up to four parameters for the
    four kinds of transitions
    
    | rate    |  meaning                         | in PAML | Implemented as |
    | Rsi     |  synonymous & transition         | k       | Rs * Ri        |
    | Rsv     |  synonymous & transversion       | 1       | Rs * Rv        |
    | Rni     |  non-synonymous & transition     | k * w   | Rn * Ri        |
    | Rnv     |  non-synonymous & transversion   | w       | Rn * Rv        |

    If fix_kappa or fix_omega are set to true, then the number of rates
    is reduced to two (kappa, not_kappa) and (omega, not_omega).
    
    
    
    """
    
    def __init__(self, **kwargs):
        TransitionsCodonML.__init__(self, **kwargs)

    def addTransitions(self, chain, codon1, codon2):
        """decides, if a transition is to be added and adds it.
        Both directions will be added.
        """
        changes = []
        for x in range(0,3):
            if codon1[x] != codon2[x]:
                changes.append( (x, codon1[x], codon2[x]) )
        
        aa1 = self.mCodonTable.forward_table[("".join(codon1)).upper()]
        aa2 = self.mCodonTable.forward_table[("".join(codon2)).upper()]
        
        if self.mFixOmega and self.mFixKappa:
            rsi, rsv, rni, rnv = "r", "r", "r", "r"
        elif self.mFixOmega:
            rsi, rsv, rni, rnv = "%(prefix)skappa", "%(prefix)snot_kappa", "%(prefix)skappa", "%(prefix)snot_kappa"
        elif self.mFixKappa:
            rsi, rsv, rni, rnv = "%(prefix)somega", "%(prefix)somega", "%(prefix)snot_omega", "%(prefix)snot_omega"
        else:
            rsi, rsv, rni, rnv = "%(prefix)sRs * %(prefix)sRi", "%(prefix)sRs * %(prefix)sRv", "%(prefix)sRn * %(prefix)sRi", "%(prefix)sRn * %(prefix)sRv"
            
            
        prefix = { 'prefix' : self.mRatePrefix }
                    
        if codon1 == codon2:
            chain.addTransition( codon1, codon2 )
        
        elif len(changes) == 1:
            change = ( changes[0][1], changes[0][2] )
            is_transition = change in (("a","g"), ("g","a"), ("t","c"), ("c","t"))
            is_synonymous = aa1 == aa2
            
            if is_synonymous:            
                if is_transition: 
                    rate = rsi % prefix
                else:
                    rate = rsv % prefix
            else:
                if is_transition:
                    rate = rni % prefix
                else:
                    rate = rnv % prefix

            chain.addTransition( codon1, codon2, 
                                 ("rate", "%s * %s" %\
                                  (rate, self.mCodonFrequencyGenerator.getPj( codon2 ) )))
            chain.addTransition( codon2, codon1, 
                                 ("rate", "%s * %s" %\
                                  (rate, self.mCodonFrequencyGenerator.getPj( codon1 ) )))

    def registerParameters(self, chain):
        """register rates of the model."""

        TransitionsCodonML.registerParameters( self, chain )
        
        ## set default rates by adding them into the dictionary
        ## mDefaultRates that is used by the RateGenerator to
        ## lookup the rates
        if self.mFixOmega and self.mFixKappa:
            rates = ("r" )
            self.mDefaultRates['%sr' % self.mRatePrefix]  = 1.0
        elif self.mFixOmega:
            rates = ("kappa" , "not_kappa")
            ## default rates are set for a kappa of 1.0 and a ds of 0.3
            self.mDefaultRates.update( { "%skappa" % self.mRatePrefix: 0.1, 
                                        "%snot_kappa" % self.mRatePrefix : 0.1 } )
        elif self.mFixKappa:
            rates = ("omega", "not_omega")
            self.mDefaultRates.update( { "%somega" % self.mRatePrefix : 1.0, 
                                        "%snot_omega" % self.mRatePrefix : 1.0 } )
        else:
            rates = ( "Rs", "Rn", "Rv", "Ri" )
            ## default rates are set for omega = 0.1 kappa = 2.0 ds = 0.1
            ## rsi=4.7951, rsv=2.3975, rni=0.4795, rnv=0.2398
            self.mDefaultRates.update( { "%sRi" % self.mRatePrefix: 40.0, 
                                         "%sRv" % self.mRatePrefix: 20.0, 
                                         "%sRs" % self.mRatePrefix: 1.0, 
                                         "%sRn" % self.mRatePrefix: 0.4, } )
            
        for r in rates:
            rate = "%s%s" % (self.mRatePrefix, r)
            chain.addParameter( (rate, self.mRateGenerator.getRate(rate) ),
                                is_explicit = True,
                                is_rate = True )
        
        
class ChainCodonML( ChainCodons ):         
    """generator for Goldman & Yang (1994) like models."""
    def __init__( self, 
                  codon_table = Bio.Data.CodonTable.standard_dna_table,
                  prefix = "",
                  codon_model = "f3x4-four",
                  fix_omega = False,
                  fix_kappa = False,
                  fix_frequencies = False,
                  shared_frequencies = False,
                  shared_rates = False,
                  codon_frequencies = None,
                  **kwargs ):
        
        if codon_model == "const":
            generator_transitions = TransitionsCodonMLConst( codon_table = codon_table,
                                                        fix_omega = fix_omega, 
                                                        fix_frequencies = fix_frequencies,
                                                        shared_frequencies = shared_frequencies,
                                                        shared_rates = shared_rates,                                                        
                                                        prefix = prefix )
        elif codon_model == "sn":
            generator_transitions = TransitionsCodonMLSn( codon_table = codon_table,
                                                          fix_omega = fix_omega,
                                                          fix_frequencies = fix_frequencies,
                                                          shared_frequencies = shared_frequencies,
                                                          shared_rates = shared_rates,                                                          
                                                          prefix = prefix )
        else:
            freqs, params = codon_model.split("-")
            if shared_frequencies:
                frequency_prefix = ""
            else:
                frequency_prefix = prefix

            if freqs == "f3x4":
                frequency_generator = CodonFrequencyGeneratorF3X4( codon_table = codon_table,
                                                                   prefix = frequency_prefix,
                                                                   fix_frequencies = fix_frequencies )
            elif freqs == "uniform":
                frequency_generator = CodonFrequencyGenerator( codon_table = codon_table,
                                                               prefix = frequency_prefix,
                                                               fix_frequencies = fix_frequencies )                
            elif freqs == "codons":
                frequency_generator = CodonFrequencyGeneratorCodons( codon_table = codon_table,
                                                                     prefix = frequency_prefix,
                                                                     fix_frequencies = fix_frequencies,
                                                                     codon_frequencies = codon_frequencies )
            else:
                raise XGram.Exceptions.UsageError ("unknown codon model %s" % codon_model)
        
            if params == "four":
                generator_transitions = TransitionsCodonMLFourParameters( codon_table = codon_table,
                                                                          fix_omega = fix_omega,
                                                                          fix_frequencies = fix_frequencies,
                                                                          shared_rates = shared_rates,
                                                                          codon_frequency_generator = frequency_generator,
                                                                          prefix = prefix )
            elif params == "fourproducts":
                generator_transitions = TransitionsCodonMLFourProducts( codon_table = codon_table,
                                                                        fix_omega = fix_omega,
                                                                        fix_frequencies = fix_frequencies,
                                                                        shared_rates = shared_rates,
                                                                        codon_frequency_generator = frequency_generator,
                                                                        prefix = prefix )
            elif params == "two":
                generator_transitions = TransitionsCodonMLTwoParameters( codon_table = codon_table,
                                                                         fix_omega = fix_omega,
                                                                         fix_frequencies = fix_frequencies,
                                                                         shared_rates = shared_rates,
                                                                         codon_frequency_generator = frequency_generator,
                                                                         prefix = prefix )
            else:
                raise XGram.Exceptions.UsageError ("unknown codon model %s" % codon_model)
            
        generator_initial_states = InitialStatesFromTransitions()
        
        ChainCodons.__init__( self, 
                              codon_table = codon_table,
                              generator_initial_states = generator_initial_states, 
                              generator_transitions = generator_transitions, 
                              **kwargs ) 
       
        self.mName = "codonml_f3x4"



def getFrequenciesPerCodonPosition( sequences):
    """get nucleotide frequencies for each codon position.
    
    returns a tuple of hashes with the frequency values.
    """
        
    frequencies = [ {} for x in range( 0, 3 ) ]
    for l in Bio.Alphabet.IUPAC.unambiguous_dna.letters:                
        for p in range( 0, 3 ):
            frequencies[p][l] = 0
            
    codon_counts = {}
    for x in Bio.Data.CodonTable.standard_dna_table.forward_table.keys():
        codon_counts[x] = 0

    for sequence in sequences:
        sequence = re.sub( "\s+", "", sequence ).upper()

        if len( sequence ) % 3:
            raise XGram.Exceptions.UsageError( "at least one sequence is not a multiple of 3." )

        for codon in [ sequence[x:x+3] for x in range( 0, len( sequence ), 3 )]:
            if "-" in codon: continue
            codon_counts[codon] += 1

    total_codons = sum(codon_counts.values() )

    for codon, count in codon_counts.items():
        for p in range(0,3):
            nuc=codon[p]
            frequencies[p][nuc] += count

    for p in range( 0, 3 ):
        f = frequencies[p]
        t = float( sum( f.values() ) )
        if t == 0:
            raise XGram.Exceptions.UsageError ("no counts for codon position %i" % p)
             
        for nuc in f.keys():
            f[nuc] /= t
    
    return frequencies

def getFrequencies( model, terminals = None ):
    """returns frequencies for codons given the model.
    
    Initial frequencies are given by parametric expressions and can
    thus not be read directly from the grammar but need to be computed
    via the parameters P0G, ....
    
    Note: this method is maybe not the most efficient way to do this and
    not the safest either ;=)
    """
    
    frequencies = {}

    if terminals:
        chain = model.mGrammar.mChains[terminals]
    else:
        if len(model.mGrammar.mChains) > 1:
            raise "more than one chain - please supply terminals."
        chain = model.mGrammar.mChains.values()[0]

    # this works by building a little code block and executing it.
    # All parameter variables are added into the name space and then
    # the expressions are added explicitely from the model.
    result = {}
    
    lines = []
    for parameter, value in model.mGrammar.getParameters().items():
        lines.append( "%s = %s" % (parameter, str(value)))

    for states, expr in chain.getInitialStates().items():
        if type(expr[1]) == StringType:
            expr = expr[1]
        else:            
            expr = " ".join(expr[1])
        
        lines.append( "result[%s] = %s" % (str(states), expr))
        
    exec("\n".join(lines))

    return result

def buildCodonML( codon_model = "f3x4-four", 
                  grammar_type = "linear",
                  num_blocks = 2,
                  explicit_extension = False,
                  fix_omega = False,
                  fix_kappa = False,
                  fix_frequencies = False,
                  shared_frequencies = False,
                  shared_rates = False,
                  annotate_terminals = None,
                  codon_frequencies = None ):
    """build a codonml-like model.
    
    fix_omega: If true, omega is set to 1. Thus the model does not distinguish
        between synonymous and non-synonymous transitions. Use this option
        to estimate the number of synonymous and non-synonymous sites.
        (The original rate parameters are per codon.)
    
    fix_kappa: If true, kappa is set to 1. Thus the model does not estimate
        the transition/transversion ratio.
     
    fix_frequencies: if True: frequencies of the codons are set to const
        and are not estimated. The corresponding rate parameters should be
        added afterwards.
    
    grammar_type can be one of "linear", "linear-blocks".
        linear: simple linear sequence with gaps
        linear-blocks: blocks of linear sequence
    
    codon_model can be one of "f3x4-two", "f3x4-four"
        f3x4-two: two parameter parameterization (a al CodonML)
            Codon frequencies are given by the nucleotide frequencies at each position.
        f3x4-four: four parameter parameterization
            (better convergence properties). 
            Codon frequencies are given by the nucleotide frequencies at each position.
        codons-four: four parameter parametrization with explicitely
            given codon frequencies. These are added as constant
            expressions into the grammar.
            
    num_blocks: the number of blocks for the "linear-blocks" grammar.
            
    shared_frequencies: the nucleotide frequencies are shared between blocks.
        The default (false) has separate nucleotide frequencies for each block.
            
    shared_rates: the rates between blocks are shared. The default (false)
        has separate rates for each block.
            
    annotate_terminals:
        The linear-blocks grammar allows to annotate the input
        alignments. Provide a hash mapping the emitted states to
        Annotations. 
        
    explicit_extension:
        use an explicit parameter for the markov chain.
    """
    
    model = XGram.Model.Model()
    grammar = XGram.Model.Grammar()
    alphabet = XGram.Model.Alphabet()
    
    if grammar_type == "linear":
        chain = XGram.Model.Chain()    
        
        generator_chain = ChainCodonML(fix_omega = fix_omega,
                                       fix_kappa = fix_kappa,
                                       fix_frequencies = fix_frequencies,
                                       codon_frequencies = codon_frequencies,
                                       codon_model = codon_model )
        generator_chain.buildGrammar( chain )
    
        grammar.addChain( chain )
        generator_grammar = GrammarLinearSequence( explicit_extension = explicit_extension )  
        
    elif grammar_type == "linear-blocks":
            
        for x in range( num_blocks ):
            chain = XGram.Model.Chain()    
            generator_chain = ChainCodonML(fix_omega = fix_omega,
                                           fix_kappa = fix_kappa,
                                           fix_frequencies = fix_frequencies,
                                           shared_frequencies = shared_frequencies,
                                           shared_rates = shared_rates,
                                           codon_model = codon_model,
                                           codon_frequencies = codon_frequencies,
                                           prefix = "B%x_" % x )
            
            generator_chain.buildGrammar( chain )
            grammar.addChain( chain )

        generator_grammar = GrammarLinearSequenceBlocks( annotate_terminals = annotate_terminals )
    else:
        raise XGram.Exceptions.UsageError( 'unknown grammar %s' % grammar_type )
        
    generator_grammar.buildGrammar( grammar )     
    grammar.setName( "codonML" )
        
    generator_alphabet = AlphabetDNA()
    generator_alphabet.buildGrammar(alphabet)
     
    model.setAlphabet( alphabet )
    model.setGrammar( grammar )
    
    return model


    