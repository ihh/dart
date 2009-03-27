import sys, time, string, re, optparse

USAGE="""module for calculating PAML-like ka/ks values with Xrate.
"""

from XGram.Generator.Prebuilt import Codons
from XGram.Model import Annotation
import XGram.Run
import Bio.Data.CodonTable

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def calculateKaKs( model, options,
                   per_site = False, 
                   terminals = None,
                   prefix = "",
                   normalize_matrix = False ):
    """calculate ks, ka and kaks from a trained grammar.
    
    Iterate over the transition matrix and collect synonymous and non-synonymous
    substitutions.
    """
    
    ## retrieve the terminal frequencies for all codons
    pi = model.evaluateTerminalFrequencies( terminals = terminals )
    
    ## retriev the rate matrix for all codons
    matrix = model.evaluateRateMatrix( terminals = terminals )

    if options.loglevel >= 4:
        options.stdlog.write("# rate matrix\n" )
        for row, col in matrix.items():
            options.stdlog.write( "# " + "".join(row), map(lambda x: "%6.4f" % x, col.values()) + "\n" )
        options.stdlog.flush()
        
    ## translate pi and the matrix to codons
    for key, value in pi.items():
        del pi[key]
        pi["".join(key).upper()] = value

    for key, value in matrix.items():
        for kkey, vvalue in value.items():
            del value[kkey]
            value["".join(kkey).upper()] = vvalue
        del matrix[key]
        matrix["".join(key).upper()] = value

    ## calculate ks
    ## ks is the sum_{i=1}^61 sum_{i!=j,aa_i=aa_j} pi_i * Q_{ij}
    codon_table = Bio.Data.CodonTable.standard_dna_table

    codons = codon_table.forward_table.keys()

    if normalize_matrix:
        ## normalize matrix
        t = 0
        for codon in codons:
            t += pi[codon] * matrix[codon][codon]
        
        t = abs(t)
        
        for codon_i in codons:
            for codon_j in codons:
                matrix[codon_i][codon_j] /= t
    # note: the matrix is symmetric, so you could
    # only do one half.
    rs, ra = 0.0, 0.0
    for codon_i in codons:
        for codon_j in codons:
            if codon_i == codon_j: continue
            if codon_table.forward_table[codon_i] == codon_table.forward_table[codon_j]:
                rs += matrix[codon_i][codon_j] * pi[codon_i]
            else:
                ra += matrix[codon_i][codon_j] * pi[codon_i]

    if not per_site:
        ## return ka/ks per codon        
        return ra, rs
    else:
        ## return ka/ks per site
        # calculate w from the model
        try:
            K0 = float(model.mGrammar.getParameter('%sK0' % prefix))
            K1 = float(model.mGrammar.getParameter('%sK1' % prefix))           
            K2 = float(model.mGrammar.getParameter('%sK2' % prefix))
            K3 = float(model.mGrammar.getParameter('%sK3' % prefix))
    
            w = (K3 + K2) / (K1 + K0)
            rs0 = rs
            ra0 = ra / w
            
            t2 = rs0 + ra0
            rs0 /= t2 
            ra0 /= t2
        
            t1 = rs + ra
            rs /= t1
            ra /= t1
            
            ks = rs / (3 * rs0)
            ka = ra / (3 * ra0)
            
            if options.loglevel >= 2:
                options.stdlog.write( "# ka=%f, ks=%f, ra=%f, rs=%f, ra0=%f, rs0=%f, w=%f, k0=%f, k1=%f, k2=%f, k3=%f, k2/k0=%f, (k3+k2)/(k1+k0)=%f\n" %\
                                    (ka, ks, ra, rs,ra0, rs0, 
                                     w, K0, K1, K2, K3,
                                     K2/K0, (K3+K2/K1+K0)))
            
            return ka, ks
        
        except KeyError:
            pass
        
        K0 = float(model.mGrammar.getParameter('%sK0' % prefix))
        K1 = float(model.mGrammar.getParameter('%sK1' % prefix))           

        w = K1
        rs0 = rs
        ra0 = ra / w
        
        t2 = rs0 + ra0
        rs0 /= t2 
        ra0 /= t2
    
        t1 = rs + ra
        rs /= t1
        ra /= t1
        
        ks = rs / (3 * rs0)
        ka = ra / (3 * ra0)
        
        if options.loglevel >= 2:
            options.stdlog.write( "# ka=%f, ks=%f, ra=%f, rs=%f, ra0=%f, rs0=%f, w=%f, k0=%f, k1=%f, t1=%f\n" %\
                                (ka, ks, ra, rs, ra0, rs0, w, K0, K1, t1 ))
        
        return ka, ks
        
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def getSequencesFromStk( filename):
    """quick parse through stk file.
    
    Only gets sequences."""
    sequences = {}
    for line in open(filename, "r"):
        if line[0] == '#': continue
        if line[:2] == '//': break
        id, s = re.split('\s+', line[:-1])
        if id not in sequences:
            sequences[id] = []
        sequences[id].append(s)
    
    for id, s in sequences.items():
        sequences[id] = "".join(s)
        
    return sequences

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def setFrequencies( model, data):
    """set frequencies in a model according to those observed in data.
    """
    
    sequences = getSequencesFromStk( data )
    frequencies = Codons.getFrequenciesPerCodonPosition( sequences.values() )

    ## build a dummy grammar to insert frequencies
    dummy_grammar = XGram.Model.Grammar()
    for x in range(0,3):
        params = []
        for a in ('A', 'C', 'G', 'T'):
            params.append( ("P%i%s" % (x, a), frequencies[x][a]) )
        dummy_grammar.addVariable( params )
    
    model.mGrammar.copyParameters( dummy_grammar, 
                                   ignore_missing = True)

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def setFixedRates( model, rates):
    """set rates in a model to specified values and fix them."""

    for x in range(len(rates)):
        k = "K%i" % x
        model.mGrammar.setParameter( k, rates[x] )
        model.mGrammar.moveVariableToConst( k )
    

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def test_codon_models( codon_models, 
                       fix_frequencies, 
                       insert_frequencies,
                       data,
                       xgram):
    """module to compare codon models.
    """
    
    t0 = time.time()

    for codon_model in codon_models:
        print "###############################################################"
        print "executing dart/codeml with codon model %s" % codon_model

        generated_model = Codons.buildCodonML(codon_model = codon_model, 
                                              fix_frequencies = fix_frequencies )

        if insert_frequencies:
            setFrequencies( generated_model, data)

        t1 = time.time()

        trained_model = xgram.train( generated_model, data )
        
        t2 = time.time()

        vka, vks = calculateKaKs( trained_model.getModel(), 
                                  per_site = True )        
        
        print "per site - implicit:\tks= %f, ka = %f, ka/ks=%f, n=%i, L=%f" %\
             (vks, vka, vka/vks, 
              trained_model.getNumIterations(),
              trained_model.getLogLikelihood() )
        # print trained_model.getIterations()
#        print "###################################"
#        print trained_model.getLog()
#        print "###################################"
                
        vka, vks = calculateKaKs( trained_model.getModel() )

        print "per codon - w variable:\tks= %f, ka = %f, ka/ks=%f, n=%i, L=%f" %\
             (vks, vka, vka/vks, 
              trained_model.getNumIterations(),
              trained_model.getLogLikelihood() )
        
        
        model_fixed = Codons.buildCodonML( codon_model = codon_model,
                                           fix_frequencies = fix_frequencies,
                                           fix_omega = True)  
    
        if insert_frequencies:
            setFrequencies( model_fixed, data )
        
        t3 = time.time()
        
        trained_model_fixed = xgram.train( model_fixed, data )
        
        t4 = time.time()
        
        fka, fks = calculateKaKs( trained_model_fixed.getModel() )    
        
        print "per codon - w fixed:\tks= %f, ka = %f, ka/ks=%f, n=%i, L=%f" % \
            (fks, fka, fka/fks, 
             trained_model_fixed.getNumIterations(),
             trained_model_fixed.getLogLikelihood() )
            
        # normalize vks, vka, fks and fka
        fks /= fks + fka
        fka /= fks + fka        
        vks /= vks + vka
        vka /= vks + vka
        ks = vks / ( 3 * fks ) 
        ka = vka / ( 3 * fka )
        
        print "per site - explicit:\tks= %f, ka = %f, ka/ks=%f" % (ks, ka, ka/ks)
        print "execution time: %i %i %i" % ( t2 - t1,  t4-t3, t4 - t1 )
    
    print 'dart/codeml finished in %i seconds.' % (time.time() - t0)
    
def test_blocks( codon_models, 
                 fix_frequencies, 
                 insert_frequencies,
                 data,
                 xgram):
    """module to best codon models: linear versus blocks.
    """
    t0 = time.time()

    num_blocks = 2

    ## build annotation
    labels = string.letters.upper()
    annotate_terminals = {}
    for x in range(num_blocks):

        annotations = []
        key = []
        
        for c in range( 0,3 ):
            t = "B%i_POS%i" % (x, c)
            key.append(t)
            annotations.append( Annotation( row = "STATE",
                                            column = t,
                                            label = labels[x % len(labels)] ))
            
        key = tuple(key)
        annotate_terminals[ key ] = annotations

    for codon_model in codon_models:
        print "###############################################################"
        print "executing dart/codeml with codon model %s" % codon_model
        
        t1 = time.time()

        
        model_free = Codons.buildCodonML(codon_model = codon_model, 
                                          grammar_type =  "linear-blocks",
                                          annotate_terminals = annotate_terminals,
                                          num_blocks = num_blocks,
                                          fix_frequencies = fix_frequencies )
        
        if insert_frequencies:
            setFrequencies( model_fixed, data )

        trained_model_free = xgram.train( model_free, data )

        print trained_model_free.getData()

        t2 = time.time()
        
        print "calculated free model in %i seconds: Likelihood = %f" % (t2 - t1, trained_model_free.getLogLikelihood())

        t1 = time.time()
        
        model_fixed = Codons.buildCodonML( codon_model = codon_model,
                                           grammar_type =  "linear-blocks",
                                           num_blocks = 2,
                                           fix_frequencies = fix_frequencies,
                                           fix_omega = True)  
    
        if insert_frequencies:
            setFrequencies( model_fixed, data )

        trained_model_fixed = xgram.train( model_fixed, data )
                
        t2 = time.time()        

        print "calculated fixed model in %i seconds: Likelihood = %f" % (t2 - t1, trained_model_fixed.getLogLikelihood())
        
        t1 = time.time()
        
        model_uniform = Codons.buildCodonML( codon_model = codon_model,
                                             grammar_type =  "linear",
                                             fix_frequencies = fix_frequencies )
    
        if insert_frequencies:
            setFrequencies( model_uniform, data )

        trained_model_uniform = xgram.train( model_uniform, data )
                
        t2 = time.time()        

        print "calculated uniform model in %i seconds: Likelihood = %f" % (t2 - t1, trained_model_uniform.getLogLikelihood())
        
        
        for id in range(num_blocks):
                
            print "## kaks for block %i" % id
            
            prefix = "B%i_" % id
            terminals = ("%sPOS0" % prefix, 
                         "%sPOS1" % prefix,
                         "%sPOS2" % prefix)
            
            vka, vks = calculateKaKs( trained_model_fixed.getModel(), 
                                      terminals = terminals,
                                      prefix = prefix,
                                      per_site = True )        
            
            print "per site - implicit:\tks= %f, ka = %f, ka/ks=%f, n=%i, L=%f" %\
                 (vks, vka, vka/vks, 
                  trained_model_free.getNumIterations(),
                  trained_model_free.getLogLikelihood() )
                 
            # print trained_model.getIterations()
    #        print "###################################"
    #        print trained_model.getLog()
    #        print "###################################"
                    
            vka, vks = calculateKaKs( trained_model_free.getModel(),
                                      terminals = terminals,
                                      prefix = prefix )
    
            print "per codon - w variable:\tks= %f, ka = %f, ka/ks=%f, n=%i, L=%f" %\
                 (vks, vka, vka/vks, 
                  trained_model_free.getNumIterations(),
                  trained_model_free.getLogLikelihood() )
            
            fka, fks = calculateKaKs( trained_model_fixed.getModel(),
                                      terminals = terminals,
                                      prefix= prefix )    
            
            print "per codon - w fixed:\tks= %f, ka = %f, ka/ks=%f, n=%i, L=%f" % \
                (fks, fka, fka/fks, 
                 trained_model_fixed.getNumIterations(),
                 trained_model_fixed.getLogLikelihood() )
                
            # normalize vks, vka, fks and fka
            fks /= fks + fka
            fka /= fks + fka        
            vks /= vks + vka
            vka /= vks + vka
            ks = vks / ( 3 * fks ) 
            ka = vka / ( 3 * fka )
            
            print "per site - explicit:\tks= %f, ka = %f, ka/ks=%f" % (ks, ka, ka/ks)
            
    print 'test finished in %i seconds.' % (time.time() - t0)

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def writeModel( grammar, section, options):
    """write a model to output file."""
    
    if section in options.write or "all" in options.write:
        outfile = open( options.output_pattern % section, "w" )
        outfile.write( "%s\n" % grammar.getGrammar())
        outfile.close()

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def calculateRates( fixed, variable, options ):
    """calculate dN and dS rates using two trained models
    with fixed and variable dN/dS ratio, respectively.
    """
    ## site-wise ka and ks
    ska, sks = calculateKaKs( variable.getModel(),
                              options, 
                              per_site = True )        

    ## codon-wise ka and ks
    vka, vks = calculateKaKs( variable.getModel(), options )
    fka, fks = calculateKaKs( fixed.getModel(), options )
                
    # normalize vks, vka, fks and fka
    n_fks = fks / (fks + fka)
    n_fka = fka / (fks + fka)
    n_vks = vks / (vks + vka)
    n_vka = vka / (vks + vka)
    ks = n_vks / ( 3 * n_fks ) 
    ka = n_vka / ( 3 * n_fka )
  
    if options.loglevel >= 1:
        options.stdlog.write( "# unit\tka\tks\tka\tks\n" )
        options.stdlog.write( "# codon\t%f\t%f\t%f\t%f\n" % (vka, vks, fka, fks ) )
        options.stdlog.write( "# site\t%f\t%f\t%f\t%f\n" % (ska, sks, ka, ks ) )        
  
    return ks, ka

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def run( codon_model, 
         data,
         xgram,
         options ):
    """run xgram on data.
    """
    
    t0 = time.time()

    if options.loglevel >= 1:
        options.stdlog.write( "# executing dart/codeml with codon model %s\n" % codon_model )

    model_variable = Codons.buildCodonML(codon_model = codon_model, 
                                         fix_frequencies = options.fix_frequencies )

    if options.insert_frequencies:
        setFrequencies( model_variable, data )

    if options.fix_rates:
        setFixedRates( model_variable, options.fix_rates ) 

    writeModel( model_variable, "input_variable", options )

    t1 = time.time()

    trained_variable = xgram.train( model_variable, data )

    writeModel( trained_variable.getModel(), "trained_variable", options )

    t2 = time.time()

    model_fixed = Codons.buildCodonML( codon_model = codon_model,
                                       fix_frequencies = options.fix_frequencies,
                                       fix_omega = True)  
    
    if options.insert_frequencies:
        setFrequencies( model_fixed, data )

    writeModel( model_fixed, "input_fixed", options )
        
    t3 = time.time()
        
    trained_fixed = xgram.train( model_fixed, data )

    writeModel( trained_fixed.getModel(), "trained_fixed", options )
    
    t4 = time.time()
        
    if options.loglevel >= 1:
        options.stdlog.write('# dart/codeml finished in %i seconds.\n' % (t4-t0)) 
        
        if options.loglevel >= 2:
            options.stdlog.write( "# execution times\t%i\t%i\t%i\t%i\t%i\n" % (t4-t0,t1-t0, t2-t1, t3-t2, t4-t3))
        
    return trained_fixed, trained_variable
                    
if __name__ ==  "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: codonml.py,v 1.7 2007/10/08 12:16:37 grepall Exp $", usage = USAGE)

    parser.add_option("-v", "--verbose", dest="loglevel", type="int", 
                      help="loglevel to use."  )

    parser.add_option("-m", "--model", dest="models", type="choice", action="append",
                      choices=("f3x4-two", "f3x4-four"),
                      help="models to use."  )

    parser.add_option("-w", "--write", dest="write", type="choice", action="append",
                      choices=("input_fixed", "trained_fixed", "input_variable", "trained_variable", "all" ),
                      help="output sections to write."  )

    parser.add_option("-o", "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern for output files." )
    
    parser.add_option("-t", "--test", dest="test", type="choice",
                      choices=("models", "blocks" ),
                      help="perform test on models or block/non-block models." )

    parser.add_option("-i", "--insert-frequencies", dest="insert_frequencies", action="store_true",
                      help="insert frequencies from input file." )

    parser.add_option("-f", "--fix-frequencies", dest="fix_frequencies", action="store_true",
                      help="set initial frequencies to const." )

    parser.add_option("-r", "--fix-rates", dest="fix_rates", type="string",
                      help="""fix rates to specified values. Note that the number of rates has to match the ones
in the model. Provide values in a comma-separated list.""")

    parser.set_defaults(
                        loglevel = 1,
                        models = [],
                        test = None,
                        insert_frequencies = False,
                        fix_frequencies = False,
                        write = [],
                        output_pattern = "%s.eg",
                        stdout = sys.stdout,
                        stdlog = sys.stdout,
                        value_format = "%6.4f",
                        fix_rates = None,
                        )

    (options, args) = parser.parse_args()
    if options.fix_rates: options.fix_rates = map( float, options.fix_rates.split(",") )
    
    xgram = XGram.XGram()

    if len(args) > 0:
        data = args[0]
    else:
        data = XGram.PATH_DATA + "/dpse_dmel.stk"
    
    if options.test:
        
        xgram.setDebug()
        data = XGram.PATH_DATA + "/dpse_dmel.stk"    
        
        if options.test=="models":
            test_codon_models( options.codon_models, 
                               options.fix_frequencies, 
                               options.insert_frequencies,
                               data,
                               xgram)
        elif do_test_blocks:
            test_blocks( options.codon_models,
                         options.fix_frequencies,
                         options.insert_frequencies,
                         data, xgram )
    
    else:
        xgram.setDebug()
        for model in options.models:
            fixed_model, variable_model = run ( model, data, xgram, options )
            
        ## write log likelihoold
        options.stdout.write('model\tniterations\tloglikelihood\n')
        options.stdout.write('trained\t%i\t%s\n' % (variable_model.getNumIterations(),
                                                    options.value_format % variable_model.getLogLikelihood()))
        options.stdout.write('fixed\t%i\t%s\n' % (fixed_model.getNumIterations(),
                                                  options.value_format % fixed_model.getLogLikelihood()))
            
        ka, ks = calculateRates( fixed_model, variable_model, options )

    