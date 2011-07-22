import sys,os 
from math import exp, log
from parseSexpr import str2sexpr
from SExpr import *
from extra_functions import float2str, get_cmd_args
allReqsFound = True
try:
    from weblogolib import *
except:
    sys.stderr.write("Warning: could not import WebLogo module.  This is needed to make sequence logo plots.  Download it from: http://weblogo.googlecode.com/ \n")
    allReqsFound = False
try:
    from rpy import *
except:
    sys.stderr.write("Warning: RPy could not be imported. \n")
    allReqsFound = False
try:
    from corebio.seq import *
except:
    sys.stderr.write("Warning: CoreBio could not be imported. \n")
    allReqsFound = False

grayColors = r.gray_colors(n=101, start=1, end=0)
def prob2color(prob, logTrans=False):
    if prob>=1:
        return grayColors[-1]
    elif logTrans:
        print prob
        if prob<0.0:
            return grayColors[ -int(100*prob / log(minPostProb)) ]
        else:
            return grayColors[ -1 ]
    else:
        return grayColors[int(prob*100)]
def getfloat(val):
    retval = None
    try:
        retval = float(val)
        return retval
    except:
        sys.stderr.write("Could not convert this to a float: " + str(val) +'\n')
        return 1e-10
        

def write_PNG(filename, weights, alphabet, type="weblogo"):
    if type == "R":
        r.png(fileName)        
        r.par(mar=[0, 0, 0, 0])
        plotDist = []
        plotChars =[]
        plotCols = []
        otherChars = ''
        otherVals = 0
        for i,val in enumerate(weights):
            if val>.01:
                plotChars.append(alphabet[i].upper()) # eventually we need to store, import the actual alphabet!!
                plotDist.append(val)
                plotCols.append(colors[i])
            else:
                otherChars += alphabet[i]
                otherVals += val
        otherChars = "other"
        plotDist.append(otherVals)
        plotChars.append(otherChars)
#         r.library("seqLogo")
#         pwm = r.makePWM(r.matrix(absDist, ncol=1,nrow=4, byrow=True))
#         r.seqLogo(pwm, ic_scale=False, yaxis=False, xaxis=False)
        r.pie(plotDist, plotChars,col=plotCols+['black'])
        r.dev_off()
    elif type=="weblogo":
        # alphabet=''.join(alphabet).upper()
        data = LogoData(counts=[weights],alphabet=unambiguous_protein_alphabet , entropy=[.5], weight=[1], entropy_interval=[[.4,.6]],length=1)
        options = LogoOptions()
        
        options.show_fineprint=False
        options.alphabet=unambiguous_protein_alphabet#''.join(alphabet).upper()
        options.show_errorbars=False
        options.creator_text=('',)
        options.logo_margin=0.01
        options.unit_name="probability"
        options.yaxis_label=''
        options.show_yaxis=True
        options.show_xaxis=False
        options.colorscheme=colorscheme.chemistry
        format = LogoFormat(data,options); 
        fout = open(fileName, 'w'); 
        png_formatter( data, format, fout); 
        fout.close()



# START
sys.stderr.write("Parsing profile SExpr...\n")
profile = SExpr(str2sexpr(''.join(sys.stdin.readlines()))[0])
profile.make_all_SExprs()

node = profile.get_value_for_tag("node")
incoming_alphabet = list("arndcqeghilkmfpstwyv".upper()) # eventually get this from the profile file
alphabet = str(unambiguous_protein_alphabet)
try:
    viterbi_path = profile.find_all_with_tag("viterbi_path")[0]
except:
    sys.stderr.write("Viterbi path not found...continuing\n")
    viterbi_path = []


outString = "digraph profile_node_"+node+"{\nedge[arrowsize=0.5];\n rankdir=LR\n"
colors = r.rainbow(20);     
sys.stderr.write("Creating PNGs for each state in profile...\n")
try:
    maxPostProb = max( [ getfloat(state.get_value_for_tag("postprob")) for state in profile.find_all_with_tag("state") \
                             if not state.get_value_for_tag("type") in ['start', 'end','wait']])
    minPostProb = min( [ getfloat(state.get_value_for_tag("postprob")) for state in profile.find_all_with_tag("state") \
                         if not state.get_value_for_tag("type") in ['start', 'end','wait']])
except:
    sys.stderr.write("Postprobs not found...continuing\n")
    maxPostProb = minPostProb = 1



statecount = 0
for state in profile.find_all_with_tag("state"):
    sys.stderr.write("Writing state number: %s\n"%statecount)
    statecount += 1
    stateName = state.get_value_for_tag("name")
    if stateName in viterbi_path:
        viterbi_state = True
    else:
        viterbi_state = False
    stateType = state.get_value_for_tag("type")
    if not stateType in ['start', 'end','wait']:
        postProb = getfloat(state.get_value_for_tag("postprob"))/maxPostProb
        #postProb = max(postProb, 0.01)
        absDist = [0.0 for i in alphabet]
        for i, val in enumerate(state.find_all_with_tag("absorb")[0][1:]):
            try:
                newCoord = alphabet.index(incoming_alphabet[i%20]) # eek!
            except:
                sys.stderr.write("Couldn't find: " + str(incoming_alphabet[i][0]) +'\n')
            absDist[newCoord] += exp(getfloat(val)) 
        all = sum(absDist)
        absDist = [i/all for i in absDist]
        # Plot it in something (R or webLogo)
        fileName = "absorb_%s.png"%stateName
        outString += 'Node_%s[image="%s"][imagescale=true][label=""]'%(stateName, fileName)
        #outString += '[width=%s][height=%s][fixedsize=true]'%(postProb, postProb)
        outString += '[fillcolor="%s"][style=filled]'%(prob2color(postProb, logTrans=False))
        if postProb < .01:
            sys.stderr.write("Making this state a 'point', for small postprob\n")
            outString += '[shape=point]'
        if viterbi_state:
            outString += '[color="red"]'
        outString += ';\n'
        write_PNG(fileName, absDist, alphabet, type="weblogo")

    else:
        if stateName != "-1":
            outString += 'Node_%s[label="%s"][color="red"][fillcolor="black"][fontcolor="white"][style=filled];\n'%(stateName, stateType.capitalize())
        else:
            outString += 'Start[color="red"][fillcolor="black"][fontcolor="white"][style=filled];\n'
    
statecount = 0
for transition in profile.find_all_with_tag("transition"):
    sys.stderr.write("Writing transition number: %s\n"%statecount)
    statecount += 1
    fromState = "Node_" + transition.get_value_for_tag("from")
    if fromState == "Node_-1":
        fromState = "Start"
    toState = "Node_" + transition.get_value_for_tag("to")
    weight = exp(getfloat(transition.get_value_for_tag("weight")))
    includeWeight = True
    if not includeWeight:
        weight =''
    else:
        weight = '%s' %getfloat('%.1g' % weight)
    outString += '%s -> %s [label="%s"]'%(fromState, toState, weight)
    if toState.replace("Node_","") in viterbi_path:
        if fromState == "Start" or viterbi_path[viterbi_path.index(toState.replace("Node_",""))-1] == fromState.replace("Node_",""):
            outString += '[color="red"][arrowsize=2]'
    outString += ';\n'



dotFileName = 'test.dot'
outputFileName = sys.argv[1]
DOTcmd = "dot -Tpdf %s > %s"%(dotFileName, outputFileName)

fh=open(dotFileName,'w') 
fh.write(outString + "}\n")
fh.close()

sys.stderr.write("DOT file written to %s, now creating pdf with the following command: \n%s\n"%(dotFileName, DOTcmd))
os.system(DOTcmd)





