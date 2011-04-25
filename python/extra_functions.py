import sys
from math import log, exp
import copy
# from Bio.Data import CodonTable
# ## BioPython codon table:
# table = CodonTable.standard_dna_table
# table.forward_table['TAG'] = 'stop'
# table.forward_table['TAA'] = 'stop'
# table.forward_table['TGA'] = 'stop'
def float2str(f, n=2):
    s = str(f)
    return s[:s.index('.')+n+1]
def unique(inList):
    out = []
    for item in inList:
        if not item in out:
            out.append(item)
    return out
def isUnique(inList):
    return len(unique(inList)) == 1

def within(a,(b,c),leftInclude=1, rightIncude=1):
    '''returns whether a is within (b,c)'''
    if leftInclude and rightInclude:
        return a>=b and a<=c
    elif leftInclude and not rightInclude:
        return a>=b and a<c
    elif not leftInclude and not rightInclude:
        return a>b and a<c
    elif not leftInclude and rightInclude:
        return a>b and a<=c
def category(point, cutoffs):
    for i,cutoff in enumerate(cutoffs):
        if i == 0 and point <= cutoff:
            return 0
        elif i == len(cutoffs)-1:
            return i
        elif point <= cutoff and point > cutoffs[i-1]:
            return i
def add_logs(x,y):
    "A fast way to add logarithms"
    if not x==0 and not y==0:
        return x+log(1+exp(y-x))
    elif x==0 and y==0:
        sys.stderr.write('problem with log values\n')
        sys.exit(1)
    elif x==0:
        return y
    elif y==0:
        return x

def add_bits(x,y):
    "A fast way to add base2 logarithms"
    if not x==0 and not y==0:
        return x+log(1+2**(y-x)/log(2))
    elif x==0 and y==0:
        sys.stderr.write('problem with log values\n')
        sys.exit(1)
    elif x==0:
        return y
    elif y==0:
        return x

def get_cmd_args():
    import sys
    argDict = {}
    argDict['options'] = []
    i = 1
    last = False
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg.startswith('-'):
            if i == len(sys.argv)-1:
                last = True
                argDict['options'].append(arg[1:])
                break
            elif sys.argv[i+1].startswith('-'):
                argDict['options'].append(arg[1:])
                i += 1
                continue
            # a new argument
            arg = arg[1:]
            i+=1
            while i < len(sys.argv):
                if not sys.argv[i].startswith('-'):
                    if argDict.has_key(arg):
                        if type(argDict[arg]) == list:
                            argDict[arg].append(sys.argv[i])
                        else:
                            argDict[arg] = [argDict[arg]]
                            argDict[arg].append(sys.argv[i])
                    else:
                        argDict[arg] = sys.argv[i]
                    i += 1
                else:
                    i -= 1
                    break

        else:
            argDict['options'].append(arg)
        i += 1
    return argDict

def parseBC(fileStream):
    out = {}
    #parses a breakpoint file into a region:count dictionary
    while 1:
        line = fileStream.readline()
        if not line:
            break
        if not len(line.split())>1:
            continue
        line = line.strip()
        line = line.strip().split()
        regionBegin = int(line[0].replace('(','').replace(',',''))
        regionEnd = int(line[1].replace(')',''))
        recCount = (float(line[2]))
        nonRecCount = float(line[3])

        if (regionBegin, regionEnd) in out.keys():
            out[(regionBegin, regionEnd)] = (out[(regionBegin, regionEnd)][0]+recCount, \
                                            out[(regionBegin, regionEnd)][1]+nonRecCount)
        else:
            out[(regionBegin,regionEnd)] = (recCount, nonRecCount)
    return out
def sigFigs(number,digits):
    number = str(number)
    return number[:(number.index('.')+digits+1)]
def coverage(fileName):
    coverageDict = {}
    fileHandle = open(fileName,'r')
    total = 0
    while 1:
        total += 1
        line = fileHandle.readline()
        if not total%100000:
            sys.stderr.write('Examining read number %s\n'%total)
        if not line:
            break
        elif line.startswith('ExAl'):
            tStart= int(line.split(',')[1])
        if tStart in coverageDict.keys():
            coverageDict[tStart] += 1
        else:
            coverageDict[tStart] = 1
    return coverageDict


def codon_synonymous(codonL1, codonL2):
    '''returns whether or not two codons code for the same AA'''
    codon1 = codon2 = ''
    for i in range(3):
        codon1+=str(codonL1[i])
        codon2+=str(codonL2[i])
    return table.forward_table[codon1] == table.forward_table[codon2]

def is_synonymous(base, position):
    start_coding=0
    try:
        pos_in_codon = (position - start_coding)%3
    except:
        print 'ERRor', pos_in_codon, position, start_coding
        sys.exit()
    start = position-pos_in_codon
    codon = list(template[start:start+3])
    newCodon = copy.deepcopy(codon)
    newCodon[pos_in_codon] = base
    return codon_synonymous(codon,newCodon)


if __name__ == '__main__':
    args = get_cmd_args()
    if 'unique' in args['options']:
        inFile = open(args['f'])
        BC_dict = parseBC(inFile)
        for key in BC_dict:
            print key,'\t',BC_dict[key][0],'\t',BC_dict[key][1]
        inFile.close()
    elif 'coverage' in args['options']:
        coverageDict = {}
        inFiles = args['f']
        if type(inFiles) != list:
            inFiles = [inFiles]
        for fileName in inFiles: 
            coverageDict.update(coverage(fileName))
        for key in coverageDict:
            print key,coverageDict[key]
    elif 'mutations' in args['options']:
        pass
def sortDict(D_in):
    D = copy.deepcopy(D_in)
    valSort = D.values()
    valSort.sort()
    out = []
    for val in valSort:
        for key in D:
            if D[key] == val:
                out.append(key)
                del D[key]
                break
    return out
