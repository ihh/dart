import re
import sys
import utilGen

class Stk(object):
    '''Stockholm alignment.'''
    def __init__(self, pStr=None, pList=None, pFile=None):
	'''Constructor.

	Input: 
	   pStr (opt) = alignment as a string
           pList (opt) = alignment as a list
	   pFile (opt) = file containing alignment'''
	# Gap characters
	self.gapchar = {'-':1, '.':1}
	# Preamble
	self.preamble = ''
	# File data
	self.gf = {}
	# Sequence name and data
	self.seqs = {}
        self.seqName = []
	# By sequence data
	self.gr = {}
	# Consensus data
	self.gc = {}
        # String representation
        self.str = ''
	if pFile: self.fromFile(pFile)
	elif pStr: self.fromStr(pStr)
        elif pList: self.fromList(pList)

    def fromStr(self, pStr):
        self.fromList(pStr.splitlines(True))

    def fromFile(self,pFile):
        stkFile = file(pFile)
        stkList = stkFile.readlines()
        self.fromList(stkList)

    def fromList(self, pList):
        self.str = ''.join(pList)
	reGr = re.compile(r'\s*\#=GR\s*(\S+)\s*(\S+)\s*(\S+)\s*$')
	reGc = re.compile(r'\s*\#=GC\s*(\S+)\s*(\S+)\s*$')
	reGf = re.compile(r'\s*\#=GF\s*(\S+)\s*(\S+)\s*$')
	rePre = re.compile(r'\s*\#')
	reSep = re.compile(r'\s*\/\/')
	reSeq = re.compile(r'\s*(\S+)\s*(\S+)\s*$')
	reNonWs = re.compile(r'\S')
	for line in pList:
	    #print line,
	    # GR lines
	    match = reGr.match(line)
	    if match:
		seqName = match.group(1)
		if seqName not in self.gr: self.gr[seqName] = []
		utilGen.append(self.gr[seqName], match.group(2), match.group(3))
		continue
	    # GC lines
	    match = reGc.match(line)
	    if match:
		utilGen.append(self.gc, match.group(1), match.group(2))
		continue
	    # GF lines
	    match = reGf.match(line)
	    if match:
		utilGen.append(self.gf, match.group(1), match.group(2))
		continue
	    # Unrecognized line starting with '#'; append to preamble
	    if rePre.match(line):
		self.preamble += line
		continue
	    # Alignment separator; exit loop
	    if reSep.match(line): 
		break
	    # Sequence line
	    match = reSeq.match(line)
	    if match:
                name = match.group(1)
                if name not in self.seqs:
                    self.seqName.append(name)
		utilGen.append(self.seqs, name, match.group(2))
		continue
	    # Unrecognized non-whitespace line
	    if reNonWs.search(line):
		print >>sys.stderr, "Ignoring line: ", line,
	# Concatenate string lists
	utilGen.join(self.seqs)
	utilGen.join(self.gc)
	utilGen.join(self.gf)
	for key in self.gr:
	    utilGen.join(self.gr[key])

    def getStr(self):
        return self.str


