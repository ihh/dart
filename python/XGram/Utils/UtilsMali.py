"""This modules provides mali, a rudimentary alignment class.

"""

import re, sys, string, time

class AlignedString:
    def __init__(self, identifier, fr, to, s):
        self.mId = identifier
        self.mFrom = fr
        self.mTo = to
        self.mString = s

    def __len__(self):
        return len(self.mString)

class Mali:
    """a rudimentary multiple alignment class.
    
    This class reads and writes alignments in various
    formats and provides a few edit operations.
    """
    
    mGapChars = ("-", ".")
    mGapPattern = re.compile("[.-]")
    # character used for gaps
    mGapChar = "-"
    # admissable gap characters
    mGapChars = ("-", ".")
    mMaskChar = "X"
    
    def __init__(self):

        self.mIdentifiers = []
        self.mMali = {}
        self.mLength = 0
        self.mAnnotations = {}
        
    def __getitem__(self, key):
        return self.mMali[key].mString

    def __len__(self):
        return len(self.mIdentifiers)


    def getLength(self):
        return self.__len__()

    def getWidth(self):
        return len(self.mMali[self.mIdentifiers[0]])
    
    def countCharacters(self, row ):
        return len(row) - len(self.mGapPattern.findall(row))
    
    def readFromFile( self, infile, format = "fasta" ):
        """read multiple alignment from file in various format."""
        
        self.mMali = {}
        self.mIdentifiers = []

        pattern_parse_ranges=re.compile("(\S+)/(\d+)-(\d+)")

        lines=  infile.readlines()

        if format not in ("stockholm"):
            # save comments 
            self.mComments = filter( lambda x: x[0] == "#", lines)
            lines = filter( lambda x: x[0] != "#", lines)
        else:
            self.mComments = []
            
        if format.lower() == "plain":

            for line in lines:
                data = line[:-1].split("\t")
                self.mIdentifiers.append( data[3] )
                self.mMali[data[3]] = AlignedString(data[3], int(data[0]) - 1, int(data[2]), data[1] )

        elif format.lower() == "fasta":
            pattern_identifier = "\S+"
            id = None
            fragments = []
            for line in lines:
                if line[0] == ">":
                    if id:
                        s = re.sub( "\s", "", string.join( fragments, ""))
                        x = pattern_parse_ranges.match( id )
                        if x:
                            id, fr, to = x.groups()
                            fr, to = int(fr) - 1, int(to)
                        else:
                            fr, to = 0, self.countCharacters( s )
                            
                        self.mIdentifiers.append(id)
                        self.mMali[id] = AlignedString(id, fr, to, s)
                        
                    id = re.search( "^(%s)" % pattern_identifier, line[1:-1]).group(0)
                    fragments = []
                    continue
                fragments.append( line[:-1] )

            if id:
                s = re.sub( "\s", "", string.join( fragments, ""))
                x = pattern_parse_ranges.match( id )
                if x:
                    id, fr, to = x.groups()
                    fr, to = int(fr) - 1, int(to)
                else:
                    fr, to = 0, self.countCharacters( s )
                    
                self.mIdentifiers.append(id)                    
                self.mMali[id] = AlignedString(id, fr, to, s)
                
        elif format.lower() == "clustal":
            ## skip header line
            del lines[0]
            fragments = {}

            ## prune lines
            lines = map( lambda x: x.strip(), lines )
            ## remove empty lines
            lines = filter( lambda x: len(x[:-1]) > 0, lines )
            
            for line in lines:
                ## remove consensus lines
                if line[0] in ("*", ":"): continue
                
                data = re.split( "\s+", line )
                if len(data) != 2:
                    raise ValueError, "parsing error in line %s" % line
                
                id, fragment = data
                if id not in fragments:
                    fragments[id] = []
                    self.mIdentifiers.append(id)
                    
                fragments[id].append( fragment )

            for id, f in fragments.items():
                s = re.sub( "\s", "", string.join( f, ""))                        
                self.mMali[id] = AlignedString(id, 0, self.countCharacters( s ), s)

        elif format.lower() == "stockholm":
            ## skip header line
            del lines[0]
            fragments = {}
            annotations = {}
            ## prune lines
            lines = map( lambda x: x.strip(), lines )
            ## remove empty lines
            lines = filter( lambda x: len(x[:-1]) > 0, lines )
            
            for line in lines:
                data = re.split( "\s+", line )
                
                if data[0] == "//": break
                
                if line[0] == '#':
                    if data[0] == "#=GC":
                        id, fragment = data[1:3]
                    else:
                        self.mComments.append( line )
                        continue
                    if id not in annotations:
                        annotations[id] = []
                    annotations[id].append( fragment )
                else:

                    if len(data) != 2:
                        raise ValueError, "parsing error in line %s" % line

                    id, fragment = data
                    
                    if id not in fragments:
                        fragments[id] = []
                        self.mIdentifiers.append(id)
                    
                    fragments[id].append( fragment )

            n = []
            for id in self.mIdentifiers:
                f = fragments[id]

                s = re.sub( "\s", "", string.join( f, ""))
                x = pattern_parse_ranges.match( id )
                if x:
                    id, fr, to = x.groups()
                    fr, to = int(fr) - 1, int(to)
                else:
                    fr, to = 0, self.countCharacters( s )

                n.append(id)                    
                self.mMali[id] = AlignedString(id, fr, to, s)
            self.mIdentifiers = n
            
            for id, f in annotations.items():
                s = re.sub( "\s", "", string.join( f, ""))
                annotations[id] = s
            self.mAnnotations = annotations
        else:
            raise "unknown alignment format %s" % format

        self.mLength = min( map( lambda x: x.mString, self.mMali.values() ) )
        
    def writeToFile( self, outfile, format = "plain"):
        """write alignment to file."""

        if format == "plain":
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]
                outfile.write("%i\t%s\t%i\t%s\n" % (
                    m.mFrom+1, m.mString, m.mTo, identifier) )
                    
        elif format == "fasta":
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]                
                outfile.write( ">%s/%i-%i\n%s\n" % (identifier, m.mFrom + 1, m.mTo, m.mString) )

        elif format == "stockholm":
            outfile.write("# STOCKHOLM 1.0\n")
            ## calculate offset:
            max_l = 0
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]
                # tab does not work as separator
                if m.mTo:
                    x = "%s/%i-%i" % (identifier, m.mFrom+1, m.mTo)
                else:
                    x = "%s" % (identier)
                max_l = max(max_l, len(x))
            for identifier in self.mAnnotations.keys():
                x = "#=GC %s" % identifier
                max_l = max(max_l, len(x))

            format = "%-" + str(max_l) + "s  %s\n"
            for identifier in self.mIdentifiers:
                m = self.mMali[identifier]
                # tab does not work as separator
                if m.mTo:
                    x = "%s/%i-%i" % (identifier, m.mFrom+1, m.mTo)
                else:
                    x = "%s" % (identier)

                outfile.write( format % (x, m.mString) )
                
            for identifier, value in self.mAnnotations.items():
                x = "#=GC %s" % identifier
                outfile.write( format % (x,value) )
            
            outfile.write("//\n")
        else:
            raise "unknown alignment format %s" % format        

    def removeUnalignedEnds( self ):
        """remove unaligned ends in the multiple alignment.

        unaligned ends correspond to lower-case characters.
        """
        pattern_start = re.compile( "^([- .a-z]+)" )
        pattern_unaligned = re.compile("[a-z]")
        
        for s in self.mMali.values():

            t0 = time.time()
            first = pattern_start.match( s.mString )
            if first:
                first = first.groups()[0]
                nchars = len( pattern_unaligned.findall( first ) )
                s.mFrom += nchars
                s.mString = self.mGapChar * len(first) + s.mString[len(first):]

            t0 = time.time()
            ## search from the back end by reversing. This is much faster than
            ## using $ from the back.
            last = pattern_start.match( s.mString[::-1] )
            if last:
                last = last.groups()[0]
                nchars = len( pattern_unaligned.findall( last ) )
                s.mTo -= nchars
                l = len(s) - len(last)
                s.mString = s.mString[:l] + self.mGapChar * l

    def removeEndGaps( self ):
        """remove end gaps.

        end gaps do not include any characters and thus
        the alignment coordinates won't change.
        """

        pattern_start_gaps = re.compile( "^([- ]+)" )
        
        min_from = self.mLength
        max_to = 0
        
        for s in self.mMali.values():

            first = pattern_start_gaps.match( s.mString )
            if first:
                first = first.groups()[0]
                min_from = min( min_from, len(first) )
                
            ## search from the back end by reversing. This is much faster than
            ## using $ from the back.
            last = pattern_start_gaps.search( s.mString[::-1] )
            if last:
                last = last.groups()[0]
                max_to = max( max_to, len(s) - len(last) )

        for s in self.mMali.values():
            s.mString = s.mString[min_from:max_to]

        self.mLength = min( map( lambda x: x.mString, self.mMali.values() ) )

    def removeGaps( self, allowed_gaps = 0 ):
        """remove gappy columns."""
            
        ngaps = [0] * self.mLength
        for s in map(lambda x: x[2], self.mMali.values()):
            for x in range(len(s)):
                if s[x] in self.mGapChars:
                    ngaps[x] += 1

        columns = []
        
        for x in range( len(ngaps)):
            if ngaps[x] <= allowed_gaps:
                columns.append(x)
                
        self.takeColumns( columns )

    def upper( self ):
        """convert all characters in mali to uppercase."""
        
        for s in self.mMali.values():
            s.mString = s.mString.upper()
            
    def lower( self ):
        """convert all characters in mali to lowercase."""
        
        for s in self.mMali.values():
            s.mString = s.mString.lower()

    def shiftAlignment( self, map_id2offset ):
        """shift alignment by offset."""

        for identifier, m in self.mMali.items():
            if identifier in map_id2offset:
                o = map_id2offset[identifier]
                m.mFrom += o
                m.mTo += o

    def markTransitions( self, map_id2transitions, mode="case" ):
        """mark transitions in the multiple alignment.

        if mode == case, then upper/lower case is used for the transitions

        Otherwise, a character given by mode is inserted.
        """

        if mode in ("case", "keep-odd", "keep-even"):
            for identifier, s in self.mMali.items():
                if identifier not in map_id2transitions: continue

                new_chars = []
                c = s.mFrom
                
                is_upper = True
                transitions = map_id2transitions[identifier]
                for x in s.mString:
                    if x in self.mGapChars:
                        pass
                    else:
                        if c in map_id2transitions[identifier]:
                            if is_upper: is_upper = False
                            else: is_upper = True
                        c += 1

                        if x in string.lowercase:
                            x = self.mMaskChar

                        if mode == "case":
                            if is_upper:
                                x = string.upper(x)
                            else:
                                x = string.lower(x)
                        elif mode == "keep-even":
                            if is_upper:
                                x = self.mGapChar
                        elif mode == "keep-odd":
                            if not is_upper:
                                x = self.mGapChar
                                
                    new_chars.append( x )

                s.mString = "".join(new_chars)
        else:
            raise "character insertion not implemented yet."

    def buildColumnMap( self, other, join_field = None ):
        """build map of columns in other to this."""

        if not join_field:
            join_field = other.mIdentifiers[0]
            
        if join_field not in other.mMali or \
           join_field not in self.mMali:
            raise "line %s not in both alignments." % (join_field)

        this_seq = self.mMali[join_field]
        other_seq = other.mMali[join_field]

        if this_seq.mFrom != other_seq.mFrom or \
           this_seq.mTo != other_seq.mTo:
            raise "residue ranges for sequence %s doe not correspond." % (join_field)

        map_this2other = []

        this_seq = this_seq.mString.upper()
        other_seq = other_seq.mString.upper()        

        # position in other
        o = 0
        for c in this_seq:
            if c in self.mGapChars:
                map_this2other.append( None )
            else:
                while other_seq[o] != c:
                    o += 1
                map_this2other.append( o )
                o += 1
        return map_this2other
                
    def copyAnnotations( self, other ):
        """copy annotations from annother mali."""

        map_this2other = self.buildColumnMap( other )
        ncols = self.getWidth()
        
        for key, annotation in other.mAnnotations.items():
            a = []
            for x in range(ncols):
                m = map_this2other[x] 
                if m !=  None:
                    a.append( annotation[m] )
                else:
                    a.append( self.mGapChar )
            
            self.mAnnotations[key] = "".join(a)
            
    def getColumns(self):
        """returns the multiple alignment columns."""
        
        ncols = self.getWidth()
        columns = [ [] for x in range(ncols) ]
        for identifier in self.mIdentifiers:
            s = self.mMali[identifier].mString
            for x in range(ncols):
                columns[x].append( s[x] )
                
        return columns
            
        
        