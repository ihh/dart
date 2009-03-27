#! /usr/bin/perl

#INPUT: stockholm file name , three letter abbreviation for the schoring scheme , binning similiarity cutoff, ouput filename for consensus files


use Stockholm;
use Stockholm::Database;

my $usage= "\nUsage: stockdb.stock SSA #.## output.stock\n
              stockdb.stock is the input file, a stockholm format file containing the multiple alignments over the same sequence\n
              SSA is a three letter abbreviation for the schoring scheme, either AMA SPS PPV or TCS\n
              #.## is the scoring similarity cutoff for different bins\n
              output.stock is the name of the output file\n
              [-h, --help] display this message\n
              [-v, --verbose] print lots of information about what we're doing\n";

if (!defined(@ARGV)){
    die $usage;
}

my $checker=1;
while (@ARGV){
    my $arg = shift @ARGV;
    if ($arg =~ /^-/ ){
	if (($arg eq "-h" ) || ($arg eq "--help")) {print "$usage"; exit;}
	elsif (($arg eq "-v") || ($arg eq "--verbose")) {$verbose=1;}
	else {die $usage;}
    }
    elsif(($arg=~ /\.stock/) && $checker){$infile=$arg; $checker=0;}
    elsif($arg=~ /AMA|SPS|PPV|TCS/){ $scoring = $arg;}
    elsif($arg=~ /\.\d\d/){ $cutoff = $arg;}
    elsif($arg=~ /\.stock/) { $outfile=$arg;}
}


#input: reference to a bin, reference to distance matrix, scheme being used
#output: array (best average score in bin , index of that stockholm object in stockArray)
sub consen_align{

    my @bin=@{$_[0]};
    my @dMatrix=@{$_[1]};
    my $scheme=$_[2];
    my $n=scalar(@bin);
    my $AVSC=0;

    #creates sum of scores array
    my @measures=();
    my %measure_key=();
    $measure_key{'scale'}=$n;

    for(my $i=0;$i<$n;$i++){
        #initialize the sum by setting it to zero
        $measures[$i]=0;

        foreach $stock (@bin){
            if ($stock<$bin[$i]){
                $measures[$i]+=$dMatrix[$stock][$bin[$i]];
            }
            elsif($stock>$bin[$i]){ #distance between 2 and 1 has to be looked up as $dmatrix[1][2]
                $measures[$i]+=$dMatrix[$bin[$i]][$stock];
            }
        }
	$measure_key{$bin[$i]}=$measures[$i];
    }

    my $bestScore=$measure_key{$bin[0]};
    my $bestIndex=0;
    #loops through, finds highest score
    for ($i=1;$i<$n;$i++){
	if ($bestScore<$measure_key{$bin[$i]}){
	    $bestScore=$measure_key{$bin[$i]};
	    $bestIndex=$i;
	}
    }
    #scale best score and return it
    return $bin[$bestIndex];
}




#PROGRAM STARTS HERE
#get database object of stockholm flatfile, and number of alignments
my $db = Stockholm::Database->from_file ($infile);
my @stockArray=@{$db};
my $n=scalar(@stockArray);
if ($verbose){print "read in stockholm database from $infile containing $n alignments\n";}

#get the score wanted based on command line argument using the following lookup table
my %scoreScheme = ( "AMA" , 0,
                    "SPS" , 1,
                    "PPV" , 2,
                    "TCS" , 3);

my $scheme=$scoreScheme{$scoring};


#create the distace matrix
my @dMatrix=();

if ($verbose) {print "using scoring scheme $scoring\n";}
if ($verbose) {print "making distance matrix...";}
for (my $i=0;$i<scalar(@stockArray);$i++){
    for (my $j=$i+1; $j<scalar(@stockArray);$j++){

 
	#outputs stockholm objects stored in array to files so they can be read by cmpalign
        $stockArray[$i]->to_file ("stock1.stock");
        $stockArray[$j]->to_file ("stock2.stock");

        #compalign outputs a string of scores deliminated by tabs, shove that into a array
        #get the score we are interested in and store it in a two by two matrix
        $scoreStr = `/nfs/dart/bin/cmpalign -s stock1.stock stock2.stock`;
        @scoreArr = split /\t/ , $scoreStr;
        $dMatrix[$i][$j]=$scoreArr[$scheme];

    }
}

`rm stock1.stock`;
`rm stock2.stock`;

if ($verbose){print "done\n;"}

#binning
my @indices=(0 .. $n-1);
my @currBin=();
my @bins=();
my $starter=0; #this is the starting inxed

if ($verbose) {print "binning...";}

while ($indices[@indices-1]!=0){       #check if last element is zero (if last element has been binned)
    push @currBin, $indices[$starter]; #add the current start element to a new bin
    $indices[$starter]=0;              #once an element is added to a bin it is replaced with zero in indices array
    $starter=0;
    foreach $index (@indices){
        if($index!=0){                 #if it doesn't equal zero (if it hasn't already been added to a bin)
            if ($dMatrix[$currBin[0]][$index]>$cutoff){ # and if the distance is more then the cutoff value
	    push @currBin , $indices[$index]; #add it to the current bin
	    $indices[$index]=0; #replace it with zero to show its been added
	}
	elsif($starter==0){ #if starter is zero then an element hasn't been found which doesn't fit in a current bin
	    $starter=$index; #mark it to be the first element in the new bin
	}
    }
}
    #add and reset
    push @bins, [@currBin];
    @currBin=();
}
if ($verbose){print "done\n";}

#output a bunch of consensus alignments
$name=$outfile;
$counter=1;
if ($verbose){print "outputting consensus alignments...";}
foreach my $binref (@bins){
    $AVSC=0;

    if(scalar(@{$binref})==1){

	$bestAlign=$stockArray[${$binref}[0]];
	$bestAlign->add_gf (AVSC , $bestAlign->get_gf("SC"));
	$bestAlign->to_file ( "$name.$counter.stock" );
    }
    else{
	consen_align($binref,\@dMatrix,$scheme);
	$bestAlign=$stockArray[consen_align($binref,\@dMatrix,$scheme)];

	foreach my $stockInd (@{$binref}){
	    $AVSC+=$stockArray[$stockInd]->get_gf("SC");
	}
	$AVSC=$AVSC/scalar(@{$binref});

	$bestAlign->add_gf (AVSC , $AVSC);
	$bestAlign->to_file ( "$name.$counter.stock" );
    }
    $counter++;
}
if($verbose){print "done\n";}





