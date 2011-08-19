#!/usr/bin/perl

use Getopt::Long;
use SequenceIterator qw(iterseq printseq revcomp);

my $usage = "";
$usage .= "$0 -- convert DNA to tokenized-codon sequence (or protein sequence)\n";
$usage .= "\n";
$usage .= "Usage: $0 [-f <reading frame>] [-revcomp] [-aa] [-decode] [filename(s)]\n";
$usage .= "\n";
$usage .= "The reading frame can be 0, 1, 2, -0, -1, -2, or a comma-separated list, or \"all\" (or equivalently \".\"); 0 is assumed by default.\n";
$usage .= "\n";

my $colWidth = 50;    # column width for output

my (%aa, %tok);
%aa = ( 'ttt'=>'F',  'tct'=>'S',  'tat'=>'Y',  'tgt'=>'C',
        'ttc'=>'F',  'tcc'=>'S',  'tac'=>'Y',  'tgc'=>'C',
        'tta'=>'L',  'tca'=>'S',  'taa'=>'!',  'tga'=>'!',
        'ttg'=>'L',  'tcg'=>'S',  'tag'=>'!',  'tgg'=>'W',
       
        'ctt'=>'L',  'cct'=>'P',  'cat'=>'H',  'cgt'=>'R',
        'ctc'=>'L',  'ccc'=>'P',  'cac'=>'H',  'cgc'=>'R',
        'cta'=>'L',  'cca'=>'P',  'caa'=>'Q',  'cga'=>'R',
        'ctg'=>'L',  'ccg'=>'P',  'cag'=>'Q',  'cgg'=>'R',
       
        'att'=>'I',  'act'=>'T',  'aat'=>'N',  'agt'=>'S',
        'atc'=>'I',  'acc'=>'T',  'aac'=>'N',  'agc'=>'S',
        'ata'=>'I',  'aca'=>'T',  'aaa'=>'K',  'aga'=>'R',
        'atg'=>'M',  'acg'=>'T',  'aag'=>'K',  'agg'=>'R',
       
        'gtt'=>'V',  'gct'=>'A',  'gat'=>'D',  'ggt'=>'G',
        'gtc'=>'V',  'gcc'=>'A',  'gac'=>'D',  'ggc'=>'G',
        'gta'=>'V',  'gca'=>'A',  'gaa'=>'E',  'gga'=>'G',
        'gtg'=>'V',  'gcg'=>'A',  'gag'=>'E',  'ggg'=>'G' );

%tok = ( 'ttt'=>'F',  'tct'=>'S',  'tat'=>'Y',  'tgt'=>'C',
	 'ttc'=>'f',  'tcc'=>'s',  'tac'=>'y',  'tgc'=>'c',
	 'tta'=>'L',  'tca'=>'5',  'taa'=>':',  'tga'=>'!',
	 'ttg'=>'l',  'tcg'=>'$',  'tag'=>';',  'tgg'=>'W',
       
	 'ctt'=>'<',  'cct'=>'P',  'cat'=>'H',  'cgt'=>'R',
	 'ctc'=>'[',  'ccc'=>'p',  'cac'=>'h',  'cgc'=>'r',
	 'cta'=>'(',  'cca'=>',',  'caa'=>'Q',  'cga'=>')',
	 'ctg'=>'{',  'ccg'=>"'",  'cag'=>'q',  'cgg'=>'}',
       
	 'att'=>'I',  'act'=>'T',  'aat'=>'N',  'agt'=>'%',
	 'atc'=>'i',  'acc'=>'t',  'aac'=>'n',  'agc'=>'#',
	 'ata'=>'|',  'aca'=>'~',  'aaa'=>'K',  'aga'=>'/',
	 'atg'=>'M',  'acg'=>'`',  'aag'=>'k',  'agg'=>'\\',
       
	 'gtt'=>'V',  'gct'=>'A',  'gat'=>'D',  'ggt'=>'G',
	 'gtc'=>'v',  'gcc'=>'a',  'gac'=>'d',  'ggc'=>'g',
	 'gta'=>'^',  'gca'=>'@',  'gaa'=>'E',  'gga'=>'_',
	 'gtg'=>'>',  'gcg'=>'&',  'gag'=>'e',  'ggg'=>'=' );

my $frame = 0;
my $revcomp = 0;
my $use_aa = 0;
my $untokenize = 0;

GetOptions ("frame=i" => \$frame,
	    "revcomp" => \$revcomp,
	    "aa"  => \$use_aa,
	    "decode" => \$untokenize) or die $usage;

my $trans_ref = $use_aa ? \%aa : \%tok;
my %untok = map (($$trans_ref{$_} => $_), keys %$trans_ref);

@ARGV = ("-") unless @ARGV;
for my $filename (@ARGV) {
    iterseq ($filename, $untokenize ? \&untokenize : \&tokenize);
}

sub tokenize {
    my ($name, $seq) = @_;
    $seq = lc $seq;
    if ($revcomp) { $seq = revcomp ($seq) }
    my $trans = "";
    for (my $pos = $frame; $pos < length($seq); $pos += 3) {
	$codon = substr ($seq, $pos, 3);
	if (exists $$trans_ref{$codon}) { $trans .= $$trans_ref{$codon} }
    }
    printseq ($name, $trans);
}

sub untokenize {
    my ($name, $seq) = @_;
    my $untrans = "n" x $frame;
    for (my $pos = 0; $pos < length($seq); ++$pos) {
	$token = substr ($seq, $pos, 1);
	if (exists $untok{$token}) { $untrans .= $untok{$token} }
    }
    printseq ($name, $untrans);
}
