#!/usr/bin/perl

use Getopt::Long;
use SequenceIterator qw(iterseq printseq revcomp);
use Stockholm;

my $usage = "";
$usage .= "$0 -- convert DNA to tokenized-codon sequence (or protein sequence)\n";
$usage .= "\n";
$usage .= "Usage: $0 [-f <frame>] [-revcomp] [-aa] [-rna] [-decode] [-truncate] [-align <alignment file>] [filename(s)]\n";
$usage .= "\n";
$usage .= "The 'frame' (i.e. reading frame) can be 0, 1, or 2.\n";
$usage .= "If '-truncate' is specified, terminal stop codons will be discarded,\n";
$usage .= "and premature stop codons will result in truncation and a warning being issued.\n";
$usage .= "\n";
$usage .= "If a (Stockholm) alignment file is specified, the output will be that alignment\n";
$usage .= "but with all sequences replaced by the correspondingly-named tokenized sequences.\n";
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
        'gtg'=>'V',  'gcg'=>'A',  'gag'=>'E',  'ggg'=>'G',

	'n'=>'z', 'nn'=>'Z', 'nnn'=>'X',
	map (($_ x 3 => $_), qw(* - ? . x)) );

# ASCII characters 33 through 126:
# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
# Characters avoided in token set:
# ( ) ; " ' (used by S-expression format)
# - * ? . (used for alignment/reconstruction/ambiguity)
# > (used by FASTA format)
%tok = ( 'ttt'=>'F',  'tct'=>'S',  'tat'=>'Y',  'tgt'=>'C',
	 'ttc'=>'f',  'tcc'=>'s',  'tac'=>'y',  'tgc'=>'c',
	 'tta'=>'L',  'tca'=>'5',  'taa'=>'0',  'tga'=>'2',
	 'ttg'=>'l',  'tcg'=>'$',  'tag'=>'1',  'tgg'=>'W',
       
	 'ctt'=>'<',  'cct'=>'P',  'cat'=>'H',  'cgt'=>'R',
	 'ctc'=>'[',  'ccc'=>'p',  'cac'=>'h',  'cgc'=>'r',
	 'cta'=>'{',  'cca'=>',',  'caa'=>'Q',  'cga'=>'+',
	 'ctg'=>'/',  'ccg'=>'8',  'cag'=>'q',  'cgg'=>'}',
       
	 'att'=>'I',  'act'=>'T',  'aat'=>'N',  'agt'=>'%',
	 'atc'=>'i',  'acc'=>'t',  'aac'=>'n',  'agc'=>'#',
	 'ata'=>'|',  'aca'=>'~',  'aaa'=>'K',  'aga'=>'@',
	 'atg'=>'M',  'acg'=>'`',  'aag'=>'k',  'agg'=>']',
       
	 'gtt'=>'V',  'gct'=>'A',  'gat'=>'D',  'ggt'=>'G',
	 'gtc'=>'v',  'gcc'=>'a',  'gac'=>'d',  'ggc'=>'g',
	 'gta'=>'^',  'gca'=>'4',  'gaa'=>'E',  'gga'=>'9',
	 'gtg'=>'7',  'gcg'=>'&',  'gag'=>'e',  'ggg'=>'6',

	 'n'=>'z', 'nn'=>'Z', 'nnn'=>'X',
	 map (($_ x 3 => $_), qw(* - ? . x)) );

my %is_stop = map (($_ => 1), qw(tag taa tga));

my $frame = 0;
my $revcomp = 0;
my $use_aa = 0;
my $is_rna = 0;
my $untokenize = 0;
my $truncate = 0;
my $align_file;

GetOptions ("frame=i" => \$frame,
	    "revcomp" => \$revcomp,
	    "aa"  => \$use_aa,
	    "rna"  => \$is_rna,
	    "truncate" => \$truncate,
	    "align=s" => \$align_file,
	    "decode" => \$untokenize) or die $usage;

if (defined($align_file) && $untokenize) { die "Can't use -decode and -align options together\n" }
if (defined($align_file) && !$truncate) { warn "Warning: -align option is best used with -truncate\n" }
my ($stock, %untranslated);
if (defined $align_file) {
    $stock = Stockholm->from_file ($align_file);
    %untranslated = map ($stock->seqdata->{$_} =~ /[^\-\*\.]/ ? ($_ => 1) : (), @{$stock->seqname});
}

my $trans_ref = $use_aa ? \%aa : \%tok;
my %untok = map (($$trans_ref{$_} => $_), keys %$trans_ref);

# Uncomment to check for duplicate tokens
for my $c (keys %tok) { die $c unless $untok{$tok{$c}} eq $c }

@ARGV = ("-") unless @ARGV;
for my $filename (@ARGV) {
    iterseq ($filename,
	     sub {
		 my ($name, $seq) = @_;
		 my $newseq = $untokenize ? untokenize($seq,$name) : tokenize($seq,$name,$trans_ref,1);
		 if (defined $align_file) {
		     if (defined ($stock->seqdata->{$name})) {
			 my $stockrow = uc $stock->seqdata->{$name};
			 my $stockseq = $stockrow;
			 $stockseq =~ s/[\-\.]//g;
			 my $aaseq = tokenize($seq,$name,\%aa,0);
			 if ($stockseq ne $aaseq) {
			     warn
				 "Translation of sequence '$name' does not match corresponding alignment row.\n",
				 " Alignment row: $stockrow\n",
				 " Alignment seq: $stockseq\n",
				 "Translated seq: $aaseq\n",
				 " Tokenized seq: $newseq\n\n";
			 } else {
			     my $newseq_pos = 0;
			     for (my $row_pos = 0; $row_pos < length($stockrow); ++$row_pos) {
				 my $row_char = substr ($stockrow, $row_pos, 1);
				 if ($row_char ne '-' && $row_char ne '.') {
				     substr ($stockrow, $row_pos, 1) = substr ($newseq, $newseq_pos++, 1);
				 }
			     }
			     $stock->seqdata->{$name} = $stockrow;
			     delete $untranslated{$name};
			 }
		     } else {
			 warn "Sequence '$name' not found in alignment; ignoring\n";
		     }
		 } else {
		     printseq ($name, $newseq);
		 }
	     });
}

if (defined $align_file) {
    for my $name (keys %untranslated) { $stock->seqdata->{$name} =~ s/[^\-\.]/*/g }
    print $stock->to_string;
    if (%untranslated) {
	warn "The following alignment rows could not be tokenized, so were replaced with wildcards:\n", join (" ", keys %untranslated), "\n";
    }
}


sub tokenize {
    my ($seq, $name, $trans_ref, $warning) = @_;
    $seq = lc $seq;
    $seq =~ s/u/t/g;  # do this even if -rna was not specified; no need to punish user
    if ($revcomp) { $seq = revcomp ($seq); $name .= " (reverse strand)" }
    my $trans = "";
  CODON: for (my $pos = $frame; $pos < length($seq); $pos += 3) {
	my $remaining_chars = length($seq) - ($pos + 3);
	$codon = substr ($seq, $pos, 3);
	if (exists $$trans_ref{$codon}) {
	    if ($truncate && $is_stop{$codon}) {
		if ($warning && $remaining_chars > 0) {
		    warn "Premature stop codon ($codon) found with $remaining_chars characters remaining while tokenizing sequence $name\n";
		}
		last CODON;
	    }
	    $trans .= $$trans_ref{$codon};
	} elsif (length($codon) == 1) { $trans .= $$trans_ref{'n'}; warn "Extra character ($codon) at end of sequence $name\n" if $warning }
	elsif (length($codon) == 2) { $trans .= $$trans_ref{'nn'}; warn "Extra characters ($codon) at end of sequence $name\n" if $warning }
	else { $trans .= $$trans_ref{'nnn'}; warn "Unrecognized codon ($codon) at position $pos of sequence $name\n" if $warning }
    }
    return $trans;
}

sub untokenize {
    my ($seq, $name) = @_;
    my $untrans = "n" x $frame;
    for (my $pos = 0; $pos < length($seq); ++$pos) {
	$token = substr ($seq, $pos, 1);
	if (exists $untok{$token}) { $untrans .= $untok{$token} }
	else { $untrans .= 'nnn'; warn "Unrecognized token ($token) at position $pos of sequence $name\n" }
    }
    if ($revcomp) { $untrans = revcomp ($untrans) }
    if ($is_rna) { $untrans =~ s/t/u/g }
    return $untrans;
}
