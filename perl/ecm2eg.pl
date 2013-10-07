#!/usr/bin/perl -w

use Getopt::Long;
use SequenceIterator qw(iterseq);

my $usage = "";
$usage .= "$0 -- convert ECM models (Kosiol, Holmes & Goldman, MBE 2007) to protpal/handalign chain format\n";
$usage .= "\n";
$usage .= "Usage: $0 [-seqfile <tokenized FASTA file>] [-pseudocount <pseudocount>] <ECM model file>\n";
$usage .= "\n";
$usage .= "Use the -s option to reparameterize codon frequencies.\n";
$usage .= "\n";

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

# this must match the token set in $DARTDIR/perl/tokenize.pl
%tok = map (($_ => join(" ",split(//,$_))), keys %aa);

my $seqfile;
my $pseudocount = 0;
GetOptions ("seqfile=s" => \$seqfile,
	    "pseudocount=f" => \$pseudocount) or die $usage;

@ARGV = ("-") unless @ARGV;
my @ecm = <>;
grep (chomp, @ecm);

my @codon = split (/\s+/, lc ("@ecm[65..68]"));
my @token = map ($tok{$_}, @codon);  # just the codons in the model (stop codons excluded)
my @aa = map ($aa{$_}, @codon);

my @all_token = sort values %tok;  # all codons, including stop codons

my @original_pi = split (/\s+/, $ecm[62]);
my @pi = @original_pi;

if (defined $seqfile) {
    my %freq = map (($_ => $pseudocount), @token);
    my $total = 0;
    iterseq ($seqfile,
	     sub {
		 my ($name, $seq) = @_;
		 for my $pos (0..length($seq)-1) {
		     ++$freq{substr($seq,$pos,1)};
		     ++$total;
		 }
		 @pi = map ($freq{$token[$_]} / $total, 0..@token-1);
	     });
}

print "(alphabet\n (name DNA)\n (token (a c g t))\n (wildcard *))\n\n";
print "(chain\n (update-policy rev)\n (terminal (X1 X2 X3))\n";
for my $i (0..@codon-1){
    print "(initial (state (";
    print join(' ',split '',$codon[$i]); 
    print ")) (prob $pi[$i]))  ;; $codon[$i] [", $aa{$codon[$i]}, "]\n"; 
}


my ($total, $original_total) = (0, 0);
for my $pass (1..2) {
    for my $j (1..60) {
	my @s_ij = split (/\s+/, $ecm[$j-1]);
	for my $i (0..$j-1) {
	    my $s_ij = $s_ij[$i];
	    if ($i != $j) {
		my ($i_cod, $j_cod) = map (lc ($codon[$_]), $i, $j);
		my $i_cod_space = join(' ',split '',$i_cod); 
		my $j_cod_space = join(' ',split '',$j_cod); 
		my ($i_tok, $j_tok) = map ($tok{$_}, $i_cod, $j_cod);
		my ($i_aa, $j_aa) = map ($aa{$_}, $i_cod, $j_cod);
		if ($pass == 1) {
		    $original_total += 2 * $original_pi[$i] * $original_pi[$j] * $s_ij;
		    $total += 2 * $pi[$i] * $pi[$j] * $s_ij;
		} else {  # $pass == 2
		    my $r_ij = $s_ij * $pi[$j] / $total;
		    my $r_ji = $s_ij * $pi[$i] / $total;
		    print " (mutate (from ($i_cod_space)) (to ($j_cod_space)) (rate $r_ij))  ;; $i_cod [$i_aa] -> $j_cod [$j_aa]\n" if $r_ij > 0;
		    print " (mutate (from ($j_cod_space)) (to ($i_cod_space)) (rate $r_ji))  ;; $j_cod [$j_aa] -> $i_cod [$i_aa]\n" if $r_ji > 0;

		}
	    }
	}
    }
}

print ")\n";

warn "Check: total rate using original codon frequencies = $original_total (should be close to 1)\n";
warn "Rate normalizer for modified codon frequencies = $total\n";

