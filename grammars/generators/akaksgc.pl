#!/usr/bin/perl -w

# flag indicating whether to use codon or nucleotide frequencies in mutation rates
# for now we use codon frequencies, as this is most similar to CodeML
#  bear in mind that using codon frequencies instead of base frequencies, but keeping rate parameters fixed,
#  effectively gives you rates that are lower (by a factor that I hazily estimate at around 48 = 3*4*4)
my $useCodonFreqs = 1;

# transition table
my %transition = qw (a g
		     c t
		     g a
		     t c);

# codon table
my %aa = qw( ttt F  tct S  tat Y  tgt C
	     ttc F  tcc S  tac Y  tgc C
	     tta L  tca S  taa X  tga X
	     ttg L  tcg S  tag X  tgg W
	     
	     ctt L  cct P  cat H  cgt R
	     ctc L  ccc P  cac H  cgc R
	     cta L  cca P  caa Q  cga R
	     ctg L  ccg P  cag Q  cgg R
	     
	     att I  act T  aat N  agt S
	     atc I  acc T  aac N  agc S
	     ata I  aca T  aaa K  aga R
	     atg M  acg T  aag K  agg R
	   
	     gtt V  gct A  gat D  ggt G
	     gtc V  gcc A  gac D  ggc G
	     gta V  gca A  gaa E  gga G
	     gtg V  gcg A  gag E  ggg G );

# alphabet
my @alph = qw(a c g t);

# usage
my $usage = "Usage: $0 <template file>\n";
die $usage if @ARGV != 1 or $ARGV[0] =~ /^-h/;;

# path to template
my $tmpl_path = shift;

# load template
open TMPL, "<$tmpl_path" or die $!;
my @tmpl = <TMPL>;
close TMPL;
my $tmpl = join "", @tmpl;

# make rates
my (@initial, @mutate);
for my $ijk (sort keys %aa) {
    my @ijk = split //, $ijk;
    my ($i, $j, $k) = @ijk;
    for my $x (@alph) {
	if ($x ne $i) { push @mutate, mutate ($ijk, "$x$j$k") }
	if ($x ne $j) { push @mutate, mutate ($ijk, "$i$x$k") }
	if ($x ne $k) { push @mutate, mutate ($ijk, "$i$j$x") }
    }
    push @initial, "(initial (state ($i $j $k)) (prob " . join(" * ",map(p_branch_var($ijk[$_],$_),0..2)) . "))  ;; $aa{$ijk}" unless $aa{$ijk} eq "X";
}


# substitute into template
my $mutPattern = '(\s+)[^\n]*==MUTATE\((\S+?)\)==';
my $nodeVar = node_var();
while ($tmpl =~ /$mutPattern/) {
    my ($spc, $arg) = ($1, $2);
    my $mutate = join ("", map ("$spc$_", @initial, @mutate));
    $mutate =~ s/$nodeVar/$arg/g;
    $tmpl =~ s/$mutPattern/$mutate/;
}

# print
print $tmpl;

# create mutate expression
sub mutate {
    my ($ijk, $xyz) = @_;
    my ($i, $j, $k) = split //, $ijk;
    my ($x, $y, $z) = split //, $xyz;
    my ($rate, $comment) = rateFunc($ijk,$xyz);
    return
	defined($rate)
	? ("(mutate (from ($i $j $k)) (to ($x $y $z)) (rate ($rate)))  ;; $comment")
	: ();
}

# get rate for a given mutation
sub rateFunc {
    my ($ijk, $xyz) = @_;

    return (undef,undef) if $aa{$ijk} eq "X" || $aa{$xyz} eq "X";

    my $nDiff = 0;
    my ($isTransition, @pVar);
    for my $pos (0..2) {
	my $ijk_pos = substr ($ijk, $pos, 1);
	my $xyz_pos = substr ($xyz, $pos, 1);
	if ($ijk_pos ne $xyz_pos) {
	    $isTransition = $transition{$ijk_pos} eq $xyz_pos;
	    push @pVar, p_branch_var ($xyz_pos,$pos);
	    ++$nDiff;
	} elsif ($useCodonFreqs) {
	    push @pVar, p_branch_var ($xyz_pos,$pos);
	}
    }
    return (undef,undef) unless $nDiff == 1;
    my $isSynonymous = $aa{$ijk} eq $aa{$xyz};
    my $rVar = branch_var ("R" . ($isSynonymous ? "s" : "n") . ($isTransition ? "i" : "v"));
    return (join (" * ", $rVar, @pVar), "$aa{$ijk} -> $aa{$xyz}");
}

# branch-specific variable name generator
sub node_var { return "NODE" }
sub branch_var { return "(. " . shift() . " _ " . node_var() . ")" }
sub p_branch_var { my ($new, $pos) = @_; return branch_var ("p$new$pos") }
