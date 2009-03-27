#!/usr/bin/perl -w

# alphabet
my @alph = qw(a c g t);

# transition & transversion tables
my %transition = qw (a g
		     c t
		     g a
		     t c);

my %transversions;
for my $a (@alph) {
    $transversions{$a} = [ grep ($_ ne $a && $_ ne $transition{$a}, @alph) ];
}

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
my @mutate;
for my $ijk (sort keys %aa) {
    my ($i, $j, $k) = split //, $ijk;
    for my $x (@alph) {
	if ($x ne $i) { push @mutate, mutate ($ijk, "$x$j$k") }
	if ($x ne $j) { push @mutate, mutate ($ijk, "$i$x$k") }
	if ($x ne $k) { push @mutate, mutate ($ijk, "$i$j$x") }
    }
    push @mutate, mutate ($ijk, $ijk);
}


# substitute into template
my $mutPattern = '(\s+)[^\n]*==MUTATE==';
if ($tmpl =~ /$mutPattern/) {
    my $spc = $1;
    $mutate = join ("", map ("$spc$_", @mutate));
    $tmpl =~ s/$mutPattern/$mutate/g;
}

# print
print $tmpl;

# create mutate expression
sub mutate {
    my ($ijk, $xyz) = @_;
    my ($i, $j, $k) = split //, $ijk;
    my ($x, $y, $z) = split //, $xyz;
    my ($rate, $comment) = rateFunc($ijk,$xyz);
    return "(mutate (from ($i $j $k)) (to ($x $y $z)) (rate ($rate)))  ;; $comment";
}

# get rate for a given mutation
sub rateFunc {
    my ($ijk, $xyz) = @_;
    my ($rFunc, $comment);
    my $k = k_branch_var();
    my $not_k = not_k_branch_var();
    my $not_w = not_w_branch_var();
    if ($ijk eq $xyz) {
	my ($selfFunc, $kFunc, $wFunc, $kwFunc) = ("", "", "");
	my @ijk = split //, $ijk;
	for my $pos (0..2) {
	    my $ijk_pos = substr ($ijk, $pos, 1);
	    my $ijk_left = join ("", @ijk[0..$pos-1]);
	    my $ijk_right = join ("", @ijk[$pos+1..2]);

	    my @nonsyn_xyz_pos = grep ($aa{$ijk} ne $aa{"${ijk_left}${_}${ijk_right}"}, @alph);
	    my @nonsyn_transitions = grep ($_ eq $transition{$ijk_pos}, @nonsyn_xyz_pos);
	    my @nonsyn_transversions = grep ($_ ne $transition{$ijk_pos}, @nonsyn_xyz_pos);

	    $selfFunc .= " + " . p_branch_var($ijk_pos,$pos);
	    $kFunc .= join ("", map (" + " . p_branch_var($_,$pos), @{$transversions{$ijk_pos}}));
	    $wFunc .= join ("", map (" + " . p_branch_var($_,$pos), @nonsyn_transitions));
	    $kwFunc .= join ("", map (" + " . p_branch_var($_,$pos), @nonsyn_transversions));
	}

	$selfFunc =~ s/^ \+ //;
	$kFunc =~ s/^ \+ //;
	$wFunc =~ s/^ \+ //;
	$kwFunc =~ s/^ \+ //;

	$rFunc = "($selfFunc + $not_k * ($kFunc)"
	    . (length($wFunc) ? " + $not_w * ($wFunc)" : "")
	    . (length($kwFunc) ? " + $k * $not_w * ($kwFunc)" : "")
	    . ")";

	$comment = "unobserved self- & rejected mutations";

    } else {
	my $nDiff = 0;
	my $isTransition;
	for my $pos (0..2) {
	    my $ijk_pos = substr ($ijk, $pos, 1);
	    my $xyz_pos = substr ($xyz, $pos, 1);
	    if ($ijk_pos ne $xyz_pos) {
		$isTransition = $transition{$ijk_pos} eq $xyz_pos;
		$rFunc = p_branch_var ($xyz_pos, $pos);
		++$nDiff;
	    }
	}
	return undef unless $nDiff == 1;
	my $isSynonymous = $aa{$ijk} eq $aa{$xyz};
	$rFunc .= " * " . k_branch_var() unless $isTransition;
	$rFunc .= " * " . w_branch_var() unless $isSynonymous;
	$comment =
	    "$aa{$ijk} -> $aa{$xyz} (" .
	    ($isSynonymous ? "" : "non") .
	    "synonymous " .
	    ($isTransition ? "transition" : "transversion") .
	    ")";

    }

    return (S_branch_var() . " * $rFunc", $comment);
}

# branch-specific variable name generator
sub branch_var { my ($var) = @_; return "(. ${var}_ NODE)" }

sub S_branch_var { return branch_var("S") }
sub p_branch_var { my ($new, $pos) = @_; return branch_var ("p$new$pos") }
sub k_branch_var { return branch_var("k") }
sub w_branch_var { return branch_var("w") }
sub not_k_branch_var { return branch_var("~k") }
sub not_w_branch_var { return branch_var("~w") }
