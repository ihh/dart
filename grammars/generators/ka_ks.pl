#!/usr/bin/perl -w

use DartSexpr;
use PhyloGram::Dna;

# config
my $ks_csv = ".1,5,2";
my $ratio_csv = ".5,1,2";
my $rate_variable_prefix = 'R';
my $params_file;

# usage
my $usage = "\nUsage: $0 [-ks <csv>] [-ratio <csv>] [-params <file>]\n";
$usage .= "where <csv> is a list of comma-separated values\n";
$usage .= "Defaults are \"-ks $ks_csv -ratio $ratio_csv\"\n\n";

# parse args
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-ks") {
	defined ($ks_csv = shift) or die $usage;
    } elsif ($arg eq "-ratio") {
	defined ($ratio_csv = shift) or die $usage;
    } elsif ($arg eq "-params") {
	defined ($params_file = shift) or die $usage;
    } elsif ($arg =~ /^-/) {
	die $usage, "Unknown option: $arg";
    } else {
	push @argv, $arg;
    }
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

# grammar
my $grammar = PhyloGram::Dna->new;
$grammar->grammar->parametric;

# range of values for Ka & Ks
my @ks = split /,/, $ks_csv;
my @ratio = split /,/, $ratio_csv;

# annotation for Ka & Ks
my $ks_annot_tag = "Ks";
my $ratio_annot_tag = "Ka/Ks";

# create table of Ka-Ks pairs
my @ka_ks;
for (my $i = 0; $i < @ks; ++$i) {
    my $ks = $ks[$i];
    my $ks_annot = length($ks) <= 3 ? $ks : chr (ord('a') + $i);
    $ks_annot = "_" x (3-length($ks_annot)) . $ks_annot;
    for (my $j = 0; $j < @ratio; ++$j) {
	my $ratio = $ratio[$j];
	my $ka = $ratio * $ks;
	$ka =~ s/^(\-?)0/$1/g;  # remove leading "0" from ka if <1
	my $ratio_annot = length($ratio) <= 3 ? $ratio : chr (ord('a') + $j);
	$ratio_annot = "_" x (3-length($ratio_annot)) . $ratio_annot;
	push @ka_ks, [$ka, $ks, $ratio_annot, $ks_annot];
    }
}

# load rate vars from file?
warn "[creating variables]\n";
my @cod = grep ($aa{$_} ne "X", keys %aa);
my %got_rate;   # hash to record which rate vars are present
if ($params_file) {
    my $params_file_sexpr = DartSexpr->from_file ($params_file);
    my $params_sexpr = $params_file_sexpr->grammar->find ('params');
    $grammar->add ($params_sexpr);
    $params_sexpr->grep_child (sub {
	my ($child) = @_;
	my $param = $child->has_tag ? $child->tag : $child->[0]->tag;
	$got_rate{$param} = 1;
    });
}

# initialise/check rate vars
my $init_rate = 0.1 / 60;   # low initial rate, helps EM convergence
my %rate;
foreach my $src_cod (@cod) {
    my $src_rate = $rate{$src_cod} = {};
    foreach my $dest_cod (@cod) {
	# only allow single-nucleotide changes
	next if $src_cod eq $dest_cod;
	next unless
	    (substr($src_cod,0,1) eq substr($dest_cod,0,1) && substr($src_cod,2,1) eq substr($dest_cod,2,1))
	    || (substr($src_cod,0,2) eq substr($dest_cod,0,2))
	    || (substr($src_cod,1,2) eq substr($dest_cod,1,2));
	# create var name
	my $var = "${rate_variable_prefix}_${src_cod}_${dest_cod}";
	$src_rate->{$dest_cod} = $var;
	# create grammar variable
	if (defined $params_file) {
	    die "Parameter $var undefined" unless $got_rate{$var};
	} else {
	    my $pgroup = $grammar->new_rate ($var, $init_rate);
	}
    }
}

# initial distribution over codons
my @init_cod = map ("p_$_", @cod);
if (!defined $params_file) {
    $grammar->new_pgroup (@init_cod);
}

# create chains
my @term;
my @prefix;
warn "[creating chains]\n";
foreach my $ka_ks (@ka_ks) {
    my ($ka, $ks, $ks_annot, $ratio_annot) = @$ka_ks;
    warn "[ka=$ka ks=$ks]\n";

    my $prefix = "ka${ka}_ks${ks}";
    my $term = [map ("${prefix}_pos$_", 1..3)];
    push @prefix, $prefix;
    push @term, $term;

    my $chain = $grammar->new_empty_chain (@$term);
    $chain->update_policy ("parametric");

    for (my $i = 0; $i < @cod; ++$i) {
	my $cod = $cod[$i];
	my $init_cod = $init_cod[$i];
	$chain->initial ($cod, $init_cod);
    }

    my %mutate_hash;
    while (my ($src_cod, $src_rate) = each %rate) {
	while (my ($dest_cod, $var) = each %$src_rate) {
	    my $mul = $aa{$src_cod} eq $aa{$dest_cod} ? $ks : $ka;
	    next unless $mul;
	    $chain->mutate ($src_cod, $dest_cod, $mul == 1 ? $var : "$mul * $var", \%mutate_hash);
	}
    }
}

# create start state
my $start = "START";
my $tprob = 1 / (@term + 1);
$grammar->add_end_transition ($start)->prob ($tprob);

# create states
warn "[creating states]\n";
for (my $i = 0; $i < @ka_ks; ++$i) {
    my $state = $prefix[$i];
    my $term = $term[$i];
    my ($ka, $ks, $ks_annot, $ratio_annot) = @{$ka_ks[$i]};
    warn "[ka=$ka ks=$ks]\n";

    my $emit = $grammar->add_emission ("@$term $state'");
    for (my $pos = 0; $pos < 3; ++$pos) {
	$emit->add (["annotate", ["row", $ks_annot_tag], ["column", $term->[$pos]], ["label", substr ($ks_annot, $pos, 1)]]);
	$emit->add (["annotate", ["row", $ratio_annot_tag], ["column", $term->[$pos]], ["label", substr ($ratio_annot, $pos, 1)]]);
    }

    $grammar->add_transition ($start, $state)->prob ($tprob);
    foreach my $dest (@prefix) {
	$grammar->add_transition ($state, $dest)->prob ($tprob);
    }
    $grammar->add_end_transition ($state)->prob ($tprob);
}

# print
print $grammar->to_string, "\n";
