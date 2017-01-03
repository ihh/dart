#!/usr/bin/env perl -w

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Stockholm::Database;
use Getopt::Long;

my $dssp3 = 0;
my $dssp3gtj = 0;
&GetOptions( "dssp3" => \$dssp3, "dssp3gtj" => \$dssp3gtj);
if (@ARGV != 1) { help(); exit(1)}
my ($dbfile) = @ARGV;

my $db = Stockholm::Database->from_file ($dbfile);
for my $stock (@$db) {
    my $cols = $stock->columns;
    my @features = grep (!exists ($stock->gc->{$_}), keys %{$stock->gr});
    for my $feature (@features) {
	my $cons = "";
	my $gr_feat = $stock->gr->{$feature};
	for (my $col = 0; $col < $cols; ++$col) {
	    my %count;
	    my $total = 0;
	    while (my ($seqname, $gr) = each %$gr_feat) {
		my $c = substr ($stock->seqdata->{$seqname}, $col, 1);
		if ($c ne "" && $c ne '.' && $c ne '-') {
		    # not a gap
		    ++$count{substr($gr,$col,1)};
		    ++$total;
		}
	    }
	    my ($max_c, $max_count) = ('.', undef);
	    for my $c (sort keys %count) {
		my $count = $count{$c};
		($max_c, $max_count) = ($c, $count) if !defined($max_count) || $count > $max_count;
	    }
            if ($feature eq "DSSP")
            {
	      if ($dssp3)
	      {
		if ($max_c =~ /[HGI]/) { $max_c = "H"}
		elsif ($max_c =~ /[BE]/) { $max_c = "E"}
		else { $max_c = "L"}
	      }
	      elsif ($dssp3gtj)
	      {
		if ($max_c eq "H") { }
		elsif ($max_c =~ /[AEP]/) { $max_c = "E"}
		elsif ($max_c =~ /[\.\-\?]/) {$max_c = "."}
		else { $max_c = "L"}
	      }
	    }
	    $cons .= $max_c;
	}
	$stock->gc->{$feature} = $cons;
    }
}

print $db->to_string;

sub help {
    print STDERR <<EOF;

$0: Majority consensus annotation for Stockholm databases

Usage: $0  [options] <Stockholm database>
Options:
   -dssp3 : map DSSP feature to 3 states
   -dssp3gtj : map DSSP feature to 3 states using GTJ method
 
EOF
}
