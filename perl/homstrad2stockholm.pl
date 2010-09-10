#!/usr/bin/perl -w

my $unwanted_prefix = "P1;";  # all homstrad protein names seem to have this at the beginning, for some reason

my $usage = "Usage: $0 <directory> [<directories...>]\n";
$usage .= "Options:\n";
$usage .= " [-annot]  keep annotation\n";
$usage .= "    [-db]  directory is top level of entire database\n";
$usage .= "            (default assumption: directory contains one family)\n";

my $get_terms = 0;
my $db = 0;
foreach my $arg (@ARGV) {
    if ($arg eq "-annot") { $get_terms = 1 }
    elsif ($arg eq "-db") { $db = 1 }
    elsif ($arg =~ /^-/) { die $usage }
    else { push @arg, $arg }
}

unless (@arg) {
    my $cwd = `pwd`;
    chomp $cwd;
    push @arg, $cwd;
}

foreach my $topdir (@arg) {
    if (!$db) {
	convert_dir ($topdir);
    } else {
	local *D;
	opendir D, $topdir;
	my @subdir = grep (!/^\./, readdir(D));
	closedir D;
	foreach my $subdir (@subdir) {
	    convert_dir ("$topdir/$subdir");
	}
    }
}

sub convert_dir
{
    my ($dir) = @_;
    my $dir_name = $dir;
    $dir_name =~ s/\/+$//;
    $dir_name =~ s/.*\///;
    my ($align_file, $annot_file, $id) = ("$dir/$dir_name.ali", "$dir/$dir_name.tem", $dir_name);
    unless (-e $align_file) {
	warn "Warning: couldn't find alignment file '$align_file'\n";
	return;
    }

    local *ALIGN;
    open ALIGN, "<$align_file";
    my (%seq, %term, %annot);
    my @name;
    my ($name, $term, $class, $family);
    my $length = 0;
    while (<ALIGN>) {
	if (/^\s*>\s*(\S+)/) {
	    $name = $1;
	    $name =~ s/^$unwanted_prefix//;
	    push @name, $name;
	    $seq{$name} = "";
	    $term{$name} = [];
	    $annot{$name} = {};
	    $length = max ($length, length $name);
	    $_ = <ALIGN>;   # discard first line after ">"
	} elsif (/\S+/ && defined $name) {
	    s/\s//g;
	    $seq{$name} .= $_;
	  } elsif (/^C; class: (.+)/)
	  {
	    $class = $1;
	    chomp($class);
	  }
	  elsif (/^C; family: (.+)/)
	  {
	    $family = $1;
	    chomp($family);
	  }
    }
    close ALIGN;

    if ($get_terms) {
	local *TERM;
	open TERM, "<$annot_file";
	while (<TERM>) {
	    if (/^\s*>\s*(\S+)/) {
		$name = $1;
		$name =~ s/^$unwanted_prefix//;
		if (!defined ($seq{$name})) {
		    undef $term;
		    next;
		}

		$term = <TERM>;
		chomp $term;
		if ($term eq "sequence") {
		    undef $term;
		    next;
		}

		$term =~ s/ /_/g;
		$term =~ s/[^A-Za-z_]//g;
		push @{$term{$name}}, $term;
		$annot{$name}->{$term} = "";
		$length = max ($length, length (tag ($name, $term)));

	    } elsif (/\S+/ && defined $term) {
		s/\s//g;
		$annot{$name}->{$term} .= $_;
	    }
	}
	close TERM;
    }

    print "# STOCKHOLM 1.0\n";
    print "#=GF ID $id\n";
    print "#=GF CLASS $class\n";
    print "#=GF FAMILY $family\n";
    foreach my $n (@name) {
	my $seq = $seq{$n};
	chop $seq;  # remove trailing "*"
	print &rowname ($n, $length), $seq, "\n";
	foreach my $t (@{$term{$n}}) {
	    my $annot = $annot{$n}->{$t};
	    chop $annot;   # remove trailing "*"
	    print &rowname (tag($n,$t), $length), $annot, "\n";
	}
    }
    print "//\n";
}

sub tag {
    my ($name, $term) = @_;
    return "#=GR $name $term";
}

sub rowname {
    my ($string, $width) = @_;
    return $string . " " x (1 + $width - length($string));
}

sub max {
    my ($a, $b) = @_;
    return $a > $b ? $a : $b;
}
