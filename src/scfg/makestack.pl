#!/usr/local/bin/perl -w

my $usage = "Usage: $0 <energy file> <stack> <hairpin> <interior> <bulge>\n";
die $usage if @ARGV != 5;
my ($energy, $stack, $hairpin, $interior, $bulge) = @ARGV;

my $A = 4;  # alphabet size
my $mul = -1000;  # energy multiplier
my $inf = "InfinityScore";

local *ENERGY;
open ENERGY, "<$energy" or die $!;

# Declarations and body of output
my @decl;
my @body;

# Array declarations
print_array_decl ($stack);
print_array_decl ($hairpin);
print_array_decl ($interior);
push @decl, "\n";

# Stacking energies
print_energy_array ($stack);

# Hairpin loop terminal mismatch energies
while (<ENERGY>) { last if /Terminal mismatch stacking energies .hairpin loops./ }
print_energy_array ($hairpin);

# Interior loop terminal mismatch energies
while (<ENERGY>) { last if /Terminal mismatch stacking energies .interior loops./ }
print_energy_array ($interior);

# Loop size destabilising energies
while (<ENERGY>) { last if /Loop destabilizing energies/ }
print_loop_vectors ($interior, $bulge, $hairpin);

close ENERGY;

# Output and exit
print @decl, "\n";
print @body;
exit;

sub print_array_decl {
    my ($name) = @_;
    push @decl, "array2d<array2d<Score> > $name ($A, $A, array2d<Score> ($A, $A, -$inf));\n";
}

sub print_vector_decl {
    my ($name) = @_;
    push @decl, "vector<Score> $name ($A, $A, array2d<Score> ($A, $A, -$inf));\n";
}

sub print_energy_array {
    my ($name) = @_;
    while (<ENERGY>) { last if /STACKING ENERGIES/ }  # regexp for recognising stacking energies block
    for (my $lastx = 0; $lastx < $A; ++$lastx) {
	while (<ENERGY>) { last if /3\' <-- 5\'/ }  # regexp for recognising data sub-block (4 rows)
	my @block = map ([], 1..$A);
	for (my $thisx = 0; $thisx < $A; ++$thisx) {
	    $_ = <ENERGY>;
	    my @f = split;
	    for (my $lasty = 0; $lasty < $A; ++$lasty) {
		push @{$block[$lasty]}, [@f[$A*$lasty..$A*$lasty+$A-1]];
	    }
	}
	for (my $lasty = 0; $lasty < $A; ++$lasty) {
	    my $newline;
	    for (my $thisx = 0; $thisx < $A; ++$thisx) {
		for (my $thisy = 0; $thisy < $A; ++$thisy) {
		    my $energy = $block[$lasty]->[$thisx]->[$thisy];
		    if ($energy ne ".") {
			push @body, $name . "($lastx,$lasty)($thisx,$thisy) = " . $mul*$energy . ";\n";
			$newline = 1;
		    }
		}
	    }
	    push @body, "\n" if $newline;
	}
    }
}

sub print_loop_vectors {
    my @name = @_;
    my $suffix = "_loop";
    my $tmp = "_tmp";
    @name = map ("$_$suffix", @name);
    my @data = map ([], @name);
    while (<ENERGY>) { last if /^-----/ }  # regexp for recognising start of loop data
    while (<ENERGY>) {
	last unless /^\s*\d/;
	my ($size, @row) = split;
	die "Missing data" if @row != @data;
	for (my $col = 0; $col < @row; ++$col) { push @{$data[$col]}, $row[$col] eq "." ? "-$inf" : $mul*$row[$col] }
    }
    for (my $n = 0; $n < @name; ++$n) {
	my $loop_size = @{$data[$n]};
	push @decl, "vector<Score> $name[$n];\n";
	push @body, "Score $name[$n]$tmp\[" . ($loop_size+1) . "\] = { -$inf, " . join(", ", @{$data[$n]}) . " };\n";
	push @body, "$name[$n] = vector<Score> ($name[$n]$tmp, $name[$n]$tmp + " . ($loop_size+1) . ");\n";
	push @body, "\n";
    }
}
