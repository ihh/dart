#!/usr/local/bin/perl

my %pstr = ( ".max"      => "{\\bf Unguided}  ",
	     ".rind"     => "{\\bf Restricted}",
	     ".rind.max" => "{\\bf Guided}    " );

my %b;
while (<>) {
    next unless /HSM(\d)(\S)(\S*)\s+(\S+)/;
    my ($C, $run, $protocol, $bits) = ($1, $2, $3, $4);
    if (exists $pstr{$protocol}) {
#	warn "protocol=$protocol C=$C run=$run bits=$bits\n";
	$b{$protocol}->{$C}->{$run} = $bits;
    }
}

my @protocol = sort keys %pstr;
foreach my $C (1..4) {
    my $cstr = $C;
    foreach my $protocol (@protocol) {
	my $pstr = $pstr{$protocol};
	print "$cstr & $pstr & ", join (" & ", map (printbits($b{$protocol}->{$C}->{$_}), qw(a b c))), " \\\\\n";
	$cstr =~ s/./ /g;
    }
}

sub printbits {
    my $bits = shift;
    return defined($bits) ? sprintf ("%.5g", $bits) : "n/a";
}
