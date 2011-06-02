#!/usr/bin/perl -w

my ($ontology, $d, $id, @is_a, @name);
while (<>) {
    if (/^ontology: (.*)/) { $ontology = $1 }
    elsif (/^\[([^\]]+)\]/) { $d = $1 }
    elsif (defined($d) && $d eq "Term") {
	if (/^id: (\S+)/) { $id = $1 }
	elsif (/^is_a: (\S+)/) { push @is_a, ["'$id", "'$1"]; }
	elsif (/^name: (.*)/) { push @name, ["'$id", "\"$1\""]; }
    }
}

print_rel ("$ontology:is_a", @is_a);
print_rel ("$ontology:name", @name);

sub print_rel {
    my ($tag, @rel) = @_;
    print "(define \%$tag\n  (\%rel ()";
    for my $rel (@rel) {
	print "\n    ((@$rel))";
    }
    print "))\n";
}
