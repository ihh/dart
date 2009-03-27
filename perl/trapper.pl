#!/usr/local/bin/perl -w

$usage = "$0: TRee ALignment Picture PlottER\n";

sub pathonly { my $path = shift; $path =~ s/(.*)\/([^\/]+)$/$1/; return $path }
sub nameonly { my $path = shift; $path =~ s/(.*)\/([^\/]+)$/$2/; return $path }
sub dosys { my $command = shift; warn $command, "\n"; system $command }

my $progname = nameonly ($0);
my $progpath = pathonly ($0);

my $usage = "\nUsage: $progname <treefile> <alignfile>\n\n";

# parse command-line options

my @argv;
while (@ARGV) {
    $arg = shift @ARGV;
    if ($arg =~ /^-/) {
	if ($arg eq "-dummy") { die "dummy option invoked!\n" }
	else { die "$usage\nUnknown option: $arg\n\n" }
    } else { push @argv, $arg }
}
die $usage unless @argv == 2;

my ($treefile, $alignfile) = @argv;

# read in tree & alignment

local *TREE;
open TREE, "<$treefile" or die "Couldn't open tree: $!";
my $tree = join ("", <TREE>);
close TREE;
$tree =~ s/\s//g;

local *ALIGN;
open ALIGN, "<$alignfile" or die "Couldn't open alignment: $!";
my @rowName;
my @rowData;
while (<ALIGN>) {
    my ($name, $data) = split;
    push @rowName, $name;
    push @rowData, $data;
}
close ALIGN;

# parse tree
my @nodes;
my $root = parse ($tree, \@nodes);

# extract nice arrays
my @parent = map (defined($_->{'parent'}) ? $_->{'parent'} : -1, @nodes);
my @lenToParent = map (defined($_->{'lenToParent'}) ? $_->{'lenToParent'} : -1, @nodes);
my @child = map (defined($_->{'child'}) ? $_->{'child'} : [], @nodes);
my @name = map (defined($nodes[$_]->{'name'}) ? $nodes[$_]->{'name'} : "\#$_", 0..@nodes-1);
my %name2node = map (($name[$_] => $_), 0..@name-1);
my @leaf = grep (defined($nodes[$_]->{'leaf'}), 0..@nodes-1);
my @internal = grep (defined($nodes[$_]->{'internal'}), 0..@nodes-1);

my @ancestors = map ([], 0..@nodes-1);
foreach my $node (0..@nodes-1) {
    my $parent = $parent[$node];
    if ($parent >= 0) {
	push @{$ancestors[$node]}, @{$ancestors[$parent]}, $parent;
    }
}

my @descendants = map ([], 0..@nodes-1);
for (my $node = @nodes-1; $node >= 0; --$node) {
    foreach my $child (@{$child[$node]}) {
	push @{$descendants[$node]}, $child, @{$descendants[$child]};
    }
}

# parse alignment
my @row2node = map (undef, @rowName);
my @node2row = map (undef, @node);
foreach my $row (0..@rowName-1) {
    my $name = $rowName[$row];
    my $node = handelNode ($name);
    if (defined $node) {
	if (defined $node2row[$node]) {
	    die "Alignment rows '", $rowName[$node2row[$node]], "' and '$name' appear to refer to the same node\n";
	}
	$row2node[$row] = $node;
	$node2row[$node] = $row;
    } else {
	die "Alignment row with name '$name' could not be found in the tree";
    }
}

# check for completeness
my @undefLeaf = grep (!defined($node2row[$_]), @leaf);
if (@undefLeaf) {
    die "The following leaf nodes were not found in the alignment: ", join(" ",map($name[_],@undefLeaf)), "\n";
}
my $gotAllNodes = grep (!defined(), @node2row) == 0;

# subroutine to find most recent common ancestor of two nodes
sub mrca {
    my ($node1, $node2) = @_;
    my %isAncOf2 = map (($_=>1), @{$ancestors[$node2]});
    while ($node1 >= 0 && !defined $isAncOf2{$node1}) {
	$node1 = $parent[$node1];
    }
    return $node1;
}

# subroutine to get node from Handel descriptor (A::B means "mrca of A & B")
sub handelNode {
    my ($desc) = @_;
    if ($desc =~ /^(\S+)::(\S+)$/) {
	my ($desc1, $desc2) = @_;
	return mrca ($name2node{$desc1}, $name2node{$desc2});
    }
    return $name2node{$desc};
}

# tree parsing subroutine; old legacy code (it works, goddammit)
sub parse {
    my ($nodeExpr, $nodesArrayRef) = @_;
    my $nodeRef;
    $nodeExpr =~ s/\s//g;
    if ($nodeExpr =~ s/^\((.*)\)$/$1/) {
	$nodeRef = { 'root'=>1, 'internal'=>1, 'child'=>[], 'len'=>[] };
	# replace top-level commas with ampersands so split command can separate out child branches
	my $sep = "&";
	if ($nodeExpr =~ /$sep/) { die "Found $sep character in $nodeExpr" }
	my ($i,$nest);
	for ($i = 0; $i < length($nodeExpr); $i++) {
	    my $c = substr ($nodeExpr, $i, 1);
	    if    ($c eq "(") { $nest++ }
	    elsif ($c eq ")") { $nest-- }
	    elsif ($c eq "," && $nest==0) { substr($nodeExpr,$i,1) = $sep }
	}
	# process each child branch
	my @expr = split /$sep/, $nodeExpr;
	foreach (@expr) {
	    /^(.*):(\-?[\d\.]+)$/ or die "Unrecognised branch:length format in $_";
	    my $len = $2;
	    my $child = parse ($1, $nodesArrayRef);
	    delete $child->{'root'};  # we know the kid isn't a root
	    $child->{'parent'} = $nodeRef;  # we know who the kid's parent is: us
	    $child->{'lenToParent'} = $len;
	    push @{$nodeRef->{'child'}}, $child;
	    push @{$nodeRef->{'len'}}, $len;
	}
	if (@{$nodeRef->{'child'}}!=2) { warn "Size of child array is ".scalar(@{$nodeRef->{'child'}}) }
    } else {
	$nodeRef = { 'root'=>1, 'leaf'=>1, 'name'=>$nodeExpr };
    }
    push @{$nodesArrayRef}, $nodeRef;
    $nodeRef->{'index'} = scalar @$nodesArrayRef;
    return $nodeRef;
}

