
=head1 NAME

    Newick.pm

    =head1 SYNOPSIS

    Lightweight Perl module encapsulating a rooted phylogenetic tree in Newick format.

    Nodes are sorted in preorder (parents before children).

    =head1 METHODS

=cut

package Newick;
use strict;

# bare-bones Newick tree class

=head2 new

    my $tree = Newick->new();

Creates an empty Newick object.

=cut

sub new {
    my ($class) = @_;
    my $self = {
	parent => [],
	node_name => [],
	branch_length => [],
    };
    bless $self, $class;
    return $self;
}

=head2 parent

    my $indexOfParent = $tree->parent->[$indexOfNode]

    Returns a reference to an array giving the indices of each node's parent.

    The root node always has index 0, and its (dummy) parent index is -1.

=head2 node_name

    my $nodeName = $tree->node_name->[$indexOfNode]

    Returns a reference to an array of node names.

=head2 branch_length

    my $branchLength = $tree->branch_length->[$indexOfNode]

    Returns a reference to an array of branch lengths.
    The N'th entry in this array is the branch length from the N'th node to its parent.

=cut

sub parent { my ($self) = @_; return $self->{'parent'} }
sub node_name { my ($self) = @_; return $self->{'node_name'} }
sub branch_length { my ($self) = @_; return $self->{'branch_length'} }

=head2 nodes

    my $nodes = $tree->nodes

    Returns the number of nodes in the tree.

=cut

sub nodes {
    my ($self) = @_;
    return @{$self->parent} + 0;
}

=head2 children

    my @kids = $tree->children ($indexOfNode)

    Returns the list of children of a given node.

=cut

sub children {
    my ($self, $node) = @_;
    return grep ($self->parent->[$_] == $node, $node+1..$self->nodes-1);
}

=head2 siblings

    my @sibs = $tree->siblings ($indexOfNode)

    Returns the list of siblings of a given node.

=cut

sub siblings {
    my ($self, $node) = @_;
    my $parent = $self->parent->[$node];
    return grep ($self->parent->[$_] == $parent && $_ != $node, $parent..$self->nodes-1);
}

=head2 add_node

    my $indexOfNewNode = $tree->add_node ($indexOfParent, $newNodeName, $branchLengthToParent)

    Appends a new node & returns the index.

    $indexOfParent should be -1 for the root node.

=cut

# add_node appends a new node & returns the index
sub add_node {
    my ($self, $parent, $node_name, $branch_length) = @_;
    my $node = $self->nodes;
    $parent = -1 unless defined $parent;
    if ($parent >= $node) {
	die "Nodes must be added in preorder (parents before children)";
    }
    push @{$self->parent}, $parent;
    push @{$self->node_name}, $node_name;
    push @{$self->branch_length}, $branch_length;
    return $node;
}

=head2 insert_node

    $tree->insert_node ($indexOfNewNode, $indexOfParent, $newNodeName, $branchLengthToParent)

    Adds a new node with a given index.

=cut

# insert_node adds a new node before $node & updates the parent array
sub insert_node {
    my ($self, $node, $parent, $node_name, $branch_length) = @_;
    if ($parent >= $node) {
	die "Nodes must be added in preorder (parents before children)";
    }
    my @map_parent = map ($_ + ($_ >= $node ? 1 : 0), @{$self->parent});
    @{$self->parent} = (@map_parent[0..$node-1],
			$parent,
			@map_parent[$node..@map_parent-1]);
    splice @{$self->node_name}, $node, 0, $node_name;
    splice @{$self->branch_length}, $node, 0, $branch_length;
}

=head2 delete_node

    $tree->delete_node ($indexOfNode);
$tree->tidy;

$tree->delete_node (@listOfNodeIndices);
$tree->tidy;

Deletes a node with a given index (or a list of nodes).

    Deleting nodes can leave the tree containing redundant nodes (i.e. internal nodes with two neighbors) or unnamed leaf nodes.
    Call $tree->tidy() after $tree->delete_node in order to "clean up" such trees.

=cut

# delete_node deletes leaf node $node & updates the parent array
sub delete_node {
    my ($self, @del) = @_;
    my %is_leaf = map (($_ => 1), $self->leaves);
    if (grep !$is_leaf{$_}, @del) {
	die "Attempt to delete internal node";
    }
    my %del = map (($_ => 1), @del);
    my @nodes_to_keep = grep (!$del{$_}, 0..$self->nodes-1);
    my %new_index = (-1 => -1,
		     map (($nodes_to_keep[$_] => $_), 0..@nodes_to_keep-1));
    for my $node (@nodes_to_keep) { my $parent = $self->parent->[$node]; if (!exists ($new_index{$parent})) { die "Node $node has parent $parent which is unmapped (while deleting nodes: @del)" } }
    @{$self->parent} = map ($new_index{$self->parent->[$_]}, @nodes_to_keep);
    @{$self->node_name} = @{$self->node_name} [@nodes_to_keep]; 
    @{$self->branch_length} = @{$self->branch_length} [@nodes_to_keep];
}

=head2 insert_node_above

    my $indexOfNewNode = $tree->insert_node_above ($indexOfNode, $newNodeName, $branchLengthToParent)

    Adds a new node by splitting the branch between a node and its parent.

    Returns the index of the new node (actually, this is redundant:
    the new node is inserted immediately before the existing one,
    so $indexOfNewNode = $indexOfNode).

=cut

# insert_node_above is like insert_node, but it splits the branch leading to $node
sub insert_node_above {
    my ($self, $node, $node_name, $branch_length) = @_;
    $self->insert_node ($node, $self->parent->[$node], $node_name, $branch_length);
    $self->parent->[$node + 1] = $node;
    return $node;
}

=head2 remove_redundant_nodes

    my $deleted = $tree->remove_redundant_nodes

    Eliminates all internal nodes of degree 2 (excluding the root).

    Returns true if any nodes were deleted.

=cut

sub remove_redundant_nodes {
    my ($self) = @_;
    my @nodes_to_delete;
    for (my $node = 1; $node < $self->nodes; ++$node) {
	my @kids = $self->children ($node);
	if (@kids == 1) {
	    $self->parent->[$kids[0]] = $self->parent->[$node];
	    $self->branch_length->[$kids[0]] += $self->branch_length->[$node];
	    push @nodes_to_delete, $node;
	}
    }
    $self->delete_node (@nodes_to_delete);
    return @nodes_to_delete > 0;
}


=head2 remove_redundant_branches

    my $deleted = $tree->remove_redundant_branches

    Eliminates all branches of length zero.

    Returns true if any branches were deleted.

=cut

sub remove_redundant_branches {
    my ($self) = @_;
    my @nodes_to_delete;
    for (my $node = 1; $node < $self->nodes; ++$node) {
	if (!defined($self->branch_length->[$node]) || $self->branch_length->[$node] == 0) {
	    my @kids = $self->children ($node);
	    for my $child (@kids) { $self->parent->[$child] = $self->parent->[$node] }
	    push @nodes_to_delete, $node;
	}
    }
    $self->delete_node (@nodes_to_delete);
    return @nodes_to_delete > 0;
}

=head2 delete_unnamed_leaf_nodes

    $tree->delete_unnamed_leaf_nodes

    Deletes all unnamed leaf nodes.

    Returns true if any nodes were deleted.

=cut

sub delete_unnamed_leaf_nodes {
    my ($self) = @_;
    my @nodes_to_delete = grep ($self->children($_)==0
				&& (!defined($self->node_name->[$_]) || length($self->node_name->[$_])==0),
				0..$self->nodes-1);
    $self->delete_node (@nodes_to_delete);
    return @nodes_to_delete > 0;
}


=head2 tidy

    $tree->tidy

    Repeatedly calls remove_redundant_nodes, remove_redundant_branches and delete_unnamed_leaf_nodes until there are no more nodes to remove.

=cut

sub tidy {
    my ($self) = @_;
    while (1) {
	my $rbn = $self->remove_redundant_branches;
	my $rrn = $self->remove_redundant_nodes;
	my $duln = $self->delete_unnamed_leaf_nodes;
	last unless $rbn || $rrn || $duln;
    }
}

=head2 subtree

    my $subtree = $tree->subtree ($node)

    Returns (as a Newick tree object) the subtree rooted at a particular node.

=cut

# subtree method
sub subtree {
    my ($self, $node) = @_;
    $node = 0 unless defined $node;
    # map from old node index to new node index
    my @new_index = map (undef, 1..$self->nodes);  # by default the old parent index is undef, indicating that it's not in the tree
    # build the subtree Newick object
    my $subtree = ref($self)->new;
    # populate the new object by scanning through the existing tree in a preorder traversal
    # can skip node 0, because it's the root, and therefore only in the subtree if $node==0
    for (my $n = 0; $n < $self->nodes; ++$n) {
	my $new_index_of_parent =  # figure out the parent index of this node in the subtree
	    $n == $node  # is this the new subtree root?
	    ? -1  # if yes, then the parent index is -1
	    : ($n == 0  # otherwise, is this the old subtree root?
	       ? undef  # if so, and it's not the new subtree root, then it has no parent in the old tree
	       : $new_index[$self->parent->[$n]]);  # look up the new parent index in our map
	if (defined $new_index_of_parent) {
	    # add a node to the subtree, and record its index
	    $new_index[$n] = $subtree->add_node ($new_index_of_parent, $self->node_name->[$n], $self->branch_length->[$n]);
	}
    }
    # return
    return $subtree;
}

=head2 ancestors

    my $ancestors = $tree->ancestors ($node)

    Returns the ancestors of a node in the tree,
    sorted in postorder (children before parents).

=cut

sub ancestors {
    my ($self, $node) = @_;
    my @anc;
    while ($node != -1) {
	push @anc, $node;
	$node = $self->parent->[$node];
    }
    return @anc;
}

=head2 descendants

    my $descendants = $tree->descendants ($node)

    Returns the descendants of a node in the tree,
    sorted in preorder (parents before children).

=cut

sub descendants {
    my ($self, $node) = @_;
    my @desc;
    my @q = ($node);
    while (@q) {
	my @c = $self->children(shift @q);
	push @desc, @c;
	push @q, @c;
    }
    return @desc;
}

=head2 lca

    my $ancestor = lca(@nodes);

Returns the lowest common ancestor of a set of nodes.

=cut

sub lca{
    my ($self, @nodes) = @_; 
    my @ancestors=(); 

    foreach my $node (@nodes){
	my @path_ancest = $self->ancestors($node);
	push @ancestors, [@path_ancest];
    }
    
    
    my $compare_path=shift @ancestors;
    my @matches=();

    foreach my $lca (@{$compare_path}){
	
	foreach my $path (@ancestors){
	    
	    foreach my $element (@{$path}){

		if ($element == $lca){
		    push @matches, $element;
		}
	    }
	}
	
	if (scalar(@matches)==scalar(@ancestors) && $matches[0]==$lca){
	    return $matches[0];
	}
	@matches=();
    }	
    
    
}

=head2 leaves

    my @leaves = $tree->leaves

    Returns the set of leaf nodes in the tree.

=cut

sub leaves {
    my ($self) = @_;
    my %is_leaf = map (($_ => 1), 0..$self->nodes - 1);
    for my $parent (@{$self->parent}) {
	delete $is_leaf{$parent} if exists $is_leaf{$parent};
    }
    return sort {$a <=> $b} keys %is_leaf;
}

=head2 distance

    my $distance = $tree->distance ($node1, $node2)

    Returns the distance (total branch length) between two nodes in the tree.

=cut

sub distance {
    my ($self, $node1, $node2) = @_;
    my @anc1 = $self->ancestors ($node1);
    my @anc2 = $self->ancestors ($node2);
#    warn "n1=$node1 a1=(@anc1) n2=$node2 a2=(@anc2)";
    my %anc1_index = map (($anc1[$_] => $_), 0..@anc1-1);
    my $i2;
    for ($i2 = 0; $i2 < @anc2 && !exists($anc1_index{$anc2[$i2]}); ++$i2) { }
    my $i1 = $i2 == @anc2 ? (@anc1-1) : $anc1_index{$anc2[$i2]};
    my @path = (@anc1[0..$i1-1],
		@anc2[0..$i2-1]);
#    warn "path=(@path)";
    my $len = 0;
    for my $node (@path) { $len += $self->branch_length->[$node] }
    return $len;
}

=head2 find_nodes

    my @nodes = $tree->find_nodes ($name)

    Returns the list of nodes with the given name.

=cut

sub find_nodes {
    my ($self, $name) = @_;
    return grep ($self->node_name->[$_] eq $name, 0..$self->nodes - 1);
}

=head2 find_node

    my $node = $tree->find_node ($name)

    Returns the unique node with the given name,
    throwing an error if there are zero or more than two such nodes.

=cut

sub find_node {
    my ($self, $name) = @_;
    my @nodes = $self->find_nodes ($name);
    die "No nodes named $name" if @nodes == 0;
    die @nodes+0, " nodes named $name" if @nodes > 1;
    return $nodes[0];
}

=head2 parse

    my $tree = Newick->parse ($newickFormatString)
    my $tree = Newick->from_string ($newickFormatString)

    Creates a new Newick object and initializes it from a Newick-format string.

    from_string is a synonym for parse.

=cut

sub from_string { return parse (@_) }

# parser ripped out of BioPerl by IH on 4/4/2007
# modified to ignore NHX extensions & deal with multiple trees on 6/22/2015
sub parse {
    my ( $class, $string ) = @_;
    my @s = $class->_split ($string);
    my @trees;
    for my $s (@s) {
	my $tree = $class->new;
	my $remainder = $s;
	my $token;
	my @tokens;
	while ( ( $token, $remainder ) = $tree->_next_token( $remainder ) ) {
	    last if ( ! defined $token || ! defined $remainder );
	    push @tokens, $token;
	}
	my $i;
	for ( $i = $#tokens; $i >= 0; $i-- ) {
	    last if $tokens[$i] eq ';';
	}
	my $root = $tree->add_node(-1);
	$tree->_parse_node_data( $root, @tokens[ 0 .. ( $i - 1 ) ] );
	$tree->_parse_clade( $tree, $root, @tokens[ 0 .. ( $i - 1 ) ] );
	push @trees, $tree;
    }
    return wantarray() ? @trees : shift(@trees);
}

=head2 from_file

    my $tree = Newick->from_file ($filename)

    Creates a new Newick object and initializes it from a Newick-format file.

=cut

# from_file
sub from_file {
    my ($class, $filename) = @_;
    local *FILE;
    open FILE, "<$filename" or die "Couldn't open $filename: $!";
    my @file = <FILE>;
    close FILE;
    return $class->parse (join ("", @file));
}

=head2 to_string

    print $tree->to_string()

    Renders a Newick tree object as a (Newick-format) string.

=cut

# to_string method (recursive)
sub to_string {
    my ($self, $node) = @_;
    $node = 0 unless defined $node;
    my @children = $self->children ($node);
    my $node_text;
    if (@children) {
	$node_text = '(' . join (',', map ($self->to_string($_), @children)) . ')';
    } else {
	$node_text = "";
    }
    my $nn = $self->node_name->[$node];
    $node_text .= $nn if defined $nn;

    my $bl = $self->branch_length->[$node];
    $node_text .= ":$bl" if defined($bl) && ($node > 0 || $bl > 0);

    $node_text .= ';' if $node == 0;

    return $node_text;
}

# original credits for parser code:

# ACKNOWLEDGEMENTS
# The author would like to thank Jason Stajich for many ideas borrowed
# from BioPerl L<http://www.bioperl.org>, and CIPRES
# L<http://www.phylo.org> and FAB* L<http://www.sfu.ca/~fabstar>
# for comments and requests.

# COPYRIGHT & LICENSE
# Copyright 2005 Rutger A. Vos, All Rights Reserved. This program is free
# software; you can redistribute it and/or modify it under the same terms as Perl
# itself.

sub _split {
    my ( $self, $string ) = @_;
    my ( $QUOTED, $COMMENTED ) = ( 0, 0 );
    my $decommented = '';
    my @trees;
  TOKEN: for my $i ( 0 .. length( $string ) ) {
      if ( ! $QUOTED && ! $COMMENTED && substr($string,$i,1) eq "'" ) {
	  $QUOTED++;
      }
      elsif ( ! $QUOTED && ! $COMMENTED && substr($string,$i,1) eq "[" ) {
	  $COMMENTED++;
	  next TOKEN;
      }
      elsif ( ! $QUOTED && $COMMENTED && substr($string,$i,1) eq "]" ) {
	  $COMMENTED--;
	  next TOKEN;
      }
      elsif ( $QUOTED && ! $COMMENTED && substr($string,$i,1) eq "'" && substr($string,$i,2) ne "''" ) {
	  $QUOTED--;
      }
      $decommented .= substr($string,$i,1) unless $COMMENTED;
      if ( ! $QUOTED && ! $COMMENTED && substr($string,$i,1) eq ';' ) {
	  push @trees, $decommented;
	  $decommented = '';
      }
  }
    return @trees;
}

sub _parse_clade {
    my ( $self, $tree, $root, @tokens ) = @_;
    my ( @clade, $depth, @remainder );
  TOKEN: for my $i ( 0 .. $#tokens ) {
      if ( $tokens[$i] eq '(' ) {
	  if ( not defined $depth ) {
	      $depth = 1;
	      next TOKEN;
	  }
	  else {
	      $depth++;
	  }
      }
      elsif ( $tokens[$i] eq ',' && $depth == 1 ) {
	  my $node = $tree->add_node ($root);
	  $self->_parse_node_data( $node, @clade );
	  $self->_parse_clade( $tree, $node, @clade );
	  @clade = ();
	  next TOKEN;
      }
      elsif ( $tokens[$i] eq ')' ) {
	  $depth--;
	  if ( $depth == 0 ) {
	      @remainder = @tokens[ ( $i + 1 ) .. $#tokens ];
	      my $node = $tree->add_node ($root);
	      $self->_parse_node_data( $node, @clade );
	      $self->_parse_clade( $tree, $node, @clade );
	      last TOKEN;
	  }
      }
      push @clade, $tokens[$i];
  }
}

sub _parse_node_data {
    my ( $self, $node, @clade ) = @_;
    my @tail;
  PARSE_TAIL: for ( my $i = $#clade; $i >= 0; $i-- ) {
      if ( $clade[$i] eq ')' ) {
	  @tail = @clade[ ( $i + 1 ) .. $#clade ];
	  last PARSE_TAIL;
      }
      elsif ( $i == 0 ) {
	  @tail = @clade;
      }
  }
    # name only
    if ( scalar @tail == 1 ) {
        $self->node_name->[$node] = $tail[0];
    }
    elsif ( scalar @tail == 2 ) {
        $self->branch_length->[$node] = $tail[-1];
    }
    elsif ( scalar @tail == 3 ) {
        $self->node_name->[$node] = $tail[0];
        $self->branch_length->[$node] = $tail[-1];
    }
}

sub _next_token {
    my ( $self, $string ) = @_;
    my $QUOTED = 0;
    my $token = '';
    my $TOKEN_DELIMITER = qr/[():,;]/;
  TOKEN: for my $i ( 0 .. length( $string ) ) {
      next TOKEN if substr($string,$i,1) =~ /\s/;
      $token .= substr($string,$i,1);
      if ( ! $QUOTED && $token =~ $TOKEN_DELIMITER ) {
	  my $length = length( $token );
	  if ( $length == 1 ) {
	      return $token, substr($string,($i+1));
	  }
	  else {
	      return substr($token,0,$length-1),substr($token,$length-1,1).substr($string,($i+1));
	  }
      }
      if ( ! $QUOTED && substr($string,$i,1) eq "'" ) {
	  $QUOTED++;
      }
      elsif ( $QUOTED && substr($string,$i,1) eq "'" && substr($string,$i,2) ne "''" ) {
	  $QUOTED--;
      }        
  }
}

1;
