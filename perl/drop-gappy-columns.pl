#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Written and maintained by Andrew Uzilov, May-August 2006.
#-------------------------------------------------------------------------------

use strict;
use Stockholm;

# alphabet of characters denoting gaps
my @gaps = ('.', '_', '-', '/', '?');

# RNA only:
# characters that annotate single-stranded (unpaired) nucleotides/columns
my @unpaired = ('.', ',', '_', '-', ':', '~');

# RNA only:
# pairs of characters that annotate base-paired nucleotides/columns
# (note that we can also use upper/lowercase letters to denote base pairs,
#  with the uppercase being 5', the lowercase being 3')
my @paired = (['(', ')'],
	      ['{', '}'],  # N.B.: the first in a pair is ALWAYS 5'!
	      ['[', ']'],
	      ['<', '>']);

for (my $upper = 65; $upper <= 90; $upper++)  # fill out the letter pairs
  {
    push @paired, [chr($upper), chr($upper + 32)];
  }

# these store arguments to the program:
#   $rna_ss = true if we are expecting an RNA secondary structure on the
#             annotation line and we are expected to use it to drop non-gappy
#             columns pairing to gappy columns
#   $threshold = stores the gap threshold
#   $in_file = path/name of input file
#   $d = true if running in debug mode
#   $width = column width of output
#   $drop_annot_gaps = true if user wants to only drop columns with gaps in the
#                      $drop_annot '#=GC' or '#=GR' annotation
#   $out_file = file to output to; undef if you want output written to STDOUT
#   $unpaired_char = character that replaces base pair symbols in the
#                    '#=GC SS_cons' line when bases paired to that column get
#                    dropped (only applicable when $rna_ss is true); undef if we
#                    just want to drop the whole column
#   $overwrite_char = if defined, do not drop columns, but instead overwrite all
#                     their GC and GR annotations with this character
#
my ($rna_ss, $threshold, $in_file, $d, $width, $drop_annot_gaps, $drop_annot,
    $out_file, $unpaired_char, $overwrite_char) =
  (0, 1, undef, 0, 80, 0, undef, undef, undef, undef);  # default values

parse_params(@ARGV);

if ($d)
  {
    # NB: yes, this will generate warnings in many cases... don't worry about it
    warn (
	  "is RNA = $rna_ss; threshold = $threshold; input file = $in_file; ",
	  "width = $width; drop gaps from #=GC = $drop_annot_gaps; ",
	  "unpaired char = $unpaired_char; overwrite char = $overwrite_char; ",
	  "output file = $out_file\n"
	 );
  }

my $stk = Stockholm->from_file($in_file);
die "ERROR: file is not flush; cannot proceed!\n" unless $stk->is_flush();

if ($out_file) {
  open (OUTFILE, ">$out_file")
    or die "ERROR: cannot open $out_file for writing ($!)\n";
}

# stores indices (1-based) of columns to be dropped
# (N.B.: do not be afraid to put duplicates in here, as we will remove them
#        prior to any duplicate-sensitive operations)
my @dropped_columns;

# build the set of columns to drop, based on what the user specified on the
# command line, which is either:
#   - columns containing gap characters in the '#=GC' and '#=GR' annotation
#     lines with the specified $drop_annot feature
# or
#   - gappy columns
# but not both.
#
if ($drop_annot_gaps)
  {
    # drop columns only from annotation lines

    foreach my $feature (keys %{$stk->{gc}})  # for each '#=GC' annotation
      {
	if ($feature eq $drop_annot)
	  {
	    my @annot = split (//, $stk->{gc}->{$feature});
	    for (my $col = 1; $col <= @annot; $col++)
	      {
		my $char = $annot[$col-1];  # convert to 0-based indexing

		# does the annotation character match any gap character?
		my $gapmap = join ('', map {$_ eq $char ? 1 : 0} @gaps);
		push (@dropped_columns, $col) unless $gapmap == 0;
	      }
	  }
      }

    foreach my $feature (keys %{$stk->{gr}})  # for each '#=GR' annotation
      {
	if ($feature eq $drop_annot)
	  {
	    foreach my $seq (keys %{$stk->{gr}->{$feature}})  # for each sequence
	      {
		my @annot = split (//, $stk->{gr}->{$feature}->{$seq});
		for (my $col = 1; $col <= @annot; $col++)
		  {
		    my $char = $annot[$col-1];  # convert to 0-based indexing

		    # does the annotation character match any gap character?
		    my $gapmap = join ('', map {$_ eq $char ? 1 : 0} @gaps);
		    push (@dropped_columns, $col) unless $gapmap == 0;
		  }
	      }
	  }
      }
  }
else
  {
    # drop columns only based on gappiness
    for (my $col = 1; $col <= $stk->columns(); $col++)  # NB: 1-based indexing!
      {
	my $gap_count = 0;
	# for each nuc in col...
        foreach my $nuc (split (//, $stk->get_column($col-1)))  # to 0-based
	  {
	    # does the nuc character match any gap character?
	    my $gapmap = join ('', map {$_ eq $nuc ? 1 : 0} @gaps);
	    $gap_count++ unless $gapmap == 0;
	  }
	warn "column $col contains $gap_count gaps\n" if $d;

	# the column gets dropped unless it is below the gap threshold
	push (@dropped_columns, $col)
	  unless ( ($gap_count / $stk->sequences()) < $threshold );
      }
  }

# remove duplicates from the list of dropped columns, so that each column is in
# the list exactly once; everything that follows may depend on there being NO
# duplicates in the dropped columns list
#
@dropped_columns = remove_duplicates(@dropped_columns);

# RNA only:
# For any base pair with a nucleotide in a dropped column (as determined by the
# annotation on the '#=GC SS_cons' line), either drop the column it is paired
# to, or (if $unpaired_char was specified), replace the annotation in the column
# it is paired to with $unpaired_char on the '#=GC SS_cons' line, so that the
# base pair nesting remains correct.
#
# Note that when we say "dropped" below, we usually mean "dropped or replaced by
# $unpaired_char".
#
if ($rna_ss)
  {
    # parse secondary structure string to a set of base pairs
    my @bp_set = parse_ss($stk->{'gc'}->{'SS_cons'});

    # keep track of dropped or replaced columns that are in base pairs
    my @bp_cols_dropped;

    # build a list of base pair-containing columns that are dropped for ANY
    # reason, even if both bases are in gappy columns already (this enables us
    # to count them)

    while (@bp_set)  # while there are still base pairs left in the set...
      {
	my $bp = shift @bp_set;  # ...take a base pair out...
	my ($nuc5, $nuc3) = ($bp->[0], $bp->[1]);

	# ...and see if any of its bases match a dropped column (this is a bit
	# inefficient algorithmically, but... N is small here)
	foreach my $col (@dropped_columns) {
	  if ( ($col == $nuc5) || ($col == $nuc3) ) {
	    push (@bp_cols_dropped, $nuc5, $nuc3);
	  }
	}
      }

    # remove duplicates resulting from our stupid (but easy to implement)
    # accounting process
    @bp_cols_dropped = remove_duplicates(@bp_cols_dropped);

    # now that duplicates are gone, number of base pairs dropped or replaced
    # (for ANY reason) should just be half of their columns
    warn "we lost ", (@bp_cols_dropped / 2), " base pair(s)\n" if $d;

    if ($unpaired_char)  # replace annotation with unpaired symbol
      {
	my @annot = split ('', $stk->{'gc'}->{'SS_cons'});
	foreach my $col (@bp_cols_dropped)
	  {
	    $annot[$col-1] = $unpaired_char;  # correct the indexing
	  }
	$stk->{'gc'}->{'SS_cons'} = join ('', @annot);
      }
    else  # drop the entire column
      {
	# merge into the dropped columns list and resolve duplication
	@dropped_columns = remove_duplicates(@dropped_columns, @bp_cols_dropped);
      }
  }

# drop the columns and output
if (scalar(@dropped_columns))
  {
    # drop the columns, or just replace the annotation?
    if ($overwrite_char)
      {
	warn "Overwriting annotation with $overwrite_char for ",
	  scalar(@dropped_columns), " columns: @dropped_columns\n"
	    if $d;

	foreach my $key (keys %{$stk->{'gc'}})  # overwrite GC annotations
	  {
	    my @string_as_array = split ('', $stk->{'gc'}->{$key});
	    foreach my $col (@dropped_columns)
	      {
		$string_as_array[$col-1] = $overwrite_char;
	      }
	    $stk->{'gc'}->{$key} = join ('', @string_as_array);
	    warn "annotation GC $key is now: ", $stk->{'gc'}->{$key}, "\n"
	      if $d;
	  }
	foreach my $key_feat (keys %{$stk->{'gr'}})  # overwrite GR annotations
	  {
	    #warn "on feature $key_feat\n";  # D!!!
	    foreach my $key_seq (keys %{$stk->{'gr'}->{$key_feat}})
	      {
		#warn "on sequence $key_seq\n";  # D!!!
		my @string_as_array =
		  split ('', $stk->{'gr'}->{$key_feat}->{$key_seq});
		foreach my $col (@dropped_columns)
		  {
		    $string_as_array[$col-1] = $overwrite_char;
		  }
		$stk->{'gr'}->{$key_feat}->{$key_seq} =
		  join ('', @string_as_array);
		warn "annotation GR $key_feat of seq $key_seq is now: ",
		  $stk->{'gr'}->{$key_feat}->{$key_seq}, "\n"
		    if $d;
	      }
	  }

	# overwriting complete... output
	if ($out_file)
	  {
	    print OUTFILE $stk->to_string($width);
	  }
	else
	  {
	    print $stk->to_string($width);
	  }
      }
    else
      {
	warn "Dropping ",
	  scalar(@dropped_columns), " columns: @dropped_columns\n"
	    if $d;

	if ($stk->drop_columns(map {$_-1} @dropped_columns))  # to 0-based
	  {
	    if ($out_file)
	      {
		print OUTFILE $stk->to_string($width);
	      }
	    else
	      {
		print $stk->to_string($width);
	      }
	  }
	else
	  {
	    die "ERROR: error dropping columns! (no output)\n";
	  }
      }
  }
else
  {
    warn "No columns dropped\n" if $d;

    if ($out_file)
      {
	print OUTFILE $stk->to_string($width);
      }
    else
      {
	print $stk->to_string($width);
      }
  }


exit 0;


#===============================================================================
#     HELPER FUNCTIONS
#===============================================================================


#
# Removes duplicates from a list (array) of numbers passed in, returns the
# sorted, filtered array.
#
sub remove_duplicates
  {
    return unless @_;  # bail if nothing is passed in, otherwise bad tings...

    # numerically sort the input
    my @sorted_input = sort { $a <=> $b } @_;

    # initialize returned value to the first item in the sorted input
    my @return_arr = shift @sorted_input;

    # to through the rest of the sorted input...
    foreach my $num (@sorted_input)
      {
	# ...and only add to returned value if not a duplicate
	unless ($num == $return_arr[$#return_arr])
	  {
	    push (@return_arr, $num);
	  }
      }

    return @return_arr;
  }


#
# Helper for converting the Stockholm '#=GC SS_cons' secondary structure string
# to a set of base pairs.  Yes, this does pseudoknots (assuming input follows
# the Rfam convention for annotating pseudoknots).
#
# The set is an array of references to 2-item arrays representing base pairs,
# storing 1-based indices of nucleotides in the pair.  The 0th item is always
# the 5' nuc, the 1st item is always the 3' nuc.
#
# If you have a nucleotide index in mind and want to figure out what it is
# paired to, if anything, you have to do a search through the array (yeah, this
# is stupid).
#
sub parse_ss {
  my @struct_string = split(//, shift);  # array representation of the
                                         # structure string
  my @bp_set;  # for returning

  # parse string, rinsing the pushdown automata massive inside the ride!

  # we're going to make multiple passes because of the possibility of
  # pseudoknotted structures; we know that Rfam labels loops in pseudoknots
  # with different symbols, so if we fix a single pair of symbols to
  # denote opening/closing (5'/3'), and do a pass for each, we should be
  # able to parse each pass using a PDA (yeah the runtime is terrible, so
  # shoot me, the input size is miniscule)
  #
  foreach my $pair (@paired)  # loop through all ways to annotate a base pair
    {
      my $nuc5prime = $pair->[0];
      my $nuc3prime = $pair->[1];

      #warn "5' = $nuc5prime, 3' = $nuc3prime\n" if $d;  # for XtrEEEm debugging

      # the PDA stack; stores, for every LEFT (5') nuc in a base pair, its
      # 1-based index
      my @stack;

      for (my $i = 1; $i <= @struct_string; $i++)
	{
	  # remember, $i is our 1-based nuc index!

	  if ($struct_string[$i-1] eq $nuc5prime)
	    {
	      push(@stack, $i);  # push "opening" symbol onto PDA stack
	    }
	  elsif ($struct_string[$i-1] eq $nuc3prime)
	    {
	      # pop "closing" symbol from PDA stack, check for correctness
	      my $matching_nuc = pop @stack;
	      die
		"ERROR: malformed SS_cons string at nuc $i:\n",
		  @struct_string, "\n"
		    unless $matching_nuc;

	      push(@bp_set, [$matching_nuc, $i]);
	    }
	}

      # if the structure string was well-formed, we should not have anything
      # left on the PDA stack
      if (@stack > 0)
	{
	  die
	    "ERROR: stack not empty after parsing SS_cons string!\n",
	      "contains: ", @stack, "\n";
	}

    }  # closes 'foreach my $pair ...'

  return @bp_set;
}


sub parse_params {
  print_usage("You must tell the program what to do.") if @_ < 1;

  # pre-filter for help options
  foreach my $opt (@_)
    {
      if ( ($opt eq '-h') or ($opt eq '--help') or ($opt eq '-?') )
	{
	  print_usage();
	  return;
	}
    }

  # check that the last option is an input file, not a parameter
  if ( ($_[$#_] =~ /^-./) || ($_[$#_ - 1] =~ /^-[anotw]/) )
    {
      print_usage("You should provide the input file as the last argument, not " .
		  "a parameter.");
    }
  else
    {
      $in_file = $_[$#_];  # get file name
    }

  print_usage("You did not specify an input file.") unless $in_file;

  # parse the rest of the options (except the last one)
  for (my $i = 0; $i < $#_; $i++)
    {
      my $opt = $_[$i];

      if ($opt eq '-a')
	{
	  @gaps = split (//, $_[++$i]);
	}
      elsif ($opt eq '-d')
	{
	  $d = 1;
	}
      elsif ($opt eq '-g')
	{
	  $drop_annot_gaps = 1;
	  $drop_annot = $_[++$i];
	}
      elsif ($opt eq '-n')
	{
	  $overwrite_char = $_[++$i];
	  print_usage ("You should only specify a single character after " .
		       "the -n argument (you specified: $overwrite_char)")
	    unless length ($overwrite_char) == 1;
	}
      elsif ($opt eq '-o')
	{
	  $out_file = $_[++$i];
	}
      elsif ($opt eq '-p')
	{
	  $unpaired_char = $_[++$i];
	  print_usage ("You should only specify a single character after " .
		       "the -p argument (you specified: $unpaired_char)")
	    unless length ($unpaired_char) == 1;
	}
      elsif ($opt eq '-r')
	{
	  $rna_ss = 1;
	}
      elsif ($opt eq '-t')
	{
	  $threshold = $_[++$i];
	  print_usage("You did not specify a threshold after the -t argument.")
	    if $threshold eq '';

	  # since threshold can be a ratio, parse it if such
	  if ($threshold =~ /^(\d+)\/(\d+)$/)
	    {
	      $threshold = $1 / $2;
	    }

	  # do correctness check
	  print_usage("Your threshold ($threshold) is not a valid threshold.")
	    unless $threshold =~ /^\d*.?\d+$/;

	  print_usage("Your threshold ($threshold) is not in the range [0, 1].")
	    unless ($threshold >=0) and ($threshold <= 1);
	}
      elsif ($opt eq '-w')
	{
	  $width = $_[++$i];
	  print_usage("You did not specify a width after the -w argument.")
	    if $width eq '';

	  print_usage ("Your width ($width) is not a valid integer width.")
	    unless $width =~ /^\d+$/;

	  print_usage("You cannot have the non-positive width $width.")
	    unless $width > 0;
	}
      else
	{
	  print_usage("Unrecognized parameter ($opt).");
	}
    }

  # sanity check
  if (defined $unpaired_char) {
    if ($rna_ss != 1) {
      print_usage ("You can only use the -p option if -r is also set.");
    }
  }
}


sub print_usage {
  my $error = shift;

  print <<END;
--------------------------------------------------------------------------------
 Drops gappy columns (columns whose gap content is at or above a specified
 threshold) from a file in Stockholm format by default, or overwrites their
 per-column annotations if the '-n' flag is specified.

 Also drops the per-column annotations (GC and GR) in the gappy columns.  If the
 '-r' (RNA secondary structure) flag is specified, non-gappy columns that base
 pair to gappy columns are also dropped to preserve the nested base pairing
 structure (base pairing is determined from the '#=GC SS_cons' line).

 Output is written to standard out by default, in Stockholm format.

 The gap alphabet is: @gaps

 Usage:
   > $0 [-a] [-d] [-g] [-h|--help|-?] [-n <char>] [-o <output file>] [-r [-p <char>] ] [-t <threshold>] [-w <width>] <Stockholm file>

 where:
   -a <char list>
       Overwrite the default gap alphabet.  Specify a list of characters, no
       whitespace, of your desired gap alphabet.

   -d
       Run in debug mode (lots of verbose output).  Default is off.

   -g <feature>
       Drop columns that contain a gap character in any '#=GC <feature>' or
       '#=GR <seqname> <feature>' (for every <seqname>) annotation.  Note that
       this overrides and ignores the gap threshold - if a <feature> annotation
       has a gap character, that column will get dropped, but no columns will be
       dropped due to their gap content.

   -h|--help|-?
       Print this usage message.

   -n <char>
       Don't drop any columns at all; instead, replace annotations for the gappy
       columns with <char> in all GC and GR lines.  When '-r' is specified,
       annotations for columns pairing to gappy columns are also thus replaced,
       unless '-p' is specified, which overrides '-n' on the '#=GC SS_cons' line
       (but nowhere else).

   -o <output file>
       Output to specified file (Stockholm format), instead of to standard out.

   -p <char>
       Works only with the '-r' option; instead of dropping non-gappy columns
       because they are base paired to gappy columns, replace the annotation in
       the '#=GC SS_cons' line with <char> (which, presumably, is a character
       representing an unpaired base), but leave the column otherwise intact.

   -r
       Specifies that input file contains RNA secondary structure annotations
       (we assume these annotations are on the lines denoted '#=GC SS_cons').
       When in debug mode (-d), number of base pairs dropped is output to to
       standard error.

       Default is we assume column annotations are independent of each other.

   -t <threshold>
       Specifies the gap threshold, which is the fraction of sequences in a
       column that contain gaps.  Columns below the threshold are kept, columns
       at or above the threshold are dropped.  Default is 1, which means that
       only columns containing 100% gaps are dropped.

       Note that you can express the threshold as a ratio of integers; e.g. for
       an alignment of 16 sequences, this:
         -t 9/16
       will specify that columns with 9 or more gaps will get dropped.

       Note that the '-g' option overrides this one.

   -w <width>
       Max width of Stockholm output, in columns (default is $width).
--------------------------------------------------------------------------------
END

  print "ERROR: $error\n" if $error;
  exit 0;
}
