#!/usr/bin/perl

use Getopt::Long;
use File::Basename;

use SequenceIterator qw(iterseq printseq revcomp);
use Stockholm;

my ($progname) = fileparse($0);

my $usage = "";
$usage .= "$progname -- convert DNA to tokenized-codon sequence (or protein sequence)\n";
$usage .= "\n";
$usage .= "Usage: $progname [-table] [-f <frame>] [-revcomp] [-aa] [-rna] [-decode] [-truncate] [-align <Stockholm alignment file>] [FASTA filename(s)]\n";
$usage .= "\n";
$usage .= "The 'frame' (i.e. reading frame) can be 0, 1, or 2.\n";
$usage .= "If '-truncate' is specified, terminal stop codons will be discarded,\n";
$usage .= "and premature stop codons will result in truncation and a warning being issued.\n";
$usage .= "\n";
$usage .= "The interpretation of the alignment file is slightly different when encoding and decoding:\n";
$usage .= "\n";
$usage .= " When encoding, the FASTA file is assumed to be DNA, and the (optional) Stockholm file is protein.\n";
$usage .= " If an alignment file is specified when tokenizing, the output will be that alignment\n";
$usage .= " with all sequences replaced by the correspondingly-named tokenized sequences.\n";
$usage .= " (The translations of the DNA sequences in the FASTA file must match the protein sequences in the Stockholm file.\n";
$usage .= "  If no FASTA sequence file is specified, the Stockholm file is assumed to be a protein-coding DNA alignment.)\n";
$usage .= "\n";
$usage .= " When decoding, the FASTA file is assumed to be tokenized, and the (optional) Stockholm file is also tokenized.\n";
$usage .= " If an alignment file is specified when decoding, the output will be that alignment\n";
$usage .= " with all sequences de-tokenized, annotated with their protein-coding translation.\n";
$usage .= " (For decoding, no FASTA file need be specified if a Stockholm file is specified, since both are tokenized.)\n";
$usage .= "\n";
$usage .= "EXAMPLES\n";
$usage .= "\n";
$usage .= "$progname -table\n";
$usage .= "...to print the table of codons, tokens and amino acids\n";
$usage .= "\n";
$usage .= "$progname DNASEQS.fasta  > TOKSEQS.fasta\n";
$usage .= "...to get a FASTA file of token sequences from the DNA sequences in DNASEQS.fasta\n";
$usage .= "\n";
$usage .= "$progname DNASEQS.fasta -aa  > PROTSEQS.fasta\n";
$usage .= "...to get a FASTA file of protein sequences from the DNA sequences in DNASEQS.fasta\n";
$usage .= "\n";
$usage .= "$progname DNASEQS.fasta -align PROTALIGN.stockholm -truncate  > TOKALIGN.stockholm\n";
$usage .= "...to get a Stockholm alignment of token sequences from the DNA sequences in DNASEQS.fasta, aligned as per the protein alignment in PROTALIGN.stockholm\n";
$usage .= "\n";
$usage .= "$progname -align DNAALIGN.stockholm  > TOKALIGN.stockholm\n";
$usage .= "...to get a Stockholm alignment of token sequences from the DNA sequences in DNAALIGN.stockholm, maintaining the alignment\n";
$usage .= "\n";
$usage .= "$progname -decode TOKSEQS.fasta  > DNASEQS.fasta\n";
$usage .= "...to get a FASTA file of DNA sequences from the token sequences in TOKSEQS.fasta\n";
$usage .= "\n";
$usage .= "$progname -decode TOKSEQS.fasta -aa  > PROTSEQS.fasta\n";
$usage .= "...to get a FASTA file of protein sequences from the token sequences in TOKSEQS.fasta\n";
$usage .= "\n";
$usage .= "$progname -decode -align TOKALIGN.stockholm  > DNAALIGN.stockholm\n";
$usage .= "...to get a Stockholm alignment of DNA sequences (annotated with protein translations) from the token alignment in TOKALIGN.stockholm\n";
$usage .= "\n";
$usage .= "$progname -decode -align TOKALIGN.stockholm -aa  > PROTALIGN.stockholm\n";
$usage .= "...to get a Stockholm alignment of protein sequences from the token alignment in TOKALIGN.stockholm\n";
$usage .= "\n";

my $gr_aa = "AA";  # 3-letter amino acid annotation
my $colWidth = 50;    # column width for output

my (%aa, %tok);
%aa = ( 'ttt'=>'F',  'tct'=>'S',  'tat'=>'Y',  'tgt'=>'C',
        'ttc'=>'F',  'tcc'=>'S',  'tac'=>'Y',  'tgc'=>'C',
        'tta'=>'L',  'tca'=>'S',  'taa'=>'!',  'tga'=>'!',
        'ttg'=>'L',  'tcg'=>'S',  'tag'=>'!',  'tgg'=>'W',
       
        'ctt'=>'L',  'cct'=>'P',  'cat'=>'H',  'cgt'=>'R',
        'ctc'=>'L',  'ccc'=>'P',  'cac'=>'H',  'cgc'=>'R',
        'cta'=>'L',  'cca'=>'P',  'caa'=>'Q',  'cga'=>'R',
        'ctg'=>'L',  'ccg'=>'P',  'cag'=>'Q',  'cgg'=>'R',
       
        'att'=>'I',  'act'=>'T',  'aat'=>'N',  'agt'=>'S',
        'atc'=>'I',  'acc'=>'T',  'aac'=>'N',  'agc'=>'S',
        'ata'=>'I',  'aca'=>'T',  'aaa'=>'K',  'aga'=>'R',
        'atg'=>'M',  'acg'=>'T',  'aag'=>'K',  'agg'=>'R',
       
        'gtt'=>'V',  'gct'=>'A',  'gat'=>'D',  'ggt'=>'G',
        'gtc'=>'V',  'gcc'=>'A',  'gac'=>'D',  'ggc'=>'G',
        'gta'=>'V',  'gca'=>'A',  'gaa'=>'E',  'gga'=>'G',
        'gtg'=>'V',  'gcg'=>'A',  'gag'=>'E',  'ggg'=>'G',

	'n'=>'z', 'nn'=>'Z', 'nnn'=>'X',
	map (($_ x 3 => $_), qw(* - ? . x)) );

my %aa_to_aa3 = ('A' => 'Ala',
		 'C' => 'Cys',
		 'D' => 'Asp',
		 'E' => 'Glu',
		 'F' => 'Phe',
		 'G' => 'Gly',
		 'H' => 'His',
		 'I' => 'Ile',
		 'K' => 'Lys',
		 'L' => 'Leu',
		 'M' => 'Met',
		 'N' => 'Asn',
		 'P' => 'Pro',
		 'Q' => 'Gln',
		 'R' => 'Arg',
		 'S' => 'Ser',
		 'T' => 'Thr',
		 'V' => 'Val',
		 'W' => 'Trp',
		 'Y' => 'Tyr',
		 map (($_ => '!!!'), qw(! 0 1 2)),
		 'z'=>'n', 'Z'=>'nn', 'X'=>'nnn',
		 map (($_ => $_ x 3), qw(* - ? . x)) );

my %aa3 = map (($_ => $aa_to_aa3{$aa{$_}}), keys %aa);

# ASCII characters 33 through 126:
# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
# Characters avoided in token set:
# ( ) ; " ' (used by S-expression format)
# - * ? . (used for alignment/reconstruction/ambiguity)
# > (used by FASTA format)
%tok = ( 'ttt'=>'F',  'tct'=>'S',  'tat'=>'Y',  'tgt'=>'C',
	 'ttc'=>'f',  'tcc'=>'s',  'tac'=>'y',  'tgc'=>'c',
	 'tta'=>'L',  'tca'=>'5',  'taa'=>'0',  'tga'=>'2',
	 'ttg'=>'l',  'tcg'=>'$',  'tag'=>'1',  'tgg'=>'W',
       
	 'ctt'=>'<',  'cct'=>'P',  'cat'=>'H',  'cgt'=>'R',
	 'ctc'=>'[',  'ccc'=>'p',  'cac'=>'h',  'cgc'=>'r',
	 'cta'=>'{',  'cca'=>',',  'caa'=>'Q',  'cga'=>'+',
	 'ctg'=>'/',  'ccg'=>'8',  'cag'=>'q',  'cgg'=>'}',
       
	 'att'=>'I',  'act'=>'T',  'aat'=>'N',  'agt'=>'%',
	 'atc'=>'i',  'acc'=>'t',  'aac'=>'n',  'agc'=>'#',
	 'ata'=>'|',  'aca'=>'~',  'aaa'=>'K',  'aga'=>'@',
	 'atg'=>'M',  'acg'=>'`',  'aag'=>'k',  'agg'=>']',
       
	 'gtt'=>'V',  'gct'=>'A',  'gat'=>'D',  'ggt'=>'G',
	 'gtc'=>'v',  'gcc'=>'a',  'gac'=>'d',  'ggc'=>'g',
	 'gta'=>'^',  'gca'=>'4',  'gaa'=>'E',  'gga'=>'9',
	 'gtg'=>'7',  'gcg'=>'&',  'gag'=>'e',  'ggg'=>'6',

	 'n'=>'z', 'nn'=>'Z', 'nnn'=>'X',
	 map (($_ x 3 => $_), qw(* - ? . x)) );

my %is_stop = map (($_ => 1), qw(tag taa tga));

my $frame = 0;
my $revcomp = 0;
my $use_aa = 0;
my $is_rna = 0;
my $untokenize = 0;
my $truncate = 0;
my $print_table = 0;
my $align_file;

GetOptions ("frame=i" => \$frame,
	    "revcomp" => \$revcomp,
	    "aa"  => \$use_aa,
	    "rna"  => \$is_rna,
	    "truncate" => \$truncate,
	    "table" => \$table,
	    "align=s" => \$align_file,
	    "decode" => \$untokenize) or die $usage;

if ($table) {
    print map (join(" ",$_,$aa3{$_},$aa{$_},$tok{$_})."\n", sort keys %tok);
    exit;
}

my ($stock, %untranslated);
if (defined($align_file) && @ARGV > 0 && !$truncate && !$untokenize) { warn "Warning: -align option is best used with -truncate\n" }
if (defined($align_file) && @ARGV > 0 && ($frame != 0 || $revcomp)) { warn "Warning: -align option is not currently compatible with -f or -revcomp options\n" }
if (defined $align_file) {
    $stock = Stockholm->from_file ($align_file);
    %untranslated = map ($stock->seqdata->{$_} =~ /[^\-\*\.]/ ? ($_ => 1) : (), @{$stock->seqname});
}

my $trans_ref = $use_aa ? \%aa : \%tok;
my $annot_ref = $use_aa ? \%aa : \%aa3;

my %untok = map (($tok{$_} => ($use_aa ? $aa{$_} : $_)), keys %tok);
my %annot = map (($tok{$_} => ($use_aa ? $aa{$_} : $aa3{$_})), keys %tok);

my $cod_size = $use_aa ? 1 : 3;

# Uncomment to check for duplicate tokens
#for my $c (keys %tok) { die $c unless $untok{$tok{$c}} eq $c }

if (defined($align_file) && @ARGV == 0) {
    for my $name (@{$stock->seqname}) {
	if (defined ($stock->seqdata->{$name})) {
	    my $stockrow = $stock->seqdata->{$name};
	    my $stockseq = $stockrow;
	    if ($untokenize) {
		$stockseq =~ s/[\-\.]//g;
		visit_seq ($name, $stockseq);
	    } else {
		visit_seq ($name, $stockseq, 1);
	    }
	}
    }
} else {
    @ARGV = ("-") unless @ARGV;
    for my $filename (@ARGV) {
	iterseq ($filename, \&visit_seq);
    }
}

if (defined $align_file) {
    for my $name (keys %untranslated) { $stock->seqdata->{$name} =~ s/[^\-\.]/*/g }
    print $stock->to_string;
    if (%untranslated) {
	warn "The following alignment rows could not be tokenized, so were replaced with wildcards:\n", join (" ", keys %untranslated), "\n";
    }
}

sub visit_seq {
    my ($name, $seq, $seq_is_alignment_row) = @_;
    $seq_is_alignment_row = 0 unless defined $seq_is_alignment_row;
    my $gr;
    my $newseq = $untokenize ? untokenize($seq,$name,\$gr) : tokenize($seq,$name,$trans_ref,1);
    if (defined $align_file) {
	if (defined ($stock->seqdata->{$name})) {
	    my $stockrow = $stock->seqdata->{$name};
	    my $stockseq = $stockrow;
	    $stockseq =~ s/[\-\.]//g;
	    if ($untokenize) {
		if ($stockseq ne $seq) {
		    warn
			"Sequence '$name' does not match corresponding alignment row.\n",
			" Alignment row: $stockrow\n",
			" Alignment seq: $stockseq\n",
			" Tokenized seq: $seq\n\n";
		} else {
		    my $newrow = "";
		    my $newgr = "";
		    my $newseq_pos = 0;
		    for (my $row_pos = 0; $row_pos < length($stockrow); ++$row_pos) {
			my $row_char = substr ($stockrow, $row_pos, 1);
			if ($row_char eq '-' || $row_char eq '.') {
			    $newrow .= $row_char x $cod_size;
			    $newgr .= $row_char x $cod_size if defined $gr;
			} else {
			    $newrow .= substr ($newseq, $newseq_pos, $cod_size);
			    $newgr .= substr ($gr, $newseq_pos, $cod_size) if defined $gr;
			    $newseq_pos += $cod_size;
			}
		    }
		    $stock->seqdata->{$name} = $newrow;
		    $stock->gr->{$gr_aa}->{$name} = $newgr if defined $gr;
		    delete $untranslated{$name};
		}
	    } elsif ($seq_is_alignment_row) {
		my $aaseq = tokenize($seq,$name,\%aa,0);
		$stock->seqdata->{$name} = $newseq;
		$stock->gr->{$gr_aa}->{$name} = $aaseq;
		delete $untranslated{$name};
	    } else {
		my $aaseq = tokenize($seq,$name,\%aa,0);
		if (uc($stockseq) ne $aaseq) {
		    warn
			"Translation of sequence '$name' does not match corresponding alignment row.\n",
			" Alignment row: $stockrow\n",
			" Alignment seq: $stockseq\n",
			"Translated seq: $aaseq\n",
			" Tokenized seq: $newseq\n\n";
		} else {
		    my $newrow = "";
		    my $newseq_pos = 0;
		    for (my $row_pos = 0; $row_pos < length($stockrow); ++$row_pos) {
			my $row_char = substr ($stockrow, $row_pos, 1);
			if ($row_char eq '-' || $row_char eq '.') {
			    $newrow .= $row_char;
			} else {
			    $newrow .= substr ($newseq, $newseq_pos++, 1);
			}
		    }
		    $stock->seqdata->{$name} = $newrow;
		    $stock->gr->{$gr_aa}->{$name} = $stockrow;
		    delete $untranslated{$name};
		}
	    }
	} else {
	    warn "Sequence '$name' not found in alignment; ignoring\n";
	}
    } else {
	printseq ($name, $newseq);
    }
}

sub tokenize {
    my ($seq, $name, $trans_ref, $warning) = @_;
    $seq = lc $seq;
    $seq =~ s/u/t/g;  # do this even if -rna was not specified; no need to punish user
    if ($revcomp) { $seq = revcomp ($seq); $name .= " (reverse strand)" }
    my $trans = "";
  CODON: for (my $pos = $frame; $pos < length($seq); $pos += 3) {
	my $remaining_chars = length($seq) - ($pos + 3);
	$codon = substr ($seq, $pos, 3);
	if (exists $$trans_ref{$codon}) {
	    if ($truncate && $is_stop{$codon}) {
		if ($warning && $remaining_chars > 0) {
		    warn "Premature stop codon ($codon) found with $remaining_chars characters remaining while tokenizing sequence $name\n";
		}
		last CODON;
	    }
	    $trans .= $$trans_ref{$codon};
	} elsif (length($codon) == 1) { $trans .= $$trans_ref{'n'}; warn "Extra character ($codon) at end of sequence $name\n" if $warning }
	elsif (length($codon) == 2) { $trans .= $$trans_ref{'nn'}; warn "Extra characters ($codon) at end of sequence $name\n" if $warning }
	else { $trans .= $$trans_ref{'nnn'}; warn "Unrecognized codon ($codon) at position $pos of sequence $name\n" if $warning }
    }
    return $trans;
}

sub untokenize {
    my ($seq, $name, $aa3_ref) = @_;
    my $untrans = "n" x $frame;
    my $annot = $untrans;
    for (my $pos = 0; $pos < length($seq); ++$pos) {
	$token = substr ($seq, $pos, 1);
	if (exists $untok{$token}) {
	    my $untok = $untok{$token};
	    my $annot_sym = $annot{$token};
	    if (defined($align_file) && length($untok) != $cod_size) {
		die "Can't cope with ($token=>$untok) tokens when decoding alignments\n";
	    }
	    if (defined($align_file) && length($annot_sym) != $cod_size) {
		die "Can't cope with ($token=>$annot_sym) tokens when decoding alignments\n";
	    }
	    $untrans .= $untok;
	    $annot .= $annot_sym;
	} else { $untrans .= 'nnn'; $annot .= '...'; warn "Unrecognized token ($token) at position $pos of sequence $name\n" }
    }
    if ($revcomp) { $untrans = revcomp ($untrans); $annot = reverse $annot }
    if ($is_rna) { $untrans =~ s/t/u/g }
    if (ref $aa3_ref) { $$aa3_ref = $annot }
    return $untrans;
}
