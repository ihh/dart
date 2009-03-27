#!/usr/local/bin/perl -w

my $fancyeisen = "fancyeisen.pl";
my $viewer = `uname -a` =~ /linux/i ? "display" : "xv";

my $table = "tmp_table";
my $annot = "tmp_annot";
my $image = "tmp_image";
END {
    unlink $table if -e $table;
    unlink $annot if -e $annot;
    unlink $image if -e $image;
}

my $usage = "Usage: $0 [-gif] [-blue|-grey] [-invert] [-nodisplay] [-eqm] [-rect] [-nointra] <HSM file> [<image file>]\n";
my @argv;
my $gif = "";
my $color = "";
my $invert = 1;
my $display = 1;
my $eqm = 0;
my $rect = "";
my $nointra = 0;
foreach my $arg (@ARGV) {
    if ($arg =~ /^-/) {
	if ($arg eq "-gif") { $gif = "-gif" }
	elsif ($arg eq "-blue") { $color = "-allblue" }
	elsif ($arg eq "-grey") { $color = "-allgrey" }
	elsif ($arg eq "-invert") { $invert = -1 }
	elsif ($arg eq "-nodisplay") { $display = 0 }
	elsif ($arg eq "-eqm") { $eqm = 1 }
	elsif ($arg eq "-rect") { $rect = "-rect" }
	elsif ($arg eq "-nointra") { $nointra = 1 }
	else { die $usage }
    } else { push @argv, $arg }
}
die $usage if @argv != 1 && @argv != 2;
my ($hsmfile, $outfile) = @argv;
$outfile = $image if !defined $outfile;

local *HSM;
open HSM, "<$hsmfile" or die $!;
my @rate;
$_=<HSM>;
my($C,$A)=split;
$_=<HSM>;
my @al = split;
$_=<HSM>;
my@p=split;
my@intra;
foreach my $c(0..$C-1){my@c;foreach my $ai(0..$A-1){$_=<HSM>;push@c,[split]}push@intra,\@c;push@rate,map(@$_,@c)}
my@inter;
foreach my $a(0..$A-1){my@c;foreach my $ci(0..$C-1){$_=<HSM>;push@c,[split]}push@inter,\@c;push@rate,map(@$_,@c)}
close HSM;

die "Bad alphabet size: $A" if $A != @al;

if ($eqm) {  # rescale equilibrium probs to get max contrast
    my $norm = 0;
    foreach my $i (0..$A*$C-1) { if ($p[$i] > $norm) { $norm = $p[$i] } }
    foreach my $i (0..$A*$C-1) { $p[$i] /= $norm }
    @rate = sort { $a <=> $b } @rate;
} else { @rate = () }   # the only use of these rates is to rescale eqm probs, so clear them if !$eqm

local *TABLE;
open TABLE, ">$table" or die $!;
print TABLE "_";
if ($eqm) { print TABLE map ("\t$_", 1..$C) }
foreach my $c (1..$C) {
    print TABLE map ($C>1 ? "\t$_$c" : "\t$_", @al);
}
print TABLE "\n";
foreach my $ai (0..$A-1) {
    print TABLE $al[$ai];
    if ($eqm) { print TABLE map ("\t".$rate[$p[$_*$A+$ai] * @rate-1], 0..$C-1) }
    foreach my $c (0..$C-1) {
	foreach my $aj (0..$A-1) {
	    if ($ai == $aj) { print TABLE "\t*" }  # null out rate to self
	    else { print TABLE "\t", $invert * $intra[$c]->[$ai]->[$aj] }
	}
    }
    print TABLE "\n";
}
if ($C > 1 && !$nointra) {
    foreach my $ci (0..$C-1) {
	print TABLE 'X';
	print TABLE $ci + 1 if $C > 1;
	if ($eqm) { print TABLE map ("\t*", 0..$C-1) }
	foreach my $cj (0..$C-1) {
	    if ($ci == $cj) { print TABLE map ("\t*", @al) }
	    else {
		foreach my $a (0..$A-1) {
		    print TABLE "\t", $invert * $inter[$a]->[$ci]->[$cj];
		}
	    }
	}
	print TABLE "\n";
    }
}
close TABLE or die $!;

local *ANNOT;
open ANNOT, ">$annot" or die $!;
print ANNOT "_";
print ANNOT map ("\t$_", @al);
print ANNOT map ("\tX$_", 1..$C) if $C > 1 && !$nointra;
print ANNOT "\n";
close ANNOT or die $!;

system "$fancyeisen $rect $gif $color -blockwidth 20 -blockheight 20 -rowtext -coltext -ignoresign -bgwhite $table $annot >$outfile";
system "$viewer $outfile" if $display;
