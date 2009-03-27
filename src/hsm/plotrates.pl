#!/usr/local/bin/perl -w

my $table = "tmp_table";
my $annot = "tmp_annot";
my $image = "tmp_image";
END {
    unlink $table if -e $table;
    unlink $annot if -e $annot;
    unlink $image if -e $image;
}

my $usage = "Usage: $0 [-gif] [-blue] <HSM file> <image file>\n";
my @argv;
my $gif = "";
my $blue = "";
foreach my $arg (@ARGV) {
    if ($arg =~ /^-/) {
	if ($arg eq "-gif") { $gif = "-gif" }
	elsif ($arg eq "-blue") { $blue = "-allblue" }
	else { die $usage }
    } else { push @argv, $arg }
}
die $usage if @argv != 1 && @argv != 2;
my ($hsmfile, $outfile) = @argv;
$outfile = $image if !defined $outfile;

local *HSM;
open HSM, "<$hsmfile" or die $!;
my @al = split (//, uc "arndcqeghilkmfpstwyv");
$_=<HSM>;
my($C,$A)=split;
$_=<HSM>;
my@p=split;
my@intra;
foreach my $c(0..$C-1){my@c;foreach my $ai(0..$A-1){$_=<HSM>;push@c,[split]}push@intra,\@c}
my@inter;
foreach my $a(0..$A-1){my@c;foreach my $ci(0..$C-1){$_=<HSM>;push@c,[split]}push@inter,\@c}
close HSM;

die "Bad alphabet size: $A" if $A != @al;

local *TABLE;
open TABLE, ">$table" or die $!;
print TABLE "_";foreach my $c(1..$C){print TABLE map("\t$_$c",@al)}print TABLE "\n";
foreach my $ci(0..$C-1){
    foreach my $ai(0..$A-1){
	print TABLE $al[$ai], $ci+1;
	foreach my $cj(0..$C-1){
	    if ($ci!=$cj) { print TABLE map("\t*",0..$ai-1),"\t",$inter[$ai]->[$ci]->[$cj],map("\t*",$ai+1..$A-1) }
	    else { print TABLE map ($_==$ai?"\t*":"\t".$intra[$ci]->[$ai]->[$_], 0..$A-1) }
	}
	print TABLE "\n";
    }
}
close TABLE or die $!;

local *ANNOT;
open ANNOT, ">$annot" or die $!;
print ANNOT "_";foreach my $c(1..$C){print ANNOT map("\t$_$c",@al)}print ANNOT "\n";
close ANNOT or die $!;

system "fancyeisen.pl $gif $blue -blockwidth 20 -blockheight 20 -rowtext -coltext -ignoresign -bgwhite $table $annot >$outfile";
system "display $outfile";
