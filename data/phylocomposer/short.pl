#!/usr/bin/perl

%s=qw(k 1 ~k 1 a 1 ~a 1 g 1 ~g 1 b 2 ~b .415 pi_T 2 pi_A 2 pi_C 2 Q_TC 2.737 Q_TA 4.322);

open FILE, "short.sxpr";

while (<FILE>) {

if(/^(;;[^\|]*\|\s*)([^\|]+?)(\s*\|\s*)[^\|]*$/){
    ($a,$p,$b)=($1,$2,$3);
    $_="";
    $o=$p;
    $p=~s/(\S+)\^2/$1 $1/;
    $l=eval(join(" + ",map($s{$_},split(/\s+/,$p))));
    $row="$a$o$b$l  [auto]";
    if (defined $printed_postprob) {
	print $row, "\n";
    } else {
	push@t,$row;
	push@l,$l;
	$t+=2**-$l;
    }
}elsif(/(^;;\s*Forward.*=[ \t]*)([^=]*)$/&&defined$t){
    for$i(0..@t-1){
	print $t[$i];
	if ($i<@l) {
	    print
		"  PostProb = ",
		2**(-$l[$i]) / $t;
	}
	print "\n";
    }
    $printed_postprob = 1;
    $_="$1 ".-log($t)/log(2)."  [auto]\n";
}elsif(defined$t){
    push@t,$_;
}

print;
}

