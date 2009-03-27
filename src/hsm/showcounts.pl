#!/usr/local/bin/perl -w
my$s="arndcqeghilkmfpstwyv";
while(<>){my@f=split;my@c=@f[9..28];my$t=0;foreach my$c(@c){$t+=$c}print$f[5];foreach my$i(0..19){my$c=$c[$i];my$res=substr($s,$i,1);print" $res ",int(1000*$c/$t)/1000}print"\n"}
