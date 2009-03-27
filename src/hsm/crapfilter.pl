#!/usr/local/bin/perl
while (<>) { s/\0/ /g; s/\/\d+-\d+//; print }
