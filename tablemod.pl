#!/usr/bin/perl -w

open(WRITE,">$tablenew.tex");
open(LOG,"$table.txt");
while ($line=<LOG>) {

$line=~s/{smcl}|{com}|{sf}|{ul off}|{res}//g;
$line=~s/\& {\\ //g;
chop($line);
$pos=index($line,">");

if ($pos==1) {
$line=~s/> //;
print WRITE ("$line");
} else {
print WRITE ("\n$line");
}

}