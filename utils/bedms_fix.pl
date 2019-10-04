#! /usr/bin/perl -w
use strict;

## File: bedms_fix.pl
## Author: Nathan D VanderKraats
##         Washington University in St. Louis, School of Medicine
## Copyright 2013 Washington University in St. Louis
## All Rights Reserved

## Copying Permission Statement
## WIMSi is free software: you can redistribute it and/or modify it under the terms 
## of the GNU General Public License as published by the Free Software Foundation, 
## either version 3 of the License, or (at your option) any later version.  WIMSi 
## is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
## without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
## PURPOSE. See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with WIMSi 
## (see gpl.txt). If not, see http://www.gnu.org/licenses/.

my $usage = "$0 <file1.bedms>\n".
	    "    Attempts to fix basic problems with a 4-column .bedms file.  Fixes:\n".
	    "      1. Makes chromosome blocks contiguous,\n".	
	    "      2. Within each chromosome, orders the entries by position 1 (column 2),\n".
	    "    Outputs new file to STDOUT.\n".
	    "";	    
die $usage unless @ARGV==1;

my %seenchr = (); my %stops = (); my %vals = ();
my ($lastchr,$lastpos,$line) = (-1,-1,0);

open FILE, $ARGV[0] or die "Error: $!\n";
print STDERR "Loading bed file.";
while( <FILE> )
{
    my ($chr, $pos1, $pos2, $val) = split /[\t\n]/;
    die "Error: Chromosome '$chr' in line '$line' doesn't start with 'chr'!\n" unless $chr =~ s/chr(.*)/$1/;
    $stops{$chr}{$pos1} = $pos2;
    $vals{$chr}{$pos1} = $val;
    print STDERR "." if( ++$line % 10000 == 0 );
}
print STDERR "\n";
close FILE;

foreach my $chr (sort mycmp keys %stops)
{
    foreach my $pos1 (sort {$a<=>$b} keys %{$stops{$chr}})
    {
	print STDOUT "chr$chr\t$pos1\t$stops{$chr}{$pos1}\t$vals{$chr}{$pos1}\n";
    }
}




sub mycmp{
    ($b=~/\D/)?( ($a=~/\D/)?($a cmp $b):-1 ):( ($a=~/\D/)?1:($a <=> $b) );
}
