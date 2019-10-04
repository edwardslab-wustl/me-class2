#! /usr/bin/perl -w
use strict;

## File: bedms_sane.pl
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
	    "    Performs basic sanity checks on a .bedms file.\n";
die $usage unless @ARGV==1;

my %seenchr = ();
my ($lastchr,$lastpos,$line) = (-1,-1,0);

open FILE, $ARGV[0] or die "Error: $!\n";
print STDERR "Processing bedms file: ";
while( <FILE> )
{
    my ($chr, $pos1, $pos2, $val) = split /[\t\n]/;
    ++$line;
    if( $chr ne $lastchr )
    {
	if( $seenchr{$chr} ){
	    die "\nError (line $line): chr=$chr has at least two separated blocks!\n";
	}
	$seenchr{$chr} = 1;

	$lastchr = $chr;
	$lastpos = -1;
	print STDERR "$chr,";
    }

    unless( $pos1 > $lastpos ){
	die "\nError (line $line): chr=$chr, pos=$pos1 not greater than previous position $lastpos !!\n";
    }

    if( $val > 1 ){
	die "\nError (line $line): Value too large!  chr=$chr, pos1=$pos1, val=$val\n";
    }
    if( $val < -1 ){
	die "\nError (line $line): Value too small!  chr=$chr, pos1=$pos1, val=$val\n";
    }
    if( $pos2 < $pos1 ){
	die "\nError (line $line): Bad pos2 value: chr=$chr, pos1=$pos1, $pos2=$pos2, $val=$val\n";
    }

    $lastpos = $pos1;
}
close FILE;

print STDERR "\nFile is okay!\n";
