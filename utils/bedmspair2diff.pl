#! /usr/bin/perl
use strict;

## File: bedmspair2diff.pl
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

my $usage = "$0 <file1.bedms> <file2.bedms>\n".
	    "     Computes the differential methylation scores from two samples.\n".
	    "     Inputs are two single-sample .bedms format files, already filtered for\n".
	    "     minimum coverage level.  Outputs a single file in .bedms format.\n".
	    "     Outputs to STDOUT:\n".
	    "          <chrom>    <chromStart>    <chromEnd>     <methdiff (score2-score1)>\n";
die $usage unless @ARGV==2;

# Width to print each site (probably 2)
my $width = 2;
my ($file1,$file2) = @ARGV;
my ($mychr,$lastchr,$line1,$line2,$tempptr);
my @score1;  my @score2;
my %cstarts2 = ();

open FILE1, $file1 or die "Error opening file '$file1': $!\n";
open FILE2, $file2 or die "Error opening file '$file2': $!\n";
print STDERR "File1: $file1\nFile2: $file2\n";

## Locate all chromosome starts in file2
$lastchr = -1;  $tempptr = tell(FILE2);
print STDERR "Reading chromosome starts for bedfile2.";
while(<FILE2>)
{
    m/^chr([^\t]+)\t/;
    $mychr = $1;
    if( $mychr ne $lastchr ){
	$cstarts2{$mychr} = $tempptr;
	$lastchr = $mychr;
	print STDERR ".";
    }
    $tempptr = tell(FILE2);
}
print STDERR "\n";
seek FILE2,0,0;


$line1 = <FILE1>;  chomp($line1);
OUTER: while(1)
{
    $line1 =~ m/^chr([^\t]+)\t/;
    $mychr = $1;
    unless( defined $cstarts2{$mychr} )
    {
	last OUTER unless($line1 = <FILE1>);  chomp($line1);
	next;
    }
    seek FILE2, $cstarts2{$mychr}, 0;
    $line2 = <FILE2>;  chomp($line2);
    print STDERR "Processing chromosome $mychr...\n";

    ## Chromosomes are now aligned; walk through both files
    @score1 = split /\t/, $line1;
    @score2 = split /\t/, $line2;

    INNER: while(1)
    {
	if( $score1[1] eq $score2[1] )
	{
	    print STDOUT $score1[0]."\t".$score1[1]."\t".($score1[1]+$width)."\t";
	    print STDOUT ($score2[3]-$score1[3])."\n";

	    ## Advance both files
	    last OUTER unless($line1 = <FILE1>);  chomp($line1);
	    last INNER unless($line2 = <FILE2>);  chomp($line2);

	    ## Termination criteria (always ensure line1 is on the first line of next chromosome
	    last INNER if( $line1 !~ m/^chr$mychr\t/ );
	    if( $line2 !~ m/^chr$mychr\t/ )
	    {
		do{
		    last OUTER unless($line1 = <FILE1>);  chomp($line1);
		}while( $line1 =~ m/^chr$mychr\t/ );
		last INNER;
	    }

	    @score1 = split /\t/, $line1;
	    @score2 = split /\t/, $line2;
	}
	elsif( $score1[1] > $score2[1] )
	{
	    unless($line2 = <FILE2>)
	    {
		do{
		    last OUTER unless($line1 = <FILE1>);  chomp($line1);
		}while( $line1 =~ m/^chr$mychr\t/ );
		last INNER;
	    }
	    chomp $line2;
	    @score2 = split /\t/, $line2;
	}
	else
	{
	    last OUTER unless($line1 = <FILE1>);  chomp $line1;
	    last INNER if( $line1 !~ m/^chr$mychr\t/ );
	    @score1 = split /\t/, $line1;
	}
    }
}

