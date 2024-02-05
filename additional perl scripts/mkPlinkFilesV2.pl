#! /usr/bin/env perl
###############################################################################################
# Arun Prasad Manoharan
# This script cuts the raw Illumina Immunochip file 
# The original file is a huge text file and we only need a few columns to generate the 
# plink input files. So we cut only the columns needed.
# The file writes out two files - .ped file which contains the sample and genotype information
# .map file which contains the SNP information - SNP name, chr and position 
###############################################################################################
use IO::File;
use strict;

die "usage: mkPlinkFilesV2.pl <illumina file> <output file> <GC score Thr>\n" if(scalar(@ARGV) < 2);
#my $scoreThr = 0.7;
my ($illumina_file, $output_file, $scoreThr) = @ARGV;
$scoreThr = 0.7 if ! defined $scoreThr;
print "$illumina_file $output_file $scoreThr\n";
my $OFP = ();
my $FP = ();

my ($afname1, $afname2, $tmpfn, $af1, $af2,$scoreName) = ("Allele1 - Forward","Allele2 - Forward","tmpfile",-1,-1,"GC Score");

#######################################################################
# This cuts out the header line from the input file
# It gets the index of Allele1, Allele2, chr pos and score fields
#######################################################################
my $header = `head $illumina_file |grep "SNP Name"`;
my @header = split /\t/,$header;
my ($chrf,$posf,$scrf) = (-1,-1,-1);
chomp @header;
foreach my $i(0..$#header){
    print "$header[$i]\n";
    $af1 = $i if($header[$i] eq $afname1);
    $af2 = $i if($header[$i] eq $afname2);
    $chrf = $i if($header[$i] eq "Chr");
    $posf = $i if($header[$i] eq "Position");
    $scrf = $i if($header[$i] eq $scoreName);

}

die "fields $afname1 or $afname2 or $scoreName not found in $illumina_file\n" if $af1 == -1 || $af2 == -1 || $scrf == -1;

#######################################################################
# Generate the SNP and genotype lists
#######################################################################
my @snphead = `head -1000 $illumina_file |sed -n '/SNP Name/,50 p' | grep -v "SNP Name" |cut -f 1`;
my @samphead = `head -1000 $illumina_file | sed -n '/SNP Name/,50 p' | grep -v "SNP Name" |cut -f 2`;
chomp @snphead;
chomp @samphead;

die "Illumina file is in SNP major format .. exiting\n" if($snphead[0] eq $snphead[1]);

$FP = new IO::File("<$illumina_file") or die "$!";
my (%chrH, %posH, @snplist, @samplelist);
my $totalGt = 0;
while(<$FP>){
    chomp; my @s=split /\t/,$_; 
    push @samplelist, $s[1] if($s[0] eq $snphead[0]);
    if ($s[1] eq $samphead[0]){
	push @snplist, $s[0];
	if($chrf != -1 && $posf != -1){
	    $chrH{$s[0]} = $s[$chrf];
	    $posH{$s[0]} = $s[$posf];
	}else{
	    $chrH{$s[0]} = $posH{$s[0]} = 0;
	}
    }
    $totalGt++;
}
close $FP;

#print "SNP $_\n" foreach(@snplist);
#print "SAM $_\n" foreach(@samplelist);

print "chrf $chrf $posf\n";
foreach (keys %chrH){
    print $_,"$chrH{$_}\n";
}
#######################################################################
# write out the map file
#######################################################################
my $mapfile = $output_file . ".map";
open OFH,">$mapfile" or die "$!";
foreach my $snp(@snplist){
    print OFH "$chrH{$snp}\t$snp\t0\t$posH{$snp}\n";
}
close OFH;

#######################################################################
# write out the ped file
#######################################################################
my $pedfile = $output_file . ".ped";
$OFP = new IO::File(">$pedfile") or die "$!";
my $cursamp = ();
my $line = ();
my $c = ();
my $noSpaceID = "";
my $dataflag = 0;
my $finalGt = 0;
$FP = new IO::File("<$illumina_file") or die "$!";
while($line = <$FP>){
    if($dataflag == 1){
	
	my @s = split(/\t/,$line);
	if($s[$scrf] < $scoreThr){
	    $s[$af1] = $s[$af2] = "0";
	}else{
	    $finalGt++;
	}
	if($s[0] eq $snplist[0]){
	    $c = 0;
	    $cursamp = $s[1];
	    print "current $cursamp\n";
	    if($snplist[$c] eq $s[0]){
		if($cursamp ne $samplelist[0]){
		    print $OFP "\n";
		}
		$s[$af1] = "0" if($s[$af1] eq "-");
		$s[$af2] = "0" if($s[$af2] eq "-");
		($noSpaceID = $s[1]) =~ s/\s+/_/g;
		print $OFP "$noSpaceID $noSpaceID 0 0 0 0 $s[$af1] $s[$af2] ";
		$c++;
	    }else{
		die "snp name mismatch expecting $s[$c], found $s[0] $c..\n";
	    }
	}else{
	    if($snplist[$c] eq $s[0] and $s[1] eq $cursamp){
		$s[$af1] = "0" if($s[$af1] eq "-");
		$s[$af2] = "0" if($s[$af2] eq "-");
		print $OFP "$s[$af1] $s[$af2] ";
		$c++;
	    }else{
		die "snp/samp name mismatch expecting $s[$c], found $s[0], samp $cursamp, $s[1]\n";
	    }
	}
    }
    $dataflag = 1 if($line =~ /Allele2/);
}

print "Total records : $totalGt\nAfter filtering base on $scoreName : $finalGt\n";

close $FP;
close $OFP;
