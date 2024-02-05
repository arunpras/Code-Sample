#! /usr/bin/env perl
use strict;

###############################################################################################
# Arun Prasad Manoharan
# This script creates the intensity files for Opticall from the raw Illumina Immunochip file
# The raw illumina file contains too many columns while opticall only requires a few columns
# Hence this script cuts the necessary columns from the input raw file to be submitted to 
# opticall genotype caller
###############################################################################################

###############
# usage
###############

die "usage: mkOpticallFiles.pl <illumina file> <output prefix>\n" if(scalar(@ARGV) != 2);


my ($illumina_file, $output_file) = @ARGV;
my ($snpcol, $samplecol, $allelecol, $xcol, $ycol) = ("SNP Name","Sample ID","SNP","X","Y");
my ($snpidx,$sampleidx,$alleleidx,$xidx,$yidx)=("","","","","");

###############
# Get field index
###############
my $header = `head $illumina_file |grep "SNP Name"`;
my @header = split /\t/,$header;
chomp @header;
foreach my $i(0..$#header){
    print "$header[$i]\n";
	$snpidx = $i if($header[$i] eq $snpcol);
    $sampleidx = $i if($header[$i] eq $samplecol);
    $alleleidx = $i if($header[$i] eq $allelecol);
    $xidx = $i if($header[$i] eq $xcol);
    $yidx = $i if($header[$i] eq $ycol);
}

print $snpidx,"\t",$sampleidx,"\t",$alleleidx,"\t",$xidx,"\t",$yidx,"\n";

##############################
# Pull out data from the file
##############################


open FH,"<$illumina_file" or die "cannot open file $!";
my (%x, %y,%allele, %snplist, %samplelist);
my $totalGt = 0;
my $linecnt=0;
foreach(<FH>){
	$linecnt++;
    if($linecnt>10){
    	chomp; my @s=split /\t/,$_; 
    	my $sample=$s[$sampleidx];
    	my $snp=$s[$snpidx];
    	my $x=$s[$xidx];
    	my $y=$s[$yidx];
    	my $allele=$s[$alleleidx];
    	$allele=~s/\[//;
    	$allele=~s/\]//;
    	$allele=~s/\///;
    	chomp($allele);
		if (length($sample)>0 && $sample!~/sample/ && $sample!~/Sample/){
			$x{$sample}{$snp} = $x;
			$y{$sample}{$snp} = $y;
			$allele{$snp} = $allele;
			$samplelist{$sample} = $sample;
			$snplist{$snp} = $snp;
		}
	}
}
close FH;

print "size of hash:  " . keys( %samplelist ) . ".\n";
while( my ($k, $v) = each %samplelist ) {
       if(length($k)>0){ print "key: $k, value: $v.\n";}
}
print "size of hash:  " . keys( %snplist ) . ".\n";
print "size of hash:  " . keys( %allele ) . ".\n";
 
 

##############################
# Print to output file
##############################
 open OUTFH,">$output_file" or die "cannot open file $!";

#print header
 my $line ="SNP\tSNP1\tALLELE";
 my $samplecnt=0;
 while (my ($samplekey, $sample) = each(%samplelist))
 {
 if (length($sample)>0 && $sample!~/sample/ && $sample!~/Sample/){
 	$line.= "\t".$sample."_A\t".$sample."_B";
 	$samplecnt++;
 	}
 }
 print "sample count:",$samplecnt,"\n";
 print OUTFH $line,"\n";

#print intensities
 my $snpcnt=0;
 while (my ($snpkey, $snp) = each(%snplist))
 {
 	$snpcnt++;
 	$line= $snp."\t".$snp."\t".$allele{$snp};
 	while (my ($samplekey, $sample) = each(%samplelist))
 	{
 		if (length($sample)>0 && $sample!~/sample/ && $sample!~/Sample/){
 			$line.= "\t".$x{$sample}{$snp}."\t".$y{$sample}{$snp};
 		}
 	}
 	print OUTFH $line,"\n";
 }
 print "snp count:",$snpcnt,"\n";

close(OUTFH);
