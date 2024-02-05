#! /usr/bin/env perl
use strict;

###############################################################################################
# Arun Prasad Manoharan
# This script cuts the raw Illumina Immunochip file 
# The original file is a huge text file and we only need a few columns to generate the 
# plink input files. So we cut only the columns needed.
###############################################################################################

###############
# usage
###############

die "usage: mkOpticallFiles.pl <illumina file> <output file>\n" if(scalar(@ARGV) != 2);
my ($illumina_file, $output_file) = @ARGV;
my @collist = ("SNPName","SampleID","Allele1-Top","Allele2-Top","Allele1-Forward","Allele2-Forward","Allele1-Design","Allele2-Design","Allele1-AB","Allele2-AB","Allele1-Plus","Allele2-Plus","Chr","Position","SNP","ILMNStrand","CustomerStrand");
#my @colidx;
my @colidx=(1,2,3,4,10,11,12,13,14,15,16,17,18,19,22,23,24);

###############
# Get field index
###############

my $header = `head $illumina_file |grep "SNP Name"`;
$header =~ s/ +//g;
$header =~ s/\t/ /g;
#print $header;

 my @header = split / /,$header;
 chomp @header;

##############################
# Pull out data from the file
##############################


open FH,"<$illumina_file" or die "cannot open file $!";
open OFH,">$output_file" or die "cannot open file $!";
print OFH join(" ", @collist),"\n";
close(OFH);
my $linecnt=0;
foreach(<FH>){
	$linecnt++;
    if($linecnt>10)
    {
    		chomp;
    		my $line=$_;
			my $cmd= "echo $line|awk '{print \$".$colidx[0];    		
			for my $i(1..$#colidx)
    		{
				$cmd.= ",\$".$colidx[$i];
  			}
			$cmd.= "}'>>$output_file";
    		system($cmd);
	}
}
close FH;
