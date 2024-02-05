#! /usr/bin/env perl
use IO::File;
use Getopt::Std;
use File::Basename;

#############################################################################################
# Arun Manoharan
# 11/30/11
# This converts opticall output to Plink format
# the usage is perl convert_opticall_to_plink.pl -i <inputfile> -o <outputfile>
#############################################################################################

getopts("i:o:");
my $inputfile="";
my $outputfile="";
if(defined $opt_i && defined $opt_o){
    $inputfile=$opt_i;
    $outputfile=$opt_o;
    }else{die("needs input parameters. Usage ex: perl convert_opticall_to_plink -i blood-results.txt -o blood\n");}

print $outputfile,"\t",$inputfile,"\n";
#############################################################################################
#To create the ped file
#line1 is snp, line2 is SNP1, line3 is Allele and genotype calls per sample start from line 4
#ped - sample, sample, 0,0,0,0,2GT for each SNP
#map- chr,SNP,0,pos
#############################################################################################
my @snp="";
my @allele="";
my @line="";
my @gt="";

#############################################################################################
# writing ped file
#############################################################################################
my $outfile=$outputfile.".ped";
print "writing ped file:",$outfile,"\n";
# The input file contains, SNP,SNP1,ALLELE and one GT for each sample
open FH,"<$inputfile" or die "cannot open file $!";
open OFH,">$outfile" or die "cannot open file $!";
my $linecnt=0;
foreach(<FH>)
{
	$line=$_;
	if($linecnt==2)
	{
		@allele=split(/,/,$line); 
#		print "no of samples:",scalar(@allele),"\n";
	}
	if($linecnt>2)
	{
		@gt=split(/,/,$line); 
		#print length(@gt),"\n";
		print OFH $gt[0],"\t",$gt[0],"\t0\t0\t0\t0";
		for($i=1;$i<@gt;$i++)
		{
			if($gt[$i]==1)
			{
				print OFH "\t",substr($allele[$i],0,1),"\t",substr($allele[$i],0,1);
			}
			if($gt[$i]==2)
			{
				print OFH "\t",substr($allele[$i],0,1),"\t",substr($allele[$i],1,1);
			}
			if($gt[$i]==3)
			{
				print OFH "\t",substr($allele[$i],1,1),"\t",substr($allele[$i],1,1);
			}
			if($gt[$i]==4)
			{
				print OFH "\t0\t0";
			}
		}
		print OFH "\n";
	}
	$linecnt++;
}
print "linecnt: $linecnt \n";
close(FH);
close(OFH);

#############################################################################################
#To create the map file
#map- chr,SNP,0,pos
#############################################################################################
my $mapfile="MW_Embark.map";
$outfile=$outputfile.".map";
my %map="";
	open FH,"<$mapfile" or die "cannot open file $!";
	foreach(<FH>)
	{
		@a=split(/\t/,$_);
		$map{$a[1]}=$_;
	}	

my $snpheader=`head -1 $inputfile`;
chomp($snpheader);
@snp=split(/,/,$snpheader);
my $headerlen=scalar(@snp);
#print $snp[0],$snp[196524],$headerlen,scalar(@snp),"\n";
#print "test:",$map{"seq-rs56191141"};

print "writing map file:",$outfile,"\n";
open OFH1,">$outfile" or die "cannot open file $!";
for($i=1;$i<$headerlen;$i++)
{
	#$line= `cat $bimfile|grep -w $snp[$i]`;
	print OFH1 $map{$snp[$i]};
}
close(OFH1);