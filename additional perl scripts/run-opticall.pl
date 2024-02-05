#! /usr/bin/env perl
use IO::File;
use Getopt::Std;
use File::Basename;

########################################################
# Arun Manoharan
# 9/29/11
# This submits the jobs which will run opticall to qsub scheduler
# the usage is perl run-opticall.pl -i <inputfile> -n <splitsize>
########################################################
#    $qoutf = "qsublog/plink.qsub";
getopts("i:n:");
my $inputfile="";
my $splitsize="";
if(defined $opt_i){
    $inputfile=$opt_i;
    $splitsize=$opt_n;
    }else{die("need input project name. Usage ex: perl run-opticall.pl -i MW_Embark -n 100\n");}

###########################################################################
# This splits the input file into 1000 SNPs each
###########################################################################

print $inputfile,"\n";
my $path= 'opticall-run2';
my $outputpath = 'opticall-run2/split';
my $infile=$path."/".$inputfile."_intensity.txt";
print $infile,"\n";
my $header = `head -1 $infile`;

open FH,"<$infile" or die "cannot open file $!";
my $linecnt=0;
my $fileidx=0;
my $outfile=$outputpath."/".$fileidx."_".$inputfile."_intensity.txt";
print $outfile;
open OFH, ">$outfile" or die "cannot open file $!";
while(<FH>){
	if($linecnt%$splitsize==0)
	{	close(OFH);
		$fileidx++;
		$outfile=$outputpath."/".$fileidx."_".$inputfile."_intensity.txt";
		open OFH, ">$outfile" or die "cannot open file $!";
		print OFH $header;
	}
	$linecnt++;
	print OFH $_; 
}
close(FH);
close(OFH);

###########################################################################
# This submits qsub jobs to run opticall for each of the individual split file
###########################################################################

my $infile=$inputfile."_intensity.txt";
 opendir (DIR,$outputpath) or die "cannot open directory";
 my @files = grep {/$infile$/} readdir DIR;
 $numfiles=scalar(@files);
 print $numfiles,"\n";
for $file(@files){

      $qoutf = "$outputpath/bsublog/$file.bsub";
      $idf = "$outputpath/bsublog/jobids.txt";
      open OFH, ">$qoutf" or die "$!";
  	$qscript = "
 
	#BSUB -J opticall
	#BSUB -q medium
	#BSUB -oo opticall.o%J
	#BSUB -eo opticall.e%J
	#BSUB -o opticall-run2/split/bsublog/
	#BSUB -e opticall-run2/split/bsublog/
	cd \$LS_SUBCWD
	";

  	print OFH $qscript;
  	close OFH;
  open OFH, ">>$qoutf" or die "$!";
  print OFH "\n"; 
  print OFH "
  # This calls the opticall 
  ~/software/opticall/opticall -nointcutoff -noblank -in $outputpath/$file -out $outputpath/$file 
  \n";
 print OFH "\n";
  	close OFH;
   open OFHID,">>$idf" or die "$!";
 	$jobid = `/opt/lsf/bin/bsub \< $qoutf`;

 	print OFHID $jobid;
 	print $jobid;
   close OFHID;
 }
# 
# 
# 
