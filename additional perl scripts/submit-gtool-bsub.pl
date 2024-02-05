#! /usr/bin/env perl
###############################################################################################
# Arun Prasad Manoharan
# This script submits jobs to the bsub scheduler in the cluster for each input file chunk
# Each individual job calls the GTOOL program which converts the impute2 results back to plink
# This is required since the number of input samples is large and running in one system will
# take a long time
###############################################################################################
use IO::File;
use Getopt::Std;
use File::Basename;


$dir="imputeRes";
opendir(DIR, $dir) or die "Can't open $dir: $!";
my @files = grep {/impute2$/} readdir DIR;
$numfiles=scalar(@files);

###################################################################
# create the individual bsub jobs for each file chunk
###################################################################
for ($i=0;$i<$numfiles;$i++){

    
    $qoutf = "$dir/bsublog/$files[$i].bsub";
    $idf = "$dir/bsublog/jobids.txt";
    open OFHID,">>$idf" or die "$!";
    open OFH, ">$qoutf" or die "$!";

print "running file", $i," out of $numfiles & file is $files[$i] \n";

@qscript = <<END_SRC;
#BSUB -J gtool
#BSUB -q medium
#BSUB -oo gtool.o%J
#BSUB -eo gtool.e%J
#BSUB -o bsublog/
#BSUB -e bsublog/
#BSUB -n 1 
#BSUB -R "rusage[mem=10]"
cd \$LS_SUBCWD
	
END_SRC
	
print OFH @qscript;
close OFH;

open OFH, ">>$qoutf" or die "$!";

print OFH "\n"; 

print OFH "
# This calls the gtool 
/gne/home/manohara/software/gtool -S --g $dir/$files[$i] --s results/EA.qced-SLE.sample --sample_id imputeRes/EA.SLE.swedish  --inclusion results/passing-snps.txt  --os results/imputed/$files[$i].SLE.swedish.sample  --og results/imputed/$files[$i].SLE.swedish.gen
/gne/home/manohara/software/gtool -S --g $dir/$files[$i] --s results/EA.qced-SLE.sample --sample_id imputeRes/EA.SLE.noswedish  --inclusion results/passing-snps.txt  --os results/imputed/$files[$i].SLE.noswedish.sample  --og results/imputed/$files[$i].SLE.noswedish.gen
\n
";

print OFH "\n";

close OFH;

###################################################################
# submit the bsub jobs
###################################################################

	$jobid = `/opt/lsf/bin/bsub \< $qoutf`;
	print OFHID $jobid;
	print $jobid;


close OFHID;
}