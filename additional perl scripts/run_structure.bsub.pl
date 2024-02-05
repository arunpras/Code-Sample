#! /usr/bin/env perl
###############################################################################################
# Arun Prasad Manoharan
# This script submits jobs to the bsub scheduler in the cluster for each input file chunk
# Each individual job calls the STRUCTURE program which performs population classification
# This is required since the number of input samples is large and running in one system will
# take a long time
###############################################################################################
use IO::File;
use Getopt::Std;
use File::Basename;


$dir=".";
opendir(DIR, $dir) or die "Can't open $dir: $!";
my @files = grep {/ped$/} readdir DIR;
$numfiles=scalar(@files);

for ($i=0;$i<$numfiles;$i++){    
    $qoutf = "$dir/bsublog/sample-split$i.bsub";
    $idf = "$dir/bsublog/jobids.txt";
    open OFHID,">>$idf" or die "$!";
    open OFH, ">$qoutf" or die "$!";

print "running file", $i," out of $numfiles & file is $files[$i] \n";

#############################################################
#This writes out the parameters required for the scheduler
#############################################################

@qscript = <<END_SRC;
#BSUB -J structure
#BSUB -q medium
#BSUB -oo structure.o%J
#BSUB -eo structure.e%J
#BSUB -o split/bsublog/
#BSUB -e split/bsublog/
cd \$LS_SUBCWD
	
END_SRC
	print OFH @qscript;
	close OFH;

open OFH, ">>$qoutf" or die "$!";

print OFH "\n"; 

#############################################################
# Generate the structure command file
#############################################################

if ($files[$i] =~/325/)
{
print "last file is:", $i,"\n";
print OFH "
cat $files[$i]|awk \'{\$1=NR;
# here the hardcoding sets the population codes for our spiked in hapmap values
if(NR <=13){\$2=0\;\$3=0\;}
else if (NR>=14 && NR <=73){\$2=1\;\$3=1\;}
else if (NR>=74 && NR <=133){\$2=2\;\$3=1\;}
else if (NR>=134 && NR <=223){\$2=3\;\$3=1\;}
\$4=\$5=\$6=\"\"\;gsub(FS\"+\",FS)
}1\' | tr \"ATCG\" \"1234\" > $files[$i].structure.in
\n
# This calls the structure 
../../structure -i $files[$i].structure.in -o $files[$i].structure.out -K 3 -m mainparams-sample-split-force-pop -e extraparams-force-pop -L 24455 -N 223 > $files[$i].k3.log
\n
";
}
else{
print OFH "
cat $files[$i]|awk \'{\$1=NR;
# here the hardcoding sets the population codes for our spiked in hapmap values
if(NR <=50){\$2=0\;\$3=0\;}
else if (NR>=51 && NR <=110){\$2=1\;\$3=1\;}
else if (NR>=111 && NR <=170){\$2=2\;\$3=1\;}
else if (NR>=171 && NR <=260){\$2=3\;\$3=1\;}
\$4=\$5=\$6=\"\"\;gsub(FS\"+\",FS)
}1\' | tr \"ATCG\" \"1234\" > $files[$i].structure.in

\n
# This calls the structure 
../../structure -i $files[$i].structure.in -o $files[$i].structure.out -K 3 -m mainparams-sample-split-force-pop -e extraparams-force-pop -L 24455 -N 260 > $files[$i].k3.log

\n
";

}

print OFH "\n";

close OFH;

#############################################################
# submit to scheduler
#############################################################
	$jobid = `/opt/lsf/bin/bsub \< $qoutf`;
	print OFHID $jobid;
	print $jobid;


close OFHID;
}