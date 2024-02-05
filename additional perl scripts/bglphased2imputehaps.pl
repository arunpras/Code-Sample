#! /usr/bin/env perl
###############################################################################################
# Arun Prasad Manoharan
# This script converts the BEAGLE phased files to  the haplotype file required by IMPUTE2 
# to perform imputation
###############################################################################################
use MyUtils;
use IO::File;
use Getopt::Std;
use File::Basename;

############################################################################
# usage
############################################################################
getopts("i:o:m:c:");
die "Usage: $0 -i <input bgl phased file> -m <input marker file> -c <chr> -o <output>" unless defined $opt_i and defined $opt_o and defined $opt_c and defined $opt_m;

$FP = new IO::File("<$opt_m") or die "$!";
while(<$FP>){
    chomp; s/\s+/ /g; my @s=split; 
    $posH{$s[0]} = $s[1];
}
close $FP;

############################################################################
# validating user input
############################################################################

if($opt_i =~ /\.gz$/){
    $FP = new IO::File("zcat $opt_i |") or die "$!";
}else{
    $FP = new IO::File("<$opt_i") or die "$!";
}
$opt_o .= ".gz" unless $opt_o =~ /\.gz/;

$OFP = new IO::File("| gzip > $opt_o") or die "$!";
print "writing genotype data to $opt_o\n";

############################################################################
# parse the input phased beagle file
# and generate the SNP marker file
############################################################################
while(<$FP>){

    chomp; s/\s+/ /g; my @s=split; 
    $ity = shift @s;
    if($ity ne "M"){
	shift @s;
	for($i=0; $i <= $#s; $i += 2){
	    push @{ $famH{$ity} }, $s[$i];
	}
	next;
    }

    $rs = shift @s;

    %cH = %aH =();
    foreach $a(@s){
	$aH{$a}++ if $a ne "?";
    }
    $n_alleles = scalar keys %aH;
    if($n_alleles != 2){
	print STDERR "$rs : ",$n_alleles," alleles found. skipping $rs\n";
	next;
    }
    unless(exists $posH{$rs}){
	print STDERR "$rs not found in $opt_m. skipping $rs\n";
	next;
    }
    @an = keys %aH;

    $A = $an[1];
    $B = $an[0];
    if($aH{$an[0]} > $aH{$an[1]}){
	$A = $an[0];
	$B = $an[1];
    }
    $cH{$B} = 1;
    $cH{$A} = 0;
    $cH{"?"} = "?";

    #write out SNP file
	print $OFP "$opt_c $rs $posH{$rs} $A $B ";
    foreach $a(@s){
	print $OFP "$cH{$a} ";
    }
    print $OFP "\n";

#    $i++;

}
close($FP);
close($OFP);

############################################################################
# write out the genotype file
############################################################################
($sampfile = $opt_o) =~ s/.gen.gz$|.gz$|.impute2$|.impute2.gz$/.sample/g;
print "writing sample data to $sampfile\n";
$OFH = new IO::File(">$sampfile") or die "$!";
@ids = @{$famH{"I"}};

print $OFH "ID_1 ID_2 missing\n0 0 0\n";
foreach $i(0..$#ids){
    $iid = $ids[$i];
    if(exists $famH{"P"}){
	$pedig = $famH{"P"}[$i];
    }else{
	$pedig = $iid;
    }
    print $OFH "$pedig $iid 0.0001\n";
}
close $OFH;
