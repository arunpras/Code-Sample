#! /usr/bin/env perl
##########################################################################################
# Written by Arun Manoharan
# This perl script combines the results of the various association results
# The input files are :
# -l: logistic regression output 
# -c: chisquared association output
# -h: hardy weinberg output
#usage: perl mkAssocResTable.pl -l <project>.logistic.assoc.logistic -c <project>.chisq.assoc -h <project>.hwe.hwe -o <project>.assoc.results.txt
##########################################################################################
use IO::File;
use Getopt::Std;

#Get input arguments
# There are two options 1: where the phenotype is discrete 2: the phenotype is continuous
getopts("l:c:h:o:");
if(defined $opt_h){
    print_plink_res_with_hardy($opt_l,$opt_c,$opt_o,$opt_h);
}else{
    print_plink_res($opt_l,$opt_c,$opt_o);
}

##########################################################
# This function reads in the hardy output file and 
# generates a hash of HWE pvalues for each SNP
##########################################################
sub read_plink_hardy{

    $input_fn = shift @_;

    $FP = new IO::File("<$input_fn") or die "$!";
    $tmp = <$FP>;
    chomp $tmp;
    $tmp =~ s/^\s+//;
    @cf = split /\s+/,$tmp;
#    $cf[3] = "RiskAllele";

    while(defined($line = <$FP>)){
	chomp $line;
	$line =~s /^\s+//;
	@f = split /\s+/,$line;
	foreach $i(0..$#cf){
	    $qres_hardyH{$f[1]}{$cf[$i]} = $f[$i];
	}
	foreach my $test(qw/ALL AFF UNAFF/){
	    if($f[2] eq $test){
		$qres_hardyH{$f[1]}{$test} = $f[5]
	    }
	}
    }
    close($FP);

    return \%qres_hardyH;
}
############

##########################################################
# This function reads in the chisquared association  
# output file and generates a hash of pvalues for each SNP
##########################################################
sub read_plink_chisq{
    $input_fn = shift @_;

    $FP = new IO::File("<$input_fn") or die "$!";
    $tmp = <$FP>;
    chomp $tmp;
    $tmp =~ s/^\s+//;
    @cf = split /\s+/,$tmp;
#    $cf[3] = "RiskAllele";

    while(defined($line = <$FP>)){
	chomp $line;
	$line =~s /^\s+//;
	@f = split /\s+/,$line;
	foreach $i(0..$#cf){
	    $qres_chisqH{$f[1]}{$cf[$i]} = $f[$i];
	}
    }
    close($FP);

    return \%qres_chisqH;
}

##########################################################
# This function reads in the logistic regression output  
# file and generates a hash of HWE pvalues for each SNP
##########################################################
sub read_plink_logistic{
    $input_fn = shift @_;

    $FP = new IO::File("<$input_fn") or die "$!";
    $tmp = <$FP>;
    chomp $tmp;
    $tmp =~ s/^\s+//;
    @cf = split /\s+/,$tmp;
#    $cf[3] = "RiskAllele";

    while(defined($line = <$FP>)){
	chomp $line;
	$line =~s /^\s+//;
	next if $line !~ /ADD/;
	@f = split /\s+/,$line;
	foreach $i(0..$#cf){
	    $qresH{$f[1]}{$cf[$i]} = $f[$i];
	}
	push @snpIds,$f[1];
    }
    close($FP);

    return (\%qresH, \@snpIds);
}

##########################################################
# This function generates the output file for discrete
# phenoype. It links the previous three hashes and for 
# each SNP writes out the columns specified in the header
##########################################################
sub print_plink_res_with_hardy{

    ($logiF, $chisqF, $outf, $hardyF) = @_;

    ($logiR, $snpIds) = read_plink_logistic("$logiF");
    %logiH = %{$logiR};
    @snpIds = @{$snpIds};
    %chisqH = %{read_plink_chisq($chisqF)};
    %hardyH = %{read_plink_hardy($hardyF)};

    $OFP = new IO::File(">$outf") or die "$!";
    @header = qw/SNP CHR BP A1 A2 F_A F_U OR L95 U95 SE P A1_g A2_g AFF_g UNAFF_g/;
    print $OFP "$_\t" foreach @header;
    print $OFP "\n";

    foreach my $snp(@snpIds){

	foreach my $fld(qw/SNP CHR BP A1/){
	    print $OFP "$logiH{$snp}{$fld}\t";
	}
	foreach my $fld(qw/A2 F_A F_U/){
	    print $OFP "$chisqH{$snp}{$fld}\t";
	}
	foreach my $fld(qw/OR L95 U95 SE P/){
	    print $OFP "$logiH{$snp}{$fld}\t";
	}
	foreach my $fld(qw/A1 A2 AFF UNAFF/){
	    print $OFP "$hardyH{$snp}{$fld}\t";
	}
	print $OFP "\n";
    }
    close $OFP;

}

##########################################################
# This function generates the output file for continuous
# phenoype. It links the previous hashes and for 
# each SNP writes out the columns specified in the header
##########################################################
sub print_plink_res{
    ($logiF, $chisqF, $outf) = @_;

    ($logiR, $snpIds) = read_plink_logistic("$logiF");
    %logiH = %{$logiR};
    @snpIds = @{$snpIds};
    %chisqH = %{read_plink_chisq($chisqF)};

    $OFP = new IO::File(">$outf") or die "$!";
    @header = qw/SNP CHR BP A1 A2 F_A F_U OR L95 U95 SE P/;
    print $OFP "$_\t" foreach @header;
    print $OFP "\n";

    foreach my $snp(@snpIds){

	foreach my $fld(qw/SNP CHR BP A1/){
	    print $OFP "$logiH{$snp}{$fld}\t";
	}
	foreach my $fld(qw/A2 F_A F_U/){
	    print $OFP "$chisqH{$snp}{$fld}\t";
	}
	foreach my $fld(qw/OR L95 U95 SE P/){
	    print $OFP "$logiH{$snp}{$fld}\t";
	}
	print $OFP "\n";
    }
    close $OFP;

}

##########################################################################################
