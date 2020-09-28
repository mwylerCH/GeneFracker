#!/usr/bin/perl

use warnings;
use strict;
use English;
use File::Temp qw/ tempdir /;
use Getopt::Long 'HelpMessage';
use File::Basename;
use Cwd;
use POSIX;
use File::Temp qw(tempfile); 


### INPUT -------------------------------------------------------------

# getwd
my $dir = getcwd . "/";

my $ALIGNMENT = $ARGV[0];
my $NRcores = $ARG[1];

# define hash with codons
my %CodonTable = ('TTT' => 'F',
	'TCT' => 'S',
	'TAT' => 'Y',
	'TGT' => 'C',
	'TTC' => 'F',
	'TCC' => 'S',
	'TAC' => 'Y',
	'TGC' => 'C',
	'TTA' => 'L',
	'TCA' => 'S',
	'TAA' => 'X',
	'TGA' => 'X',
	'TTG' => 'L',
	'TCG' => 'S',
	'TAG' => 'X',
	'TGG' => 'W',
	'CTT' => 'L',
	'CCT' => 'P',
	'CAT' => 'H',
	'CGT' => 'R',
	'CTC' => 'L',
	'CCC' => 'P',
	'CAC' => 'H',
	'CGC' => 'R',
	'CTA' => 'L',
	'CCA' => 'P',
	'CAA' => 'Q',
	'CGA' => 'R',
	'CTG' => 'L',
	'CCG' => 'P',
	'CAG' => 'Q',
	'CGG' => 'R',
	'ATT' => 'I',
	'ACT' => 'T',
	'AAT' => 'N',
	'AGT' => 'S',
	'ATC' => 'I',
	'ACC' => 'T',
	'AAC' => 'N',
	'AGC' => 'S',
	'ATA' => 'I',
	'ACA' => 'T',
	'AAA' => 'K',
	'AGA' => 'R',
	'ATG' => 'M',
	'ACG' => 'T',
	'AAG' => 'K',
	'AGG' => 'R',
	'GTT' => 'V',
	'GCT' => 'A',
	'GAT' => 'D',
	'GGT' => 'G',
	'GTC' => 'V',
	'GCC' => 'A',
	'GAC' => 'D',
	'GGC' => 'G',
	'GTA' => 'V',
	'GCA' => 'A',
	'GAA' => 'E',
	'GGA' => 'G',
	'GTG' => 'V',
	'GCG' => 'A',
	'GAG' => 'E',
	'GGG' => 'G',
	'GTN' => 'V',
   	'TCN' => 'S',
    	'CCN' => 'P',
    	'ACN' => 'T',
    	'GCN' => 'A',
    	'GGN' => 'G',
	'CGN' => 'R');

### Define SUBROUTINES -------------------------------------------------------------

# converts DNA hash to AA hash
sub DNAtoAA {
	# input need to be a hash containing DNA fasta
	my %DNAhash = @_ ;
	my %AAhash ;
	# translate each sequence 
	foreach my $header (keys %DNAhash){
		my $DNAsequence = $DNAhash{$header};
		# get length of the sequence
		my $GeneLen = length $DNAsequence;
		# take one triplet after the other
  		for (my $i=0; $i <= $GeneLen; $i += 3){
        		my $triplet = substr($DNAsequence, $i, 3); # 3 nucleotides
        		if(exists($CodonTable{"$triplet"})) { # remove last nucleotides
        		    $AAhash{"$header"} .= $CodonTable{"$triplet"};
        		    }
        	}
	}
	return %AAhash ;
}

# check presence of Pfam dataset and download in case
sub checkPfam{
	unless (-e 'Pfam-A.hmm') {
		print "No Pfam-A.hmm file detected. \nTry to download and press.\n";
		system 'rm -f Pfam-A.hmm.h3p; rm -f Pfam-A.hmm.h3m ; rm -f Pfam-A.hmm.h3f ; rm -f Pfam-A.hmm.h3i';
		system 'wget --quiet ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.gz -O Pfam-A.hmm.gz';
		if (-e 'Pfam-A.hmm.gz') {
			system 'gunzip Pfam-A.hmm.gz';
			system 'hmmpress Pfam-A.hmm';
		}else{
			print "Unable to donwload the Pfam database.\n";
			print "Check the internet connection or copy manually the DB in the current folder.\n";
		exit;
		}
	}
}


### WORKER ALIGNMENT -------------------------------------------------------------

# test if alignment is fasta format
open(FH, '<', $ALIGNMENT) or die "Can't open '$ALIGNMENT': $!";
my $ALIGNMENTfirstLine = <FH>;
if ($ALIGNMENTfirstLine !~ /^>/) {
	print "ERROR: '--precise' requires a fasta alignment.";
	exit;
}
close(FH);


# read in alignment
my %ALNseqs;
open(FH, $ALIGNMENT) or die "Can't open '$ALIGNMENT': $!";
	## read alignment inside an hash
	my $ALNheader = '';
	while (my $aln = <FH>){
	    chomp $aln;
	    $aln =~ s/^\s*$//; # remove empty lines
	    if ($aln =~ m/^>(.*)$/){
		# remove sequence description from header
		$ALNheader = (split ' ', $1)[0];
	    } else {
		# remove '-'
		$aln =~ s/-//g;
		$ALNseqs{"$ALNheader"} .= $aln;
	    }
	}
close(FH);

# convert to AA
%ALNseqs = &DNAtoAA (%ALNseqs);

# print out AA sequence to tmp file
my ($fh, $TEMPaaFasta) = tempfile ();
open($fh, '>', $TEMPaaFasta) or die $!;
while ( my ($key, $value) = each %ALNseqs ) {
	print $fh ">$key\n$value\n";
}
close($fh);


#### run Pfam scan
# test presence of Pfam DB
&checkPfam;

# output hammer scan in tmp file
my $TEMPhammerOut_Alignment = tempdir( DIR => $dir, CLEANUP => 1 ); 
system "hmmsearch --tblout $TEMPhammerOut_Alignment/hammerOut -E 1e-5 --cpu 1 Pfam-A.hmm $TEMPaaFasta >/dev/null 2>&1";

# read in domains for each gene
my %HammerOutput;
my $FILEhammerout =  $TEMPhammerOut_Alignment . '/hammerOut';
open(FH, $FILEhammerout);
while (<FH>){
	if ($_ !~ '#'){
	        chomp;
	        my @HMMrow = (split ' ', $_);
	        my $GeneName = $HMMrow[0];
	        # add all domain in one row, comma separated
	        if(exists($HammerOutput{"$GeneName"})) {
			$HammerOutput{"$GeneName"} = $HammerOutput{"$GeneName"} . ',' . $HMMrow[2];
        	}else{
            		$HammerOutput{"$GeneName"} = $HMMrow[2];
            	}
	}
}
close(FH);

# extract all domains found into list
my @domains = values %HammerOutput;

# count number of domains over all sequences of the alignment
my %count = ();
foreach my $element (@domains) {
    $count{$element}++;
}


# most common domain
my $MainDomain = (reverse (sort {$count{$a} <=> $count{$b}} keys %count))[0] ;
print $MainDomain;

