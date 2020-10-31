#!/usr/bin/perl
# Domain plotting for candidates CDS
# Wyler Michele, October 31th 2020

use warnings;
use strict;
use English;
use File::Temp qw/ tempdir /;
use Getopt::Long 'HelpMessage';
use File::Basename;
use Cwd;
use POSIX;


### INPUT -------------------------------------------------------------

# getwd
my $dir = getcwd . "/";

# test input on the command line
GetOptions(
    'Fasta=s' => \my $FASTA,
    'outPrefix=s' => \my $PREFIX,
    'help'     =>   sub { HelpMessage(0) },
) or HelpMessage(1);

# check presence of mandatory inputs
HelpMessage(1) unless $FASTA;

=head1 NAME

Annotation tester

=head1 SYNOPSIS

  --Fasta, -F            Genes of interest in fasta format (required)
  --outPrefix, -o        Prefix of the output files, allows redirection (optional)
  
Help
  --help, -h             Print this help
  
=head1 VERSION

0.01

=cut

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
		$DNAsequence = uc($DNAsequence);
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
			print "ERROR: Unable to donwload the Pfam database.\n";
			print "Check the internet connection or copy manually the DB in the current folder.\n";
		exit;
		}
	}
}


### WORKER DEPENDENCY and INPUTS -------------------------------------------------------------

# test if user want redirection of output
if (! defined $PREFIX){
    $PREFIX = $dir;
}

# test if hmm is installed
if ( (`hmmsearch`)[1] !~ /^Usage/ ) {
	print "hmmer seems not installed.\n";
	print "Try 'sudo apt-get install hmmer'.\n";
}

# test if Pfam files are available
&checkPfam;

### Read in Genes
## read input sequence inside an HASH (geneID "\t" fasta)

# unzip if needed
if ($FASTA =~ /.gz$/) {
    open(IN, "gunzip -c $FASTA |") or die "can't open $FASTA";
}
else {
    open(IN, $FASTA) or die "can't open $FASTA";
}
## read input sequence inside an hash

my %seqs = ();
my $header = '';

while (my $line = <IN>){
    chomp $line;
    $line =~ s/^\s*$//;
    if ($line =~ m/^>(.*)$/){
        # remove description from header (after space)
        $header = (split ' ', $1)[0];
    } else {
        $seqs{"$header"} .= $line;
    }
}
close (IN);

## change frames
my %AA6frame;
my $Header;
my $Sequence;

# forward frames
foreach my $FRAME (1, 2, 3){
	while ( ($Header, $Sequence) = each %seqs) {
		my $FrameSeq = substr($Sequence, $FRAME-1);
		my $FrameHeader = $Header . '_' . $FRAME;
		#print ">$FrameHeader\n$FrameSeq\n";
		$AA6frame{"$FrameHeader"} = $FrameSeq;
	}
}
# reverse frames
foreach my $FRAME (1, 2, 3){
	while ( ($Header, $Sequence) = each %seqs) {
		# reverse
		$Sequence = reverse $Sequence;
		# complement
		$Sequence =~ tr/ACGTacgt/TGCAtgca/;

		my $FrameSeq = substr($Sequence, $FRAME-1);
		my $FrameHeader = $Header . '_-' . $FRAME;
		#print ">$FrameHeader\n$FrameSeq\n";
		$AA6frame{"$FrameHeader"} = $FrameSeq;
	}
}

%AA6frame = &DNAtoAA (%AA6frame) ;

## print out all into temp folder
my $TEMPfolderOutput = tempdir( DIR => $dir, CLEANUP => 1 ); 

while ( ($Header, $Sequence) = each %AA6frame) {
	my $FILE = $TEMPfolderOutput . '/fastaAA.fa' ;
        open(FH, '>>', $FILE) or die $!;
	print FH ">$Header\n$Sequence\n";
        close(FH);
}

## Run domain Scan
my $HMMOUT = $TEMPfolderOutput . '/candidate.hmm';
system "hmmsearch --domtblout $HMMOUT -E 1e-5 --cpu 1 Pfam-A.hmm $TEMPfolderOutput/fastaAA.fa";

# parse HMM results
# fgrep -v '#' candidate.hmm | awk '{print $1,$3,$4,$19,$20}'
my @HMMtable;
my $HMM_header = "geneName\tgeneLength\tDomain\tStart\tEnd";
push(@HMMtable, $HMM_header);
 
open(IN, $HMMOUT ) or die "can't open $HMMOUT";
	while(<IN>){
		# don't keep lines with comments
		if ($_ !~ m/^#/){
		$_ =~ s/ +/\t/g;
		# make one array for each row
		my @HMMrow = (split ' ', $_);
		my $ParseHMMrow = "$HMMrow[0]\t$HMMrow[2]\t$HMMrow[3]\t$HMMrow[17]\t$HMMrow[18]";
		push(@HMMtable, $ParseHMMrow);
		}
	}
close(IN);

# print out coordinates
my $DomCoordinates = $TEMPfolderOutput . '/DomainCoordinates.txt';
foreach (@HMMtable) {
	open(FH, '>>', $DomCoordinates) or die $!;
	print FH "$_\n";
	close(FH);
}

# plot 
my $SCRIPTPATH = dirname($0);
my $PlotName = getcwd() . "/". $PREFIX . "DomainPlot.pdf";


system "Rscript $SCRIPTPATH/subScriptPlotting $DomCoordinates $PlotName";





