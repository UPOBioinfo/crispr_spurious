#!/usr/bin/perl
use strict;

#usage: perl crispr_spurious_pfp_multifasta_gff3.pl trembl_bacteria.fasta_def

#open downloaded ncbi file

my $file = $ARGV[0];

open (IN, "$file");

my $id;
my $fasta;
my $gff3;
my $gen;
my $pos1;
my $pos2;

#Create all fasta files

while (my $line=<IN>){
	chomp $line;
	if ($line=~/^>(.+)/){
		$id=$line;
		$id=~s/>lcl\|//g;
		$id=~s/\.(.+)//g;
		$fasta=$id.".fasta";
		open (OUT, ">>$fasta");
		print OUT "$line\n";
	} if ($line=~/^[A-Z](.+)/){
		print OUT "$line\n";
	} if ($line=~/^\n/){
		close OUT;
	}
}


my @list= `ls *.fasta`;

#Create all gff3 files

foreach my $var(@list){
	chomp $var;
	$gff3=$var;
	$gff3=~s/fasta/gff3/g;
	#print "$gff3\n";
	open (IN, "$var");
	open (OUT, ">$gff3");
	my $id_2=$var;
	$id_2=~s/.fasta//g;
	print OUT "$id_2\tGenbank\tregion\t1\t10000000\t\.\t+\t\.\n"; #Minimum and Maximum length
	while (my $line=<IN>){
		chomp $line;
		if ($line=~/^>(.+)/){
			$id=$line;
			$id=~s/>lcl\|//g;
			$id=~s/\.(.+)//g;
			$gen=$line;
			$gen=~s/(.+)locus_tag=//g;
			$gen=~s/(.+)gene=//g;
			$gen=~s/](.+)//g;
			$pos1=$line;
			$pos1=~s/(.+)location=complement\(//g;
			$pos1=~s/(.+)location=//g;
			$pos1=~s/join\(//g;
			$pos1=~s/\>//g;
			$pos1=~s/\<//g;
			$pos2=$pos1;
			$pos1=~s/\.\.(.+)//g;
			$pos2=~s/(.+)\.\.//g;
			$pos2=~s/\](.+)//g;
			$pos2=~s/\)(.+)//g;
			$pos2=~s/\)//g;
			print OUT "$id\tGenbank\tCDS\t$pos1\t$pos2\t\.\t+\t\.\tgen=$gen;locus_tag=$gen;\n";
		}
	}
	close IN;
	close OUT;
}
