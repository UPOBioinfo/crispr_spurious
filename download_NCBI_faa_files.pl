#!/usr/bin/perl
use strict;
use LWP::Simple;
use LWP::UserAgent;

#usage: perl download_NCBI_faa_files.pl

#Save all identifiers files
my @files=`ls x*`;

foreach my $var(@files){
	chomp $var;
	open (IN, "$var");

	my $OUT=$var.".fasta";
	open (OUT, ">$OUT");
	print "Downloading $var\n";

	my 	@acc_array="";

	while (my $linea=<IN>){
		chomp $linea;
		push @acc_array,$linea;
	}

	#append [accn] field to each accession
	for (my $i=0; $i < @acc_array; $i++) {
	   $acc_array[$i] .= "[accn]";
	}

	#join the accessions with OR
	my $query = join('+OR+',@acc_array);

	#assemble the esearch URL
	my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base . "esearch.fcgi?db=nuccore&term=$query&usehistory=y";

	#post the esearch URL
	my $output = get($url);

	#parse WebEnv and QueryKey
	my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

	#assemble the efetch URL (fasta_cds_aa)
	$url = $base . "efetch.fcgi?db=nuccore&query_key=$key&WebEnv=$web";
	$url .= "&rettype=fasta_cds_aa&retmode=text";

	#post the efetch URL
	my $fasta = get($url);
	print OUT "$fasta";

	print "Download complete $var\n";
}

