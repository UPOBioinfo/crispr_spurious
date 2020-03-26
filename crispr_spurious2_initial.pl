#!/usr/bin/perl
use strict;

##  crispr_spurious2_initial scanner
##  Created by Pablo MIER on November 2019
##  Copyright © 2019 Pablo MIER. All rights reserved.

my $minimumTR = 7;			my $maximumTR = 20;			#Minimum length of TR					#Maximum length of TR
my $minimumSP = 7;			my $maximumSP = 20;			#Minimum length of spacer			#Maximum length of spacer
#my $minimumUN = 2;			my $maximumUN = 4;			#Minimum number of repeats		#Maximum number of repeats
		
my $file = $ARGV[0];			#First argument -> Input file
my $output = $ARGV[1];		#Second argument -> Output file

open (OUT,">>$output"); print OUT "ID\tNrRepeatedUnits\tRepeatedUnit\tUnitLength\tNrDiffAA\n";			#Print headline in output file 		
#my $mul = Bio::SeqIO->new( -file  => "<$file" , '-format' => 'fasta');															#Input into multifasta file
#my $ind = 0;

open (IN, "$file");

my %hash;
my $id;
my $seq;

while (my $line=<IN>){
	chomp $line;
	if ($line=~/^>(.+)/){
		$id=$line;
		$seq="";
	} else {
		$seq.=$line;
		$hash{$id}=$seq;
	}
}



foreach my $aa (keys %hash) {  																																#Go over the multifasta file, each sequence at a time
	my $prot_name = $aa;	my $prot_seq = $hash{$aa};		my $length = length($prot_seq);			#PROTEIN INFORMATION
	#$ind++; print "$ind\n";
	my $ya = '';																			#To account for considered units in current protein
	for (my $in = 0; $in <= $length; $in++) {					#Go over the complete length of the sequence
		my $extra = 0; my $cc = 0;											#Auxiliary variables
		while ($cc >= 0) {															#Exit variable
			my $substr = substr($prot_seq,$in,$minimumTR+$extra);					#Get seq from current position (in) till position in + /$minimumTR/ (minimum length of TR) + /$extra/ (additional aa when TR is found)
			if 		((length($substr) < $minimumTR) || (length($substr) > $maximumTR)) 	{		$in = $length+1; 	$cc = -1;	}			#If current substring is shorter or larger than allowed, finish the search
			else {																				#The substring has an allowed length 
				my $count = () = $prot_seq =~ /$substr/g;		#Check how many times the current substring is present in the full sequence
				if ($count > 1) {														#If it is part of a TR, more than one; extend the TR
					$cc++;		$extra++;												#Number of times a TR has been extended
				} else {	$cc = -1; 	}											#Not repeated if not part of a TR, or the TR already finished; exit the loop
			}	
			my $exitloop = $in+$minimumTR+$extra+1;				#End coordinate of next iteration of the loop
			if ($exitloop>$length) { $cc = -1;	} 				#Next iteration would not be possible because the sequence would be finished; end the loop
		}
		if ($extra > 0) {																								#There is at least one repetition
			my $substr2 = substr($prot_seq,$in,$minimumTR+$extra-1);			#TR from current position till the last coordinate minus 1
			my $number = () = $prot_seq =~ /$substr2/g;										#Check how many times the TR is present in the full sequence
			$in = $in + $extra;																						#Advance the complete repeat and re-start the search in the end of the current repeat
			next if ((length($substr2) < $minimumTR) || (length($substr2) > $maximumTR));		#Check that the TR has an allowed length
			my $number2 = () = $ya =~ /$substr2/g;												#Check if the repeated unit has already been found
			if ($number2 == 0) {																					#If the repeated unit is new, check that it is interspaced correctly
				my @check = split(/$substr2/,$prot_seq);										#Cut the complete sequence by the repeated unit
				my $num_check = $#check;																		#Number of spacers separated by the repeated unit, plus beginning (first) and end (last) of the sequence
				my $incorrect = 0;																					#Variable to count spacers too short or too long
				for (my $aux = 1; $aux < $num_check; $aux++) {							#Go over all spacers (not the first and the last fragment, they are not spacers)
					my $length_spacer = length($check[$aux]);									#Current spacer length
					if (($length_spacer < $minimumSP) || ($length_spacer > $maximumSP)) {		$incorrect++; }		#If spacer is too short or too long, count it
				}
				my $ll = length($substr2);																	#Unit length
				if ($incorrect == 0) {																			#If there are no incorrect spacers, take the repeated unit as true result

				#Calculate number of different amino acids in repeated unit
				my $count = 0; my $diff = 0;
				$count = () = $substr2 =~ /A/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /C/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /D/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /E/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /F/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /G/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /H/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /I/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /K/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /L/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /M/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /N/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /P/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /Q/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /R/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /S/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /T/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /V/g; if ($count > 0) { $diff++;} $count = 0;
				$count = () = $substr2 =~ /W/g; if ($count > 0) { $diff++;} $count = 0;	$count = () = $substr2 =~ /Y/g; if ($count > 0) { $diff++;} $count = 0;

					print OUT "$prot_name\t$number\t$substr2\t$ll\t$diff\n";					#AC		NrRepeatedUnits		RepeatedUnit		UnitLength	
			}	}
			$ya = $ya."$substr2/";																				#Add considered unit to the found ones for the current protein
}	}	}
exit;
