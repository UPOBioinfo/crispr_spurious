#!/usr/bin/perl
use strict;

#usage: perl crispr_spurious_pfp.pl uniprot_database candidates maximun_distances crispr_database minimun_identity minimun_coverage
#usage: perl crispr_spurious_pfp.pl uniprot_sprot_archaea.dat uniprot_sprot_archaea.tsv 15000 /home/alex/Escritorio/Proyecto_CRISPR/script_CRISPR/smp/ 25 70

my $file_uniprot = $ARGV[0]; #uniprot database

if ($file_uniprot=~/(.+)\.dat/){
	&parse_uniprot($file_uniprot);
	$file_uniprot=~s/.dat/_RM.txt/g;
}

my $file_result_i = $ARGV [1]; #candidates file

my $inter= $ARGV [2]; #maximum distance

my $path_db= $ARGV[3]; #path to crispr database

my $ID = $ARGV[4] || 25; # identity

my $SC = $ARGV[5] || 70; # subject (CDD Domain) coverage

####Candidates selection

print "Parsed database and candidates selection\n";

my @list_embl= &select_candidates($file_uniprot, $file_result_i);

####Download EMBL files

&download_file(@list_embl);

print "Download EMBL files completed\n";

####Identifiers files creation

&create_sub_file($file_uniprot, $file_result_i, $inter);

####Multifasta builder

&getfasta($file_uniprot, $file_result_i);

print "Build Multifasta completed\n";

####Build crispr database

&build_crisprdb($path_db);

####Blast search against crispr database

&casfinder($file_uniprot, $file_result_i);

print "Casfinder completed\n";

&filterRPSBlast($file_uniprot, $file_result_i, $ID, $SC);




sub parse_uniprot{

#In order to select specific fields

	my $uniprot=shift;

	my $table_resume=$uniprot;

	$table_resume=~s/.dat/_RM.txt/g;

	open (UNIP, ">$table_resume");

	open (IN, "$uniprot");

	my $iden;
	my $date;
	my $locus;
	my $embl;
	my $pfam;
	my $seq;
	my $phylum;

	while (my $line=<IN>){
		chomp $line;
		if ($line=~/^ID\s{3}(.+)/){
			$iden = $line;
			$iden =~ s/ID\s{3}//g;
			$iden =~ s/\s{1,}(.+)//g;
			$seq="";
			$date="";
			$locus="";
			$embl="";
			$phylum="";
			$pfam="";
		}
		if ($line=~/DT(.+)integrated(.+)into(.+)/){
			$date = $line;
			$date =~ s/,(.+)//g;
			$date =~ s/(.+)-//g;
		}
		if ($line=~/DT(.+)sequence version 1/){
			$date = $line;
			$date =~ s/(.+)-//g;
			$date =~ s/,(.+)//g;
		}
		if (($line=~/^OC\s{3}Archaea(.+)/) or ($line=~/^OC\s{3}Bacteria(.+)/)){
			$phylum=$line;
			$phylum=~s/OC\s{3}//g;
			$phylum=~s/;(.+)//g;
		}
		if ($line=~/^GN\s{3}(.+)/){
			$locus .= $line;
		}
		if ($line=~/DR\s{3}EMBL(.+)/){
			$embl = $line;
			$embl =~ s/DR(.+)EMBL;\s//g;
			$embl =~ s/;(.+)//g;
		}
		if ($line=~/DR(.+)Pfam(.+)/){
			$pfam = $line;
			$pfam =~ s/DR(.+)Pfam;\s//g;
			$pfam =~ s/;(.+)//g;
		}
		if ($line=~/^\/\//){
			$line=~s/\s//g;
			$locus =~ s/GN(.+)OrderedLocusNames=//g;
			$locus =~ s/GN(.+)ORFNames=//g;
			$locus =~ s/GN(.+)Name=//g;
			$locus =~ s/;(.+)//g;
			$locus =~ s/;//g;
			$locus =~ s/,(.+)//g;
			$locus =~ s/\s(.+)//g;
			if ($phylum=~/Archaea/ or $phylum=~/Bacteria/){
				print UNIP "$iden\t$locus\t$embl\t$date\t$pfam\t$phylum\n";
			}
		}
	}
	close IN;
	close OUT;

}

sub dat_extract{

#To create a multi dimensional hash
	
	my $parsed_uniprot=shift;

	open (IN, "$parsed_uniprot");

	my $iden;
	my $date;
	my $locus;
	my $embl;
	my $pfam;
	my %hash_filter;

	while (my $line=<IN>){
		chomp $line;
		my @list= split (/\t/, $line);
		$iden=$list[0];
		$locus=$list[1];
		$embl=$list[2];
		$date=$list[3];
		$pfam=$list[4];
		$hash_filter{$iden}{"PFAM"}=$pfam;
		$hash_filter{$iden}{"LOCUS"}=$locus;
		$hash_filter{$iden}{"EMBL"}=$embl;
		$hash_filter{$iden}{"DATE"}=$date;
	}
	return %hash_filter;
	close IN;

}

sub extract_id_result{

	my $result= shift;

	open (IN, "$result");

	my $id;
	my $repeat;
	my %hash_id;
	my @lista_id;

	while (my $line=<IN>){
		chomp $line;
		next if ($line=~/^ID(.+)NrRepeatedUnits(.+)/);
		my @lista= split (/\t/, $line);
		$id=$lista[0];
		$repeat=$lista[3];
		push @lista_id, $id;
		$hash_id{$id}{"REPEAT"}=$repeat
	}
	return %hash_id;
	close IN;
}

sub select_candidates{

	my @EMBL_list;
	my $uniprot=shift;
	my $candidates=shift;
	my $table=$candidates;
	$table=~s/.tsv/_resume01.txt/g;
	open (OUT, ">$table");
	my %db_filter = &dat_extract($uniprot);
	my %result = &extract_id_result($candidates);

	print STDERR "TSV resume\n";
	foreach my $var (keys %result){
		print OUT "$var\t$db_filter{$var}{LOCUS}\t$db_filter{$var}{EMBL}\t$db_filter{$var}{DATE}\t$db_filter{$var}{PFAM}\n";
		push @EMBL_list, $db_filter{$var}{EMBL};
	}
	return @EMBL_list;
	close OUT;

}

sub download_file{

	#In order to download fasta and gff3 files from NCBI (máx 300 queries)

	my @list=@_;

	open (OUT, ">downloaded_files01.txt");
	my $yes="n";
	foreach my $var (@list){
		my $fasta=$var.".fasta";
		my $gff3=$var.".gff3";
		if (-e $fasta) {
			print "$fasta found\n";
			next;
		}
		`wget -O $gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=$var"`;
		`wget -O $fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta_cds_aa&id=$var"`;
		print "downloaded $fasta\n";
	}
	foreach my $var2 (@list){
		my $fasta=$var2.".fasta";
		my $gff3=$var2.".gff3";
		if (-e $fasta) {
			$yes="y";
			print OUT "$var2\t$yes\n";
		} else {
			print OUT "$var2\t$yes\n";
		}	
	}
}

sub downloaded_files{

	open (IN, "downloaded_files01.txt");
	my %hash_download;
	while (my $line=<IN>){
		chomp $line;
		my @list=split (/\t/, $line);
		$hash_download{$list[0]}=$list[1];
	}
	`rm downloaded_files01.txt`;
	return %hash_download;

}

sub extract_locus{

	my $uniprot=shift;
	my $candidates=shift;
	my %db_filter = &dat_extract($uniprot);
	my %result = &extract_id_result($candidates);
	my %hash_locus;
	foreach my $var (keys %result){
		$hash_locus{$db_filter{$var}{LOCUS}}=$db_filter{$var}{EMBL};
	}
	return %hash_locus;

}

sub extract_coord{

	open (OUT, ">gene_location01.txt");
	my $uniprot=shift;
	my $candidates=shift;
	my %locus= &extract_locus($uniprot, $candidates);
	my %hash_coor;
	foreach my $var (keys %locus){
		my $gff3=$locus{$var}.".gff3";
		open (IN, "$gff3");
		my $locate="n";
		my @lista;
		my $init;
		my $end;
		while (my $linea=<IN>){
			chomp $linea;
			if ($linea=~/\t(.+)region\t(.+)/){
				@lista= split (/\t/, $linea);
				$init=$lista[3];
				$end=$lista[4];
			}
			if ($linea=~/(.+)CDS(.+)$var(.+)/ or $linea=~/(.+)pseudogene(.+)$var(.+)/){
				my @list= split (/\t/, $linea);
				$hash_coor{$var}{"SEQ_S"}=$list[3];
				$hash_coor{$var}{"SEQ_E"}=$list[4];
				$hash_coor{$var}{"ID_GEM"}=$gff3;
				$hash_coor{$var}{"GEM_S"}=$init;
				$hash_coor{$var}{"GEM_E"}=$end;
				$locate="y";
			}
		}
		print OUT "$var\t$locate\n";
	close IN
	}
	return %hash_coor;

}

sub locate_gene{

	open (IN, "gene_location01.txt");
	my %hash_locate;
	while (my $line=<IN>){
		chomp $line;
		my @list=split (/\t/, $line);
		$hash_locate{$list[0]}=$list[1];
	}
	`rm gene_location01.txt`;
	return %hash_locate;

}

sub create_sub_file{

	#To create a identifiers list which are within the distance window

	my $uniprot=shift;
	my $candidates=shift;
	my $intervale= shift;
	my %hash_coor= &extract_coord($uniprot, $candidates);
	my $inter_init;
	my $inter_init_c;
	my $inter_end;
	my $inter_end_c;
	my $complete;
	foreach my $var (keys %hash_coor){
		my $gff3=$hash_coor{$var}{ID_GEM};
		$inter_init=$hash_coor{$var}{SEQ_S} - $intervale;
		$inter_end=$hash_coor{$var}{SEQ_E} + $intervale;
		$inter_init_c=$inter_init;
		$inter_end_c=$inter_end; 
		if ($hash_coor{$var}{GEM_E}>500000){
			$complete="y";
		} elsif ($hash_coor{$var}{GEM_E}<500000){
			$complete="n";
		}
		if ($inter_init<0){
			$inter_init_c= $hash_coor{$var}{GEM_E} + $inter_init;
			$inter_init=1;
		}
		my $id_locus="id_".$var;
		open (OUT2, ">>distance_matrix01.txt");
		open (OUT, ">$id_locus.txt");
		open (IN, "$gff3");
		while (my $linea=<IN>){
			chomp $linea;
			my @list=split (/\t/, $linea);
			if ($list[2]=~/CDS/){
				my $pos_i= $list[3];
				my $pos_e= $list[4];
				my $locus_name= $list[8];
				my $distance;
				my $distance1;
				my $distance2;
				my $mean_distance;
				my $pos_max;
				my $point_m=($hash_coor{$var}{SEQ_S}+$hash_coor{$var}{SEQ_E})/2;
				$locus_name=~s/(.+)locus_tag=//g;
				$locus_name=~s/(.+)gene=//g;
				$locus_name=~s/;(.+)//g;
				$locus_name=~s/;//g;
				if ($inter_init !=1 and $inter_init_c<$pos_i and $inter_end_c>$pos_e){
					print OUT "$locus_name\n";
					$distance1= $point_m-$pos_i;
					$distance1=~s/-//g;
					$distance2= $point_m-$pos_e;
					$distance2=~s/-//g;
					$mean_distance=($distance1+$distance2)/2;
					if ($pos_e<$hash_coor{$var}{SEQ_S}){
						$distance=$hash_coor{$var}{SEQ_S}-$pos_e;
					} if ($pos_i>$hash_coor{$var}{SEQ_E}){
						$distance=$pos_i-$hash_coor{$var}{SEQ_E};
					}
					$distance=~s/-//g;
					print OUT2 "$var\t$locus_name\t$distance\t$complete\n";
				} if ($inter_init ==1 and $inter_init_c<$pos_i){
					print OUT "$locus_name\n";
					$pos_max=$hash_coor{$var}{GEM_E}+$hash_coor{$var}{SEQ_S};
					$distance1=$pos_max-$pos_i;
					$distance2=$pos_max-$pos_e;
					$mean_distance=($distance1+$distance2)/2;
					$distance=$pos_max-$pos_e;
					$distance=~s/-//g;
					print OUT2 "$var\t$locus_name\t$distance\t$complete\n";
				} if ($inter_init ==1 and $inter_end_c>$pos_e){
					print OUT "$locus_name\n";
					$distance1= $point_m-$pos_i;
					$distance1=~s/-//g;
					$distance2= $point_m-$pos_e;
					$distance2=~s/-//g;
					$mean_distance=($distance1+$distance2)/2;
					if ($pos_e<$hash_coor{$var}{SEQ_S}){
						$distance=$hash_coor{$var}{SEQ_S}-$pos_e;
					} if ($pos_i>$hash_coor{$var}{SEQ_E}){
						$distance=$pos_i-$hash_coor{$var}{SEQ_E};
					}
					$distance=~s/-//g;
					print OUT2 "$var\t$locus_name\t$distance\t$complete\n";
				}
			}
		}
		close IN;
		close OUT;
	}
}

sub distance_matrix{
	
	open (IN, "distance_matrix01.txt");
	my %hash_distance;
	while (my $line=<IN>){
		chomp $line;
		my @list=split (/\t/, $line);
		$hash_distance{$list[0]}{$list[1]}{"DISTANCE"}=$list[2];
		$hash_distance{$list[0]}{"SEQUENCE"}{"COMPLETE"}=$list[3];
	}
	`rm distance_matrix01.txt`;
	return %hash_distance;
}

sub getfasta{

	#To create each multifasta file

	my $uniprot=shift;
	my $candidates=shift;
	my %locus= &extract_locus($uniprot, $candidates);

	foreach my $var (keys %locus){
		my $id="id_".$var.".txt";
		my $fasta=$var.".fasta";
		my @list_id;

		open (OUT, ">$fasta");
		open (IN, "$id");

		while (my $linea=<IN>){
			chomp $linea;
			push @list_id, $linea;
		}

		my $ident;
		my $seq;
		my %hash_seq;

		my $EMBL_fasta=$locus{$var}.".fasta";

		open (IN2, "$EMBL_fasta");

		while (my $linea=<IN2>){
			chomp $linea;
			if ($linea=~/^>(.+)/){
				$ident= $linea;
				$ident=~s/(.+)locus_tag=//g;
				$ident=~s/(.+)gene=//g;
				$ident=~s/](.+)//g;
				$ident=~s/;//g;
				$seq="";
		} else {
				$seq.=$linea;
				$hash_seq{$ident}=$seq;
			}
		}

		foreach my $var2 (@list_id){
			print OUT ">$var2\n$hash_seq{$var2}\n";
		}
	}
}

sub build_crisprdb{

	#To build crispr database

	my $path= shift;

	`ls $path*.smp > cdd_crispr.pn`;
	`makeprofiledb -in cdd_crispr.pn -title cdd_crispr -dbtype 'rps'`;

	open (OUTFILE, ">mapping_smp.tsv");

	my @smp = `ls $path`;
	foreach my $file (@smp) {
		chomp $file;
		my $file_path=$path.$file;
  		open (INFILE, "$file_path");
		my $tag_id;
 		 while (my $linea=<INFILE>) {
    		chomp;
    		if ($linea=~/^(.+)tag(.+)id(.+)/) {
				$tag_id=$linea;
				$tag_id=~s/(.+)tag(.+)id\s//g;
				$tag_id=~s/\n//g;
      			print OUTFILE "$tag_id\t";
   			 } elsif ($linea=~/(.+)title(.+)/) {
				my @list= split (/, /, $linea);
    	  		print OUTFILE "$list[1]\n";
   	 		}
  		}
 	 	close INFILE;
	}
	close OUTFILE;
}

sub casfinder{

	#To execute rpsblast against crispr database

	my $uniprot=shift;
	my $candidates=shift;
	my %db_filter = &dat_extract($uniprot);
	my %result = &extract_id_result($candidates);
	foreach my $var (keys %result){
		my $fasta=$db_filter{$var}{LOCUS}.".fasta";
		my $output_blast=$db_filter{$var}{LOCUS}.".blast";
		`rpsblast -query $fasta -db cdd_crispr.pn -evalue 1e-05 -outfmt '6 qseqid sseqid pident qcovs qcovhsp length qlen slen evalue qstart qend sstart send' -max_target_seqs 1 -out $output_blast`;
	}
}

sub filterRPSBlast{

	#To filter results

	my $uniprot=shift;
	my $candidates=shift;
	my $ID=shift;
	my $SC=shift;
	my $yes="n";
	my %db_filter = &dat_extract($uniprot);
	my %result = &extract_id_result($candidates);
	my %hash_distance = &distance_matrix();
	my %hash_download = &downloaded_files();
	my %hash_locate = &locate_gene();
	open (OUT, ">final_candidates.txt");
	print OUT "#UniProt_ID\tLocus_ID\tGenome_ID\tYear\tPfam\tnumber_of_Cas\tgenome\t5000kb\tfound_gene\tCas_genes\tCas_locus_list\tmin_distance\n";
	foreach my $var (keys %result){
		my $output_blast=$db_filter{$var}{LOCUS}.".blast";
		my $cont=0;
		my %smp;
		open SMP, "./mapping_smp.tsv";
		while (<SMP>) {
  			chomp;
			my ($smp, $name) = split/\t/;
  			$smp{$smp} = $name;
		}
		close SMP;

		my $list_Cas="";
		my $list_locus="";
		my @distance="";
		my $min=0;
		my $complete;
		open INFILE, $output_blast;
		while (<INFILE>) {
  			chomp;
  			my ($qseqid, $sseqid, $pident, $qcovs, $qcovhsp, $length, $qlen, $slen, $evalue, $qstart, $qend, $sstart, $send) = split/\t/;
			my ($n) = (split/\|/, $sseqid)[2];
  			my $scov = ($length / $slen) * 100;
			if ($scov >= $SC && $pident >= $ID){
				$list_locus.=$qseqid.",";
				$list_Cas.=$smp{$n}.",";
				if ($qseqid=~/$db_filter{$var}{LOCUS}/){
					next;
				}
				push @distance, "$hash_distance{$db_filter{$var}{LOCUS}}{$qseqid}{DISTANCE}";
				$cont++;
			}
		}
		chop $list_locus;
		chop $list_Cas;
		@distance = sort { $a <=> $b} @distance;
		$min= $distance[1];
		if ($hash_download{$db_filter{$var}{EMBL}}){
			$yes=$hash_download{$db_filter{$var}{EMBL}};
		};
		print OUT "$var\t$db_filter{$var}{LOCUS}\t$db_filter{$var}{EMBL}\t$db_filter{$var}{DATE}\t$db_filter{$var}{PFAM}\t$cont\t$yes\t$hash_distance{$db_filter{$var}{LOCUS}}{SEQUENCE}{COMPLETE}\t$hash_locate{$db_filter{$var}{LOCUS}}\t";
		print OUT  "$list_Cas\t$list_locus\t$min\n";
		close INFILE;
	}
}
