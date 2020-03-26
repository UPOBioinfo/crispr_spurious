#!/usr/bin/perl
# crispr_spurious1_initial.pl
# AJPerez, feb 2019

use strict;

# Inputs and variables
my $FASTA  = $ARGV[0] || "./uniprot_sprot_archaea.fasta";
my $CRISPR = $ARGV[1] || "./20190618_dr_34.fasta";
my $ONE = 0 || $ARGV[2]; # 1 = one only result per protein
my $minSP = 7;
my $maxSP = 20;
my $OUTFILE = $FASTA;
$OUTFILE =~ s/\.fasta/\.tsv/;
$| = 1;

# Gather genome sequence
my %p;
my $id;
open fasta, $FASTA || die "Error! Problem opening fasta file\n";
while (<fasta>) {
  chomp;

  if (/^>(.+)/) {
    $id = $1;
  } else {
    $p{$id} .= $_;
  }
}
close fasta;

# Run through repeats
my $n = 0;
my %trip;
open repeats, $CRISPR || die "Error! Problem opening crispr db fasta\n";
while (<repeats>) {
  chomp;

  next if /^>/;
  next if /[^ACGT]/;
  $n++;
  my $sseq = $_;
  for (my $cadena = 0; $cadena <= 1; $cadena++) {
    my @aa;
    my $strand = "+";
    if ($cadena == 1) {
      $sseq =~ tr/ACGT/TGCA/;
      $sseq = reverse $sseq;
      $strand = "-";
      $n++;
    }
    #print "s$n\t";
    for (my $pos = 0; $pos < 3; $pos++) {
      my $pos_6 = $pos;
      $pos_6 += 3 if $strand eq "-";
      my $s = substr($sseq, $pos);
      my $aa = &translate($s);
      next if $aa =~ /X/;
      push @aa, $aa;
    }
    
    # Check tri-peptides from the same strand
    next if @aa == ();
    my $trip = join "|", sort @aa;
    push @{$trip{$trip}}, $sseq;
  }
}
close repeats;

# Run though proteins
my $n = 0;
open FILE, ">$OUTFILE" || die "Error! Problem creating $OUTFILE\n";
protein: for my $id (keys %p) {
  $n++;
  my $pseq = $p{$id};
  repeat: foreach my $t (keys %trip) {
    my (@sp) = split/$t/, $pseq;

    next unless ($#sp >= 2);
    for (my $i = 1; $i <= $#sp - 1; $i++) {
      my $l = length $sp[$i];
      next repeat if $l < $minSP || $l > $maxSP;
    }
    print FILE "$id\t";
    print FILE join ",", @{$trip{$t}};
    print FILE "\t$t\t$#sp\t$pseq\n";
    last if $ONE == 1;
  }
  print "$n proteins analysed...\n" if $n =~ /000$/;
}
close FILE;
  
print "$n proteins analysed...\n";
print "$OUTFILE file created";

exit;

##############
# SUBRUTINES #
##############
sub translate () {
  my ($s) = @_;

  # Código genético
  my (%cg) = (
    'TCA' => 'S', # Serina
    'TCC' => 'S', # Serina
    'TCG' => 'S', # Serina
    'TCT' => 'S', # Serina
    'TCN' => 'S', # Serina
    'TTC' => 'F', # Fenilalanina
    'TTT' => 'F', # Fenilalanina
    'TTA' => 'L', # Leucina
    'TTG' => 'L', # Leucina
    'TAC' => 'Y', # Tirosina
    'TAT' => 'Y', # Tirosina
    'TAA' => 'X', # Stop
    'TAG' => 'X', # Stop
    'TGC' => 'C', # Cisteina
    'TGT' => 'C', # Cisteina
    'TGA' => 'X', # Stop
    'TGG' => 'W', # Triptofano
    'CTA' => 'L', # Leucina
    'CTC' => 'L', # Leucina
    'CTG' => 'L', # Leucina
    'CTT' => 'L', # Leucina
    'CTN' => 'L', # Leucina
    'CCA' => 'P', # Prolina
    'CCC' => 'P', # Prolina
    'CCG' => 'P', # Prolina
    'CCT' => 'P', # Prolina
    'CCN' => 'P', # Prolina
    'CAC' => 'H', # Histidina
    'CAT' => 'H', # Histidina
    'CAA' => 'Q', # Glutamina
    'CAG' => 'Q', # Glutamina
    'CGA' => 'R', # Arginina
    'CGC' => 'R', # Arginina
    'CGG' => 'R', # Arginina
    'CGT' => 'R', # Arginina
    'CGN' => 'R', # Arginina
    'ATA' => 'I', # Isoleucina
    'ATC' => 'I', # Isoleucina
    'ATT' => 'I', # Isoleucina
    'ATG' => 'M', # Methionina
    'ACA' => 'T', # Treonina
    'ACC' => 'T', # Treonina
    'ACG' => 'T', # Treonina
    'ACT' => 'T', # Treonina
    'ACN' => 'T', # Treonina
    'AAC' => 'N', # Asparagina
    'AAT' => 'N', # Asparagina
    'AAA' => 'K', # Lisina
    'AAG' => 'K', # Lisina
    'AGC' => 'S', # Serina
    'AGT' => 'S', # Serina
    'AGA' => 'R', # Arginina
    'AGG' => 'R', # Arginina
    'GTA' => 'V', # Valina
    'GTC' => 'V', # Valina
    'GTG' => 'V', # Valina
    'GTT' => 'V', # Valina
    'GTN' => 'V', # Valina
    'GCA' => 'A', # Alanina
    'GCC' => 'A', # Alanina
    'GCG' => 'A', # Alanina
    'GCT' => 'A', # Alanina
    'GCN' => 'A', # Alanina
    'GAC' => 'D', # Acido Aspartico
    'GAT' => 'D', # Acido Aspartico
    'GAA' => 'E', # Acido Glutamico
    'GAG' => 'E', # Acido Glutamico
    'GGA' => 'G', # Glicina
    'GGC' => 'G', # Glicina
    'GGG' => 'G', # Glicina
    'GGT' => 'G', # Glicina
    'GGN' => 'G', # Glicina
  );

  my $aa;
  for (my $x = 0; $x <= length($s)-1; $x += 3) {
    my $tri = substr($s, $x, 3);
    
    my $tri_len = length($tri);
    last if ($tri_len < 3);
    if (!$cg{$tri}) {
      if ($tri =~ /N/) {
        $aa .= "X";
      } else {
        if ($tri =~ /[^ACGT]/) {
          die "Error! with codon $tri\n";
        }
      }
    } else { # if exist
      $aa .= $cg{$tri};
    }
  }

  return $aa;
}

