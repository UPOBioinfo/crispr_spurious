# Spurious proteins from CRISPR sequences

Scripts used in the article "CRISPR sequences contaminate public databases with spurious proteins containing spaced repeats"

## Introduction
These scripts can be used to search for spurious protein sequences originating from CRISPR sequences. We used two different approaches which are descripted below. 

## First approach (search for already-annotated repeats)
It consisted in searching for translations of repeat sequences from the **CRISPRCasdb** database separated by putative spacers.

`Script: crispr_spurious1_initial.pl proteins.fasta crisprcasdb.fasta`

You can use with another argument (1) to obtain one only line per peptide dataset.

## Second approach (search for peptide repeats)
It consisted in searching for amino acid repeats separated by putative spacers directly in the protein sequences of **UniProtKB** database.

`Script: crispr_spurious2_initial.pl proteins.fasta output.tsv`

## Search for cas genes (discovery of Putative False Proteins)
Finally, the initial candidates from the two approaches are mapped to their corresponding genomic sequences, and cas genes were searched within 15 kb around the candidate.

`Script: crispr_spurious_pfp.pl proteins.dat initial_candidates.tsv path`

You can choose thresholds for both cas domain coverage and identity. The path is a folder where the cas domain profiles are stored (please see [crispr repository GitHub](https://github.com/UPOBioinfo/crispr) by UPOBioinfo for more details).

When the number of sequences is higher than 500 the following scripts, which use the **NCBI API**, should be used before:
`Script: download_NCBI_faa_files.pl`

`Script: crispr_spurious_pfp_multifasta_gff3.pl proteins.fasta`

## Files from the article
Database|First approach (initial)|First approach (final)|Second approach (initial)|Second approach (final)
---|---|---|---|---
sprot_archaea|[sprot_archaea1.1](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/initial/uniprot_sprot_archaea.tsv)|[sprot_archaea1.2](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/pfp_searching/sprot_archaea.tsv)|[sprot_archaea2.1](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/initial/uniprot_sprot_archaea.tsv)|[sprot_archaea2.2](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/pfp_searching/sprot_archaea.tsv)
sprot_bacteria|[sprot_bacteria1.1](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/initial/uniprot_sprot_bacteria.tsv)|[sprot_bacteria1.2](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/pfp_searching/sprot_bacteria.tsv)|[sprot_bacteria2.1](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/initial/uniprot_sprot_bacteria.tsv)|[sprot_bacteria2.2](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/pfp_searching/sprot_bacteria.tsv)
trembl_archaea|[trembl_archaea1.1](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/initial/uniprot_trembl_archaea.tsv)|[trembl_archaea1.2](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/pfp_searching/trembl_archaea.tsv)|[trembl_archaea2.1](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/initial/uniprot_trembl_archaea.tsv)|[trembl_archaea2.2](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/pfp_searching/trembl_archaea.tsv)
trembl_bacteria|[trembl_bacteria1.1](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/initial/uniprot_trembl_bacteria.tsv)|[trembl_bacteria1.2](http://www.bioinfocabd.upo.es/crispr_spurious/first_approach/pfp_searching/trembl_bacteria.tsv)|[trembl_bacteria2.1](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/initial/uniprot_trembl_bacteria.tsv)|[trembl_bacteria2.2](http://www.bioinfocabd.upo.es/crispr_spurious/second_approach/pfp_searching/trembl_bacteria.tsv)
