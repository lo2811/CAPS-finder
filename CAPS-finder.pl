#!/usr/bin/env perl
# Mike Covington
# created: 2014-03-03
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Data::Printer;

my $id1 = 'R500';
my $id2 = 'IMB211';

my $enzymes = restriction_enzymes();
my $sites   = restriction_sites($enzymes);
my $snps    = import_snps( $id1, $id2 );
my $seqs    = get_sequences( $id1, $id2 );
my $matches = marker_enzymes( $sites, $seqs );

say $_ for @$matches;

sub restriction_enzymes {
    return {
        BamHI   => 'G/GATCC',
        DpnI    => 'GA/TC',
        DpnII   => '/GATC',
        EcoRI   => 'G/AATTC',
        EcoRV   => 'GAT/ATC',
        HindIII => 'A/AGCTT',
    };
}

sub restriction_sites {
    my $enzymes = shift;
    my %sites;
    for ( keys $enzymes ) {
        my $site = $$enzymes{$_};
        $site =~ s|/||;
        push @{ $sites{$site} }, $_;
    }
    return \%sites;
}

sub import_snps {
    my ( $id1, $id2 ) = @_;

    my %snps;
    <DATA>;
    while (<DATA>) {
        next if /(?:INS)|(?:del)/;
        my ( $chr, $pos, $ref, $alt, $alt_geno ) = split /\t/;
        my $ref_geno = $alt_geno eq $id2 ? $id1 : $id2;
        $snps{$chr}{$pos}{$ref_geno} = $ref;
        $snps{$chr}{$pos}{$alt_geno} = $alt;
    }

    return \%snps;
}

sub get_sequences {
    my ( $id1, $id2 ) = @_;
    return {
        $id1 => 'cttctctagaggctttctcgcgaataagcGAATTCacgtgaagcgtcgatagccttcgat',
        $id2 => 'cttctctagaggctttctcgcgaataagcGAtTTCacgtgaagcgtcgatagccttcgat',
    }
}

sub marker_enzymes {
    my ( $sites, $seqs ) = @_;

    my %diffs;
    for my $site ( keys $sites ) {
        my $count = 0;
        $count += $$seqs{$id1} =~ /$site/ig;
        $count -= $$seqs{$id2} =~ /$site/ig;
        $diffs{$site} = $count;
    }

    my @matching_sites = grep { $diffs{$_} > 0 } keys %diffs;
    my @matching_enzymes;
    push @matching_enzymes, @{$$sites{$_}} for @matching_sites;

    return \@matching_enzymes;
}

__DATA__
chr	pos	ref_base	snp_base	genotype	insert_position	SNP_CLASS
A01	116746	T	C	IMB211	NA	SNP
A01	669437	C	T	R500	NA	SNP
A01	782857	T	A	R500	NA	SNP
A01	783874	T	C	R500	NA	SNP
A01	1148697	G	C	IMB211	NA	SNP
A01	1172157	C	A	IMB211	NA	SNP
A01	1224938	C	G	R500	NA	SNP
A01	1309920	A	G	R500	NA	SNP
A01	1309921	C	T	R500	NA	SNP
A01	1512215	T	G	IMB211	NA	SNP
A01	1843646	A	G	R500	NA	SNP
A01	1843647	T	G	R500	NA	SNP
A01	1884572	G	A	R500	NA	SNP
A01	2409594	T	C	IMB211	NA	SNP
A01	2413313	A	G	R500	NA	SNP
A01	2414553	C	T	IMB211	NA	SNP
A01	2430867	T	C	R500	NA	SNP
A01	2437915	C	A	R500	NA	SNP
A01	2461616	A	T	IMB211	NA	SNP
A01	2461634	G	A	IMB211	NA	SNP
A01	2461636	C	T	IMB211	NA	SNP
A01	2461647	A	G	IMB211	NA	SNP
A01	2461958	G	A	IMB211	NA	SNP
A01	2461994	A	G	IMB211	NA	SNP
A01	2462063	T	G	IMB211	NA	SNP
A01	2462111	A	G	IMB211	NA	SNP
A01	2462168	G	A	IMB211	NA	SNP
A01	2462216	A	G	IMB211	NA	SNP
A01	2462222	A	G	IMB211	NA	SNP
A01	2462228	A	T	IMB211	NA	SNP
A01	2462264	C	T	IMB211	NA	SNP
A01	2462278	G	A	IMB211	NA	SNP
A01	2462294	C	T	IMB211	NA	SNP
A01	2462309	T	G	IMB211	NA	SNP
A01	2462402	G	A	IMB211	NA	SNP
A01	2622870	C	T	IMB211	NA	SNP
A01	2622894	INS	A	IMB211	1	SNP
A01	2622894	INS	C	IMB211	2	SNP
A01	2622894	INS	T	IMB211	3	SNP
A01	2622936	T	C	IMB211	NA	SNP
A01	3770936	C	T	R500	NA	SNP
A01	3771291	G	C	R500	NA	SNP
A01	3771578	G	del	IMB211	NA	SNP
A01	3771580	G	del	IMB211	NA	SNP
A01	3771598	C	T	R500	NA	SNP
A01	3771640	C	T	R500	NA	SNP
A01	3772879	C	G	R500	NA	SNP
A01	3772891	C	A	R500	NA	SNP
A01	3772975	G	A	IMB211	NA	SNP
A01	3772978	C	T	IMB211	NA	SNP
A01	3773071	G	T	R500	NA	SNP
