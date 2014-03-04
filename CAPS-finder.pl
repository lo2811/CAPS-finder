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

my @snp_files = 'sample-files/polyDB.A05.nr';
# my @snp_files = @ARGV;

my $id1 = 'R500';
my $id2 = 'IMB211';

my $enzymes = restriction_enzymes();
my $sites   = restriction_sites($enzymes);
my $snps    = import_snps( \@snp_files, $id1, $id2 );
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
    my ( $snp_files, $id1, $id2 ) = @_;

    my %snps;
    for my $file (@$snp_files) {
        open my $snp_fh, "<", $file;
        <$snp_fh>;
        while (<$snp_fh>) {
            next if /(?:INS)|(?:del)/;
            my ( $chr, $pos, $ref, $alt, $alt_geno ) = split /\t/;
            my $ref_geno = $alt_geno eq $id2 ? $id1 : $id2;
            $snps{$chr}{$pos}{$ref_geno} = $ref;
            $snps{$chr}{$pos}{$alt_geno} = $alt;
        }
        close $snp_fh;
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
