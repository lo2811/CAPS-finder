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

my $fa = 'sample-files/B.rapa_genome_sequence_0830.fa';
my $region = 'A05:8000001-9000000';
# my $region = '';

my $id1 = 'R500';
my $id2 = 'IMB211';

my $enzymes = restriction_enzymes();
my $sites   = restriction_sites($enzymes);
my $snps    = import_snps( \@snp_files, $id1, $id2, $region );
find_caps_markers( $snps, $sites, $id1, $id2, $fa );

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
    my ( $snp_files, $id1, $id2, $region ) = @_;
    my ( $roi_chr, $roi_start, $roi_end )
        = $region =~ /^([^:]+):?(\d*)-?(\d*)/;

    my %snps;
    for my $file (@$snp_files) {
        open my $snp_fh, "<", $file;
        <$snp_fh>;
        while (<$snp_fh>) {
            next if /(?:INS)|(?:del)/;
            my ( $chr, $pos, $ref, $alt, $alt_geno ) = split /\t/;
            next if defined $roi_chr && $chr ne $roi_chr;
            next if defined $roi_start && $pos < $roi_start;
            next if defined $roi_end   && $pos > $roi_end;

            my $ref_geno = $alt_geno eq $id2 ? $id1 : $id2;
            $snps{$chr}{$pos}{$ref_geno} = $ref;
            $snps{$chr}{$pos}{$alt_geno} = $alt;
        }
        close $snp_fh;
    }

    return \%snps;
}

sub find_caps_markers {
    my ( $snps, $sites, $id1, $id2, $fa ) = @_;
    for my $chr ( sort keys $snps ) {
        my $chr_seq = get_chr_seq( $fa, $chr );
        for my $pos ( sort { $a <=> $b } keys $$snps{$chr} ) {
            my $seqs = get_sequences( \$chr_seq, $chr, $pos, $snps, $id1, $id2 );
            my $matches = marker_enzymes( $sites, $seqs );
            say join "\t", $chr, $pos, join ",", @$matches if @$matches;
        }
    }
}

sub get_chr_seq {
    my ( $fa, $chr ) = @_;

    my $cmd = "samtools faidx $fa $chr";
    my ( $fa_id, @seq ) = `$cmd`;
    chomp @seq;

    return join "", @seq;
}

sub get_sequences {
    my ( $chr_seq, $chr, $pos, $snps, $id1, $id2 ) = @_;

    my $flank  = 10;
    my $offset = $pos - ( $flank + 1 );
    my $length = 2 * $flank + 1;

    my $seq = substr $$chr_seq, $offset, $length;
    my $pre_snp = substr $seq, 0, $flank;
    my $post_snp = substr $seq, -$flank;

    return {
        $id1 => $pre_snp . $$snps{$chr}{$pos}{$id1} . $post_snp,
        $id2 => $pre_snp . $$snps{$chr}{$pos}{$id2} . $post_snp,
    };
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
