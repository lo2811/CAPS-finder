#!/usr/bin/env perl
# Mike Covington
# created: 2014-03-12
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Data::Printer;

my $caps_file = "CAPS.sample";

my $flank = 100;

my %caps;
open my $caps_fh, "<", $caps_file;
my ( $id1, $id2 );
while (<$caps_fh>) {
    if ($. == 1) {
        ( $id1, $id2 ) = (split)[ 3 .. 4 ];
        next;
    }
    my ( $chr, $pos, $enzyme_list, $seq1, $seq2 ) = split;
    my @enzymes = split /,/, $enzyme_list;
    $caps{$chr}{$pos} = {
        'enzymes' => \@enzymes,
        'seq1'    => $seq1,
        'seq2'    => $seq2,
    };
}
close $caps_fh;

my $primer3_path = "/Users/mfc/Downloads/release-2.3.6";
my $size_range = "90-120 120-150 150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000";
my %pcr;
for my $chr ( sort keys %caps ) {
    for my $pos ( sort { $a <=> $b } keys $caps{$chr} ) {
        my $seq1 = $caps{$chr}{$pos}{seq1};
        my $seq2 = $caps{$chr}{$pos}{seq2};

        my @snp_positions = get_snp_positions( $seq1 );
        my $excluded = join " ",
            map {"$_,1"} grep { $_ != $flank } @snp_positions;
        $caps{$chr}{$pos}{excluded} = $excluded;

        my $primer3_in = <<EOF;
SEQUENCE_ID=${chr}_$pos
SEQUENCE_TEMPLATE=$seq1
PRIMER_TASK=generic
SEQUENCE_TARGET=$flank,1
SEQUENCE_EXCLUDED_REGION=$excluded
PRIMER_PRODUCT_SIZE_RANGE=$size_range
PRIMER_NUM_RETURN=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$primer3_path/primer3_config/
=
EOF
        chomp $primer3_in;

        my $results = `echo '$primer3_in' | $primer3_path/primer3_core`;

        my ($lt_primer)    = $results =~ /PRIMER_LEFT_0_SEQUENCE=([a-z]+)/;
        my ($rt_primer)    = $results =~ /PRIMER_RIGHT_0_SEQUENCE=([a-z]+)/;
        my ($lt_primer_tm) = $results =~ /PRIMER_LEFT_0_TM=(\d+\.\d+)/;
        my ($rt_primer_tm) = $results =~ /PRIMER_RIGHT_0_TM=(\d+\.\d+)/;
        my ($pcr_size)     = $results =~ /PRIMER_PAIR_0_PRODUCT_SIZE=(\d+)/;
        my ( $lt_pos, $lt_len ) = $results =~ /PRIMER_LEFT_0=(\d+),(\d+)/;
        my ( $rt_pos, $rt_len ) = $results =~ /PRIMER_RIGHT_0=(\d+),(\d+)/;

        my $amplicon1 = substr $seq1, $lt_pos, $pcr_size;
        my $amplicon2 = substr $seq2, $lt_pos, $pcr_size;

        $pcr{$chr}{$pos} = {
            'lt_primer'    => $lt_primer,
            'rt_primer'    => $rt_primer,
            'lt_primer_tm' => $lt_primer_tm,
            'rt_primer_tm' => $rt_primer_tm,
            'pcr_size'     => $pcr_size,
            'amplicon'     => {
                $id1 => $amplicon1,
                $id2 => $amplicon2,
                }
            }
    }
}

for my $chr ( sort keys %pcr ) {
    for my $pos ( sort { $a <=> $b } keys $pcr{$chr} ) {
        my $amplicon1 = $pcr{$chr}{$pos}{amplicon}{$id1};
        my $amplicon2 = $pcr{$chr}{$pos}{amplicon}{$id2};

        my $site = "ttt/aaa";

        $pcr{$chr}{$pos}{digest_lengths}{$id1}
            = get_fragment_lengths( $amplicon1, $site );
        $pcr{$chr}{$pos}{digest_lengths}{$id2}
            = get_fragment_lengths( $amplicon2, $site );
    }
}

p %pcr;
# p %caps;

sub get_snp_positions {
    my $seq = shift;
    my @positions;
    while ($seq =~ /[ACGT]/g) {
        push @positions, $-[0];
    }
    return @positions;
}

sub get_fragment_lengths {
    my ( $seq, $site ) = @_;

    my ( $lt_site, $rt_site ) = split /\//, $site;

    my @fragments = sort { $b <=> $a }
        map { length $_ } split /(?<=$lt_site)(?=$rt_site)/i, $seq;

    return \@fragments;
}
