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

my $id1 = 'WT';
my $id2 = 'mutant';

my $enzymes = restriction_enzymes();
my $sites   = restriction_sites($enzymes);
my $seqs    = get_sequences( $id1, $id2 );

sub restriction_enzymes {
    return {
        BamHI   => 'G/GATCC',
        EcoRI   => 'G/AATTC',
        EcoRV   => 'GAT/ATC',
        HindIII => 'A/AGCTT',
    };
}

sub restriction_sites {
    my $enzymes = shift;
    my %sites = map { s|/||; $_ => 1 } values $enzymes;
    return [ keys %sites ];

sub get_sequences {
    my ( $id1, $id2 ) = @_;
    return {
        $id1 => 'cttctctagaggctttctcgcgaataagcGAATTCacgtgaagcgtcgatagccttcgat',
        $id2 => 'cttctctagaggctttctcgcgaataagcGAtTTCacgtgaagcgtcgatagccttcgat',
    }
}
}
