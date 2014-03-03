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

my $enzymes = restriction_enzymes();
my $sites   = restriction_sites($enzymes);

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
}
