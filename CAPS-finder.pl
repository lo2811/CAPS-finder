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

sub get_sequences {
    my ( $id1, $id2 ) = @_;
    return {
        $id1 => 'cttctctagaggctttctcgcgaataagcGAATTCacgtgaagcgtcgatagccttcgat',
        $id2 => 'cttctctagaggctttctcgcgaataagcGAtTTCacgtgaagcgtcgatagccttcgat',
    }
}
}
