#!/usr/bin/env perl
# Mike Covington
# created: 2014-03-04
#
# Description: Parse strider format restriction enzyme files for CAPS-finder.pl
#
use strict;
use warnings;
use autodie;
use feature 'say';

my $strider_file = $ARGV[0];    # e.g., striderc.403
open my $strider_in_fh,  "<", $strider_file;
open my $strider_out_fh, ">", "$strider_file.parsed";
while (<$strider_in_fh>) {
    next if $. <= 10;

    # Ignore asymmetric cutters for now
    next if /^#/;

    my ( $enzyme, $site ) = split /,/;

    # Ignore sites with ambiguous bases for now
    next if $site =~ /[^ACGT\/]/i;

    say $strider_out_fh join ",", $enzyme, $site;
}
close $strider_in_fh;
close $strider_out_fh;

exit;
