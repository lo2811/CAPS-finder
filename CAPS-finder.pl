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
use Getopt::Long;

my $current_version = '0.0.0';

my ( $id1, $id2, $fa, $region, $outdir ) = cli_options($current_version);
my @snp_files = @ARGV;
my $enzymes = restriction_enzymes();
my $sites   = restriction_sites($enzymes);
my $snps    = import_snps( \@snp_files, $id1, $id2, $region );
my $caps    = find_caps_markers( $snps, $sites, $id1, $id2, $fa );
output_caps_markers( $caps, $outdir, $id1, $id2, $region );

sub cli_options {
    my $current_version = shift;
    my ( $id1, $id2, $fa, $help, $version );
    my $region  = '';
    my $outdir  = '.';
    my $options = GetOptions(
        "id1=s"    => \$id1,
        "id2=s"    => \$id2,
        "fa=s"     => \$fa,
        "region=s" => \$region,
        "outdir=s" => \$outdir,
        "help"     => \$help,
        "version"  => \$version,
    );

    die "$0 $current_version\n" if $version;
    usage() if $help;
    usage() unless defined $id1 && defined $id2 && defined $fa;

    return $id1, $id2, $fa, $region, $outdir;
}

sub usage {
    die <<EOF;

USAGE
  $0 --id1 GENOTYPE --id2 GENOTYPE --fa FASTA [--region CHR:POS-POS] SNPFILE(S)

DESCRIPTION
  Finds potential sites for CAPS markers.

OPTIONS
  -h, --help                  Print this help message
  -v, --version               Print version number
  --id1          GENOTYPE     ID for genotype 1
  --id1          GENOTYPE     ID for genotype 2
  -f, --fa       FASTA        Path to FASTA reference file
  -r, --region   CHR:POS-POS  Chromosome and/or region of interest
  -o, --outdir   DIR          Directory to save output file
                               (Default is current directory)

EOF
}

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

    my %caps;
    for my $chr ( sort keys $snps ) {
        my $chr_seq = get_chr_seq( $fa, $chr );
        for my $pos ( sort { $a <=> $b } keys $$snps{$chr} ) {
            my $seqs
                = get_sequences( \$chr_seq, $chr, $pos, $snps, $id1, $id2 );
            my $matches = marker_enzymes( $sites, $seqs );
            if (@$matches) {
                $caps{$chr}{$pos}{enzymes} = join ",", @$matches;
                $caps{$chr}{$pos}{seqs} = $seqs;
            }
        }
    }

    return \%caps;
}

sub get_chr_seq {
    my ( $fa, $chr ) = @_;

    my $cmd = "samtools faidx $fa $chr";
    my ( $fa_id, @seq ) = `$cmd`;
    chomp @seq;

    die "Problem getting sequences for $chr.\n" unless @seq;
    return join "", @seq;
}

sub get_sequences {
    my ( $chr_seq, $chr, $pos, $snps, $id1, $id2 ) = @_;

    my $flank  = 10;
    my $offset = $pos - ( $flank + 1 );
    my $length = 2 * $flank + 1;

    my $seq = substr $$chr_seq, $offset, $length;
    $seq =~ tr/A-Z/a-z/;
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

sub output_caps_markers {
    my ( $caps, $outdir, $id1, $id2, $region ) = @_;

    $region =~ s/:/_/;
    my $output = "$outdir/CAPS.$id1-$id2";
    $output .= ".$region" if $region;

    open my $caps_fh, ">", $output;
    say $caps_fh join "\t", 'chr', 'pos', 'enzymes', $id1, $id2;
    for my $chr ( sort keys $caps ) {
        for my $pos ( sort { $a <=> $b } keys $$caps{$chr} ) {
            my $enzymes = $$caps{$chr}{$pos}{enzymes};
            my $seq1 = $$caps{$chr}{$pos}{seqs}{$id1};
            my $seq2 = $$caps{$chr}{$pos}{seqs}{$id2};
            say $caps_fh join "\t", $chr, $pos, $enzymes, $seq1, $seq2;
        }
    }
    close $caps_fh;
}
