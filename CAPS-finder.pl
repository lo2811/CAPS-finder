#!/usr/bin/env perl
# Mike Covington
# created: 2014-03-03
#
# Description: Find sites of potential CAPS markers.
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use File::Path 'make_path';
use Log::Reproducible;
reproduce();

# TODO: Find non-palindromic sites
# TODO: Find sites for enzymes with ambiguous bases
# TODO: Replace flank length with value related to min fragment length diff resolvable on gel
# TODO: Output restriction sites (next to enzyme names or in summary of CAPS marker enzymes)
# TODO: [longer term]/[maybe] Incorporate primer design
# TODO: Deal with $primer3_path and primer3 parameters
# TODO: Validate primer3 version

my $current_version = '0.3.1';

my $primer3_path = "/Users/mfc/Downloads/release-2.3.6";

my ( $id1, $id2, $fa, $region, $outdir, $flank, $multi_cut, $verbose )
    = cli_options($current_version);
my @snp_files = @ARGV;
my $enzymes   = restriction_enzymes();
my $sites     = restriction_sites($enzymes);
my $snps      = import_snps( \@snp_files, $id1, $id2, $region, $verbose );
my $caps
    = find_caps_markers( $snps, $sites, $id1, $id2, $fa, $flank, $multi_cut,
    $verbose );
make_path($outdir);
my $primers
    = design_primers( $caps, $id1, $id2, $flank, $primer3_path, $verbose );
digest_amplicons( $caps, $primers, $enzymes, $id1, $id2, $verbose );
output_caps_markers( $caps, $outdir, $id1, $id2, $region );
output_primers( $primers, $outdir, $enzymes, $id1, $id2, $region, $verbose );
exit;

sub cli_options {
    my $current_version = shift;
    my ( $id1, $id2, $fa, $help, $version );
    my $region  = '';
    my $outdir  = '.';
    my $flank   = 100;
    my $options = GetOptions(
        "id1=s"     => \$id1,
        "id2=s"     => \$id2,
        "fa=s"      => \$fa,
        "region=s"  => \$region,
        "outdir=s"  => \$outdir,
        "flank=i"   => \$flank,
        "multi_cut" => \$multi_cut,
        "verbose"   => \$verbose,
        "help"      => \$help,
        "version"   => \$version,
    );

    die "$0 $current_version\n" if $version;
    usage() if $help;
    usage() unless defined $id1 && defined $id2 && defined $fa;

    return $id1, $id2, $fa, $region, $outdir, $flank, $multi_cut, $verbose;
}

sub usage {
    die <<EOF;

USAGE
  $0 --id1 GENOTYPE --id2 GENOTYPE --fa FASTA [--region CHR:POS-POS] SNPFILE(S)

DESCRIPTION
  Find sites of potential CAPS markers.

OPTIONS
  -h, --help                   Print this help message
  --version                    Print version number
  --id1           GENOTYPE     ID for genotype 1
  --id1           GENOTYPE     ID for genotype 2
  --fa            FASTA        Path to FASTA reference file
  -r, --region    CHR:POS-POS  Chromosome and/or region of interest
  -o, --outdir    [.]          Directory to save output file
  --flank         [10]         Length of SNP-flanking sequence to use
  -m, --multi_cut              Allow multiple cuts per amplicon
  --verbose                    Verbose mode

EOF
}

sub restriction_enzymes {
    my %enzymes;
    while (<DATA>) {
        chomp;
        my ( $name, $site ) = split /,/;
        $enzymes{$name} = $site;
    }
    return \%enzymes;
}

sub restriction_sites {
    my $enzymes = shift;
    my %sites;
    for ( sort keys $enzymes ) {
        my $site = $$enzymes{$_};
        $site =~ s|/||;
        push @{ $sites{$site} }, $_;
    }
    return \%sites;
}

sub import_snps {
    my ( $snp_files, $id1, $id2, $region, $verbose ) = @_;
    my ( $roi_chr, $roi_start, $roi_end )
        = $region =~ /^([^:]+):?(\d*)-?(\d*)/;

    say join " ", "Importing SNPs from", scalar @$snp_files, "files"
        if $verbose;

    my %snps;
    for my $file (@$snp_files) {
        open my $snp_fh, "<", $file;
        <$snp_fh>;
        while (<$snp_fh>) {
            next if /(?:INS)|(?:del)/;
            my ( $chr, $pos, $ref, $alt, $alt_geno ) = split /\t/;
            next if defined $roi_chr && $chr ne $roi_chr;
            next if $roi_start       && $pos < $roi_start;
            next if $roi_end         && $pos > $roi_end;

            my $ref_geno = $alt_geno eq $id2 ? $id1 : $id2;
            $snps{$chr}{$pos}{$ref_geno} = $ref;
            $snps{$chr}{$pos}{$alt_geno} = $alt;
        }
        close $snp_fh;
    }

    return \%snps;
}

sub find_caps_markers {
    my ( $snps, $sites, $id1, $id2, $fa, $flank, $multi_cut, $verbose ) = @_;

    say "Finding CAPS markers on chromosome:" if $verbose;
    my %caps;
    for my $chr ( sort keys $snps ) {
        say "  $chr" if $verbose;
        my ( $chr_seq1, $chr_seq2 )
            = get_chr_seq( $fa, $chr, $snps, $id1, $id2 );
        for my $pos ( sort { $a <=> $b } keys $$snps{$chr} ) {
            my $seqs = get_sequences( \$chr_seq1, \$chr_seq2, $pos, $flank );
            my $matches = marker_enzymes( $sites, $seqs, $flank, $multi_cut );
            if (@$matches) {
                $caps{$chr}{$pos}{enzymes} = $matches;
                $caps{$chr}{$pos}{seqs} = $seqs;
            }
        }
    }

    return \%caps;
}

sub get_chr_seq {
    my ( $fa, $chr, $snps, $id1, $id2 ) = @_;

    my $cmd = "samtools faidx $fa $chr";
    my ( $fa_id, @seq ) = `$cmd`;
    chomp @seq;

    die "Problem getting sequences for $chr.\n" unless @seq;
    my $seq1 = join "", @seq;
    $seq1 =~ tr/A-Z/a-z/;
    my $seq2 = $seq1;

    for my $pos ( sort { $a <=> $b } keys $$snps{$chr} ) {
        substr $seq1, $pos - 1, 1, $$snps{$chr}{$pos}{$id1};
        substr $seq2, $pos - 1, 1, $$snps{$chr}{$pos}{$id2};
    }

    return $seq1, $seq2;
}

sub get_sequences {
    my ( $chr_seq1, $chr_seq2, $pos, $flank ) = @_;

    my $offset = $pos - ( $flank + 1 );
    my $length = 2 * $flank + 1;

    my $seq1 = substr $$chr_seq1, $offset, $length;
    my $seq2 = substr $$chr_seq2, $offset, $length;

    return {
        $id1 => $seq1,
        $id2 => $seq2,
    };
}

sub marker_enzymes {
    my ( $sites, $seqs, $flank, $multi_cut ) = @_;

    my %diffs;
    for my $site ( keys $sites ) {

        my $max = $flank;
        my $min = $max - ( 1 + length $site );
        next
            if !$multi_cut
            && $$seqs{$id1} =~ /$site/i
            && $$seqs{$id2} =~ /$site/i;
        my $count = 0;
        $count += $$seqs{$id1} =~ /^[ACGT]{$min,$max}$site[ACGT]{$min,$max}$/i;
        $count -= $$seqs{$id2} =~ /^[ACGT]{$min,$max}$site[ACGT]{$min,$max}$/i;
        $diffs{$site} = $count;
    }

    my @matching_sites = grep { $diffs{$_} > 0 } keys %diffs;
    my @matching_enzymes;
    push @matching_enzymes, @{$$sites{$_}} for @matching_sites;

    return \@matching_enzymes;
}

sub design_primers {
    my ( $caps, $id1, $id2, $flank, $primer3_path, $verbose ) = @_;

    say "Designing primers" if $verbose;
    my $primer3_parameters_path
        = write_primer3_parameters( $caps, $id1, $id2, $flank, $region,
        $outdir );
    my $primer3_out = run_primer3( $primer3_parameters_path, $primer3_path );
    my $primers = parse_primer3_results( $primer3_out, $caps );

    return $primers;
}

sub write_primer3_parameters {
    my ( $caps, $id1, $id2, $flank, $region, $outdir ) = @_;
    my $size_range
        = "90-120 120-150 150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000";

    my $primer3_parameters;
    for my $chr ( sort keys $caps ) {
        for my $pos ( sort { $a <=> $b } keys $$caps{$chr} ) {
            my $seq1 = $$caps{$chr}{$pos}{seqs}{$id1};
            my $seq2 = $$caps{$chr}{$pos}{seqs}{$id2};
            my $excluded = exclude_snps( $seq1, $flank );

            $primer3_parameters
                .= primer3_input_record( $chr, $pos, $seq1, $flank, $excluded,
                $size_range, $primer3_path );
        }
    }
    chomp $primer3_parameters;

    $region =~ s/:/_/;
    my $output = "$outdir/primer3_parameters.$id1-$id2";
    $output .= ".$region" if $region;

    open my $primer3_fh, ">", $output;
    print $primer3_fh $primer3_parameters;
    close $primer3_fh;

    return $output;
}

sub exclude_snps {
    my ( $seq1, $flank ) = @_;

    my @snp_positions;
    while ( $seq1 =~ /[ACGT]/g ) {    # Only SNPs are upper-case
        push @snp_positions, $-[0];
    }

    my $excluded = join " ",
        map {"$_,1"} grep { $_ != $flank } @snp_positions;
    return $excluded;
}

sub primer3_input_record {
    my ( $chr, $pos, $seq, $flank, $excluded, $size_range, $primer3_path )
        = @_;

    return <<EOF;
SEQUENCE_ID=${chr}_$pos
SEQUENCE_TEMPLATE=$seq
PRIMER_TASK=generic
SEQUENCE_TARGET=$flank,1
SEQUENCE_EXCLUDED_REGION=$excluded
PRIMER_PRODUCT_SIZE_RANGE=$size_range
PRIMER_NUM_RETURN=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$primer3_path/primer3_config/
=
EOF
}

sub run_primer3 {
    my ( $primer3_parameters_path, $primer3_path ) = @_;

    my $primer3_out = `$primer3_path/primer3_core $primer3_parameters_path`;
    return $primer3_out;
}

sub parse_primer3_results {
    my ( $primer3_out, $caps ) = @_;

    my $old_input_rec_sep = $/;
    $/ = '^=$';
    my @primer3_results = split /\n=\n/, $primer3_out;
    $/ = $old_input_rec_sep;

    my %primers;
    for my $result (@primer3_results) {

        my ($success) = $result =~ /PRIMER_PAIR_NUM_RETURNED=(\d+)/;
        next if $success == 0;

        my ($lt_primer)    = $result =~ /PRIMER_LEFT_0_SEQUENCE=([a-z]+)/;
        my ($rt_primer)    = $result =~ /PRIMER_RIGHT_0_SEQUENCE=([a-z]+)/;
        my ($lt_primer_tm) = $result =~ /PRIMER_LEFT_0_TM=(\d+\.\d+)/;
        my ($rt_primer_tm) = $result =~ /PRIMER_RIGHT_0_TM=(\d+\.\d+)/;
        my ($pcr_size)     = $result =~ /PRIMER_PAIR_0_PRODUCT_SIZE=(\d+)/;
        my ( $lt_pos, $lt_len ) = $result =~ /PRIMER_LEFT_0=(\d+),(\d+)/;
        my ( $rt_pos, $rt_len ) = $result =~ /PRIMER_RIGHT_0=(\d+),(\d+)/;

        my ( $chr, $pos ) = $result =~ /SEQUENCE_ID=(.+)_(\d+)/;
        my $seq1      = $$caps{$chr}{$pos}{seqs}{$id1};
        my $seq2      = $$caps{$chr}{$pos}{seqs}{$id2};
        my $amplicon1 = substr $seq1, $lt_pos, $pcr_size;
        my $amplicon2 = substr $seq2, $lt_pos, $pcr_size;

        $primers{$chr}{$pos} = {
            'lt_primer'    => $lt_primer,
            'rt_primer'    => $rt_primer,
            'lt_primer_tm' => $lt_primer_tm,
            'rt_primer_tm' => $rt_primer_tm,
            'pcr_size'     => $pcr_size,
            'amplicon'     => {
                $id1 => $amplicon1,
                $id2 => $amplicon2,
            }
        };

    }

    return \%primers;
}

sub digest_amplicons {
    my ( $caps, $primers, $enzymes, $id1, $id2, $verbose ) = @_;

    say "Performing virtual digests" if $verbose;
    for my $chr ( sort keys $primers ) {
        for my $pos ( sort { $a <=> $b } keys $$primers{$chr} ) {
            my $amplicon1 = $$primers{$chr}{$pos}{amplicon}{$id1};
            my $amplicon2 = $$primers{$chr}{$pos}{amplicon}{$id2};

            for my $enzyme ( @{ $$caps{$chr}{$pos}{enzymes} } ) {
                my $site = $$enzymes{$enzyme};

                $$primers{$chr}{$pos}{digest_lengths}{$enzyme}{$id1}
                    = get_fragment_lengths( $amplicon1, $site );
                $$primers{$chr}{$pos}{digest_lengths}{$enzyme}{$id2}
                    = get_fragment_lengths( $amplicon2, $site );
            }
        }
    }
}

sub get_fragment_lengths {
    my ( $seq, $site ) = @_;

    my ( $lt_site, $rt_site ) = split /\//, $site;
    my @fragments = sort { $b <=> $a }
        map { length $_ } split /(?<=$lt_site)(?=$rt_site)/i, $seq;

    return \@fragments;
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
            my $enzymes = join ",", @{ $$caps{$chr}{$pos}{enzymes} };
            my $seq1    = $$caps{$chr}{$pos}{seqs}{$id1};
            my $seq2    = $$caps{$chr}{$pos}{seqs}{$id2};
            say $caps_fh join "\t", $chr, $pos, $enzymes, $seq1, $seq2;
        }
    }
    close $caps_fh;
}

sub output_primers {
    my ( $primers, $outdir, $enzymes, $id1, $id2, $region ) = @_;

    $region =~ s/:/_/;
    my $output = "$outdir/primers.$id1-$id2";
    $output .= ".$region" if $region;

    say "Writing CAPS markers to file: $output";
    open my $primers_fh, ">", $output;
    say $primers_fh join "\t",
        'chr', 'pos', 'digest(s)', 'fwd_primer,rev_primer', 'fwd_Tm,rev_Tm',
        "${id1}_amplicon,${id2}_amplicon", 'length';

    for my $chr ( sort keys $primers ) {
        for my $pos ( sort { $a <=> $b } keys $$primers{$chr} ) {

            my $pcr_size    = $$primers{$chr}{$pos}{pcr_size};

            my $lt_primer      = $$primers{$chr}{$pos}{lt_primer};
            my $lt_primer_tm   = $$primers{$chr}{$pos}{lt_primer_tm};
            my $rt_primer      = $$primers{$chr}{$pos}{rt_primer};
            my $rt_primer_tm   = $$primers{$chr}{$pos}{rt_primer_tm};
            my $primer_info    = "$lt_primer,$rt_primer";
            my $primer_tm_info = "$lt_primer_tm,$rt_primer_tm";

            my $amplicon1 = $$primers{$chr}{$pos}{amplicon}{$id1};
            my $amplicon2 = $$primers{$chr}{$pos}{amplicon}{$id2};
            my $amplicons = "$amplicon1,$amplicon2";

            my @enzyme_info;
            for my $enzyme ( sort keys $$primers{$chr}{$pos}{digest_lengths} ) {
                my $site = $$enzymes{$enzyme};

                my $frags1 = join "+",
                    @{ $$primers{$chr}{$pos}{digest_lengths}{$enzyme}{$id1} };
                my $frags2 = join "+",
                    @{ $$primers{$chr}{$pos}{digest_lengths}{$enzyme}{$id2} };

                push @enzyme_info, "$enzyme($site):$frags1,$frags2";
            }
            my $enzyme_summary = join "|", @enzyme_info;

            say $primers_fh join "\t", $chr, $pos, $enzyme_summary,
                $primer_info, $primer_tm_info, $amplicons, $pcr_size;
        }
    }
    close $primers_fh;
}

__DATA__
AatII,gacgt/c
AbsI,cc/tcgagg
Acc65I,g/gtacc
AclI,aa/cgtt
AfeI,agc/gct
AflII,c/ttaag
AgeI,a/ccggt
AluI,ag/ct
AoxI,/ggcc
ApaI,gggcc/c
ApaLI,g/tgcac
AscI,gg/cgcgcc
AseI,at/taat
Asi256I,g/atc
AsiSI,gcgat/cgc
AvrII,c/ctagg
BamHI,g/gatcc
BclI,t/gatca
BfaI,c/tag
BglII,a/gatct
BmtI,gctag/c
BsiWI,c/gtacg
BspEI,t/ccgga
BspHI,t/catga
BsrGI,t/gtaca
BssHII,g/cgcgc
BstBI,tt/cgaa
BstKTI,gat/c
BstUI,cg/cg
BstZ17I,gta/tac
ChaI,gatc/
ClaI,at/cgat
CviAII,c/atg
CviQI,g/tac
DpnI,ga/tc
DraI,ttt/aaa
EagI,c/ggccg
EcoRI,g/aattc
EcoRV,gat/atc
Eco53kI,gag/ctc
EsaBC3I,tc/ga
FatI,/catg
FseI,ggccgg/cc
FspI,tgc/gca
GlaI,gc/gc
HaeIII,gg/cc
HhaI,gcg/c
HinP1I,g/cgc
HindIII,a/agctt
HpaI,gtt/aac
HpaII,c/cgg
HpyCH4IV,a/cgt
HpyCH4V,tg/ca
KasI,g/gcgcc
KpnI,ggtac/c
KroI,g/ccggc
MauBI,cg/cgcgcg
MboI,/gatc
McaTI,gcgc/gc
MfeI,c/aattg
MluI,a/cgcgt
MluCI,/aatt
MreI,cg/ccggcg
MscI,tgg/cca
MseI,t/taa
NaeI,gcc/ggc
NarI,gg/cgcc
NcoI,c/catgg
NdeI,ca/tatg
NgoMIV,g/ccggc
NheI,g/ctagc
NlaIII,catg/
NotI,gc/ggccgc
NruI,tcg/cga
NsiI,atgca/t
PabI,gta/c
PacI,ttaat/taa
PciI,a/catgt
PluTI,ggcgc/c
PmeI,gttt/aaac
PmlI,cac/gtg
Ppu10I,a/tgcat
PsiI,tta/taa
PspOMI,g/ggccc
PstI,ctgca/g
PvuI,cgat/cg
PvuII,cag/ctg
RsaI,gt/ac
SacI,gagct/c
SacII,ccgc/gg
SalI,g/tcgac
SbfI,cctgca/gg
ScaI,agt/act
SciI,ctc/gag
SelI,/cgcg
SfoI,ggc/gcc
SgrDI,cg/tcgacg
SmaI,ccc/ggg
SnaBI,tac/gta
SpeI,a/ctagt
SphI,gcatg/c
SrfI,gccc/gggc
SspI,aat/att
Sth302II,cc/gg
StuI,agg/cct
SwaI,attt/aaat
TaiI,acgt/
TaqI,t/cga
XbaI,t/ctaga
XhoI,c/tcgag
XmaI,c/ccggg
ZraI,gac/gtc
