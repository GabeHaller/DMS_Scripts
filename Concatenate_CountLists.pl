#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;
use Cwd 'abs_path';

# ---------------------------
# Usage:
#   perl make_matrix_from_DNA_XisStop.pl --dir <folder> --fasta <dna.fa> [--out matrix.tsv] [--frame 0|1|2]
# Notes:
#   - X represents STOP (no '*' used).
#   - 21 possible AAs per position: A C D E F G H I K L M N P Q R S T V W Y X
# ---------------------------

my ($dir, $fasta, $out, $frame) = ("", "", "AA_change_matrix.tsv", 0);
GetOptions(
    "dir=s"   => \$dir,
    "fasta=s" => \$fasta,
    "out=s"   => \$out,
    "frame=i" => \$frame,
) or die "Error in options\n";

die "Usage: $0 --dir <folder_with_CountList.txt> --fasta <dna.fa> [--out matrix.tsv] [--frame 0|1|2]\n"
  unless $dir && $fasta && $frame =~ /^[012]$/;

# 20 canonical AAs + X (STOP)
my @AA = qw(A C D E F G H I K L M N P Q R S T V W Y X);
my %AA_ok = map { $_ => 1 } @AA;

# ---------- Genetic code (standard code 1) ----------
# Map STOP codons to 'X' (not '*')
my %gencode = (
  'TTT'=>'F','TTC'=>'F','TTA'=>'L','TTG'=>'L',
  'CTT'=>'L','CTC'=>'L','CTA'=>'L','CTG'=>'L',
  'ATT'=>'I','ATC'=>'I','ATA'=>'I','ATG'=>'M',
  'GTT'=>'V','GTC'=>'V','GTA'=>'V','GTG'=>'V',

  'TCT'=>'S','TCC'=>'S','TCA'=>'S','TCG'=>'S',
  'CCT'=>'P','CCC'=>'P','CCA'=>'P','CCG'=>'P',
  'ACT'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T',
  'GCT'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A',

  'TAT'=>'Y','TAC'=>'Y','TAA'=>'X','TAG'=>'X',  # STOP => X
  'CAT'=>'H','CAC'=>'H','CAA'=>'Q','CAG'=>'Q',
  'AAT'=>'N','AAC'=>'N','AAA'=>'K','AAG'=>'K',
  'GAT'=>'D','GAC'=>'D','GAA'=>'E','GAG'=>'E',

  'TGT'=>'C','TGC'=>'C','TGA'=>'X','TGG'=>'W',  # STOP => X
  'CGT'=>'R','CGC'=>'R','CGA'=>'R','CGG'=>'R',
  'AGT'=>'S','AGC'=>'S','AGA'=>'R','AGG'=>'R',
  'GGT'=>'G','GGC'=>'G','GGA'=>'G','GGG'=>'G',
);

# ---------- Read DNA FASTA ----------
open my $FA, "<", $fasta;
my $dna = "";
while (my $line = <$FA>) {
    next if $line =~ /^>/;
    chomp $line;
    $line =~ s/\s+//g;
    $dna .= $line;
}
close $FA;

$dna = uc($dna);
$dna =~ tr/U/T/;
die "Error: empty DNA sequence in $fasta\n" unless length $dna;

# ---------- Translate DNA (forward strand, selected frame) ----------
sub translate_dna {
    my ($seq, $frm) = @_;
    my $prot = "";
    for (my $i = $frm; $i + 2 < length($seq); $i += 3) {
        my $codon = substr($seq, $i, 3);
        $codon =~ s/[^ACGT]/N/g;
        my $aa = exists $gencode{$codon} ? $gencode{$codon} : 'X';  # unknown -> X (treated as stop)
        $prot .= $aa;
    }
    return $prot;
}

my $prot = translate_dna($dna, $frame);
my @refAA = split //, $prot;
my @valid_pos = (1 .. scalar @refAA);  # include all positions (including X/stop)

# ---------- Find input files ----------
opendir my $DH, $dir;
my @files = grep { /\.CountList\.txt$/i && -f "$dir/$_" } readdir($DH);
closedir $DH;
die "No *CountList.txt files found in $dir\n" unless @files;

my @paths    = map { abs_path("$dir/$_") } @files;
my @colnames = @files;

# ---------- Parse each CountList into change => count ----------
my %data;
for my $fi (0 .. $#paths) {
    my $path = $paths[$fi];
    open my $IN, "<", $path;
    while (my $line = <$IN>) {
        next if $line =~ /^\s*$/;
        chomp $line;
        my ($label, $count) = split /\s+/, $line, 2;
        next unless defined $label && defined $count;

        # Accept labels like A228A or X228L (X denotes stop)
        if ($label =~ /^([A-Z])(\d+)([A-Z])$/) {
            my ($refAA_l, $pos, $altAA) = (uc($1), $2, uc($3));
            next unless $AA_ok{$refAA_l} && $AA_ok{$altAA};
            if ($pos >= 1 && $pos <= @refAA) {
                # Optional: ensure label's ref matches translated ref; comment out to ignore mismatches
                # next unless $refAA[$pos-1] eq $refAA_l;
                $data{$fi}{"$refAA_l$pos$altAA"} = $count + 0;
            }
        }
    }
    close $IN;
}

# ---------- Generate all possible changes (21 per position, including X) ----------
my @all_changes;
for my $pos (@valid_pos) {
    my $r = $refAA[$pos-1];
    $r = 'X' unless $AA_ok{$r};  # coerce any unexpected symbol to X (stop)
    for my $alt (@AA) {
        push @all_changes, ($r.$pos.$alt);
    }
}

# ---------- Write output matrix ----------
open my $OUT, ">", $out;
print $OUT join("\t", "AA_change", @colnames), "\n";
for my $chg (@all_changes) {
    my @row = ($chg);
    for my $fi (0 .. $#paths) {
        my $val = exists $data{$fi}{$chg} ? $data{$fi}{$chg} : 0;
        push @row, $val;
    }
    print $OUT join("\t", @row), "\n";
}
close $OUT;

print "Wrote matrix to $out (frame=$frame; X=STOP). Protein length = ", scalar(@refAA), " AA positions.\n";