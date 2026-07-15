#!/usr/bin/perl
use strict;
use warnings;

my ($fasta_file, $sam_file) = @ARGV;
die "Usage: $0 <ref.fasta> <alignments.sam>\n" unless @ARGV == 2;

(my $prefix = $sam_file) =~ s/\.[^\.]+$//;
my $table_file       = $prefix . ".CountTable.txt";
my $list_file        = $prefix . ".CountList.txt";
my $reads_file       = $prefix . ".ReadVariants.txt";
my $codon_list_file  = $prefix . ".CountListByCodon.txt";   # NEW: AA change broken down by codon

# Codon table (X = stop)
my %codon_table = (
    TTT=>'F', TTC=>'F', TTA=>'L', TTG=>'L', CTT=>'L', CTC=>'L', CTA=>'L', CTG=>'L',
    ATT=>'I', ATC=>'I', ATA=>'I', ATG=>'M', GTT=>'V', GTC=>'V', GTA=>'V', GTG=>'V',
    TCT=>'S', TCC=>'S', TCA=>'S', TCG=>'S', CCT=>'P', CCC=>'P', CCA=>'P', CCG=>'P',
    ACT=>'T', ACC=>'T', ACA=>'T', ACG=>'T', GCT=>'A', GCC=>'A', GCA=>'A', GCG=>'A',
    TAT=>'Y', TAC=>'Y', TAA=>'X', TAG=>'X', CAT=>'H', CAC=>'H', CAA=>'Q', CAG=>'Q',
    AAT=>'N', AAC=>'N', AAA=>'K', AAG=>'K', GAT=>'D', GAC=>'D', GAA=>'E', GAG=>'E',
    TGT=>'C', TGC=>'C', TGA=>'X', TGG=>'W', CGT=>'R', CGC=>'R', CGA=>'R', CGG=>'R',
    AGT=>'S', AGC=>'S', AGA=>'R', AGG=>'R', GGT=>'G', GGC=>'G', GGA=>'G', GGG=>'G'
);

# Read reference cDNA
open my $fa, '<', $fasta_file or die $!;
my $refseq = '';
while (<$fa>) {
    chomp;
    next if /^>/;
    $refseq .= uc($_);
}
close $fa;

# Translate reference
sub translate {
    my ($cdna) = @_;
    my $protein = '';
    for (my $i = 0; $i < length($cdna) - 2; $i += 3) {
        my $codon = substr($cdna, $i, 3);
        $protein .= $codon_table{$codon} // 'X';
    }
    return $protein;
}
my $ref_protein = translate($refseq);
my $prot_len = length($ref_protein);

# Initialize counts
my (%counts, %flat_counts, %codon_counts, %read_variants);
for my $aa (values %codon_table) {
    for my $pos (1 .. $prot_len) {
        $counts{$aa}{$pos} = 0;
    }
}

my ($no_mutation_count, $single_mutation_count,$multi_mutation_count) = (0, 0, 0);

# Parse CIGAR string to get aligned portion excluding soft/hard clips
sub get_aligned_bases {
    my ($cigar, $seq) = @_;
    my $result = '';
    my $pos = 0;

    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($len, $op) = ($1, $2);
        if ($op eq 'M' || $op eq '=' || $op eq 'X') {
            $result .= substr($seq, $pos, $len);
            $pos += $len;
        } elsif ($op eq 'I' || $op eq 'S') {
            $pos += $len;
        } elsif ($op eq 'D' || $op eq 'N' || $op eq 'H' || $op eq 'P') {
            # Do nothing to read sequence (consumes ref only or padding)
        }
    }
    return $result;
}

# Process SAM
open my $sam, '<', $sam_file or die $!;
while (<$sam>) {
    next if /^@/;
    my @fields = split /\t/;
    my ($qname, $pos, $cigar, $raw_seq) = @fields[0,3,5,9];
    $pos--;

    my $aligned_seq = get_aligned_bases($cigar, uc($raw_seq));
    next if length($aligned_seq) <= 12;

    # Trim first and last 6 bases
    $aligned_seq = substr($aligned_seq, 6, length($aligned_seq) - 12);
    $pos += 6;

    next if $pos < 0 or $pos + length($aligned_seq) > length($refseq);

    my $modseq = $refseq;
    substr($modseq, $pos, length($aligned_seq)) = $aligned_seq;

    my $start_codon = int($pos / 3);
    my $end_codon   = int(($pos + length($aligned_seq) - 1) / 3);

    # Track AA changes with codon detail
    my %aa_changes; # aa_pos => [refAA, altAA, refCodon, altCodon]
    for my $i ($start_codon .. $end_codon) {
        my $start = $i * 3;
        next if $start + 3 > length($refseq);
        my $ref_codon  = substr($refseq, $start, 3);
        my $read_codon = substr($modseq, $start, 3);
        next unless length($ref_codon) == 3 && length($read_codon) == 3;

        my $ref_aa = $codon_table{$ref_codon} // 'X';
        my $alt_aa = $codon_table{$read_codon} // 'X';
        my $aa_pos = $i + 1;

        if ($ref_codon ne $read_codon) {
            $aa_changes{$aa_pos} = [$ref_aa, $alt_aa, $ref_codon, $read_codon];
        }
    }

    my $change_count = scalar keys %aa_changes;

    if ($change_count == 0) {
        $no_mutation_count++;
    } elsif ($change_count == 1) {
        my ($apos) = keys %aa_changes;
	$single_mutation_count++;
        my ($ref_aa, $alt_aa, $ref_codon, $alt_codon) = @{ $aa_changes{$apos} };

        # Increment existing per-AA matrix
        $counts{$alt_aa}{$apos}++ if exists $counts{$alt_aa}{$apos};

        # Flat AA-only change (A228T)
        my $aa_change = "$ref_aa$apos$alt_aa";
        $flat_counts{$aa_change}++;

        # NEW: AA broken down by codon (e.g., A228T (GCA>ACA))
        my $aa_codon_change = "$aa_change ($ref_codon>$alt_codon)";
        $codon_counts{$aa_codon_change}++;

        # Read-level annotation with codon detail
        $read_variants{$qname} = "$aa_change($ref_codon>$alt_codon)";
    } else {
        $multi_mutation_count++;
        my @ann;
        foreach my $apos (sort { $a <=> $b } keys %aa_changes) {
            my ($ref_aa, $alt_aa, $ref_codon, $alt_codon) = @{ $aa_changes{$apos} };
            push @ann, "$ref_aa$apos$alt_aa($ref_codon>$alt_codon)";
        }
        $read_variants{$qname} = join(",", @ann);
    }
}
close $sam;

# Output .CountTable.txt (same as before)
open my $out1, '>', $table_file or die $!;
print $out1 join("\t", "AltAA", map { substr($ref_protein, $_ - 1, 1) . $_ } 1 .. $prot_len), "\n";
foreach my $alt_aa (sort keys %counts) {
    print $out1 $alt_aa;
    for my $pos (1 .. $prot_len) {
        print $out1 "\t", $counts{$alt_aa}{$pos} // 0;
    }
    print $out1 "\n";
}
close $out1;

# Output .CountList.txt (AA-only)
open my $out2, '>', $list_file or die $!;
foreach my $mut (sort keys %flat_counts) {
    print $out2 "$mut\t$flat_counts{$mut}\n";
}
close $out2;

# NEW: Output .CountListByCodon.txt (AA broken down by codon)
open my $out4, '>', $codon_list_file or die $!;
foreach my $mut (sort keys %codon_counts) {
    print $out4 "$mut\t$codon_counts{$mut}\n";
}
close $out4;

# Output .ReadVariants.txt (now includes codon detail per change)
open my $out3, '>', $reads_file or die $!;
foreach my $read (sort keys %read_variants) {
    print $out3 "$read\t$read_variants{$read}\n";
}
close $out3;

# Summary
print "Total reads with NO amino acid changes: $no_mutation_count\n";
print "Total reads with exactly 1 amino acid position change (included in .ReadVariants): $single_mutation_count\n";
print "Total reads with >1 amino acid position change (included in .ReadVariants): $multi_mutation_count\n";
print "Matrix written to: $table_file\n";
print "List written to:   $list_file\n";
print "By-codon list to:  $codon_list_file\n";
print "Read variants to:  $reads_file\n";
