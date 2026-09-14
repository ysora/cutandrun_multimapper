#!/usr/bin/perl
use strict;
use warnings;
use Bio::Cigar;

if (defined $ENV{SEED}) { srand($ENV{SEED}); }

my $MAX_ALTS = 5000;

while (my $line = <STDIN>) {
    chomp $line;

    next if $line =~ /^@/;

    my @f = split /\t/, $line;
    next unless @f >= 11; 

    my ($qname, $flag, $chr, $pos, $mapq, $cigar) = @f[0..5];

    next if ($flag & 0x4);
    next if ($flag & 0x100);

    my ($x0, $xa) = (undef, undef);
    for my $tag (@f[11..$#f]) {
        $x0 = $1 if $tag =~ /^X0:i:(\d+)/;
        $xa = $1 if $tag =~ /^XA:Z:(.*)/;
    }
    $x0 ||= 1; 

    my $start  = $pos - 1;                     
    my $strand = ($flag & 16) ? '-' : '+';    
    my $len    = Bio::Cigar->new($cigar)->reference_length;
    my $end    = $start + $len;

    # list of candidates ([chr,start,end,qname,weight,strand])
    my @cands = ( [ $chr, $start, $end, $qname, 1, $strand ] );

    my $n_alt = ($x0 > 1) ? ($x0 - 1) : 0;
    $n_alt = $MAX_ALTS if $n_alt > $MAX_ALTS;

    if (defined $xa && $n_alt > 0) {
        my @xa_items = grep { $_ ne '' } split /;/, $xa;

        $n_alt = @xa_items < $n_alt ? scalar(@xa_items) : $n_alt;

        for my $i (0 .. $n_alt - 1) {
            my $xa_item = $xa_items[$i];
            my ($alt_chr, $alt_pos_strand, $alt_cigar) = (split /,/, $xa_item)[0..2];
            next unless defined $alt_chr && defined $alt_pos_strand && defined $alt_cigar;

            my $alt_strand    = substr($alt_pos_strand, 0, 1);   
            my $alt_start_1b  = substr($alt_pos_strand, 1);      
            next unless $alt_start_1b =~ /^\d+$/;
            my $alt_start     = $alt_start_1b - 1;               

            my $alt_len       = Bio::Cigar->new($alt_cigar)->reference_length;
            my $alt_end       = $alt_start + $alt_len;

            push @cands, [ $alt_chr, $alt_start, $alt_end, $qname, 1, $alt_strand ];
        }
    }

    my $chosen = $cands[ int(rand(@cands)) ];
    print join("\t", @$chosen), "\n";
}
