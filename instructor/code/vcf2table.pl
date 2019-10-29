#! /usr/bin/perl -w

use strict;
use POSIX qw(floor);
#use Data::Dumper;


my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;

my $usage = "\n";
$usage .=   "Usage: $progname GFF [options]\n";
$usage .=   "Converts the vcf (variant calling format) data on STDIN to an \n";
$usage .=   "annotated variant table and writes to STDOUT. GFF file is used as \n";
$usage .=   "the source of annotations. If the reference genome fasta file is \n";
$usage .=   "provided (-g), variants in coding regions are translated into amino \n";
$usage .=   "acid changes.\n";
$usage .=   "       [-g FNA] Fasta file with the reference genome sequence.\n";
$usage .=   "\n";


my %AA = (
	'AAA' => 'Lys', 'AAC' => 'Asn', 'AAG' => 'Lys', 'AAT' => 'Asn',
	'ACA' => 'Thr', 'ACC' => 'Thr', 'ACG' => 'Thr', 'ACT' => 'Thr',
	'AGA' => 'Arg', 'AGC' => 'Ser', 'AGG' => 'Arg', 'AGT' => 'Ser',
	'ATA' => 'Ile', 'ATC' => 'Ile', 'ATG' => 'Met', 'ATT' => 'Ile',

	'CAA' => 'Gln', 'CAC' => 'His', 'CAG' => 'Gln', 'CAT' => 'His',
	'CCA' => 'Pro', 'CCC' => 'Pro', 'CCG' => 'Pro', 'CCT' => 'Pro',
	'CGA' => 'Arg', 'CGC' => 'Arg', 'CGG' => 'Arg', 'CGT' => 'Arg',
	'CTA' => 'Leu', 'CTC' => 'Leu', 'CTG' => 'Leu', 'CTT' => 'Leu',

	'GAA' => 'Glu', 'GAC' => 'Asp', 'GAG' => 'Glu', 'GAT' => 'Asp',
	'GCA' => 'Ala', 'GCC' => 'Ala', 'GCG' => 'Ala', 'GCT' => 'Ala',
	'GGA' => 'Gly', 'GGC' => 'Gly', 'GGG' => 'Gly', 'GGT' => 'Gly',
	'GTA' => 'Val', 'GTC' => 'Val', 'GTG' => 'Val', 'GTT' => 'Val',

	'TAA' => 'AMB', 'TAC' => 'Tyr', 'TAG' => 'OCH', 'TAT' => 'Tyr',
	'TCA' => 'Ser', 'TCC' => 'Ser', 'TCG' => 'Ser', 'TCT' => 'Ser',
	'TGA' => 'OPL', 'TGC' => 'Cys', 'TGG' => 'Trp', 'TGT' => 'Cys',
	'TTA' => 'Leu', 'TTC' => 'Phe', 'TTG' => 'Leu', 'TTT' => 'Phe'
);

my ($gff,$fna);


while (@ARGV) {
  my $arg = shift;
  if ($arg eq '-h' or $arg eq '-help') {
		die "$usage";
	} elsif ($arg eq '-g' or $arg eq '-genome') {
		defined ($fna=shift) or die "$usage";
	} else {
		$gff=$arg;
	}
}
die "$usage" unless defined $gff and -f $gff;


# read in the reference genome file, if provided
my %genomes=();
my $cid;
if (defined $fna) {
	open my $fhg, "$fna" or die "$!";
	while (<$fhg>) {
		chomp;
		next if /^\s*$/;
		
		if (/^>/) {
			my @a = split /\s/;
			my @b = split /\|/, $a[0];
			$cid = pop @b;
			$cid = pop(@b) if $cid eq '';
			$cid =~ s/^>//;
			$genomes{$cid} = '';
		} else {
			$genomes{$cid} = $genomes{$cid} . "$_";
		}
	}
	close $fhg;
}
#print Dumper(\%genomes);
#die;


my %table=();

# read in the annotations
open my $fh, "$gff" or die "$!";
while (<$fh>) {
	chomp;
	next if /^#/;
	
	my ($chr, $src, $type, $start, $stop, $dot1, $strand, $dot2, $info) = split /\t/;
	next if lc $type eq 'gene' or lc $type eq 'source' or lc $type eq 'region' or lc $type eq 'exon';
	
	$table{$chr}={} unless defined $table{$chr};

	if ($stop<$start) {
		my $swap=$start;
		$start=$stop;
		$stop =$swap;
	}
	
	my ($gid, $pid, $prod)=('', '', '');
	
	$gid = $1 if $info =~ /note=(.+?);/i;
	$gid = "GI:$1" if $info =~ /db_xref=GI;:(.+?) /i;
	$gid = $1 if $info =~ /locus_tag=(.+?);/i;
	$gid = $1 if $info =~ /gene=(.+?);/i;
	$gid =~ s/"//g;
	
	$pid = $1 if $info =~ /protein_id=(.+?);/i;
	$pid = $gid if $pid eq '';

	$prod = $1 if $info =~ /product=(.+?)"?;/i;
	$prod = $pid if $prod eq '';
	$prod =~ s/"//g;
	
	
	$table{$chr}->{$start}={
		'type'   => $type,
		'start'  => $start,
		'stop'   => $stop,
		'strand' => $strand,
		'id'     => $gid,
		'prodid' => $pid,
		'prod'   => $prod,
		'info'   => $info
	};
	
}
close $fh;


# read in the variants
my $incr=1;
while (<>) {
	chomp;
	next if /^\s*$/ or /^#/;
	
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = split /\t/;
	if ($chr =~ /\|/) {
		my @a=split(/\|/, $chr);
		my $x=pop @a;
		$x = pop(@a) if $x eq '';
		$chr=$x;
	}
	
	if ($id eq '.') {
		$id='GV' . sprintf("%05d", $incr);
	}
	$incr++;
	my $depth='';
	$depth=$1 if $info =~ /DP=(.+?);/i;
	
	my $vtype='SNP';
	$vtype='INDEL' if $info =~ /INDEL;/i;
	
	$table{$chr}={} unless defined $table{$chr};
	$pos = $pos+0.1 if defined $table{$chr}->{$pos};
	
	$table{$chr}->{$pos}={
		'vid'     => $id,
		'type'    => 'variant',
		'vtype'   => $vtype,
		'start'   => $pos,
		'stop'    => $pos + length($ref)-1,
		'ref'     => $ref,
		'alt'     => $alt,
		'quality' => $qual,
		'depth'   => $depth,
		'info'    => $info
	};
	
}

# assemble and print the table
my @hdr=qw(
	variant_id
	chr
	feature_type
	feature_name
	feature_acc
	product
	feature_start_in_ref
	feature_stop_in_ref
	variant_position_in_ref
	variant_type
	variant_context
	variant_placement
	variation_nt
	variation_aa
	feature_sequence
	variant_sequence
	variant_quality
	variant_depth
	variant_info_full
	feature_info_full
);
print "#" . join("\t", @hdr) . "\n";

#print Dumper(\%table);
#die;

foreach my $chr (keys %table) {
	my @positions = sort {$a <=> $b} keys %{$table{$chr}};
	for (my $i=0; $i<scalar(@positions); $i++) {
		
		# find the variant feature and its start position
		my $pos=$positions[$i];
		next unless $table{$chr}->{$pos}->{'type'} eq 'variant';
		
		my $variant = $table{$chr}->{$pos};
		my $vstart = $variant->{'start'};
		$vstart =~ s/\.\d+//gi;
		
		#print STDERR "$vstart\n";
		
		# get the 5' (non-variant) feature, if one exists
		my $j = $i;
		my $feature5;
		while (!defined $feature5 and $j>0) {
			$j--;
			$feature5 = $table{$chr}->{$positions[$j]} unless $table{$chr}->{$positions[$j]}->{'type'} eq 'variant';
		}

		# get the 3' (non-variant) feature, if one exists
		$j=$i;
		my $feature3;
		while (!defined $feature3 and $j<scalar(@positions)-1) {
			$j++;
			$feature3 = $table{$chr}->{$positions[$j]} unless $table{$chr}->{$positions[$j]}->{'type'} eq 'variant';
		}
		
		
		# add characteristics of the feature that contains the variant
		# if the variant falls in an IGR, this is a combination of 5' and 3' features
		my %featureV = (
		'type'   => '',
		'id'     => '',
		'prodid' => '',
		'prod'   => '',
		'start'  => '',
		'stop'   => '',
		'info'   => ''
		);
		
		my $variant_nt  ='';
		my $variant_aa  ='';
		my $variant_loc ='';
		my $variant_plac='';
		
		# variant is in the 3' IGR of a defined 5' feature, or the 5' IGR of a defined 3' feature
		if (defined $feature5) {
			# 5' feature is defined

			$variant_nt   = "$variant->{ref}" . ($vstart - $feature5->{'stop'});
			$variant_nt  .= " -> $variant->{alt}" if $variant->{'vtype'} eq 'SNP';
			$variant_plac = $vstart - $feature5->{'stop'};

			if ($vstart > $feature5->{'stop'}) {
				# variant is downstream of the 5' feature stop
				# since it also occurs before the start of any 3' feature, it is in the IGR
				
				$variant_loc = 'IGR';
				$variant_plac .= ' down';
				
				# record the 5' feature details
				$featureV{'type'}   = "$feature5->{type} | ";
				$featureV{'id'}     = "$feature5->{id} | ";
				$featureV{'prodid'} = "$feature5->{prodid} | ";
				$featureV{'prod'}   = "$feature5->{prod} | ";
				$featureV{'start'}  = "$feature5->{start} | ";
				$featureV{'stop'}   = "$feature5->{stop} | ";
				$featureV{'info'}   = "$feature5->{info} | ";

				if (defined $feature3) {
					# variant is in between two defined features
					# we record details for both features
					$variant_plac .= ' | ' . ($feature3->{'start'} - $variant->{'stop'}) . ' up';
					$featureV{'type'}   .= "$feature3->{type}";
					$featureV{'id'}     .= "$feature3->{id}";
					$featureV{'prodid'} .= "$feature3->{prodid}";
					$featureV{'prod'}   .= "$feature3->{prod}";
					$featureV{'start'}  .= "$feature3->{start}";
					$featureV{'stop'}   .= "$feature3->{stop}";
					$featureV{'info'}   .= "$feature3->{info}";
				} else {
					# there is no defined 3' feature
					# we can't record any info, so just say none
					$featureV{'type'}   .= "none";
					$featureV{'id'}     .= "none";
					$featureV{'prodid'} .= "none";
					$featureV{'prod'}   .= "none";
					$featureV{'start'}  .= "none";
					$featureV{'stop'}   .= "none";
					$featureV{'info'}   .= "none";
				}

			} else {
				# 5' feature is defined, but variant position is between the start and stop of the 5' feature
				# in other words, the variant is within the 5' feature

				$featureV{'type'}   = "$feature5->{type}";
				$featureV{'id'}     = "$feature5->{id}";
				$featureV{'prodid'} = "$feature5->{prodid}";
				$featureV{'prod'}   = "$feature5->{prod}";
				$featureV{'start'}  = "$feature5->{start}";
				$featureV{'stop'}   = "$feature5->{stop}";
				$featureV{'info'}   = "$feature5->{info}";

				if ($feature5->{'type'} eq 'repeat_region') {
					$variant_loc  = 'RR';
					$variant_plac = "$variant->{ref}" . ($vstart - $feature5->{'start'});
					$variant_nt   = "$variant_plac";
					$variant_nt  .= " -> $variant->{alt}" if $variant->{'vtype'} eq 'SNP';
					
				} else {
					$variant_loc = 'CDS';

					$variant_nt = (($vstart - $feature5->{'start'} + 3) % 3) + 1;
					my $codon_ref    = substr $genomes{$chr}, ($vstart - $variant_nt), 3;
					my $aa_ref_num   = floor(($vstart - $feature5->{'start'})/3)+1;
					my $aa_ref       = $AA{$codon_ref};
			
					my $pos = (($vstart - $feature5->{'start'} + 3) % 3) + 1;
					if ($variant->{'vtype'} eq 'SNP' and defined $fna) {
				
						$variant_plac = "$aa_ref$aa_ref_num";
						$variant_aa   = "$aa_ref$aa_ref_num";

						my $i=0;
						foreach my $alt_nuc (split /\s*,\s*/, $variant->{'alt'}) {
							$variant_aa .= ' | ' if $i>0;
							$variant_nt .= ' | ' if $i>0;
							$i++;
				
							my $codon_alt = substr $genomes{$chr}, ($vstart - $pos), 3;
							substr($codon_alt, $pos-1, 1, $alt_nuc);
							$variant_nt .= "$codon_ref -> $codon_alt";
				
							my $aa_alt=$AA{$codon_alt};
					
							if ("$aa_ref" eq "$aa_alt") {
								$variant_aa.= ' (silent)';
							} else {
								$variant_aa.= " -> $aa_alt";
							}
						}
				
					} else {
						$variant_plac = "$aa_ref$aa_ref_num";
						$variant_aa   = "$aa_ref$aa_ref_num (N$variant_nt)";
						$variant_nt   = "$codon_ref (N$variant_nt)";
					}
				}
			}
		} else {
			# 5' feature is undefined
			$featureV{'type'}   = "none | ";
			$featureV{'id'}     = "none | ";
			$featureV{'prodid'} = "none | ";
			$featureV{'prod'}   = "none | ";
			$featureV{'start'}  = "none | ";
			$featureV{'stop'}   = "none | ";
			$featureV{'info'}   = "none | ";
			if (defined $feature3) {
				$featureV{'type'}   .= "$feature3->{type}";
				$featureV{'id'}     .= "$feature3->{id}";
				$featureV{'prodid'} .= "$feature3->{prodid}";
				$featureV{'prod'}   .= "$feature3->{prod}";
				$featureV{'start'}  .= "$feature3->{start}";
				$featureV{'stop'}   .= "$feature3->{stop}";
				$featureV{'info'}   .= "$feature3->{info}";
			} else {
				$featureV{'type'}   .= "none";
				$featureV{'id'}     .= "none";
				$featureV{'prodid'} .= "none";
				$featureV{'prod'}   .= "none";
				$featureV{'start'}  .= "none";
				$featureV{'stop'}   .= "none";
				$featureV{'info'}   .= "none";
			}
		}
		
		
		print "$variant->{vid}\t";
		print "$chr\t";
		print "$featureV{type}\t";
		print "$featureV{id}\t";
		print "$featureV{prodid}\t";
		print "$featureV{prod}\t";
		print "$featureV{start}\t";
		print "$featureV{stop}\t";
		print "$vstart\t";
		print "$variant->{vtype}\t";
		print "$variant_loc\t";
		print "$variant_plac\t";
		print "$variant_nt\t";
		print "$variant_aa\t";
		print "$variant->{ref}\t";
		print "$variant->{alt}\t";
		print "$variant->{quality}\t";
		print "$variant->{depth}\t";
		print "$variant->{info}\t";
		print "$featureV{info}\n";
		
	}
	
}


exit;


