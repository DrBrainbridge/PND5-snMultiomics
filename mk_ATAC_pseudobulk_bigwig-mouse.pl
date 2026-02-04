#!/usr/bin/perl

## Requirements:
##  - List of cells from the cluster/genotype/whatever of interest.  {$cellListFile}
##  - BAM file from CellRanger.  {$inputBAM}
##  - ~.chromSizes file  {$chromSizes}
##  - Tools needed = samtools, bedtools, and bedGraphToBigWig.  {see lines 30-32}
##
## Steps:
##  - Pull out reads from cells of interest. Filter to keep only properly paired and ignore PCR/optical duplicates.
##  - Replace SAM header and filter to keep hits to canonical chrom only.
##  - Apply coordinate-sort to BAM file. (output #1)
##  - Apply queryname-sort, then convert to PE BED format. Strip out unnecessary columns to make fragment BED format. Sort by coordinate.
##  - Convert to bedGraph then bigWig format. (output #2)
##  - Get usable fragment count ... might want this if planning to depth-normalize coverage tracks downstream. (output #3)
##
## Important note!
##  - This version of the script is hard-coded to assume mouse data, uses only chr1-19,X,Y. If you want something different, modify @chrlist accordingly.

## User must define the following:
$cellListFile = $ARGV[0]; # txt file with list of cells, one cell ID per line
$inputBAM = $ARGV[1];     # expects path to atac_possorted_bam.bam file from CellRanger
$outroot = $ARGV[2];      # some identifier for naming output files
$workingDir = $ARGV[3];   # tmp folder; this script will not delete it afterwards, but can be scrapped by user if everything seemed to work
$chromSizes = $ARGV[4];   # mm10.chromSizes
unless (-e $cellListFile && -e $inputBAM && -e $chromSizes) { print "ERROR: One or more expected input files not found; exiting.\n"; exit; }
if (-d $workingDir) { print "ERROR: Intended tmp folder \'$workingDir\' already exists.\n"; }  # disallow to avoid overwriting existing files with new intermediates

## Tools:
$samtools_path = "fixme";
$bedtools_path = "fixme";
$bg2bw_tool_path = "fixme";

## Genome reference info:
@chrlist = (); foreach $c (1..19) { push @chrlist, "chr$c"; } push @chrlist, "chrX","chrY";   ## this is the chrom sort order to be used
%allowedChrom = map {$_ => 1} @chrlist;

## These are the (final) outputs to be generated:
$outputBAM = "$outroot.pseudobulkATAC.bam";
$outputBW = "$outroot.pseudobulkATAC.bigWig";
$outputFragCt = "$outroot.pseudobulkATAC.usable_fragments.txt";
if (-e $outputBAM || -e $outputBW || -e $outputFragCt) { print "ERROR: One or more intended output files already exists; exiting.\n"; exit; }

## Intermediate files; can be deleted at the end if everything worked properly...
system "mkdir $workingDir";
$samfile1 = "$workingDir/$outroot.tmp1.sam";
$samfile2 = "$workingDir/$outroot.tmp2.sam";
$bamfileQuerySort = "$workingDir/$outroot.qsort.bam";
$bedfilePE = "$workingDir/$outroot.PE.bed";
$bedfileFrag = "$workingDir/$outroot.frag.unsorted.bed";
$bedfileFragSorted = "$workingDir/$outroot.frag.coordsort.bed";
$bgfile = "$workingDir/$outroot.bedGraph";


%cells = (); open(IN, "$cellListFile"); while (<IN>) { chomp $_; $cells{$_} = 1; } close(IN);
open(TMPSAM1, ">$samfile1");
open(IN, "$samtools_path view -h -f2 -F1024 $inputBAM |");
while (<IN>) {
  if ($_ =~ /^\@/) { print TMPSAM1 "$_"; next; }
  foreach $tag (split/\t/, $_) {
    next unless ($tag =~ /^CB:Z:/);
    @ar = split/\:/, $tag;
    if ($cells{$ar[2]} == 1) { print TMPSAM1 "$_"; }
    last;
  }
}
close(IN); close(TMPSAM1);
%cs = (); open(CS, "$chromSizes"); while (<CS>) { chomp $_; ($c, $s) = split/\t/, $_; $cs{$c} = $s; } close(CS);
open(TMPSAM2, ">$samfile2");
print TMPSAM2 "\@HD\tSO:unsorted\n";
foreach $chr (@chrlist) { print TMPSAM2 "\@SQ\tSN:$chr\tLN:$cs{$chr}\n"; }
open(IN, "$samtools_path view -S --no-PG $samfile1 |");
while (<IN>) {
  @ar = split/\t/, $_;
  if ($allowedChrom{$ar[2]} == 1) { print TMPSAM2 "$_"; }
}
close(IN); close(TMPSAM2);
system "$samtools_path sort --no-PG $samfile2 > $outputBAM";
system "$samtools_path sort -n --no-PG $outputBAM > $bamfileQuerySort";
system "$bedtools_path bamtobed -i $bamfileQuerySort -bedpe > $bedfilePE";
system "cut -f1,2,6,7 $bedfilePE > $bedfileFrag";
foreach $chr (@chrlist) {
  system "grep -a -w ^$chr $bedfileFrag | sort -k2,2n -k3,3n >> $bedfileFragSorted";
}
system "$bedtools_path genomecov -bg -g $chromSizes -i $bedfileFragSorted > $bgfile";
system "$bg2bw_tool_path $bgfile $chromSizes $outputBW";

$x = `cat $bedfileFragSorted | wc -l`; $x =~ s/\s//g;
open(OUT, ">$outputFragCt"); print OUT "$x\n"; close(OUT);
$bedfileFragSorted = "$workingDir/$outroot.frag.coordsort.bed";


print "\nThe bedGraph file is at $bgfile ... might be worth keeping if you intend to apply a depth normalization to the coverage tracks in post-processing.\n";
print "Otherwise, if everything seems OK, you can delete the contents of $workingDir/.\n\n";

