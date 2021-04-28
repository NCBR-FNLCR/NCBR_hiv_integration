#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';
use List::Util qw( min max );

#INPUT
my $table=$ARGV[0];
my $out= $ARGV[1];
my $sample = $ARGV[2];
my @line=();
my @molbar=();
my $umi='';
my $nextumi='null';
my @chr=();
my @start=();
my @stop=();
my @barcode=();
my @mapq=();
my @strand=();
my $code='';
my @mapqs=();
my @nextmolbar=();
my $nextcode='';
my $tag='';
my $maxmapq='';

open C, ">$out";

open B, "<$table";
while (<B>){
	chomp;
	@line = split;
	push @chr, $line[0];
	push @start, $line[1];
	push @stop, $line[2];
	push @barcode, $line[3];
	push @mapq, $line[4];
	push @strand, $line[5];
}

my $a=0;
my @uniques = ();
my @bases = ();
my @uniquesbases = ();
my $diffs=0;

###Find the unique barcodes

for ($a=0; $a<@chr; $a++) {
	@molbar = split '_', $barcode[$a];
	if ($strand[$a] eq '+') {
		$code = $chr[$a] . '_' . $start[$a];
	}
	else {
		if ($strand[$a] eq '-') {
			$code = $chr[$a] . '_' . $stop[$a];
		}
	}

	my $z=0;
	my $status='null';

	for ($z=0; $z<@uniques; $z++) {
		my @cats = split '_', $uniques[$z];
		my $location = $cats[0] . '_' . $cats[1];
		if ($code eq $location) {
			@bases = split '', $molbar[0];
			@uniquesbases = split '', $cats[2];
			my $m=0;
			$diffs=0;
			for ($m=0; $m<@bases; $m++) {
				if ($bases[$m] ne $uniquesbases[$m]) {
					$diffs++;
				}
			}
			if ($diffs < 2) {
				$status='skip';
			}
		}
	}
	if ($status ne 'skip') {
		push @uniques, ($code . '_' . $molbar[0]);
	}
}

###Collapse insertions that are duplicate barcodes

my $x=0;
my @chrb=();
my @startb=();
my @stopb=();
my @tagb=();
my @maxmapqb=();
my @strandb=();

for ($x=0; $x<@uniques; $x++) {
	my @cats = split '_', $uniques[$x];
	my $location = $cats[0] . '_' . $cats[1];

	my @chra=();
	my @starta=();
	my @stopa=();
	my @barcodea=();
	my @mapqa=();
	my @stranda=();

	my $y=0;

	for ($y=0; $y<@chr; $y++) {
		@molbar = split '_', $barcode[$y];
		if ($strand[$y] eq '+') {
			$code = $chr[$y] . '_' . $start[$y];
		}
		else {
			if ($strand[$y] eq '-') {
				$code = $chr[$y] . '_' . $stop[$y];
			}
		}
		if ($location eq $code) {
			@bases = split '', $molbar[0];
			@uniquesbases = split '', $cats[2];
			my $n=0;
			$diffs=0;
			for ($n=0; $n<@bases; $n++) {
				if ($bases[$n] ne $uniquesbases[$n]) {
					$diffs++;
				}
			}
			if ($diffs < 2) {
				push @chra, $chr[$y];
				push @starta, $start[$y];
				push @stopa, $stop[$y];
				push @barcodea, $barcode[$y];
				push @mapqa, $mapq[$y];
				push @stranda, $strand[$y];
			}
		}
	}


	my $d=0;
	my $count=0;
	my @mapqs=();
	my @molbar=();

	for ($d=0; $d<@chra; $d++) {
		@molbar = split '_', $barcodea[$d];
		$count+=$molbar[2];
		push @mapqs, $mapqa[$d];
	}

	@molbar = split '_', $barcodea[0];
	my $tag = $molbar[0] . '_molbar_' . $count;
	my $maxmapq = max @mapqs;
#	print C "$chr[$a]\t$start[$a]\t$stop[$a]\t$tag\t$maxmapq\t$strand[$a]\t$sample\n";
	push @chrb, $chra[0];
	push @startb, $starta[0];
	push @stopb, $stopa[0];
	push @tagb, $tag;
	push @maxmapqb, $maxmapq;
	push @strandb, $stranda[0];
}

my $b=0;
my @mapqsa=();
my @molbara=();
my $reps='';

for ($b=0; $b<@chrb; $b++) {
	push @molbara, $tagb[$b];
#	$count+=$molbar[2];
	push @mapqsa, $maxmapqb[$b];
#	$umi = substr $molbar[0], 0, 8;
	if ($strandb[$b] eq '+') {
		$code = $chrb[$b] . '_' . $startb[$b];
	}
	else {
		if ($strandb[$b] eq '-') {
			$code = $chrb[$b] . '_' . $stopb[$b];
		}
	}
#	@nextmolbar = split '_', $barcode[$a+1];
#	$nextumi = substr $nextmolbar[0], 0, 8;
	if ($strandb[$b] eq '+') {
		$nextcode = $chrb[$b+1] . '_' . $startb[$b+1];
	}
	else {
		if ($strandb[$b] eq '-') {
			$nextcode = $chrb[$b+1] . '_' . $stopb[$b+1];
		}
	}
#	print "$umi vs $nextumi\n";
	if ($code ne $nextcode) {
		$maxmapq = max @mapqsa;
		$reps = "$sample" . ':' . @molbara;
		print C "$chrb[$b]\t$startb[$b]\t$stopb[$b]\t$reps\t";
		my $c=0;
		for ($c=0; $c<@molbara; $c++) {
			print C "$molbara[$c]";
			if (($c+1)<@molbara) {
				print C ',';
			}
		}
		print C "\t$maxmapq\t$strandb[$b]\t$sample\n";
		@mapqsa=();
		@molbara=();
	}
}

close B;
close C;