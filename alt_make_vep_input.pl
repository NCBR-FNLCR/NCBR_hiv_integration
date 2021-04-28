#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';
use List::Util qw( min max );
use List::MoreUtils qw(uniq);

#INPUT
my $table=$ARGV[0];
my $out= $ARGV[1];
my @line=();
my @molbar=();
my $umi='';
my $nextumi='null';
my $count=0;
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
my @sample=();
my @samples=();
my @barcodes=();
my @reps=();
my @allreps=();

open C, ">$out";

open B, "<$table";
while (<B>){
	chomp;
	@line = split;
	push @chr, $line[0];
	push @start, $line[1];
	push @stop, $line[2];
	push @reps, $line[3];
	push @barcode, $line[4];
	push @mapq, $line[5];
	push @strand, $line[6];
	push @sample, $line[7];
}

my $a=0;

for ($a=0; $a<@chr; $a++) {
	push @mapqs, $mapq[$a];
	push @samples, $sample[$a];
	push @barcodes, $barcode[$a];
	push @allreps, $reps[$a];
	if ($strand[$a] eq '+') {
		$code = $chr[$a] . '_' . $start[$a] . '_' . $strand[$a];
		$nextcode = $chr[$a+1] . '_' . $start[$a+1] . '_' . $strand[$a+1];
	}
	else {
		$code = $chr[$a] . '_' . $stop[$a] . '_' . $strand[$a];
		$nextcode = $chr[$a+1] . '_' . $stop[$a+1] . '_' . $strand[$a+1];
	}
#	print "$umi vs $nextumi\n";
	if ($code ne $nextcode) {
		print C "$chr[$a]\t$start[$a]\t$stop[$a]\tINS\t$strand[$a]\t$code\t";
		my @uniqsam = uniq @samples;
		print C join(',', @allreps), "\t";
		print C join(',', @uniqsam), "\t";
		print C scalar @uniqsam. "\t";
		print C join(',', @mapqs), "\t";
		print C join(',', @barcodes), "\n";
#		print join(',', @samples), "\n";
		$count=0;
		@mapqs=();
		@samples=();
		@barcodes=();
		@allreps=();
	}
}

close B;
close C;