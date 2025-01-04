#!/usr/bin/perl
#
# Copyleft O. Sotolongo-Grau (osotolongo@fundacioace.org)
#
# Script para convertir los fastq en CRAM
#
use strict;
use warnings;
use SLURMACE qw(send2slurm);
use File::Find::Rule;
use Cwd;
use File::Temp qw( :mktemp tempdir);
use Data::Dump qw(dump);
my $cfile;
my $outdir;
my %wesconf;
my $workdir = getcwd;
my $debug = 0;
my $init;
my $tmpdir =  $ENV{TMPDIR};
############################################################
# Variables con los PATHS. Cambiar aqui lo que sea necesario
#############################################################
#my $src_dir =  '/nas/Genomica/01-Data/02-WXS/01-Raw.data/202211_WES_PSP-DEGESCO/fqdata/';
my $ref_dir = '/nas/Genomica/01-Data/00-Reference_files/02-GRCh38/00_Bundle/';
my $ref_name = 'Homo_sapiens_assembly38';
my $tmp_shit = $ENV{TMPDIR};
#my $search_pattern = '_1.fq.gz';
#my $alt_pattern = '_2.fq.gz';
my $ref_fa = $ref_dir.'/'.$ref_name.'.fasta';
#my $platform = 'ILLUMINA';
#my $libraries = 'KAPPA_TE';
#################################################################
#################################################################
#################################################################
my $samtools = '/nas/software/samtools/bin/samtools';
my $gatk = 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /greebo:/greebo /nas/usr/local/opt/gatk4.simg gatk --java-options "-Djava.io.tmpdir='.$ENV{'TMPDIR'}.' -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"';
my $bwa = '/nas/usr/local/bin/bwa mem -t 4 -M';

@ARGV = ("-h") unless @ARGV;
while (@ARGV and $ARGV[0] =~ /^-/) {
	$_ = shift;
	last if /^--$/;
	if (/^-c/) { $cfile = shift; chomp($cfile);}
	if (/^-g/) { $debug = 1;}
	if (/^-i/) { $init = shift; chomp($init);}
}
die "Should supply init data file\n" unless $init;
open IDF, "<$init";
while (<IDF>){
	if (/^#.*/ or /^\s*$/) { next; }
	my ($n, $v) = /(\S*)\s*=\s*(\S*)/;
	$wesconf{$n} = $v;
}
close IDF;
$wesconf{outdir} = $workdir.'/output' unless $wesconf{outdir};
mkdir $wesconf{outdir} unless -d $wesconf{outdir};
my $slurmdir = $wesconf{outdir}.'/slurm';
mkdir $slurmdir unless -d $slurmdir;

my %ptask = (cpus => 8, job_name => 'wes', time => '8:0:0', mem_per_cpu => '4G', debug => $debug);

die "No such directory mate\n" unless -d $wesconf{src_dir};
my @content = find(file => 'name' => "*$wesconf{search_pattern}", in => $wesconf{src_dir});
my %pollos = map {/.*\/(\w+?)$wesconf{search_pattern}$/; $1 => $_} @content;
my @cuts;
if ($cfile) {
	open IDF, "<$cfile" or die "No such input file!\n";
	@cuts = <IDF>;
	chomp @cuts;
	close IDF;
}


foreach my $shit (sort keys %pollos){
	my $go = 0;
	if ($cfile) {
		if (grep {/$shit/} @cuts) {$go = 1;}
	}else{
		$go = 1;
	}
	if (-f $pollos{$shit} and $go) {
		$ptask{job_name} = 'compact_data';
		$ptask{filename} = $slurmdir.'/'.$shit.'.sh';
		$ptask{output} = $slurmdir.'/'.$shit.'.out';
		(my $another = $pollos{$shit}) =~ s/$wesconf{search_pattern}/$wesconf{alt_pattern}/;
		my $rg = '"@RG\\tID:'.$shit.'\\tPL:'.$wesconf{platform}.'\\tLB:'.$wesconf{libraries}.'\\tSM:'.$shit.'"';
		$ptask{command} = $bwa.' -R '.$rg.' '.$ref_fa.' '.$pollos{$shit}.' '.$another.' | '.$gatk.' SortSam -I /dev/stdin -O '.$tmpdir.'/'.$shit.'_sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true'." --TMP_DIR $tmp_shit\n";
		$ptask{command}.= $samtools.' view -@ 8 -T '.$ref_fa.' -C -o '.$wesconf{outdir}.'/'.$shit.'.cram '.$tmpdir.'/'.$shit.'_sorted.bam'."\n";
		$ptask{command}.= 'rm '.$tmpdir.'/'.$shit.'_sorted.bam';
		send2slurm(\%ptask);
	}
}

unless ($debug) {
	my %warn = (filename => $slurmdir.'/compact_end.sh', output => $slurmdir.'/compact_end.out', job_name => 'compact_data', mailtype => 'END', dependency => 'singleton');
	send2slurm(\%warn);
}
