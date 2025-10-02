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
use File::Slurp qw(read_file);
use Cwd;
use File::Temp qw( :mktemp tempdir);
use FindBin; 
use lib "$FindBin::Bin";
use wxsInit;
use String::Random;
use Data::Dump qw(dump);
my $cfile;
my $outdir;
my $workdir = getcwd;
my $debug = 0;
my $init;
my $tmpdir =  $ENV{TMPDIR};
############################################################
# Variables con los PATHS. Cambiar aqui lo que sea necesario
#############################################################
my %dpaths = data_paths();
my $tmp_shit = $ENV{TMPDIR};
my $ref_fa = $dpaths{ref_dir}.'/'.$dpaths{ref_name}.'.fasta';
#################################################################
#################################################################
#################################################################
my %epaths = exec_paths();

@ARGV = ("-h") unless @ARGV;
while (@ARGV and $ARGV[0] =~ /^-/) {
	$_ = shift;
	last if /^--$/;
	if (/^-c/) { $cfile = shift; chomp($cfile);}
	if (/^-g/) { $debug = 1;}
	if (/^-i/) { $init = shift; chomp($init);}
}
die "Should supply init data file\n" unless $init;
my %wesconf = init_conf($init);
$wesconf{outdir} = $workdir.'/output' unless $wesconf{outdir};
mkdir $wesconf{outdir} unless -d $wesconf{outdir};
my $slurmdir = $wesconf{outdir}.'/slurm';
mkdir $slurmdir unless -d $slurmdir;

my %ptask = (cpus => 8, job_name => 'wes', time => '8:0:0', mem_per_cpu => '4G', debug => $debug);

my $debug_file = $wesconf{outdir}.'/.debug_report_'.getLoggingTime().'.log';
open STDOUT, ">$debug_file" or die "Can't redirect stdout";
open STDERR, ">&STDOUT" or die "Can't dup stdout";
my %sname;
if (exists($wesconf{samples})  and $wesconf{samples}){
	%sname = map { /(.*),(.*)/; $2 =>$1 } read_file $wesconf{samples};
}

die "No such directory mate\n" unless -d $wesconf{src_dir};
my $chick_rule = File::Find::Rule->new;
my @chicks = $chick_rule->directory->mindepth(1)->maxdepth(1)->in($wesconf{src_dir});
foreach my $chick_path (@chicks){
	my @content = find(file => 'name' => "*$wesconf{search_pattern}", in => $chick_path);
	my %pollos = map {/.*\/(\w+?)$wesconf{search_pattern}$/; $1 => $_} @content;
	my @chick_group;
	my @chick_tmp;
	my @jobs;
	my ($chick_name) = $chick_path =~ /.*\/(.*?)$/;
	delete $ptask{dependency};
	foreach my $shit (sort keys %pollos){
		if (-f $pollos{$shit}) {
			my $ultra = String::Random->new;
			my $ultra_shit = $shit.$ultra->randpattern("CCCCCC");
			my $sample;
			if (exists($sname{$chick_name}) and $sname{$chick_name}){
				$sample = $sname{$chick_name};
			}else{
				($sample = $shit) =~  s/_\d+$//;
			}
			$ptask{job_name} = 'create_bam_'.$ultra_shit;
			$ptask{filename} = $slurmdir.'/'.$ultra_shit.'.sh';
			$ptask{output} = $slurmdir.'/'.$ultra_shit.'.out';
			(my $another = $pollos{$shit}) =~ s/$wesconf{search_pattern}/$wesconf{alt_pattern}/;
			my $rg = '"@RG\\tID:'.$shit.'\\tPL:'.$wesconf{platform}.'\\tLB:'.$wesconf{libraries}.'\\tSM:'.$sample.'"';
			$ptask{command} = $epaths{bwa}.' -R '.$rg.' '.$ref_fa.' '.$pollos{$shit}.' '.$another.' | '.$epaths{gatk}.' SortSam -I /dev/stdin -O '.$tmpdir.'/'.$ultra_shit.'_sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true'." --TMP_DIR $tmp_shit\n";
			$ptask{command}.= $epaths{samtools}.' calmd -bAr '.$tmpdir.'/'.$ultra_shit.'_sorted.bam '.$ref_fa.' > '.$tmpdir.'/'.$ultra_shit.'_calmd.bam'."\n";
			push @chick_group, $tmpdir.'/'.$ultra_shit.'_calmd.bam';
			push @chick_tmp, $tmpdir.'/'.$ultra_shit.'_sorted.bam';
			my $job_id = send2slurm(\%ptask); 
			push @jobs, $job_id;
		}
	}
	my $rname = String::Random->new;
	my $chick = $rname->randpattern("CCCCCCCC");
	print "$chick,$chick_path\n";
	$ptask{job_name} = 'compact_data';
	$ptask{filename} = $slurmdir.'/'.$chick.'.sh';
	$ptask{output} = $slurmdir.'/'.$chick.'.out';
	my $chl = join ' ',@chick_group;
	my $chtmp = join ' ',@chick_tmp;
	$ptask{command} = $epaths{samtools}.' merge '.$tmpdir.'/'.$chick.'.bam '.$chl."\n";
	#$ptask{command}.= $epaths{gatk}.' MarkDuplicates -I '.$tmpdir.'/'.$chick.'.bam -O '.$tmpdir.'/'.$chick.'_rmdup.bam --METRICS_FILE '.$wesconf{outdir}.'/'.$chick.'_metrics.txt --QUIET TRUE --MAX_RECORDS_IN_RAM 2000000 --ASSUME_SORTED TRUE --CREATE_INDEX TRUE --TMP_DIR '.$tmp_shit."\n";
	#$ptask{command}.= $epaths{samtools}.' view -@ 8 -T '.$ref_fa.' -C -o '.$wesconf{outdir}.'/'.$chick.'.cram '.$tmpdir.'/'.$chick.'_rmdup.bam'."\n";
	$ptask{command}.= $epaths{samtools}.' view -@ 8 -T '.$ref_fa.' -C -o '.$wesconf{outdir}.'/'.($sname{$chick_name}?$sname{$chick_name}:$chick).'.cram '.$tmpdir.'/'.$chick.'.bam'."\n";
	$ptask{command}.= $epaths{gatk}.' ValidateSamFile -I '.$wesconf{outdir}.'/'.($sname{$chick_name}?$sname{$chick_name}:$chick).'.cram -O '.$wesconf{outdir}.'/'.($sname{$chick_name}?$sname{$chick_name}:$chick).'.validation --MODE SUMMARY -R '.$ref_fa."\n";
	#$ptask{command}.= 'rm '.$tmpdir.'/'.$chick.'.bam '.$chl.' '.$tmpdir.'/'.$chick.'_rmdup.bam';
	$ptask{command}.= 'rm '.$tmpdir.'/'.$chick.'.bam '.$chl.' '.$chtmp;
	$ptask{dependency} = 'afterok:'.join(',afterok:',@jobs);
	send2slurm(\%ptask);
}
unless ($debug) {
	my %warn = (filename => $slurmdir.'/compact_end.sh', output => $slurmdir.'/compact_end.out', job_name => 'compact_data', mailtype => 'END', dependency => 'singleton');
	send2slurm(\%warn);
}
