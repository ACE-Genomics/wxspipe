#!/usr/bin/perl
#
# Copyleft 2025 O. Sotolongo <osotolongo@fundacioace.org>
#
# This is a simple way to initializer variables for WXS scripts
#
use strict; use warnings;
package wxsInit;
require Exporter;
our @ISA= qw( Exporter );
our @EXPORT_OK = qw (exec_paths data_paths init_conf);
our @EXPORT = qw(exec_paths data_paths init_conf);
our $VERSION = 0.1;

=head1 wxsInit

Just a way to initialize variable for all scripts. You should edit paths to reference library and executables inside the functions. Yep, maybe not the simplest way but the easiest one to implement in my mind

=over

=item exec_paths

Send executable paths to the scripts this way. Directly edit the paths here

=cut

sub exec_paths {
	my %epaths = ('fastqc' => '/nas/usr/local/bin/fastqc',
		'bwa' => '/nas/usr/local/bin/bwa mem -t 4 -M',
		'samtools' => '/nas/software/samtools/bin/samtools',
		'freemix' => 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /greebo:/greebo /nas/usr/local/opt/singularity/freemix.simg VerifyBamID --SVDPrefix /scripts/1000g.phase3.10k.b38.exome.vcf.gz.dat --NumThread 8 --max-depth 1000 --DisableSanityCheck',
		'gatk' => 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /greebo:/greebo /nas/usr/local/opt/gatk4.simg gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"',
	);
	return %epaths;
}

=item data_paths

Send reference data paths to the scripts this way. Directly edit the paths here

=cut 

sub data_paths {
	my %dpaths = ('ref_dir' => '/nas/Genomica/01-Data/00-Reference_files/02-GRCh38/00_Bundle/',
	'ref_name' => 'Homo_sapiens_assembly38',
	'known1' => 'Homo_sapiens_assembly38.known_indels.vcf.gz',
	'known2' => 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
	'dbsnp' => 'Homo_sapiens_assembly38.dbsnp138.vcf',
	'hcsnps' => '1000G_phase1.snps.high_confidence.hg38.vcf.gz',
	'omni' => '1000G_omni2.5.hg38.vcf.gz',
	'hapmap' => 'hapmap_3.3.hg38.vcf.gz',
	);
	return %dpaths;
}

=item init_conf

Read init file and export a hash of variables. 

Usage:
		
	my %wesconf = init_conf('project.init');

=cut

sub init_conf{
	my $ifile = shift;
	my %conf;
	open IDF, "<$ifile";
	while (<IDF>){
		if (/^#.*/ or /^\s*$/) { next; }
		my ($n, $v) = /(\S*)\s*=\s*(\S*)/;
		$conf{$n} = $v;
	}
	close IDF;
	return %conf;
}

=back
