#!/usr/bin/perl 
#
# Copyleft 2024 O. Sotolongo <osotolongo@fundacioace.org>  
#
# This script is intended for build the gVCF files from Fasta
use strict; 
use warnings; 
use SLURMACE; 
use File::Find::Rule;
use File::Basename;
use Cwd;
use Data::Dump qw(dump); 
############################################# 
# See: 
#   - For WES pipeline: http://detritus.fundacioace.com/wiki/doku.php?id=genetica:wes 
#   - For execution into SLURM: https://github.com/asqwerty666/acenip/blob/main/doc/SLURMACE.md 
############################################# 
#
#
# Data Paths 
# 
my $ref_dir = '/nas/Genomica/01-Data/00-Reference_files/02-GRCh38/00_Bundle/'; 
my $ref_name = 'Homo_sapiens_assembly38'; 
my $ref_fa = $ref_dir.'/'.$ref_name.'.fasta'; 
my $tmp_shit = $ENV{TMPDIR} || '/ruby/'.$ENV{USER}.'/tmp/'; 
my $known1 = 'Homo_sapiens_assembly38.known_indels.vcf.gz';
my $known2 = 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz';
my $dbsnp = 'Homo_sapiens_assembly38.dbsnp138.vcf';
my $hcsnps = '1000G_phase1.snps.high_confidence.hg38.vcf.gz';
# 
#
# Executable Paths 
#
#
my $fastqc = '/nas/usr/local/bin/fastqc'; 
my $bwa = '/nas/usr/local/bin/bwa mem -t 4 -M'; 
my $samtools = '/nas/software/samtools/bin/samtools'; 
my $verifyBamID = '/nas/usr/local/bin/verifyBamID'; 
my $freemix = 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/singularity/freemix.simg VerifyBamID --SVDPrefix /scripts/1000g.phase3.10k.b38.exome.vcf.gz.dat --NumThread 8 --max-depth 1000 --DisableSanityCheck';
my $gatk = 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/gatk4.simg gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"'; 
my $snpEff = 'java -Xmx8g -jar /nas/software/snpEff/snpEff.jar';
# 
#
# 
# Get CLI inputs 
#
#
my $cfile; 
my $init;
my %wesconf;
my $debug = 0; 
my $test = 0; 
while (@ARGV and $ARGV[0] =~ /^-/) {
	$_ = shift;         
	last if /^--$/;         
	if (/^-c/) { $cfile = shift; chomp($cfile);}         
	if (/^-i/) { $init = shift; chomp($init);}         
	if (/^-g/) { $debug = 1;}         
	if (/^-t/) { $test = 1;}         
} 
die "Should supply init data file\n" unless $init;
open IDF, "<$init";
while (<IDF>){         
	if (/^#.*/ or /^\s*$/) { next; }         
	my ($n, $v) = /(\S*)\s*=\s*(\S*)/;         
	$wesconf{$n} = $v; 
} 
close IDF;
my $workdir = getcwd;
$wesconf{outdir} = $workdir.'/output' unless $wesconf{outdir};
mkdir $wesconf{outdir} unless -d $wesconf{outdir};
my $slurmdir = $wesconf{outdir}.'/slurm';
mkdir $slurmdir unless -d $slurmdir;
# Do you want to process just a subset? Read the supplied list of subjects 
my @plist; 
if ($cfile and -f $cfile) {         
	open my $handle, "<$cfile";         
	chomp (@plist = <$handle>);         
	close $handle; 
}
my %ptask = (cpus => 8, time => '24:0:0', mem_per_cpu => '4G', debug => $test);
die "No such directory mate\n" unless -d $wesconf{src_dir};
my @content = find(file => 'name' => qr/$wesconf{search_pattern}/, in => $wesconf{src_dir});
# Clean the arrays 
@content = grep {!/.*$wesconf{cleaner}.*/} @content;
my %pollos = map {/.*\/(\w+?)$wesconf{search_pattern}$/; $1 => $_} @content;
my @jobs;
foreach my $pollo (sort keys %pollos){
	my @ljobids;
	my $go = 0;
	if ($cfile) {
		if (grep {/$pollo/} @plist) {$go = 1;}
	}else{
		$go = 1;
	}
	if (-f $pollos{$pollo} and $go){
		my $rdir = "$wesconf{outdir}/$pollo/results";
		my $tdir = "$wesconf{outdir}/$pollo/tmp";
		my $pofile = basename($pollos{$pollo});
		$ptask{job_name} = $pollo.'_sortIndex';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_sortIndex.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_sortIndex.out';
		$ptask{command} = "mkdir -p $rdir; mkdir -p $tdir; cp $pollos{$pollo} $tdir/$pofile\n";
		$ptask{command}.= "$gatk SortSam -I $tdir/$pofile -O $tdir/$pollo"."_sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true --TMP_DIR $tmp_shit\n";
		if(exists($ptask{'dependency'})){ delete($ptask{'dependency'}) };
		my $jid = send2slurm(\%ptask);
		# MarkDuplicates
		$ptask{job_name} = $pollo.'_markDuplicates';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_markDuplicates.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_markDuplicates.out';
		$ptask{command} = "$gatk MarkDuplicates -I $tdir/$pollo"."_sorted.bam -O $tdir/$pollo"."_rmdups.bam --METRICS_FILE $rdir/$pollo"."_metrics.txt --QUIET TRUE --MAX_RECORDS_IN_RAM 2000000 --ASSUME_SORTED TRUE --CREATE_INDEX TRUE --TMP_DIR $tmp_shit\n";
		$ptask{dependency} = "afterok:$jid";
		$jid = send2slurm(\%ptask);
		# VerifyBamID (freemix)
		$ptask{job_name} = $pollo.'_verifyBamID';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_verifyBamID.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_verifyBamID.out';
		$ptask{command} = "$freemix --BamFile $tdir/$pollo"."_rmdups.bam --Reference $ref_fa --Output $rdir/$pollo".".vbid2\n";
		$ptask{dependency} = "afterok:$jid";
		my $jid0 = send2slurm(\%ptask);
		push @ljobids, $jid0;
		# BaseRecalibrator
		$ptask{job_name} = $pollo.'_baseRecalibrator';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_baseRecalibrator.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_baseRecalibrator.out';
		my $unions = $wesconf{panel_dir}.'/'.$wesconf{unions};
		$ptask{command} = "$gatk BaseRecalibrator -I $tdir/$pollo"."_rmdups.bam -R $ref_fa -L $unions --known-sites $ref_dir/$known1 --known-sites $ref_dir/$known2 --known-sites $ref_dir/$dbsnp -O $rdir/$pollo"."_recal_data.table\n";
	      	$ptask{dependency} = "afterok:$jid";
		$jid = send2slurm(\%ptask);
		# ApplyBQSR, depende de BaseRecalibrator
		$ptask{job_name} = $pollo.'_applyBQSR';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_applyBQSR.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_applyBQSR.out';
		$ptask{command} = "$gatk ApplyBQSR -R $ref_fa -I $tdir/$pollo"."_rmdups.bam -L $unions -bqsr-recal-file $rdir/$pollo"."_recal_data.table -O $rdir/$pollo"."_recal.bam\n";
		$ptask{dependency} = "afterok:$jid";
		$jid0 = send2slurm(\%ptask);
		# AnalyzeCovariates, depende de BaseRecalibrator
		$ptask{job_name} = $pollo.'_analyzeCovariates';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_analyzeCovariates.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_analyzeCovariates.out';
		$ptask{command} = "$gatk AnalyzeCovariates -bqsr $rdir/$pollo"."_recal_data.table --plots $rdir/$pollo"."_AnalyzeCovariates.pdf\n";
		$ptask{dependency} = "afterok:$jid";
		my $jid1 = send2slurm(\%ptask);
		push @ljobids, $jid1;
		# CollectWGSMetrics, depende de ApplyBQSR
		$ptask{cpus} = 4;
		$ptask{job_name} = $pollo.'_collectRawMetrics';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_collectRawMetrics.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_collectRawMetrics.out';
		$ptask{command} = "$gatk  DepthOfCoverage -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_raw_wgs_metrics.txt -R $ref_fa -L $unions --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 40 --summary-coverage-threshold 50 --summary-coverage-threshold 60 --summary-coverage-threshold 70 --summary-coverage-threshold 80 --summary-coverage-threshold 90 --summary-coverage-threshold 100 --omit-depth-output-at-each-base true --omit-interval-statistics true --omit-locus-table true\n";
		$ptask{dependency} = "afterok:$jid0";
		$jid1 = send2slurm(\%ptask);
		push @ljobids, $jid1;
		$ptask{job_name} = $pollo.'_collectWgsMetrics';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_collectWgsMetrics.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_collectWgsMetrics.out';
		$ptask{command} = "$gatk  DepthOfCoverage -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_wgs_metrics.txt -R $ref_fa -L $unions --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 40 --summary-coverage-threshold 50 --summary-coverage-threshold 60 --summary-coverage-threshold 70 --summary-coverage-threshold 80 --summary-coverage-threshold 90 --summary-coverage-threshold 100 --omit-depth-output-at-each-base true --omit-interval-statistics true --omit-locus-table true --min-base-quality 20 -RF MappingQualityReadFilter --minimum-mapping-quality 20\n";
		$ptask{dependency} = "afterok:$jid0";
		$jid1 = send2slurm(\%ptask);
		push @ljobids, $jid1;
		$ptask{job_name} = $pollo.'_collectPaddedMetrics';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_collectPaddedMetrics.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_collectPaddedMetrics.out';
		$ptask{command} = "$gatk  DepthOfCoverage -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_padded_wgs_metrics.txt -R $ref_fa -L $unions --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 40 --summary-coverage-threshold 50 --summary-coverage-threshold 60 --summary-coverage-threshold 70 --summary-coverage-threshold 80 --summary-coverage-threshold 90 --summary-coverage-threshold 100 --omit-depth-output-at-each-base true --omit-interval-statistics true --omit-locus-table true -ip 100 --min-base-quality 20 -RF MappingQualityReadFilter --minimum-mapping-quality 20\n";
		$ptask{dependency} = "afterok:$jid0";
		$jid1 = send2slurm(\%ptask);
		push @ljobids, $jid1;
		# HaplotypeCaller, depende de ApplyBQSR
		$ptask{cpus} = 8;
		$ptask{job_name} = $pollo.'_haplotypeCaller';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_haplotypeCaller.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_haplotypeCaller.out';
		$ptask{command} = "$gatk HaplotypeCaller -R $ref_fa -L $unions -I $rdir/$pollo"."_recal.bam -G StandardAnnotation -G AS_StandardAnnotation -ERC GVCF --dbsnp $ref_dir/$dbsnp -O $rdir/$pollo"."_raw.snps.indels.g.vcf.gz\n";
		$ptask{command}.= "$gatk VariantEval -R $ref_fa -L $unions -D $ref_dir/$hcsnps -O $rdir/$pollo"."_eval.gatkreport --eval $rdir/$pollo"."_raw.snps.indels.g.vcf.gz\n";
		$ptask{dependency} = "afterok:$jid0";
		$jid1 = send2slurm(\%ptask);
		push @ljobids, $jid1;
		# Clean tmp files for subject
		$ptask{job_name} = $pollo.'_closeSubject';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_closeSubject.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_closeSubject.out';
		$ptask{dependency} = 'afterok:'.join(',afterok:', @ljobids);
		$ptask{command} = $debug?":\n":"rm -rf $tdir\n";
		my $fjob = send2slurm(\%ptask);
		push @jobs, $fjob;
	}
}
unless ($test) {
	my %wtask = (cpus => 1, job_name => 'end_wes', filename => $slurmdir.'/end.sh', output => $slurmdir.'/end.out', dependency => 'afterok:'.join(',afterok:',@jobs));
	send2slurm(\%wtask);
}
