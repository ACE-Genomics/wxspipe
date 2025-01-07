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
use FindBin;
use lib "$FindBin::Bin";
use wxsInit; 
############################################# 
# See: 
#   - For WES pipeline: http://detritus.fundacioace.com/wiki/doku.php?id=genetica:wes 
#   - For execution into SLURM: https://github.com/asqwerty666/acenip/blob/main/doc/SLURMACE.md 
############################################# 
#
# Data Paths 
# 
my %dpaths = data_paths();
my $ref_dir = $dpaths{ref_dir}; 
my $ref_name = $dpaths{ref_name}; 
my $ref_fa = $ref_dir.'/'.$ref_name.'.fasta'; 
my $tmp_shit = $ENV{TMPDIR}; 
my $known1 = $dpaths{known1};
my $known2 = $dpaths{known2};
my $dbsnp = $dpaths{dbsnp};
my $hcsnps = $dpaths{hcsnps};
# 
#
# Executable Paths 
#
#
my %epaths = exec_paths();
my $fastqc = $epaths{fastqc}; 
my $bwa = $epaths{bwa}; 
my $samtools = $epaths{samtools}; 
my $freemix = $epaths{freemix};
my $gatk = $epaths{gatk};
# 
#
# 
# Get CLI inputs 
#
#
my $cfile; 
my $debug = 0; 
my $test = 0;
my $mode = 'wes';
my $init;
while (@ARGV and $ARGV[0] =~ /^-/) {
	$_ = shift;         
	last if /^--$/;         
	if (/^-c/) { $cfile = shift; chomp($cfile);}         
	if (/^-i/) { $init = shift; chomp($init);}         
	if (/^-g/) { $debug = 1;}         
	if (/^-t/) { $test = 1;}
	if (/^-m/) { $mode = shift; chomp($mode);}
} 
die "Should supply init data file\n" unless $init;
my %wesconf = init_conf($init);
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
my %ptask = (cpus => 8, time => ($mode  eq 'wgs')?'72:0:0':'24:0:0', mem_per_cpu => '4G', debug => $test);
die "No such directory mate\n" unless -d $wesconf{src_dir};
my @content = find(file => 'name' => "*$wesconf{search_pattern}", in => $wesconf{src_dir});
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
		my $qdir = "$wesconf{outdir}/$pollo/qc";
		my $rdir = "$wesconf{outdir}/$pollo/results";
		my $tdir = "$wesconf{outdir}/$pollo/tmp";
		(my $another = $pollos{$pollo}) =~ s/$wesconf{search_pattern}/$wesconf{alt_pattern}/;
		# FASTQC	
		$ptask{job_name} = $pollo.'_fastQc';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_fastQc.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_fastQc.out';
		$ptask{command} = "mkdir -p $qdir\n";
		$ptask{command}.= "$fastqc -o $qdir $pollos{$pollo} $another\n";
		if(exists($ptask{'dependency'})){ delete($ptask{'dependency'}) };
		my $jid = send2slurm(\%ptask);
		push @ljobids, $jid;
		# Create SAM
		$ptask{job_name} = $pollo.'_createSam';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_createSam.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_createSam.out';
		my $rg = '"@RG\\tID:'.$pollo.'\\tPL:'.$wesconf{platform}.'\\tLB:'.$wesconf{libraries}.'\\tSM:'.$pollo.'"';
		$ptask{command} = "mkdir -p $tdir; mkdir -p $rdir\n";
		$ptask{command} .= "$bwa -R ".$rg." $ref_fa $pollos{$pollo} $another | $gatk SortSam -I /dev/stdin -O $tdir/$pollo"."_sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true  --TMP_DIR $tmp_shit\n";
		$jid = send2slurm(\%ptask);
		# MarkDuplicates
		$ptask{job_name} = $pollo.'_markDuplicates';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_markDuplicates.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_markDuplicates.out';
		$ptask{command} = "$gatk MarkDuplicates -I $tdir/$pollo"."_sorted.bam -O $tdir/$pollo"."_rmdups.bam --METRICS_FILE $rdir/$pollo"."_metrics.txt --QUIET TRUE --MAX_RECORDS_IN_RAM 2000000 --ASSUME_SORTED TRUE --CREATE_INDEX TRUE  --TMP_DIR $tmp_shit\n";
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
		my $unions = ($mode eq 'wgs')?'':$wesconf{panel_dir}.'/'.$wesconf{unions};
		$ptask{command} = "$gatk BaseRecalibrator -I $tdir/$pollo"."_rmdups.bam -R $ref_fa".(($mode eq 'wgs')?' ':" -L $unions ")."--known-sites $ref_dir/$known1 --known-sites $ref_dir/$known2 --known-sites $ref_dir/$dbsnp -O $rdir/$pollo"."_recal_data.table\n";
	      	$ptask{dependency} = "afterok:$jid";
		$jid = send2slurm(\%ptask);
		# ApplyBQSR, depende de BaseRecalibrator
		$ptask{job_name} = $pollo.'_applyBQSR';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_applyBQSR.sh';
		$ptask{output} = $slurmdir.'/'.$pollo.'_applyBQSR.out';
		$ptask{command} = "$gatk ApplyBQSR -R $ref_fa -I $tdir/$pollo"."_rmdups.bam".(($mode eq 'wgs')?' ':" -L $unions ")."-bqsr-recal-file $rdir/$pollo"."_recal_data.table -O $rdir/$pollo"."_recal.bam\n";
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
		unless ($mode eq 'wgs'){
			# CollectWGSMetrics, depende de ApplyBQSR
			$ptask{cpus} = 4;
			$ptask{job_name} = $pollo.'_collectRawMetrics';
			$ptask{filename} = $slurmdir.'/'.$pollo.'_collectRawMetrics.sh';
			$ptask{output} = $slurmdir.'/'.$pollo.'_collectRawMetrics.out';
			$ptask{command} = "$gatk  DepthOfCoverage -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_raw_wgs_metrics.txt -R $ref_fa".(($mode eq 'wgs')?' ':" -L $unions ")."--summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 40 --summary-coverage-threshold 50 --summary-coverage-threshold 60 --summary-coverage-threshold 70 --summary-coverage-threshold 80 --summary-coverage-threshold 90 --summary-coverage-threshold 100 --omit-depth-output-at-each-base true --omit-interval-statistics true --omit-locus-table true\n";
			$ptask{dependency} = "afterok:$jid0";
			$jid1 = send2slurm(\%ptask);
			push @ljobids, $jid1;
			$ptask{job_name} = $pollo.'_collectWgsMetrics';
			$ptask{filename} = $slurmdir.'/'.$pollo.'_collectWgsMetrics.sh';
			$ptask{output} = $slurmdir.'/'.$pollo.'_collectWgsMetrics.out';
			$ptask{command} = "$gatk  DepthOfCoverage -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_wgs_metrics.txt -R $ref_fa".(($mode eq 'wgs')?' ':" -L $unions ")."--summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 40 --summary-coverage-threshold 50 --summary-coverage-threshold 60 --summary-coverage-threshold 70 --summary-coverage-threshold 80 --summary-coverage-threshold 90 --summary-coverage-threshold 100 --omit-depth-output-at-each-base true --omit-interval-statistics true --omit-locus-table true --min-base-quality 20 -RF MappingQualityReadFilter --minimum-mapping-quality 20\n";
			$ptask{dependency} = "afterok:$jid0";
			$jid1 = send2slurm(\%ptask);
			push @ljobids, $jid1;
			$ptask{job_name} = $pollo.'_collectPaddedMetrics';
			$ptask{filename} = $slurmdir.'/'.$pollo.'_collectPaddedMetrics.sh';
			$ptask{output} = $slurmdir.'/'.$pollo.'_collectPaddedMetrics.out';
			$ptask{command} = "$gatk  DepthOfCoverage -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_padded_wgs_metrics.txt -R $ref_fa".(($mode eq 'wgs')?' ':" -L $unions ")."--summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 30 --summary-coverage-threshold 40 --summary-coverage-threshold 50 --summary-coverage-threshold 60 --summary-coverage-threshold 70 --summary-coverage-threshold 80 --summary-coverage-threshold 90 --summary-coverage-threshold 100 --omit-depth-output-at-each-base true --omit-interval-statistics true --omit-locus-table true -ip 100 --min-base-quality 20 -RF MappingQualityReadFilter --minimum-mapping-quality 20\n";
			$ptask{dependency} = "afterok:$jid0";
			$jid1 = send2slurm(\%ptask);
			push @ljobids, $jid1;
		}else{
			$ptask{cpus} = 4;
			$ptask{job_name} = $pollo.'_collectWgsMetrics';
			$ptask{filename} = $slurmdir.'/'.$pollo.'_collectWgsMetrics.sh';
			$ptask{output} = $slurmdir.'/'.$pollo.'_collectWgsMetrics.out';
			$ptask{command} = "$gatk CollectWgsMetrics -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_wgs_metrics.txt -R $ref_fa\n";
			$ptask{dependency} = "afterok:$jid0";
			$jid1 = send2slurm(\%ptask);
			push @ljobids, $jid1;
			$ptask{job_name} = $pollo.'_collectRawMetrics';
			$ptask{filename} = $slurmdir.'/'.$pollo.'_collectRawMetrics.sh';
			$ptask{output} = $slurmdir.'/'.$pollo.'_collectRawMetrics.out';
			$ptask{command} = "$gatk CollectRawWgsMetrics -I $rdir/$pollo"."_recal.bam -O $rdir/$pollo"."_raw_wgs_metrics.txt -R $ref_fa\n";
			$ptask{dependency} = "afterok:$jid0";
			$jid1 = send2slurm(\%ptask);
			push @ljobids, $jid1;
		}
		# HaplotypeCaller, depende de ApplyBQSR
		$ptask{cpus} = 8;
		$ptask{job_name} = $pollo.'_haplotypeCaller';
		$ptask{filename} = $slurmdir.'/'.$pollo.'_haplotypeCaller.sh'; 
		$ptask{output} = $slurmdir.'/'.$pollo.'_haplotypeCaller.out'; 
		$ptask{command} = "$gatk HaplotypeCaller -R $ref_fa".(($mode eq 'wgs')?' ':" -L $unions ")."-I $rdir/$pollo"."_recal.bam -G StandardAnnotation -G AS_StandardAnnotation -ERC GVCF --dbsnp $ref_dir/$dbsnp -O $rdir/$pollo"."_raw.snps.indels.g.vcf.gz\n"; 
		$ptask{command}.= "$gatk VariantEval -R $ref_fa".(($mode eq 'wgs')?' ':" -L $unions ")."-D $ref_dir/$hcsnps -O $rdir/$pollo"."_eval.gatkreport --eval $rdir/$pollo"."_raw.snps.indels.g.vcf.gz\n"; 
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
