#!/usr/bin/perl 
#
# Copyleft 2024 O. Sotolongo <osotolongo@fundacioace.org>  
#
# This script is intended for making the joint call with the gVCF files
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
# Edit info here
#
my $src_dir = '/nas/Genomica/01-Data/02-WXS/02-Processed/202201_TS_EOAD-DEGESCO/02-gVCF/';
my $panel_dir = '/nas/Genomica/01-Data/00-Reference_files/02-GRCh38/202201_NGS_TS_EOAD-DEGESCO/panel/';my $baits = $panel_dir.'/EOAD_hg38_bait.interval_list';
my $targets = $panel_dir.'/EOAD_hg38_target.interval_list';
my $unions = $panel_dir.'/EOAD_hg38_union.interval_list';
my $outdir = '/nas/osotolongo/wes1/output';
my $search_pattern = '_raw.snps.indels.g.vcf.gz';
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
my $omni = '1000G_omni2.5.hg38.vcf.gz';
my $hapmap = 'hapmap_3.3.hg38.vcf.gz';
my @chrs = (1 .. 22, "X", "Y");
#
# Executable Paths
#
my $gatk = 'singularity run --cleanenv -B /nas:/nas -B /ruby:/ruby -B /the_dysk:/the_dysk /nas/usr/local/opt/gatk4.simg gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16G"';
# 
# 
# Get CLI inputs
#
#
my $cfile; 
my $debug = 0; 
my $test = 0; 
while (@ARGV and $ARGV[0] =~ /^-/) {         
	$_ = shift;         
	last if /^--$/;         
	if (/^-c/) { $cfile = shift; chomp($cfile);}         
	if (/^-o/) { $outdir = shift; chomp($outdir);}         
	if (/^-g/) { $debug = 1;}         
	if (/^-t/) { $test = 1;}         
	if (/^-s/) { $src_dir = shift; chomp($src_dir);} 
} 
my $workdir = getcwd; 
$outdir = $workdir.'/output' unless $outdir; 
mkdir $outdir unless -d $outdir; 
my $slurmdir = $outdir.'/slurm'; 
mkdir $slurmdir unless -d $slurmdir; 
# Do you want to process just a subset? Read the supplied list of subjects  
my @plist; 
if ($cfile and -f $cfile) {         
	open my $handle, "<$cfile";         
	chomp (@plist = <$handle>);         
	close $handle; 
}
my @content = find(file => 'name' => "*$search_pattern", in => $src_dir);
my %pollos = map {/.*\/(\w+?)$search_pattern$/; $1 => $_} @content;
my $hencoop = "$outdir/subjects.list";
open LDF, ">$hencoop" or die "Could not create the subject list\n";
foreach my $pollo (sort keys %pollos){
	print LDF "$pollo\t$pollos{$pollo}\n" if grep {/$pollo/} @plist;
}
close LDF;
my %ptask = (cpus => 16, time => '72:0:0', mem_per_cpu => '4G', debug => $test);
my @jobs;
my @tmp_joint;
my $resss_snp = '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 '.$ref_dir.'/'.$hapmap + ' --resource:omni,known=false,training=true,truth=false,prior=12.0 '.$ref_dir.'/'.$omni.' --resource:1000G,known=false,training=true,truth=false,prior=10.0 '.$ref_dir.'/'.$hcsnps.' --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 '.$ref_dir.'/'.$dbsnp;
my $resss_indel = '--resource:mills,known=false,training=true,truth=true,prior=12.0 '.$ref_dir.'/'.$known2.' --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 '.$ref_dir.'/'.$dbsnp;
my $troptions_snp = '-an QD -an ReadPosRankSum -an FS -an SOR -an MQ -an MQRankSum -mode SNP --max-gaussians 2 --trust-all-polymorphic -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.0 -tranche 90.0';
my $troptions_indel = '-an QD -an ReadPosRankSum -an FS -an SOR -an MQ -an MQRankSum -mode INDEL --max-gaussians 4 --trust-all-polymorphic -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.0 -tranche 90.0';
foreach my $chr (@chrs){
	$ptask{job_name} = 'GenoypeGVCFs_'.$chr;
	$ptask{filename} = $slurmdir.'/'.$chr.'_GenoypeGVCFs_'.$chr.'.sh';
	$ptask{output} = $slurmdir.'/'.$chr.'_GenoypeGVCFs_'.$chr.'.out';
	my $dbdir = "$outdir/wes.chr$chr.db";
	$ptask{command} = "$gatk  GenomicsDBImport --genomicsdb-workspace-path $dbdir --batch-size 50 --sample-name-map $hencoop  -L chr $chr  --reader-threads 8\n";
	$ptask{command}.= "$gatk  GenotypeGVCFs -R $ref_fa  -V gendb://$dbdir -G StandardAnnotation -G AS_StandardAnnotation -O $outdir/wes.chr$chr.snps.indels.g.vcf.gz\n";
	my $jid = send2slurm(\%ptask);
	push @jobs, $jid;
	push @tmp_joint, "$outdir/wes.chr$chr.snps.indels.g.vcf.gz";
}
my %gtask = (cpus => 12, time => '72:0:0', mem_per_cpu => '4G', debug => $test, job_name => 'gather', filename => "$slurmdir/gather.sh", 'output' => "$slurmdir/gather.out", dependency => 'afterok:'.join(',afterok:', @jobs));
my $ppool = join(' -I ', @tmp_joint);
$gtask{command} = "$gatk GatherVcfs -I $ppool -O $outdir/wes_joint_chr_norec.vcf.gz\n";
$gtask{command}.= "$gatk IndexFeatureFile -I $outdir/wes_joint_chr_norec.vcf.gz\n";
$gtask{command}.= "$gatk VariantRecalibrator -AS -R $ref_fa -V $outdir/wes_joint_chr_norec.vcf.gz $resss_snp $troptions_snp -O $outdir/wes_joint_chr.snps.recal  --tranches-file $outdir/wes_joint_chr.snps.recalibrate.tranches --rscript-file $outdir/wes_joint_chr.snps.recalibrate.plots.R\n";
$gtask{command}.= "$gatk VariantRecalibrator -AS -R $ref_fa -V $outdir/wes_joint_chr_norec.vcf.gz $resss_indel $troptions_indel -O $outdir/wes_joint_chr.indels.recal --tranches-file $outdir/wes_joint_chr.indels.recalibrate.tranches --rscript-file $outdir/wes_joint_chr.indels.recalibrate.plots.R\n";
$gtask{command}.= "$gatk  ApplyVQSR -R $ref_fa -V $outdir/wes_joint_chr_norec.vcf.gz  -mode SNP --truth-sensitivity-filter-level 99.7 --recal-file $outdir/wes_joint_chr.snps.recal  --tranches-file $outdir/wes_joint_chr.snps.recalibrate.tranches -O $outdir/wes_joint_chr.snps.g_recalibrated.vcf.gz\n";
$gtask{command}.= "$gatk  ApplyVQSR -R $ref_fa -V $outdir/wes_joint_chr.snps.g_recalibrated.vcf.gz -mode INDEL --truth-sensitivity-filter-level 99.7 --recal-file $outdir/wes_joint_chr.indels.recal --tranches-file $outdir/wes_joint_chr.indels.recalibrate.tranches -O $outdir/wes_joint_chr.snps.indels.g_recalibrated.vcf.gz\n";
$gtask{mailtype} = 'FAIL,TIME_LIMIT,STAGE_OUT,END';
send2slurm(\%gtask);
