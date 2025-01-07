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
my $ref_fa = $dpaths{ref_dir}.'/'.$dpaths{ref_name}.'.fasta';
my $tmp_shit = $ENV{TMPDIR};
my @chrs = (1 .. 22, "X", "Y");
#
# Executable Paths
#
my %epaths = exec_paths();
# 
# 
# Get CLI inputs
#
#
my $cfile; 
my $debug = 0; 
my $test = 0; 
my $init; 
my $mode = 'wes';
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
my @content = find(file => 'name' => "*$wesconf{search_pattern}", in => $wesconf{src_dir});
my %pollos = map {/.*\/(\w+?)$wesconf{search_pattern}$/; $1 => $_} @content;
my $hencoop = "$wesconf{outdir}/subjects.list";
open LDF, ">$hencoop" or die "Could not create the subject list\n";
foreach my $pollo (sort keys %pollos){
	print LDF "$pollo\t$pollos{$pollo}\n" if not @plist or grep {/$pollo/} @plist;
}
close LDF;
my %ptask = (cpus => 16, time => '72:0:0', mem_per_cpu => '4G', debug => $test);
my @jobs;
my @tmp_joint;
my $resss_snp = '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 '.$dpaths{ref_dir}.'/'.$dpaths{hapmap}.' --resource:omni,known=false,training=true,truth=false,prior=12.0 '.$dpaths{ref_dir}.'/'.$dpaths{omni}.' --resource:1000G,known=false,training=true,truth=false,prior=10.0 '.$dpaths{ref_dir}.'/'.$dpaths{hcsnps}.' --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 '.$dpaths{ref_dir}.'/'.$dpaths{dbsnp};
my $resss_indel = '--resource:mills,known=false,training=true,truth=true,prior=12.0 '.$dpaths{ref_dir}.'/'.$dpaths{known2}.' --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 '.$dpaths{ref_dir}.'/'.$dpaths{dbsnp};
my $troptions_snp = '-an QD -an ReadPosRankSum -an FS -an SOR -an MQ -an MQRankSum -mode SNP --max-gaussians 2 --trust-all-polymorphic -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.0 -tranche 90.0';
my $troptions_indel = '-an QD -an ReadPosRankSum -an FS -an SOR -an MQ -an MQRankSum -mode INDEL --max-gaussians 4 --trust-all-polymorphic -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.0 -tranche 90.0';
foreach my $chr (@chrs){
	$ptask{job_name} = 'GenoypeGVCFs_'.$chr;
	$ptask{filename} = $slurmdir.'/'.$chr.'_GenoypeGVCFs_'.$chr.'.sh';
	$ptask{output} = $slurmdir.'/'.$chr.'_GenoypeGVCFs_'.$chr.'.out';
	my $dbdir = "$wesconf{outdir}/wes.chr$chr.db";
	$ptask{command} = "$epaths{gatk} GenomicsDBImport --genomicsdb-workspace-path $dbdir --batch-size 50 --sample-name-map $hencoop  -L chr$chr  --reader-threads 8\n";
	$ptask{command}.= "$epaths{gatk}  GenotypeGVCFs -R $ref_fa  -V gendb://$dbdir -G StandardAnnotation -G AS_StandardAnnotation -O $wesconf{outdir}/wes.chr$chr.snps.indels.g.vcf.gz\n";
	my $jid = send2slurm(\%ptask);
	push @jobs, $jid;
	push @tmp_joint, "$wesconf{outdir}/wes.chr$chr.snps.indels.g.vcf.gz";
}
my %gtask = (cpus => 12, time => '72:0:0', mem_per_cpu => '4G', debug => $test, job_name => 'gather', filename => "$slurmdir/gather.sh", 'output' => "$slurmdir/gather.out", dependency => 'afterok:'.join(',afterok:', @jobs));
my $ppool = join(' -I ', @tmp_joint);
$gtask{command} = "$epaths{gatk} GatherVcfs -I $ppool -O $wesconf{outdir}/wes_joint_chr_norec.vcf.gz\n";
$gtask{command}.= "$epaths{gatk} IndexFeatureFile -I $wesconf{outdir}/wes_joint_chr_norec.vcf.gz\n";
unless ($mode eq 'wgs'){
	$gtask{command}.= "$epaths{gatk} VariantRecalibrator -AS -R $ref_fa -V $wesconf{outdir}/wes_joint_chr_norec.vcf.gz $resss_snp $troptions_snp -O $wesconf{outdir}/wes_joint_chr.snps.recal  --tranches-file $wesconf{outdir}/wes_joint_chr.snps.recalibrate.tranches --rscript-file $wesconf{outdir}/wes_joint_chr.snps.recalibrate.plots.R\n";
	$gtask{command}.= "$epaths{gatk} VariantRecalibrator -AS -R $ref_fa -V $wesconf{outdir}/wes_joint_chr_norec.vcf.gz $resss_indel $troptions_indel -O $wesconf{outdir}/wes_joint_chr.indels.recal --tranches-file $wesconf{outdir}/wes_joint_chr.indels.recalibrate.tranches --rscript-file $wesconf{outdir}/wes_joint_chr.indels.recalibrate.plots.R\n";
	$gtask{command}.= "$epaths{gatk}  ApplyVQSR -R $ref_fa -V $wesconf{outdir}/wes_joint_chr_norec.vcf.gz  -mode SNP --truth-sensitivity-filter-level 99.7 --recal-file $wesconf{outdir}/wes_joint_chr.snps.recal  --tranches-file $wesconf{outdir}/wes_joint_chr.snps.recalibrate.tranches -O $wesconf{outdir}/wes_joint_chr.snps.g_recalibrated.vcf.gz\n";
	$gtask{command}.= "$epaths{gatk}  ApplyVQSR -R $ref_fa -V $wesconf{outdir}/wes_joint_chr.snps.g_recalibrated.vcf.gz -mode INDEL --truth-sensitivity-filter-level 99.7 --recal-file $wesconf{outdir}/wes_joint_chr.indels.recal --tranches-file $wesconf{outdir}/wes_joint_chr.indels.recalibrate.tranches -O $wesconf{outdir}/wes_joint_chr.snps.indels.g_recalibrated.vcf.gz\n";
}
$gtask{mailtype} = 'FAIL,TIME_LIMIT,STAGE_OUT,END';
send2slurm(\%gtask);
