#!/usr/bin/perl
#
# Copyleft 2024 O. Sotolongo <osotolongo@fundacioace.org> 
#
# This script is intended for manipulate the compressed cram (or bam) from b37 to hg38 
#
use strict;
use warnings;
use SLURMACE;
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
#
# Executable Paths
#
my %epaths = exec_paths(); 
#
# Get CLI inputs
#
my $cfile;
my $debug = 0;
my $test = 0;
my $init; 
while (@ARGV and $ARGV[0] =~ /^-/) {
	$_ = shift;
	last if /^--$/;
	if (/^-c/) { $cfile = shift; chomp($cfile);}
	if (/^-i/) { $init = shift; chomp($init);}
	if (/^-g/) { $debug = 1;}
	if (/^-t/) { $test = 1;}
}
my %wesconf = init_conf($init);
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
# Now read the subdirs in the source dir. I will assume each one counts as a different subject
# This can change, so I should the logic according to the organization of the samples
opendir my $dh, $wesconf{src_dir} or die "Could not open directory: $!";
my @mlist = grep {-d "$wesconf{src_dir}/$_" && ! /^\.{1,2}$/} readdir ($dh); 
# Now, a have the list of the samples. If there is any supplied list of subjkects I am going 
# to choose only those from the existing ones
if ($cfile) {
      @mlist = grep {/\S/} map { $a = $_; (grep /^$a$/, @mlist)?$a:''; } @plist;
}
# Lets process now
my %cdata = (cpus => 4, time => '24:0:0', mem_per_cpu => '4G', debug => $test);
my @jobs;
foreach my $pollo (@mlist){
	my $tmpdir = $tmp_shit.'/'.$pollo;
	my @inqs = split "\n", qx"$epaths{samtools} view -H $wesconf{src_dir}/$pollo/$pollo.bam";
	my %params;
	foreach my $inq (@inqs){
	       if( $inq	=~ /^\@RG.*/){
		       %params = ($inq =~ /^\@RG\t(\w+):(\w+)\t(\w+):(\w+)\t(\w+):(\w+)\t(\w+):(\w+)$/);
	       }
	}
	$cdata{job_name} = $pollo.'_RevertSam';
	$cdata{filename} = $slurmdir.'/'.$pollo.'_RevertSam.sh';
	$cdata{output} =  $slurmdir.'/'.$pollo.'_RevertSam.out';
	$cdata{command} = "mkdir -p $tmpdir\n";
	$cdata{command}.= "cp $wesconf{src_dir}/$pollo/$pollo.b* $tmpdir/\n";
	$cdata{command}.= "$epaths{gatk} RevertSam -I $tmpdir/$pollo.bam -O $tmpdir/$pollo"."_u.bam -R $ref_fa\n";
	$cdata{command}.= "$epaths{gatk} SortSam -I $tmpdir/$pollo"."_u.bam -O $tmpdir/$pollo"."_us.bam -SORT_ORDER queryname --TMP_DIR $tmp_shit\n";
	$cdata{command}.= "$epaths{gatk} SamToFastq -I $tmpdir/$pollo.bam -FASTQ $tmpdir/$pollo"."_a.fastq -F2 $tmpdir/$pollo"."_b.fastq -R $ref_fa\n";
	$cdata{command}.= "$epaths{bwa} -R \"\@RG\\tID:$params{ID}\\tPL:$params{PL}\\tLB:$params{LB}\\tSM:$params{SM}\" $ref_fa $tmpdir/$pollo"."_a.fastq $tmpdir/$pollo"."_b.fastq > $tmpdir/$pollo"."_alignment.sam\n";
	$cdata{command}.= "$epaths{gatk} SortSam -I $tmpdir/$pollo"."_alignment.sam -O $tmpdir/$pollo"."_salignment.sam  -SORT_ORDER queryname --TMP_DIR $tmp_shit\n";
	$cdata{command}.= "$epaths{gatk} MergeBamAlignment -ALIGNED $tmpdir/$pollo"."_salignment.sam -UNMAPPED $tmpdir/$pollo"."_us.bam -R $ref_fa -O $tmpdir/$pollo"."_merged.bam\n";
	$cdata{command}.= "$epaths{gatk} AddOrReplaceReadGroups -I  $tmpdir/$pollo"."_merged.bam -O $wesconf{outdir}/$pollo.sam -ID $params{ID} -PL $params{PL} -LB $params{LB} -PU $params{ID} -SM $params{SM}\n";
	$cdata{command}.= "rm -rf $tmpdir\n" unless $debug;
	my $jid = send2slurm(\%cdata);
	push @jobs, $jid;
}
unless ($test) {
	my %wtask = (cpus => 1, job_name => 'end_wes', filename => $slurmdir.'/end.sh', output => $slurmdir.'/end.out', dependency => 'afterok:'.join(',afterok:',@jobs));
	send2slurm(\%wtask);
}

