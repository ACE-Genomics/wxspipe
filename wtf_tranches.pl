#!/usr/bin/perl
use strict;
use warnings;
use SLURMACE qw(send2slurm wait4jobs);
use Cwd;
use FindBin; 
use lib "$FindBin::Bin";
use wxsInit;
my $init;
while (@ARGV and $ARGV[0] =~ /^-/) {
	$_ = shift;
	last if /^--$/;
	if (/^-i/) { $init = shift; chomp($init);}
}
die "Should supply init data file\n" unless $init;
my %wesconf = init_conf($init);
my $base_dir = $wesconf{src_dir};
die "Should supply source data directory in init file\n" unless $base_dir and -d $base_dir;
my $wdir =  $wesconf{outdir} || cwd();
my %epaths = exec_paths();
my $GATK = $epaths{gatk};
my %dpaths = data_paths();
my $ref_dir = $dpaths{ref_dir};
my $ref_name = $dpaths{ref_name};
my $ref = $ref_dir.'/'.$ref_name.'.fasta';
my @tranches = (100.0, 99.95, 99.9, 99.5, 99.0, 97.0, 96.0, 95.0, 94.0, 93.5, 93.0, 92.0, 91.0, 90.0);
my $sdir = $wdir.'/slurm';
my $tables_dir = $wdir.'/tables';
my $ofile = $tables_dir.'/wes_joint_alltranches.table';
mkdir $sdir unless -d $sdir;
mkdir $tables_dir unless -d $tables_dir;
my @jobs;
my %trjob = ('cpus' => 4, 'time' => '8:0:0', 'mem_per_cpu' => '4G', 'job_name' => 'fucktranches');
foreach my $tranche (@tranches) {
	(my $str_tranche = $tranche) =~ s/\.//;
	$trjob{'filename'} = $sdir.'/'.$str_tranche.'.sh';
	$trjob{'output'} = $sdir.'/'.$str_tranche.'.out';
	$trjob{'command'} = "$GATK ApplyVQSR -R $ref -V $base_dir/wes_joint_chr_norec.vcf.gz -mode SNP --truth-sensitivity-filter-level $tranche --recal-file $base_dir/wes_joint_chr.snps.recal --tranches-file $base_dir/wes_joint_chr.snps.recalibrate.tranches -O $tables_dir/wes_joint_snps.$str_tranche.vcf.gz\n";
	$trjob{'command'} .= "$GATK VariantsToTable  -R $ref -V $tables_dir/wes_joint_snps.$str_tranche.vcf.gz -F CHROM -F POS -F DP -F FILTER -F MQ -F QD -F FS -F SOR -F MQRankSum -F ReadPosRankSum -O $tables_dir/wes_joint.$str_tranche.table";
	my $jobid = send2slurm(\%trjob);
	push @jobs, $jobid;
}
wait4jobs(@jobs);
open ODF, ">$ofile";
my $written = 0;
foreach my $tranche (@tranches) {
	(my $str_tranche = $tranche) =~ s/\.//;
	my $isfirst = 1;
	open IDF, "$tables_dir/wes_joint.$str_tranche.table";
	while (<IDF>) {
		chomp;
		if($isfirst) {
			print ODF "$_\tTRANCHE\n" unless $written;
			$isfirst = 0;
			$written = 1;
		}else{
			print ODF "$_\t$tranche\n";
		}
	}
	close IDF;
}
close ODF;
print "\n\nShit done!!!!!\n\n";

#my %warn = ('job_name' => 'fucktranches', 'filename' => $sdir.'/end.sh', 'mailtype' => 'END', 'dependency' => 'singleton', 'output' =>  $sdir.'/end.out');
#send2slurm(\%warn);
