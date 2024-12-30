#!/usr/bin/perl
use strict;
use warnings;
use File::Temp qw(:mktemp);
use Spreadsheet::Write;
use Statistics::Descriptive;
use Data::Dump qw(dump);
my %wesconf;
my $init;
while (@ARGV and $ARGV[0] =~ /^-/) {
	$_ = shift;
	last if /^--$/;
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
my $ipath = $wesconf{src_dir}.'/';
die "Should supply source data directory in init file\n" unless $ipath and -d $ipath;
############ Paths #######################################
my $wgs_suffix = '_wgs_metrics.txt';
my $raw_suffix = '_raw_wgs_metrics.txt';
my $padded_suffix = '_padded_wgs_metrics.txt';
my $eval_suffix = '_eval.gatkreport';
my $freemix_suffix = '.vbid2.selfSM';
############ Busca las muestras  ######################### 
my @idirs = glob( $ipath.'*' ); 
my @pollos; 
foreach my $idir (@idirs) {         
	if ( -d $idir and -d $idir.'/results') {                 
		my ($pollo) = $idir =~ /$ipath\/*(.*)$/;                 
		push @pollos, $pollo;         
	} 
} 
my %evals;
################ Lee los reports #####################################
foreach my $pollo (@pollos) {
	my $ifile = $ipath.$pollo.'/results/'.$pollo.$raw_suffix.'.sample_summary';
	open IDF, "<$ifile";
	while (<IDF>){
		#sample_id,total,mean,granular_third_quartile,granular_median,granular_first_quartile,%_bases_above_10,%_bases_above_15,%_bases_above_20,%_bases_above_30,%_bases_above_40,%_bases_above_50,%_bases_above_60,%_bases_above_70,%_bases_above_80,%_bases_above_90,%_bases_above_100
		#22D28227621,7340059881,167.12,203,167,130,99.1,99.0,98.8,98.5,98.1,97.8,97.3,96.4,94.9,92.7,89.6
		if (/^$pollo.*/){
			my ($coverage, $pct10, $pct20, $pct30, $pct40, $pct50, $pct60, $pct70, $pct80, $pct90, $pct100) = /^$pollo,\d+,(\d+\.\d+),\d+,\d+,\d+,(\d+\.\d+),\d+\.\d+,(\d+\.\d+),(\d+\.\d+),(\d+\.\d+),(\d+\.\d+),(\d+\.\d+),(\d+\.\d+),(\d+\.\d+),(\d+\.\d+),(\d+\.\d+)$/;
			$evals{$pollo}{'RawCoverage'} = $coverage;
			$evals{$pollo}{'PCT_10x'} = $pct10;
			$evals{$pollo}{'PCT_20x'} = $pct20;
			$evals{$pollo}{'PCT_30x'} = $pct30;
			$evals{$pollo}{'PCT_40x'} = $pct40;
			$evals{$pollo}{'PCT_50x'} = $pct50;
			$evals{$pollo}{'PCT_60x'} = $pct60;
			$evals{$pollo}{'PCT_70x'} = $pct70;
			$evals{$pollo}{'PCT_80x'} = $pct80;
			$evals{$pollo}{'PCT_90x'} = $pct90;
			$evals{$pollo}{'PCT_100x'} = $pct100;
		}
	}
	close IDF;
	$ifile = $ipath.$pollo.'/results/'.$pollo.$wgs_suffix.'.sample_summary';
	open IDF, "<$ifile";
	while (<IDF>){
		if (/^$pollo.*/){
			my ($rcoverage) = /^$pollo,\d+,(\d+\.\d+),.*/;
		 	$evals{$pollo}{'MeanCoverage'} = $rcoverage;
		}
	}
	close IDF;

	$ifile = $ipath.$pollo.'/results/'.$pollo.$padded_suffix.'.sample_summary';
	open IDF, "<$ifile";
	while (<IDF>){
		if (/^$pollo.*/){
			my ($rcoverage) = /^$pollo,\d+,(\d+\.\d+),.*/;
		 	$evals{$pollo}{'PaddedCoverage'} = $rcoverage;
		}
	}
	close IDF;

	$ifile  = $ipath.$pollo.'/results/'.$pollo.$eval_suffix;
	open IDF, "<$ifile";
	while (<IDF>){
		if (/^CountVariants\s+dbsnp/){
		#CountVariants  dbsnp             eval              none            all            43921240       455666    423084         32582   0.00074183     1348.00000000  29298     10          862        1222
			my ($Novelty, $nSNPs, $nInsertions, $nDeletions) = /CountVariants\s+dbsnp\s+eval\s+none\s+(\w+)\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+)\s+\d+\s+(\d+)\s+(\d+)\s+.*/;
			$evals{$pollo}{'nSNPs'}{$Novelty} = $nSNPs;
			$evals{$pollo}{'nInsertions'}{$Novelty} = $nInsertions;
			$evals{$pollo}{'nDeletions'}{$Novelty} = $nDeletions;
		}
		#MetricsCollection  CompFeatureInput  EvalFeatureInput  JexlExpression  Novelty  concordantRate  nSNPs  nSNPloci  nIndels  nIndelLoci  indelRatio  indelRatioLociBased  tiTvRatio 
		#MetricsCollection  dbsnp             eval              none            all               99.94  29322     29298     3499        3259        0.73                 0.71       2.76 
		if (/^MetricsCollection\s+dbsnp/){
			my ($Novelty, $tiTvRatio) = /MetricsCollection\s+dbsnp\s+eval\s+none\s+(\w+)\s+\d+\.\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(?:\d+\.\d+|NA)\s+\d+\.\d+\s+(\d+\.\d+)/;

			$evals{$pollo}{'tiTvRatio'}{$Novelty} = $tiTvRatio;
		}

	}
	close IDF;
	$ifile = $ipath.$pollo.'/results/'.$pollo.$freemix_suffix;
	open IDF, "<$ifile";
	while (<IDF>){
		if (/^\w+/){
			#SEQ_ID RG CHIP_ID #SNPS #READS AVG_DP FREEMIX FREELK1 FREELK0 FREE_RH FREE_RA CHIPMIX CHIPLK1 CHIPLK0 CHIP_RH CHIP_RA DPREF RDPHET RDPALT 
			#22D28227616 NA NA 9892 1120289 117.913 0.019543 -206006 -217696 NA NA NA NA NA NA NA NA NA NA 
			my ($freemix) = /^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+.*/;
			$evals{$pollo}{'Freemix'} = $freemix;
		}
	}
	close IDF;
}
##################################### Escribe archivos finales en varios formatos ###################################
my $ofile0 = 'report_all.csv';
open ODF, ">$ofile0";
print ODF "Sample,RawCoverage,MeanCoverage,PaddedCoverage,PCT_10x,PCT_20x,PCT_30x,PCT_40x,PCT_50x,PCT_60x,PCT_70x,PCT_80x,PCT_90x,PCT_100x,dbSNP_nSNPs_all,dbSNP_nSNPs_known,dbSNP_nSNPs_novel,dbSNP_nInsertions_all,dbSNP_nInsertions_known,dbSNP_nInsertions_novel,dbSNP_nDeletions_all,dbSNP_nDeletions_known,dbSNP_nDeletions_novel,tiTvRatio_all,tiTvRatio_known,tiTvRatio_novel,Freemix\n";
foreach my $pollo (@pollos) {
	print ODF "$pollo,$evals{$pollo}{'RawCoverage'},";
	print ODF "$evals{$pollo}{'MeanCoverage'},";
	print ODF "$evals{$pollo}{'PaddedCoverage'},";
	print ODF "$evals{$pollo}{'PCT_10x'},";
	print ODF "$evals{$pollo}{'PCT_20x'},";
	print ODF "$evals{$pollo}{'PCT_30x'},";
	print ODF "$evals{$pollo}{'PCT_40x'},";
	print ODF "$evals{$pollo}{'PCT_50x'},";
	print ODF "$evals{$pollo}{'PCT_60x'},";
	print ODF "$evals{$pollo}{'PCT_70x'},";
	print ODF "$evals{$pollo}{'PCT_80x'},";
	print ODF "$evals{$pollo}{'PCT_90x'},";
	print ODF "$evals{$pollo}{'PCT_100x'},";
	foreach my $dbsnp ('nSNPs', 'nInsertions', 'nDeletions', 'tiTvRatio') {
		foreach my $isnovelty ('all', 'known', 'novel') {
			print ODF "$evals{$pollo}{$dbsnp}{$isnovelty},";
		}
	}
	print ODF "$evals{$pollo}{'Freemix'}\n";
}
close ODF;

my $ofile = mktemp($ENV{TMPDIR}.'/dbsnp.XXXXXXX');

open ODF, ">$ofile";
print ODF "Sample,nSNPs,nInsertions,nDeletions,tiTvRatio,Novelty\n";
foreach my $pollo (@pollos) {
	foreach my $isnovelty ('all', 'known', 'novel') {
		print ODF "$pollo,";
		foreach my $dbsnp ('nSNPs', 'nInsertions', 'nDeletions', 'tiTvRatio') {
			 print ODF "$evals{$pollo}{$dbsnp}{$isnovelty},";
		}
		print ODF "$isnovelty\n"; 
	 }
}
close ODF;

my $ofile1 = 'report_all.xlsx';
my $workbook = Spreadsheet::Write->new(file => $ofile1, sheet => 'DATA');
$workbook->addrow(['Sample','RawCoverage','MeanCoverage','PaddedCoverage','PCT_10x','PCT_20x','PCT_30x','PCT_40x','PCT_50x','PCT_60x','PCT_70x','PCT_80x','PCT_90x','PCT_100x','dbSNP_nSNPs_all','dbSNP_nSNPs_known','dbSNP_nSNPs_novel','dbSNP_nInsertions_all','dbSNP_nInsertions_known','dbSNP_nInsertions_novel','dbSNP_nDeletions_all','dbSNP_nDeletions_known','dbSNP_nDeletions_novel','tiTvRatio_all','tiTvRatio_known','tiTvRatio_novel','Freemix']);
foreach my $pollo (@pollos) {
	my @datarow;
	push @datarow, $pollo;
	push @datarow, $evals{$pollo}{'RawCoverage'};
	push @datarow, $evals{$pollo}{'MeanCoverage'};
	push @datarow, $evals{$pollo}{'PaddedCoverage'};
	push @datarow, $evals{$pollo}{'PCT_10x'};
	push @datarow, $evals{$pollo}{'PCT_20x'};
	push @datarow, $evals{$pollo}{'PCT_30x'};
	push @datarow, $evals{$pollo}{'PCT_40x'};
	push @datarow, $evals{$pollo}{'PCT_50x'};
	push @datarow, $evals{$pollo}{'PCT_60x'};
	push @datarow, $evals{$pollo}{'PCT_70x'};
	push @datarow, $evals{$pollo}{'PCT_80x'};
	push @datarow, $evals{$pollo}{'PCT_90x'};
	push @datarow, $evals{$pollo}{'PCT_100x'};
	foreach my $dbsnp ('nSNPs', 'nInsertions', 'nDeletions', 'tiTvRatio') {
		foreach my $isnovelty ('all', 'known', 'novel') {
			push @datarow, $evals{$pollo}{$dbsnp}{$isnovelty};
		}
	}
	push @datarow, $evals{$pollo}{'Freemix'};
	$workbook->addrow(\@datarow);
}
$workbook->close();

my $cfile = mktemp($ENV{TMPDIR}.'/coverage.XXXXXXX');
open ODF, ">$cfile";
print ODF "Sample,Type,Coverage\n";
foreach my $pollo (@pollos) {
	foreach my $kind ('RawCoverage', 'MeanCoverage', 'PaddedCoverage'){
		print ODF "$pollo,$kind,$evals{$pollo}{$kind}\n";
	}
}
close ODF;
my $rplt = mktemp($ENV{TMPDIR}.'/plotter.XXXXXXX');
open RSF, ">$rplt";
	print RSF "library(ggplot2)\n";
	print RSF "dbsnp <- read.csv(\"$ofile\")\n";
	print RSF "ggplot(dbsnp, aes(y=tiTvRatio, fill=Novelty)) + geom_boxplot() + theme(axis.title.y=element_blank(), legend.position = \"bottom\") + ggtitle(\"TiTvRatios\")-> p\n";
	print RSF "ggsave(p, filename = \"titvratios.png\", device = png, type = \"cairo\", dpi=300, width = 4, height = 3, units = \"in\")\n";
	print RSF "qc <- read.csv(\"$ofile0\")\n";
	print RSF "ggplot(qc, aes(y=Freemix)) + geom_boxplot(color=\"blue\", outlier.color=\"navy\") +  theme(axis.title.y=element_blank())  + ggtitle(\"Freemix\")-> p\n";
	print RSF "ggsave(p, filename = \"freemix.png\", device = png, type = \"cairo\", dpi=300, width = 4, height = 3, units = \"in\")\n";
	print RSF "cvg <- read.csv(\"$cfile\")\n";
	print RSF "ggplot(cvg, aes(y=Coverage, fill=Type)) + geom_boxplot() +  theme(axis.title.y=element_blank(), legend.position = \"bottom\")  + ggtitle(\"Coverage\") -> p\n";
	print RSF "ggsave(p, filename = \"coverage.png\", device = png, type = \"cairo\", dpi=300, width = 4, height = 3, units = \"in\")\n";
close RSF;
system("Rscript $rplt");
unlink $rplt;
system("montage titvratios.png freemix.png coverage.png -tile 3x1 -geometry +0+0 report.png");
unlink 'titvratios.png', 'freemix.png', 'coverage.png';
#my $rmark = mktemp($ENV{TMPDIR}.'/rmd.XXXXXXX');
my $rmark = 'qc_report.rmd';
open RSF, ">$rmark";
	print RSF "---\n";
	print RSF "title: QC Report\n";
	print RSF "author: O. Sotolongo-Grau\n";
	print RSF "output: word_document\n";
	print RSF "---\n";
	print RSF "# QC summary statistics\n";
        print RSF "| | Mean | STDDev | Min | Max | Thresholds | Failed |\n";
	print RSF "|:---|---:|---:|---:|---:|---:|---:|\n";
	my @rc; my @frm; my @titv;
	my $frc = 0; my $ffrm = 0; my $ftitv = 0;
	foreach my $pollo (@pollos) {
		push @rc, $evals{$pollo}{'RawCoverage'};
		$frc++ if ($evals{$pollo}{'RawCoverage'} < 100.0);
		push @frm, $evals{$pollo}{'Freemix'};
		$ffrm++ if ($evals{$pollo}{'Freemix'} gt 0.05);
		push @titv, $evals{$pollo}{'tiTvRatio'}{'all'};
		$ftitv++ if (($evals{$pollo}{'tiTvRatio'}{'all'} lt 3.0) or ($evals{$pollo}{'tiTvRatio'}{'all'} gt 3.3));
	}
	my $src = Statistics::Descriptive::Full->new();
	$src->add_data(@rc);	
	print RSF "|Raw Coverage|".sprintf("%.2f",$src->mean())."|".sprintf("%.2f",$src->standard_deviation())."|".sprintf("%.2f",$src->min())."|".sprintf("%.2f",$src->max())."|<100|".sprintf("%d",$frc)."|\n";
	my $srm = Statistics::Descriptive::Full->new();
	$srm->add_data(@frm);
	print RSF "|Freemix|".sprintf("%.3f",$srm->mean())."|".sprintf("%.3f",$srm->standard_deviation())."|".sprintf("%.3f",$srm->min())."|".sprintf("%.3f",$srm->max())."|>0.05|".sprintf("%d",$ffrm)."|\n";
	my $stt = Statistics::Descriptive::Full->new();
	$stt->add_data(@titv);
	print RSF "|dbSNP TiTv Ratio|".sprintf("%.2f",$stt->mean())."|".sprintf("%.2f",$stt->standard_deviation())."|".sprintf("%.2f",$stt->min())."|".sprintf("%.2f",$stt->max())."|<3.0 Or >3.3|".sprintf("%d",$ftitv)."|\n";	
	print RSF "## Descriptive Plots\n";
	print RSF '```{r warning = FALSE}'."\n";
	print RSF "library(ggplot2)\n";
	print RSF "dbsnp <- read.csv(\"$ofile\")\n";
	print RSF "ggplot(dbsnp, aes(y=tiTvRatio, fill=Novelty)) + geom_boxplot() + theme(axis.title.y=element_blank(), legend.position = \"bottom\") + ggtitle(\"TiTvRatios\")-> p\n";
	print RSF "print(p)\n\n";
	print RSF "qc <- read.csv(\"$ofile0\")\n";
	print RSF "ggplot(qc, aes(y=Freemix)) + geom_boxplot(color=\"blue\", outlier.color=\"navy\") +  theme(axis.title.y=element_blank())  + ggtitle(\"Freemix\")-> p\n";
	print RSF "print(p)\n\n";
	print RSF "cvg <- read.csv(\"$cfile\")\n";
	print RSF "ggplot(cvg, aes(y=Coverage, fill=Type)) + geom_boxplot() +  theme(axis.title.y=element_blank(), legend.position = \"bottom\")  + ggtitle(\"Coverage\") -> p\n";
	print RSF "print(p)\n\n";
	print RSF '```'."\n\n";
        print RSF "# Samples showing signs of contamination\n\n";
	print RSF "| Sample | RawCoverage | Freemix | TiTvRatio | PCT_20x | PCT_30x | Comments |\n";
	print RSF "|:---|---:|---:|---:|---:|---:|:---|\n";
	foreach my $pollo (@pollos) {
		if ($evals{$pollo}{'Freemix'} gt 0.05){
			print RSF "|$pollo|".sprintf("%.2f",$evals{$pollo}{'RawCoverage'})."|".sprintf("%.2f",$evals{$pollo}{'Freemix'})."|".sprintf("%.2f",$evals{$pollo}{'tiTvRatio'}{'all'})."|".sprintf("%.3f",$evals{$pollo}{'PCT_20x'})."|".sprintf("%.3f",$evals{$pollo}{'PCT_30x'})."| |\n";
		}
	}
	print RSF "# Samples with low coverage\n\n";
	print RSF "| Sample | RawCoverage | Freemix | TiTvRatio | PCT_20x | PCT_30x | Comments |\n";
	print RSF "|:---|---:|---:|---:|---:|---:|:---|\n";
	foreach my $pollo (@pollos) {
		if ($evals{$pollo}{'RawCoverage'} < 100.0){
			 print RSF "|$pollo|".sprintf("%.2f",$evals{$pollo}{'RawCoverage'})."|".sprintf("%.2f",$evals{$pollo}{'Freemix'})."|".sprintf("%.2f",$evals{$pollo}{'tiTvRatio'}{'all'})."|".sprintf("%.3f",$evals{$pollo}{'PCT_20x'})."|".sprintf("%.3f",$evals{$pollo}{'PCT_30x'})."| |\n";
		 }
	 }
close RSF;

$rplt = mktemp($ENV{TMPDIR}.'/render.XXXXXXX');
open RSF, ">$rplt";
	print RSF "library(\"rmarkdown\")\n";
	print RSF "render(\"$rmark\", output_file = \"report_tmp.docx\", output_dir = getwd())\n";
close RSF;
system("Rscript $rplt");
unlink $rplt;
#####################################################################################
# freaking hack
# R produce un docx correcto pero que no se puede editar en el office365
# Asi que hacemos que libreoffice lo convierta en un nuevo docx
# que es identico pero que ahora si se puede abrir sin problemas con el 365
# gracias Microsoft por la compatibilidad entre tus propios formatos!
#####################################################################################
system("unoconv -d document -f docx -o report_data report_tmp.docx");
unlink "report_tmp.docx";
my $hff = "contamination.csv";
open ODF, ">$hff";
print ODF "Sample,RawCoverage,Freemix,TiTvRatio,PCT_20x,PCT_30x,Comments\n";
foreach my $pollo (@pollos) {
	if ($evals{$pollo}{'Freemix'} gt 0.05){
		print ODF "$pollo,$evals{$pollo}{'RawCoverage'},$evals{$pollo}{'Freemix'},$evals{$pollo}{'tiTvRatio'}{'all'},$evals{$pollo}{'PCT_20x'},$evals{$pollo}{'PCT_30x'},\n";
	}
}
close ODF;
$ofile = "contamination.xlsx";
$workbook = Spreadsheet::Write->new(file => $ofile, sheet => 'Samples contamination');
$workbook->addrow(['Sample','RawCoverage','Freemix','TiTvRatio','PCT_20x','PCT_30x','Comments']);
foreach my $pollo (@pollos) {
		$workbook->addrow([$pollo,$evals{$pollo}{'RawCoverage'},$evals{$pollo}{'Freemix'},$evals{$pollo}{'tiTvRatio'}{'all'},$evals{$pollo}{'PCT_20x'},$evals{$pollo}{'PCT_30x'},'']);
}
$workbook->close();
