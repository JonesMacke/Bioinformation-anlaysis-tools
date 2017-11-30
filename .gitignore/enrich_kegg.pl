#!/usr/bin/perl -w
use strict;
use Text::NSP::Measures::2D::Fisher::right;
use Statistics::ChisqIndep;
use Excel::Writer::XLSX;
use Statistics::R;
if(@ARGV<2){
	print STDERR "perl $0 pathway_annotation dif_protein_list\n";
	exit;
}
my $pathway_annotation = shift;
my $dif_protein_list = shift;
my $signi = shift;
my $min_number = shift;
unless($signi){
	$signi=0;
}
$signi-=1;
unless($min_number){
	$min_number=1.5;
}
open FILE,$dif_protein_list;
my %sing;
my %quant;
my $dif_header=<FILE>;
chomp $dif_header;
my @dif_header=split /\t/,$dif_header;
my %quant1;
while(<FILE>){
	chomp;
	my @a  = split /\t/;
	$quant{$a[0]}=$_;
	if($a[3]){
		$quant1{$a[0]}{$a[3]}{$a[1]}=$_;
		$quant1{$a[0]}{"all"}{$a[1]}=$_;
		$sing{$a[0]}{$a[3]}=1;
		$sing{$a[0]}{all}=1;
	}
}
open FILE,$pathway_annotation;
<FILE>;
my %name;
my %dif;
my %dif_n;
my %all_n;
my %all;
my $pathway_tag;
while(<FILE>){
	chomp;
	s/'s//g;
	my @a = split /\t/;
	if($quant{$a[0]} and $a[3]){
		my @b = split /; /,$a[3];
		my $k = $a[0];
		foreach my $pathway_tag(@b){
			$pathway_tag=~s/ - [^\-]*$//;
			$all_n{$pathway_tag}{$k}=1;
			$all{$k}=1;
			if($sing{$k}){
				foreach my $kn(keys %{$sing{$k}}){
					$dif_n{$kn}{$pathway_tag}{$k}=1;
					$dif{$kn}{$k}=1;
				}
			}
		}
	}
}
my $workbook = Excel::Writer::XLSX->new('KEGG_pathway_enrichment.xlsx');
foreach my $k(keys %dif){
	KEGG($k,$k,$workbook);
}
sub KEGG{
	my $sigl_tag = shift;
	my $file_tag = shift;
	my $workbook = shift;
	my $worksheet= $workbook->add_worksheet("$sigl_tag");
	my $worksheet2=$workbook->add_worksheet("$sigl_tag Pathway2Protein");
	my $format_h = $workbook->add_format(font=>'calibri',size=>10,bold=>1,text_wrap=>1,align=>'center',valign=>'vcenter');
	my $format   = $workbook->add_format(font=>'calibri',size=>10);
	my $chi = new Statistics::ChisqIndep;
	my %out_line;
	my %out;
	my %rank;
	my $mm=0;
	my $mm1=0;
	my $mm2=0;
	my $npp=keys %all;
	my $np1=keys %{$dif{$sigl_tag}};
	open OUT,">${file_tag}_protein_pathway_enrichment.xls";
	print OUT "KEGG pathway\tMapping\tBackground\tAll Mapping\tAll Background\tFold enrichment\tp value\tRelated proteins\n";
	my @temp=("KEGG pathway","Mapping","Background","All Mapping","All Background","Fold enrichment","Fisher's exact test p value","-log10(p value)","Related proteins");
	my @temp2=("KEGG pathway","Fisher's exact test p value",@dif_header);
	$worksheet->write_row(0,0,\@temp,$format_h);
	$worksheet2->write_row(0,0,\@temp2,$format_h);
	foreach my $k(keys %{$dif_n{$sigl_tag}}){
		my $n11 = keys %{$dif_n{$sigl_tag}{$k}};
		my @n11 = keys %{$dif_n{$sigl_tag}{$k}};
		my $n1p = keys %{$all_n{$k}};
		if($n11>1){
			my $p_value = calculateStatistic(n11=>$n11+$signi,n1p=>$n1p,np1=>$np1+$signi,npp=>$npp);
			if($p_value==0){
				#	$p_value=1e-100;
				my @load_data = ([$n11,$np1],[$n1p,$npp]);
				$chi->load_data(\@load_data);
				$p_value=$chi->{p_value};
			}
			my $fold = sprintf("%.2f",$n11/$n1p/$np1*$npp);
			my $log_fold = sprintf("%.2f",-log($p_value)/log(10));
			if($fold>$min_number){
				$out_line{$k}="$k\t$n11\t$n1p\t$np1\t$npp\t$fold\t$p_value\t@n11";
				$out{$k}="$k\t$n11\t$n1p\t$np1\t$npp\t$fold\t$p_value\t$log_fold\t@n11";
				foreach my $ttn(@n11){
					foreach my $nk(keys %{$quant1{$ttn}{$sigl_tag}}){
						my @temp2=split /\t/,"$k\t$p_value\t$quant1{$ttn}{$sigl_tag}{$nk}";
						$mm2++;
						$worksheet2->write_row($mm2,0,\@temp2,$format);
					}
				}
				$rank{$k}=$p_value;
			}
		}
	}
	foreach my $k(sort{$rank{$a}<=>$rank{$b}} keys %rank){
		print OUT "$out_line{$k}\n";
		$mm1++;
		if($rank{$k}<0.05 and $mm<25){
			$mm++;
			my @temp = split /\t/,$out{$k};
			$worksheet->write_row($mm,0,\@temp,$format);
		}
	}
	my $chart = $workbook->add_chart(type=>'bar',embedded=>1);
	$chart->add_series(
		categories=>[$sigl_tag,1,$mm,0,0],
		values=>[$sigl_tag,1,$mm,7,7],
		data_labels=>{value=>1},
		fill=>{color=>'red'},
		gap=>50,
	);
	$chart->set_legend(none=>1);
	$chart->set_x_axis(name=>"-log10(Fisher' exact test p value)",major_gridlines=>{visible=>0});
	$chart->set_y_axis(reverse=>1);
	$chart->set_style(10);
	$chart->set_size(width=>720,height=>100+30*$mm);
	$worksheet->insert_chart('J2',$chart);
	print STDERR "$file_tag\n";
#	my $R=Statistics::R->new();
my $cmds=<<EOF;
library(ggplot2)
data=read.table(file="${file_tag}_protein_pathway_enrichment.xls",sep="\\t",header=T)
data
subdata=subset(data,p.value<0.05)
n=nrow(subdata)*0.3+1
if(n<3.5){n=3.5}
pdf(file="${file_tag}_protein_pathway_enrichment.pdf",height=n)
ggplot()+geom_point(data=subdata,aes(x=Fold.enrichment,y=KEGG.pathway,size=Mapping,fill=p.value),shape=21,colour="black")+scale_size_area(max_size=8)+scale_fill_gradientn(colours=c("red","yellow","green"))+theme(panel.background=element_rect(fill="white",colour="black"),panel.grid.major=element_line(color="gray"),panel.grid.minor=element_line(color="gray"),axis.text=element_text(color="black"),axis.title=element_text(color="black",face="bold"))
dev.off()
EOF
	open SOUT,">cmd.R";
	print SOUT "$cmds\n";
	`Rscript cmd.R`;
	unlink "cmd.R";
#	$R->run($cmds);
#	$R->stop;
}
