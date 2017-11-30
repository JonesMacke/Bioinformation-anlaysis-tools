#!/usr/bin/perl -w
use strict;
use Text::NSP::Measures::2D::Fisher::right;
use Statistics::ChisqIndep;
use Excel::Writer::XLSX;
my $input = shift;
my $diff  = shift;
my $signal = shift;
my $min_number = shift;
unless($signal){
	$signal=0;
}
$signal-=1;
unless($min_number){
	$min_number=1.5;
}
my %diff;
my %quant;
open FILE,$diff or die "Can't open file $diff~\n";
my $diff_header=<FILE>;
my @diff_header=split /\t/,$diff_header;
my %quant1;
while(<FILE>){
	chomp;
	my @a = split /\t/;
	$quant{$a[0]}=$_;
	if($a[3]){
		$quant1{$a[0]}{$a[3]}{$a[1]}=$_;
		$quant1{$a[0]}{"all"}{$a[1]}=$_;
		$diff{$a[0]}{$a[3]}=1;
		$diff{$a[0]}{all}=1;
	}
}
my %np1;
my %n11;
my %npp;
my %n1p;
my %s_n11;
my %s_np1;
my %hash;
open FILE,$input or die "Can't open file $input~\n";
while(<FILE>){
	chomp;
	my @a = split /\t/;
	$hash{$a[1]}{$a[2]}{"level"}=$a[3];
	$hash{$a[1]}{$a[2]}{"name"}=$a[4];
	if($quant{$a[0]}){
		$npp{$a[1]}{$a[0]}=1;
		$n1p{$a[1]}{$a[2]}{$a[0]}=1;
		if($diff{$a[0]}){
			foreach my $k(keys %{$diff{$a[0]}}){
				$s_n11{$k}{$a[1]}{$a[2]}{$a[0]}=1;
				$s_np1{$k}{$a[1]}{$a[0]}=1;
			}
		}
	}
}
my $workbook = Excel::Writer::XLSX->new('GO_enrichment.xlsx');
foreach my $k(keys %s_n11){
	Enrich($k,$workbook);
}
sub Enrich{
	my $tag = shift;
	my $workbook = shift;
	my $worksheet=$workbook->add_worksheet("$tag");
	my $worksheet1=$workbook->add_worksheet("$tag GO2Protein");
	my $format_h = $workbook->add_format(font=>'calibri',size=>10,bold=>1,text_wrap=>1,align=>'center',valign=>'vcenter');
	my $format   = $workbook->add_format(font=>'calibri',size=>10);
	my $chi = new Statistics::ChisqIndep;
	open OUT,">${tag}_protein_GO_enrichment.xls";
	print OUT "GO Terms Level 1\tGO Terms ID\tGO Terms Level\tGO Terms Description\tMapping\tBackground\tAll Mapping\tAll Background\tFold Enrichment\tFisher's exact Pvalue\tRelated Protein\n";
	my @raw_string = ("GO Terms Level 1","GO Terms Description","GO Terms Level","GO Terms ID","Mapping","Background","All Mapping","All Background","Fold Enrichment","Fisher' exact test P value","-log10(p value)","Related Protein");
	my @raw_string1 = ("GO Terms Level 1","GO Terms ID","GO Terms Level","GO Terms Description","Fisher's exact test P value",@diff_header);
	$worksheet->write_row(0,0,\@raw_string,$format_h);
	$worksheet1->write_row(0,0,\@raw_string1,$format_h);
	my $row_count=0;
	my $row_count2=0;
	my $row_count1=0;
	my %go_count;
	$go_count{"Cellular Component"}=8;
	$go_count{"Biological Process"}=14;
	$go_count{"Molecular Function"}=8;
	foreach my $k1(keys %{$s_n11{$tag}}){
		my $row_signl;
		my %rank;
		my %out;
		my %out1;
		foreach my $k2(keys %{$s_n11{$tag}{$k1}}){
			my $n11 = keys %{$s_n11{$tag}{$k1}{$k2}};
			my @n11 = keys %{$s_n11{$tag}{$k1}{$k2}};
			my $n1p = keys %{$n1p{$k1}{$k2}};
			my $np1 = keys %{$s_np1{$tag}{$k1}};
			my $npp = keys %{$npp{$k1}};
			if($hash{$k1}{$k2}{"level"}>1 and $n11>1){
				my $pvalue = calculateStatistic(n11=>$n11+$signal,n1p=>$n1p,np1=>$np1+$signal,npp=>$npp);
				if($pvalue==0){
					#	$pvalue=1e-100;
					my @load_data = ([$n11,$np1],[$n1p,$npp]);
					$chi->load_data(\@load_data);
					$pvalue=$chi->{p_value};
					if($pvalue==0){
						$pvalue=1e-175;
					}
				}
				my $fold = sprintf("%.2f",$n11/$np1/$n1p*$npp);
				my $log_pvalue = sprintf("%.2f",-log($pvalue)/log(10));
				if($fold>$min_number){
					$rank{$k2}=$pvalue;
					$out{$k2}="$k1\t$k2\t$hash{$k1}{$k2}{level}\t$hash{$k1}{$k2}{name}\t$n11\t$n1p\t$np1\t$npp\t$fold\t$pvalue\t@n11\n";
					$out1{$k2}="$k1\t$hash{$k1}{$k2}{name}\t$hash{$k1}{$k2}{level}\t$k2\t$n11\t$n1p\t$np1\t$npp\t$fold\t$pvalue\t$log_pvalue\t@n11";
					foreach my $kn(@n11){
						foreach my $nk(keys %{$quant1{$kn}{$tag}}){
							my $tem_s="$k1\t$k2\t$hash{$k1}{$k2}{level}\t$hash{$k1}{$k2}{name}\t$pvalue\t$quant1{$kn}{$tag}{$nk}\n";
							my @temp1=split /\t/,$tem_s;
							$row_count1++;
							$worksheet1->write_row($row_count1,0,\@temp1,$format);
						}
					}
				}
			}
		}
		my $mm=0;
		foreach my $k(sort{$rank{$a}<=>$rank{$b}}keys %rank){
			print OUT $out{$k};
			if($rank{$k}<0.05 and $mm<$go_count{$k1}){
				$row_count++;
				$mm++;
				my @temp = split /\t/,$out1{$k};
				if($row_signl){
					shift @temp;
					$worksheet->write_row($row_count,1,\@temp,$format);
				}
				else{
					$row_signl=1;
					$worksheet->write_row($row_count,0,\@temp,$format);
				}
			}
		}
	}
	my $chart = $workbook->add_chart(type=>'bar',embedded=>1);
	$chart->add_series(
		categories=>[$tag,1,$row_count,0,1],
		values=>[$tag,1,$row_count,10,10],
		data_labels=>{value=>1},
		fill=>{color=>'red'},
		gap=>50,
	);
	$chart->set_legend(none=>1);
	$chart->set_x_axis(name=>"-log10(Fisher' exact test p value)",major_gridlines=>{visible=>0});
	$chart->set_y_axis(reverse=>1);
	$chart->set_style(10);
	$chart->set_size(width=>720,height=>100+$row_count*30);
	$worksheet->insert_chart('M2',$chart);
}
