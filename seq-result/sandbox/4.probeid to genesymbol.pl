open IN,"probeid_genesymbol.txt";
while(<IN>){
	chomp;
	@line = split "\t";
	$hash{$line[0]} = $line[1];
}
close IN;

open DIFF,"probeid.FDR0.05.exprs.txt";
open OUT,">genesymbol.FDR0.05.exprs.txt";
open INPUT,">all_genes.txt";
while(<DIFF>){
	chomp;
	if($. == 1){
		print OUT "$_\n";
		next;
	}
	@line = split "\t";
	if(defined($hash{$line[0]})){
	$line[0] = $hash{$line[0]};
	print OUT "$line[0]";
	print INPUT "$line[0]\n";
	for(1..@line-1){
		print OUT "\t$line[$_]";
	}
	print OUT "\n";
	}else{
		next;
	}
}
close DIFF;
close OUT;
close INPUT;		

open IN,"diff.probeid.FDR0.05.txt";
open OUT,">diff.genesymbol.FDR0.05.txt";
while(<IN>){
	chomp;
	if($. == 1){
		print OUT "$_\n";
		next;
	}
	@line = split "\t";
	if(defined $hash{$line[0]}){
	$line[0] = $hash{$line[0]};
	print OUT "$line[0]";
	for(1..@line-1){
		print OUT "\t$line[$_]";
	}
	print OUT "\n";
	}else{
		next;
	}
}
close IN;
close OUT;
