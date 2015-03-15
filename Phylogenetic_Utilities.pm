# Perl often used functions
#
# POD documentation 

=head1 NAME

Phylogenetic_Utilites - DESCRIPTION of Object

=head1 SYNOPSIS

This is a workhorse Perl Module for basic phylogenetic analysis.

=head1 DESCRIPTION

=head1 FEEDBACK

=head1 AUTHOR - Corey M. Hudson

Email coreymhudson-AT-gmail-DOT-com

=cut

# Let the code begin...

package PhylogeneticUtilities;
use strict;
use Statistics::Descriptive;

=description
Step 1: Read list of genes and create protein fasta file
=cut

=head1 FUNCTION
parse_fasta
=head1 DESCRIPTION
Takes a fasta file and returns a hash reference in the form (genename => sequence).
The hash references needs to be dereferenced in the main program e.g.
	$hash_ref = Sequence_functions::parse_fasta($filename);
	%hash = %$hash_ref;
=head1 USAGE
$hash_ref = &parse_fasta(file.fasta)
=cut

sub parse_fasta {
	my $filename = $_[0];
	open(FILE, "<", $filename);
	my %fasta_hash = ();
	my $i = 0;
	my $taxa = "";
	my $sequence = "";
	while ( ! eof(FILE)){
		my $line = readline(FILE);
		chomp $line;
		if($line =~ m/^>/){
			if($i > 0){
				$fasta_hash{$taxa} = $sequence;
			}
			my @line_array = split(/\|/, $line);
			$line = $line_array[0];
			$line = substr($line, 1);
			$taxa = $line;
			$sequence = "";
			$i++;
		}
		else{
			$sequence = $sequence.$line;
		}
	}
	$fasta_hash{$taxa} = $sequence;
	
	return (\%fasta_hash);
}

=head1 FUNCTION
search_gene_fasta
=head1 DESCRIPTION
Takes a fasta file and a list of genes returns a hash reference in the form (found genename => sequence).
The hash references needs to be dereferenced in the main program e.g.
	$hash_ref = Sequence_functions::parse_fasta($filename);
	%hash = %$hash_ref;
=head1 USAGE
found_genes_hash_reference = &search_genes_fasta(db_file, gene_list)
=cut

sub search_genes_fasta {
	my $listname = $_[0];
	open(FILE, "<", $listname);
	my @file = <FILE>;
	close FILE;
	my %gene_list = ();
	foreach my $line (@file){
		chomp $line;
		my $gene = ">".$line;
		$gene_list{$line} = 1;
	}
	my $fasta_file = $_[1];
	my $fasta_hash_ref = &parse_fasta($fasta_file);
	my %fasta_hash = %$fasta_hash_ref;
	my %found_hash = ();
	foreach my $key (keys %gene_list){
		$found_hash{$key} = $fasta_hash{$key};
	}
	return (\%found_hash);
}

sub muscle_alignment {
	#usage alignment_file = &muscle_alignment(input, @params);
	my $input = $_[0];
	my $output = $input.".aln";
	my $input_parameter = "-in $input";
	my $output_parameter = "-out $output";
	my $run = `muscle $input_parameter $output_parameter -quiet`;
	return $output;
}

sub extract_gene_name {
	#usage gene = &extract_gene_name($key);
	my $key = $_[0];
	chomp $key;
	my $gene;
	if($key =~ m/^>/){
		$gene = substr $key, 1;
	}
	else{
		$gene = $key;
	}
	return $gene;
}

sub create_phylip {
	#usage file.phy = &create_phylip(file.fasta)
	my $filename = $_[0];
	my $sequence_hash_ref = &parse_fasta($filename);
	my %sequence_hash = %$sequence_hash_ref;
	my $taxa = scalar keys %sequence_hash;
	my $i = 0;
	my $outfile = $filename.".phy";
	open(FILE, ">>", $outfile);
	my $sequence_length;
	my @sequence_array;
	foreach my $key (keys %sequence_hash){
		if($i == 0){
			@sequence_array = split(//, $sequence_hash{$key});
			$sequence_length = @sequence_array;
			print FILE "$taxa $sequence_length\n";
		}
		$taxa = &extract_gene_name($key);
		print FILE "$taxa     $sequence_hash{$key}\n";
		$i++;
	}
	close FILE;
	return $outfile;
}

sub codon_alignment {
	#usage protein_alignment.fasta.codon = &codon_alignment(protein_alignment, gene_fasta)
	my $protein_hash_ref = &parse_fasta($_[0]);
	my $outfile = $_[0].".codon";
	my $nucleotide_hash_ref = &parse_fasta($_[1]);
	my %protein_hash = %$protein_hash_ref;
	my %nucleotide_hash = %$nucleotide_hash_ref;
	my $gene_name;
	my $protein_sequence;
	my $cdna_sequence;
	my @protein_array;
	my @cdna_array;
	my $char;
	my $i;
	open(FILEO, ">>", $outfile);
	foreach my $key (keys %protein_hash){
		$gene_name = $key;
		chomp $gene_name;
		$protein_sequence = $protein_hash{$gene_name};
		$cdna_sequence = $nucleotide_hash{$gene_name};
		chomp $protein_sequence;
		chomp $cdna_sequence;
		@protein_array = split(//, $protein_sequence);
		@cdna_array = split(//, $cdna_sequence);
		print FILEO ">$gene_name\n";
		$i = 0;
		foreach $char (@protein_array){
			if($char eq "-"){
				print FILEO "---";
			}
			else{
				print FILEO $cdna_array[$i];
				$i++;
				print FILEO $cdna_array[$i];
				$i++;
				print FILEO $cdna_array[$i];
				$i++;
			}
		}
		print FILEO "\n";
	}
	close FILEO;
	return $outfile;
}

sub bootstrap {
	my $random_number = int(rand(100000));
	my $tree = $_[0];
	my $alignment = $_[1];
	my $name = $_[2];
	my $best_tree = $name;
	my $algorithm = "-b $random_number -\#10";
	my $run = "raxmlHPC $algorithm -s $alignment -t $tree -n $name -m GTRGAMMA";
	system($run);
	my $trees = "RAxML_bootstrap.".$name;
	my $outpart = $name."_part";
	my $run2 = "raxmlHPC -f b -m GTRGAMMA -s $alignment -z $trees -t $tree -n $outpart"; 
	system($run2);
	return $best_tree;
}

sub raxml_tree {
	#usage tree_name = &raxml_tree(input);
	my $input = $_[0];
	my $output = $input."_raxmltree";
	my	$input_parameter = "-s $input";
	my	$tree_name_parameter = "-n $output";
	my	$sub_mat_parameter = "-m GTRGAMMA";
	my	$algorithm_parameter = "-f d";
	my $run = `./raxmlHPC $input_parameter $tree_name_parameter $sub_mat_parameter $algorithm_parameter`;
	my $treename = "RAxML_result.".$output;
	return $treename;
}

sub assign_flux {
	print "gene\tgenefamily\tduplications\tlosses\tspecies\tmax\n";
	my $flux_values = $_[0];
	my $genefamilies;
	#$genefamilies = $genefamilies."*.trees";
	my @array = glob "scripts/Notung-2.6/gene_trees/*.tree";
	my $line;
	my @taxa;
	my @line_array;
	my %hash;
	foreach my $filename (@array){
		open(FILEI, "<", $filename);
		$line = readline(FILEI);
		@taxa = split(/[,:()-.]/, $line);
		foreach (@taxa){
			if ($_ =~ m/^AT/){
				@line_array = split(/[\/.]/, $filename);
				$hash{$_} = $line_array[4];
			}
		}
		close FILEI;
	}
	open(FILE1, "<", "notung_output.txt");
	my @file1 = <FILE1>;
	close FILE1;
	my %gf;
	my $gf_name;
	foreach $gf_name (@file1){
		chomp $gf_name;
		@line_array = split(/\t/, $gf_name);
		$gf{$line_array[0]} = $gf_name;
	}
	open(FILE, "<", $flux_values);
	my @file = <FILE>;
	close FILE;
	my %flux_hash;
	my ($line, @line_array, $gene, $flux, $stat, $value, $mean, $sd, $range, $min, $max, $context);
	foreach $line (@file){
		@line_array = split(/\t/, $line);
		$flux = $line_array[1];
		$gene = $line_array[2];
		$context = $line_array[3];
		if ($context =~ m/root/){
			push @{ $flux_hash{$gene}}, $flux;
		}
	}
	my $gf_name_grab;
	foreach my $gene_name (keys %flux_hash){
		if (exists $hash{$gene_name}){
			$gf_name_grab = $hash{$gene_name};
			print "$gene_name\t$gf{$gf_name_grab}\t";
			$stat = Statistics::Descriptive::Full->new();
			foreach $value  (@{$flux_hash{$gene_name}}){
				$stat->add_data(abs($value));
			}
			$min = $stat->min();
			$max = $stat->max();
			$range = $stat->sample_range();
			print "$max\n";
		}
	}
	
	return 0;
}
;1
