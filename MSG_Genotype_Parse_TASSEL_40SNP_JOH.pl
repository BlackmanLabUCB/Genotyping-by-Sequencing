#!/usr/bin/perl -w
use warnings;
use strict;

#input genotype file(s) in hapmap format from TASSEL and specify output filename root
my @filelist = @ARGV or die "no filenames provided\n usage: perl MSG_Genotype_Parse_TASSEL_40SNP_JOH.pl hapmap_file(s) outfile";
my %genotype_hash;
my $genotype_hash_ref =\%genotype_hash;
my $output_file = pop @filelist;
my $output_file2 = $output_file . ".calls.txt";
my $output_file3 = $output_file . ".qtl_input.csv";
$output_file = $output_file . ".txt";
my $total_window_number = 0;
my $names_ref;

foreach my $hapmap (@filelist) {
	my ($names_reference, $window_number) = create_genotype_hash($hapmap, $genotype_hash_ref);
	$total_window_number += $window_number;
	$names_ref = $names_reference;
}
foreach my $hapmap (@filelist) {
	print_genotype_count_hash($genotype_hash_ref, $output_file, $names_ref, $total_window_number);
	print_genotype_code_hash($genotype_hash_ref, $output_file2, $output_file3, $names_ref, $total_window_number);
}


sub create_genotype_hash {
	my $gtype_file = shift;
	my $gtype_hash_ref = shift;
	open (GTYPE, "$gtype_file") or die "Couldn't open genotype file";
	my $line = 0;
	my $window_counter = 1;
	my $SNP_count = 0;
	my $position;
	my $chromosome;
	my $open_window_start;
	my @name_array;	my $name_array_ref = \@name_array;
	
	foreach my $SNP_entry (<GTYPE>) {
		chomp $SNP_entry; $SNP_entry =~ s/\s/\t/ig;
		if ($line == 0) {
			@name_array = split ("\t", $SNP_entry);
			splice (@name_array,0,11);
			$line = 1;
		} else { 
			my @marker_array = split ("\t", $SNP_entry);
			my $marker = shift(@marker_array);
			if ($marker =~ /S(\d+)_(\d+)/) {$chromosome = $1; $position = $2;}
			my $bases = shift(@marker_array);
			my $major_allele; my $minor_allele;
			my $cons_base; my $SNP_base;
			if ($bases =~ /([AGCT])\/([AGCT])/) {
				$major_allele = $1; $minor_allele = $2;
			} else {
			next;
			}
			splice (@marker_array,0,9);
			if ($marker_array[0] eq $marker_array[1]) {
				next;
			} elsif (($major_allele eq $marker_array[0]) && ($minor_allele eq $marker_array[1])) {
				$cons_base = $major_allele;  $SNP_base = $minor_allele;
			} elsif (($major_allele eq $marker_array[1]) && ($minor_allele eq $marker_array[0])) {
				$SNP_base = $major_allele;  $cons_base = $minor_allele;
			} else {
				next;
			}
			my $marker_array_ref = \@marker_array;
			if ($SNP_count != 0 && $SNP_count < 40) {
				tabulate_SNPs($name_array_ref, $marker_array_ref, $open_window_start, $gtype_hash_ref, $cons_base, $SNP_base, $chromosome);
				${$gtype_hash_ref}{$chromosome}{$open_window_start}{'end'} = $position if $SNP_count == 39;
				$SNP_count++;
			} else {
				$window_counter++;
				$open_window_start = $position;
				$SNP_count = 0;
				initiate_window($name_array_ref, $open_window_start, $gtype_hash_ref, $chromosome);
				tabulate_SNPs($name_array_ref, $marker_array_ref, $open_window_start, $gtype_hash_ref, $cons_base, $SNP_base, $chromosome);
				$SNP_count++;
			}
		}
			
	}
	${$gtype_hash_ref}{$chromosome}{$open_window_start}{'end'} = $position;
	close GTYPE;
	return ($name_array_ref, $window_counter);
}

sub initiate_window {
	my $name_aref = shift;
	my $window_start = shift;
	my $genotype_href = shift;
	my $chr = shift;
	foreach my $sample (@{$name_aref}) {
		${$genotype_href}{$chr}{$window_start}{$sample}{'A'} = 0;
		${$genotype_href}{$chr}{$window_start}{$sample}{'B'} = 0;
		${$genotype_href}{$chr}{$window_start}{$sample}{'H'} = 0;
		${$genotype_href}{$chr}{$window_start}{$sample}{'N'} = 0;
	}
}

sub tabulate_SNPs {
# Add A, B, H, and Ns to genotype counts in windows nested in individual.
#${$href}{$name}{$window}{A, B, H, or N}{count}
	my $name_aref = shift;
	my $marker_aref = shift;
	my $window = shift;
	my $genotype_href = shift;
	my $consbase = shift; my $altbase = shift;
	my $chr = shift;
	my $cell_count = scalar (@{$name_aref});
	die "uneven numbers in name and genotype arrays" if ($cell_count != scalar (@{$marker_aref}));
	
	for (my $cell = 0; $cell < $cell_count; $cell++) {
		my $sample = ${$name_aref}[$cell];
		if (${$marker_aref}[$cell] eq $consbase) {
			${$genotype_href}{$chr}{$window}{$sample}{'A'}++;
		} elsif (${$marker_aref}[$cell] eq $altbase) {
			${$genotype_href}{$chr}{$window}{$sample}{'B'}++;
		} elsif (${$marker_aref}[$cell] eq 'N') {
			${$genotype_href}{$chr}{$window}{$sample}{'N'}++;
		} elsif (${$marker_aref}[$cell] =~ /[TRYSWKM]/) {
			${$genotype_href}{$chr}{$window}{$sample}{'H'}++
		}
	}
}

sub print_genotype_count_hash {
	my $genotype_href = shift;
	my $output_filename = shift;
	my $name_aref = shift;
	my $windows = shift;
	open (OUT, ">$output_filename");
	print OUT "Start\tEnd\tGenotype\t";
	print OUT join("\t",@{$name_aref}) . "\n";
	foreach my $chr (sort {$a <=> $b} keys %{$genotype_href}) {
		foreach my $window (sort {$a <=> $b} keys %{${$genotype_href}{$chr}}) {
			print OUT "$window\t${$genotype_href}{$chr}{$window}{'end'}\tA\t";
			foreach my $name (@{$name_aref}) {
				print OUT "${$genotype_href}{$chr}{$window}{$name}{'A'}\t";
			}
			print OUT "\n\t\tB\t";
			foreach my $name (@{$name_aref}) {
				print OUT "${$genotype_href}{$chr}{$window}{$name}{'B'}\t";
			}
			print OUT "\n\t\tH\t";
			foreach my $name (@{$name_aref}) {
				print OUT "${$genotype_href}{$chr}{$window}{$name}{'H'}\t";
			}
			print OUT "\n\t\tN\t";
			foreach my $name (@{$name_aref}) {
				print OUT "${$genotype_href}{$chr}{$window}{$name}{'N'}\t";
			}
			print OUT "\n";
		}
	}
	close OUT;
}

sub print_genotype_code_hash {
	my $genotype_href = shift;
	my $output_filename = shift;
	my $qtl_filename = shift;
	my $name_aref = shift;
	my $windows = shift;
	open (OUT2, ">$output_filename");
	open (OUT3, ">$qtl_filename");
	print OUT2 "Chr\tStart\tEnd\tLength\t";
	print OUT2 join("\t",@{$name_aref}) . "\n";
	print OUT3 "id,,";
	print OUT3 join(",",@{$name_aref}) . "\n";
	print OUT3 "sex,,";
	for (my $k = 0; $k < (scalar @{$name_aref}); $k++) {
		print OUT3 "f,";
	} 
	print OUT3 "\npgm,,";
	for (my $k = 0; $k < (scalar @{$name_aref}); $k++) {
		print OUT3 "0,";
	} 
	print OUT3 "\n";
	foreach my $chr (sort {$a <=> $b} keys %{$genotype_href}) {
		foreach my $window (sort {$a <=> $b} keys %{${$genotype_href}{$chr}}) {
			print OUT2 "$chr\t$window\t${$genotype_href}{$chr}{$window}{'end'}\t" . (${$genotype_href}{$chr}{$window}{'end'} - $window) . "\t";
			print OUT3 "$window,$chr,";
			foreach my $name (@{$name_aref}) {
				my $total_all = ${$genotype_href}{$chr}{$window}{$name}{'A'} + ${$genotype_href}{$chr}{$window}{$name}{'B'} + ${$genotype_href}{$chr}{$window}{$name}{'H'};
				if ($total_all < 7) {print OUT2 "N\t"; print OUT3 "-,";next;}
				my $larger;
				if (${$genotype_href}{$chr}{$window}{$name}{'A'} > ${$genotype_href}{$chr}{$window}{$name}{'B'}) {$larger = ${$genotype_href}{$chr}{$window}{$name}{'A'};} else {$larger = ${$genotype_href}{$chr}{$window}{$name}{'B'};}
				if ((${$genotype_href}{$chr}{$window}{$name}{'A'} > 0) && (${$genotype_href}{$chr}{$window}{$name}{'B'} == 0) && (${$genotype_href}{$chr}{$window}{$name}{'H'} == 0)) {
					print OUT2 "A $total_all\t"; print OUT3 "A,";
				} elsif ((${$genotype_href}{$chr}{$window}{$name}{'A'} == 0) && (${$genotype_href}{$chr}{$window}{$name}{'B'} > 0) && (${$genotype_href}{$chr}{$window}{$name}{'H'} == 0)) {
					print OUT2 "B $total_all\t"; print OUT3 "B,";
				} elsif (((${$genotype_href}{$chr}{$window}{$name}{'A'} > 0) && (${$genotype_href}{$chr}{$window}{$name}{'B'} > 0)) && (${$genotype_href}{$chr}{$window}{$name}{'H'} > $larger / 5)) {
					print OUT2 "H $total_all\t"; print OUT3 "H,";
				} elsif (((${$genotype_href}{$chr}{$window}{$name}{'A'} > 0) && (${$genotype_href}{$chr}{$window}{$name}{'B'} > 0)) && (${$genotype_href}{$chr}{$window}{$name}{'H'} == 0)) {
					my $Aobs = ${$genotype_href}{$chr}{$window}{$name}{'A'}; my $Bobs = ${$genotype_href}{$chr}{$window}{$name}{'B'};  my $total = $Aobs+$Bobs;
					my $Aexp = 0.5 * $total; my $Bexp = 0.5 * $total;
					my $chi_square = (($Aobs - $Aexp)**2 / $Aexp) + (($Bobs - $Bexp)**2 / $Bexp);
					if (($chi_square < 5.02 && $total > 20) || ($chi_square < 3.84 && $total > 10 && $total <= 20) || ($chi_square < 2.71 && $total <= 10)) {
						print OUT2 "H $Aobs $Bobs $total $chi_square\t"; print OUT3 "H,";
					} elsif (($chi_square < 9.14 && $total > 20) || ($chi_square < 6.63 && $total > 10 && $total <= 20) || ($chi_square < 3.84 && $total <= 10)) {
						print OUT2 "N $Aobs $Bobs $total $chi_square\t"; print OUT3 "-,";
					} elsif ($Aobs > $Bobs) {
						print OUT2 "A $Aobs $Bobs $total $chi_square\t"; print OUT3 "A,";
					} elsif ($Aobs < $Bobs) {
						print OUT2 "B $Aobs $Bobs $total $chi_square\t"; print OUT3 "B,";
					}
				} elsif (${$genotype_href}{$chr}{$window}{$name}{'H'} > 0) {
					if ((${$genotype_href}{$chr}{$window}{$name}{'A'} > ${$genotype_href}{$chr}{$window}{$name}{'B'}) && (${$genotype_href}{$chr}{$window}{$name}{'A'} / 5 >= (${$genotype_href}{$chr}{$window}{$name}{'H'} + ${$genotype_href}{$chr}{$window}{$name}{'B'}))) {
						print OUT2 "A ${$genotype_href}{$chr}{$window}{$name}{'A'} $total_all\t"; print OUT3 "A,";
					} elsif ((${$genotype_href}{$chr}{$window}{$name}{'A'} < ${$genotype_href}{$chr}{$window}{$name}{'B'}) && (${$genotype_href}{$chr}{$window}{$name}{'B'} / 5 >= (${$genotype_href}{$chr}{$window}{$name}{'H'} + ${$genotype_href}{$chr}{$window}{$name}{'A'}))) {
						print OUT2 "B ${$genotype_href}{$chr}{$window}{$name}{'B'} $total_all\t"; print OUT3 "B,";
					} else {
						print OUT2 "H ${$genotype_href}{$chr}{$window}{$name}{'H'} $total_all\t"; print OUT3 "H,";
					}
				} else {
					print OUT2 "N\t"; print OUT3 "-,";
				}
			}
			print OUT2 "\n"; print OUT3 "\n";
		}
	}
	close OUT2;
	close OUT3;
}
