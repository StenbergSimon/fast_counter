#!/usr/bin/env perl 

use warnings;
use strict;
use Bio::SeqIO;
use List::Util qw( min max );

#Rewritten version

my $options = options(\@ARGV);
my @options = @{$options};
my $file = $options[0];
my $format = $options[1];
my $length = $options[2];
my $genomesize = $options[3];
my $num_seq = 0;
my @lengths;
my $sum = 0;
my $cov = 0;
my $gccount = 0;

my $inseq = Bio::SeqIO->new(
							-file => "<$file",
							-format => "$format",
							);
while( my $seqref = $inseq->next_seq){
		
		my $seq = $seqref->seq();
		my @array = split(//,$seq);
		foreach (@array){
					if($_=~m/[GC]{1}/){
						++$gccount;
						}
				}
		push (@lengths, length($seq));
		
				
		
}

@lengths = sort {$a <=> $b} @lengths;

if($length == 1){
	my $pos = 1;
	open(FH,">lengths.tsv");
	foreach(@lengths){
		print FH $pos, "\t", $_, "\n";
		++$pos;
		}
		
	}
	
foreach (@lengths){
	$sum = $sum + $_;

}
	
my $min = min @lengths;
my $max = max @lengths;


print "\n-----------------------------------------------------	
		\nfast_counter.pl results:
		\nNumber of Sequences: ",scalar(@lengths),
		"\nTotal number of bases: $sum",
		"\nAverage sequence length: ",sprintf('%.3f',average(\@lengths)),
		"\nMedian sequence lentgh: ", median(\@lengths),
		"\nShortest sequence $min",
		"\nLongest sequence: $max",
		"\nGC Percentage: ", sprintf('%.3f',(($gccount/$sum)*100));
		
unless ($genomesize == 0){
	$cov = $sum / $genomesize;
	print "\nGenome coverage: ", sprintf('%.3f',$cov), "\n";
	}
print "\n-----------------------------------------------------	\n";

###### SUBROUTINES ######

sub help{
	die "\n------------------------------------------------------------
		\n \t fast_counter.pl
		\n Usage: fast_counter.pl file.fastx [options]\n
		\n Options:
		\n -fa\t#Input is in FASTA format (Default FASTQ)
		\n -lo\t#Print lengths into a tab-separated list where this script is run, output is lengths.tsv  (will overwrite without question)
		\n -g [int]\t#Input genomesize in bp to caculate the coverage - total length of sequences divided by the genomesize\n";
	}

sub options{ 
	my $format;
	my $length;
	my $genomesize = 0;
	my $i2 = 0;
	my @options = @{$_[0]};
	my @sender;
	
	foreach(@options){
	

		if($_ eq "-g"){
						$genomesize = $options[($i2+1)];
								if($genomesize =~ m/\D/){
								die("\nError: Illegal genomesize, has to be numerical\n", help());
								}
				
						}
		elsif($_ eq "-fa"){
						 $format = "fasta";
						}
					
		
		elsif(($_ eq "-h") or ($_ eq "--help")){
						
						die(help());
					
						}
		elsif($_ eq "-lo"){
						$length = 1;
						
						}
		$i2++;		
		
	
		}

		$file = $options[0];
		
		unless($file =~ m/.fasta$/ || $file =~ m/.fa/ || $file =~ m/.fastq/){
			die("\n Illegal file format." , help());
			}
	unless (defined $format){
		print "\nAssuming fastq...\n";
		$format = "fastq";
		}
	unless (defined $length){
		$length = 0;
		}
	
	

	# Check that all arguments submitted are valid

	foreach(@ARGV){
		unless($_ =~ m/-h/ || $_ =~ m/--help/ || $_ =~ m/-lo/ || $_ =~ m/-fa/ || $_ =~ m/$genomesize/ || $_ =~ m/$file/ || $_ =~ m/-g/ ){
			print "\n Error: Invalid arguments\n ";
			die (help());
			}
		}
	
	push(@sender, $file);	
	push(@sender, $format);
	push(@sender, $length);
	push(@sender, $genomesize);
	

	return (\@sender);
}
	

	
sub median{

	my @array_unsorted = @{$_[0]};

	my @array = sort {$a <=> $b} @array_unsorted;
	my $pos = ((scalar(@array)/2) - 1);
	
	unless(defined($array[$pos]) && defined($array[($pos+1)])){
	
		return "0";
		}
	
	else{
	
		 if( (scalar(@array)) % 2){

			 
			 my $median = $array[($pos+0.5)];
			 return $median;
	
	
		 } else {
	
			 my $pos = ((scalar(@array)/2) - 1);
	
			 my $tot = ($array[$pos] + $array[($pos+1)]);
			 my $median = $tot/2;
	
			 return $median;
	
			 }
	}
}

sub average{

	my @array = @{$_[0]};
	my $sum = 0;
	my $mean = 0;
	
	foreach (@array){
	
		$sum = $sum + $_ ; 
	
	}
	
	unless(scalar(@array) == 0){
	
		$mean = $sum / (scalar(@array));
		return $mean;
		}	
	else{
		

	return 0;
	}
}



