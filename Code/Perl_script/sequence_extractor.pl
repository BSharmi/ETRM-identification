#!/usr/local/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $help = "";
my $ref_genome = "";
my $chromosome = "";
my $position = 0;
my $sequence_length = 0;

my $commend_line = GetOptions(
       'h|help' => \$help,
       'r|ref=s' => \$ref_genome,
       'c|chromosome=s' => \$chromosome,
       'p|position=i' => \$position,
       'l|length=i' => \$sequence_length);

unless($commend_line)
{
	die "please respecify option(s)...\n";
}

if($help || !-f $ref_genome || $ref_genome !~ /\.fa/)
{
	print "usage:\n";
	print "perl sequence_extractor.pl [-h/--help] <-r/--ref mm10.fa> <-c/--chromosome chr1> <-p/--position 34242122> <-l/--length 2> <loci>\n";
	exit(1);
}

warn "\n" . "=" x 50 . "\n";
warn "Parsing reference genome at $ref_genome\n";
warn "Please wait...\n";

if($chromosome ne "" && $position > 0 && $sequence_length > 0) {
  &create_revised_ref_genome($ref_genome, $chromosome, $position, $sequence_length);
}
elsif(@ARGV > 0) {
  &batch_extractor($ref_genome, $ARGV[0]) if($ARGV[0] =~ /\.txt/ || $ARGV[0] =~ /\.bed/);
}

sub create_revised_ref_genome
{
	my $ref_genome = shift;
	my $chromosome = shift;
	my $position = shift;
	my $sequence_length = shift;
	my $chrom;
	my $len = "";
	my $my_seq;
	
	my $subscribe = 0;
	my $rows = 0;
	
	open(FA, "<$ref_genome") or die "cannot open the given ref genome: $ref_genome\n";
	while(<FA>)
	{
		if(/^>/)
		{
			# to avoid some errors of chromosome identifier for dna chromosome descriprion.
			$chrom = (split(/\s/))[0]; 
			$chrom =~ s/^>//;
			
			# if chromosome is the mode '1 dna:chromosome chromosome:GRCm38:....', to transform into 'chr1'.
			$chrom = "chr$chrom" if $chrom =~ /^\d+$/; 
			chomp $chrom;
						
			$rows = 0;
			$subscribe = 0;
			
			next;
		}
		
		if(uc($chrom) eq uc($chromosome))
		{
			my $tmp = $_; chomp $tmp;
			next if $tmp eq "";
			
			$len = length($tmp) unless $len ne "";
			
			$rows ++;
			my $start = $len * ($rows - 1) + 1;
			my $end = $len == length($tmp) ? $len * $rows : ($len * ($rows - 1) + length($tmp));

			if($position >= $start && $position <= $end)
			{
				my $relative_coord = $position % $len;
				
				if($len - $relative_coord + 1 >= $sequence_length)
				{
				  print "The sequence specified:\n$chromosome\t$position\t$sequence_length\n", uc(substr($tmp, $relative_coord - 1, $sequence_length)), "\n";
				  last;
				}
				else
				{
					$my_seq = substr($tmp, $relative_coord - 1, $len - $relative_coord + 1);
					while(length($my_seq) != $sequence_length)
					{
						$tmp = <FA>; chomp $tmp; 
						next if $tmp eq "";
						
						if(length($my_seq) + length($tmp) < $sequence_length)
						{
							$my_seq .= $tmp;
						}
						else
						{
							$my_seq .= substr($tmp, 0, $sequence_length - length($my_seq));
						}
					}
					
					print "The sequence specified:\n$chromosome\t$position\t$sequence_length\n", uc($my_seq), "\n";
					last;
				}
			}
		}
	}
	close(FA);
}

sub batch_extractor
{
	my $ref_genome = shift;
  my $bedfile = shift; # format: chrom start end 

	my $chrom;
	my $len;
  my %region;
  my @array;
	my $seq;
  
  open(IN, "<$bedfile") or die "cannot open $bedfile \n";
  my $head = <IN>; chomp $head;
  while(<IN>) {
    chomp;
    next if $_ eq "";
    
    @array = split("\t");
    
    push(@{$region{$array[0]}}, [$array[1]-1, $array[2]]);
  }
  close(IN);
  
  $bedfile =~ s/.bed//;
	$bedfile =~ s/.txt//;
  open(OUT, ">$bedfile\.withSquence.fa") or die "Cannot create $bedfile\.withSquence.fa\n";
  #print OUT "$head\tsequence\n";
	
  $/ = ">";
	open(FA, "<$ref_genome") or die "cannot open the given ref genome: $ref_genome\n";
	while(<FA>)
	{
		chomp;
    next if $_ eq "";
    
    @array = split("\n");
    
    # to avoid some errors of chromosome identifier for dna chromosome descriprion.
		$chrom = (split(/\s/, $array[0]))[0]; 	
	  # if chromosome is the mode '1 dna:chromosome chromosome:GRCm38:....', to transform into 'chr1'.
		$chrom = "chr$chrom" if $chrom =~ /^\d+$/; 
    
    shift @array;
    if(exists $region{$chrom}) {
      $len = length($array[0]);
      foreach my $el (@{$region{$chrom}}) {
        if(int($$el[1] / $len) >= int($$el[0] / $len) + 2) {
          my $seq = substr($array[int($$el[0] / $len)], $$el[0] % $len, $len - $$el[0] % $len).join("", (@array)[(int($$el[0] / $len) + 1)..(int($$el[1] / $len) - 1)]).substr($array[int($$el[1] / $len)], 0, $$el[1] % $len);
          print OUT ">$chrom:".(1 + $$el[0])."-$$el[1]\n$seq\n";
        }
        elsif(int($$el[1] / $len) >= int($$el[0] / $len) + 1) {
          my $seq = substr($array[int($$el[0] / $len)], $$el[0] % $len, $len - $$el[0] % $len).substr($array[int($$el[1] / $len)], 0, $$el[1] % $len);
          print OUT ">$chrom:".(1 + $$el[0])."-$$el[1]\n$seq\n";
        }
        else {
          my $seq = substr($array[int($$el[0] / $len)], $$el[0] % $len, $$el[1] - $$el[0]);
          print OUT ">$chrom:".(1 + $$el[0])."-$$el[1]\n$seq\n";
        }
      }
    }
    else {
      next;
    }
	}
	close(FA);
}