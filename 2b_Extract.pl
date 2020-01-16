#!/usr/bin/env perl

# -- program name
print "-"x60, "\n";
print "2b_Extract.pl v 3.0\n";
print "Created 02 Mar 2010 E Meyer\n";
print "Last modified 07 Mar07 2013\n";
print "-"x60, "\n";

# -- program description and required arguments
unless ($#ARGV == 2)
        {print "Counts the number of occurences of a type IIb recognition site\n";
	print "or sequence motif in a set of DNA sequences.\n";
        print "Please note -- as written, this only works for symmetrical sites.\n";
        print "Output:\t a fasta file of those sites, named by position.\n";
        print "Usage:\t script input site output\n";
        print "Arguments:\n";
        print "\t input\t a fasta file containing the sequences to be searched \n";
	print "\t site\t recognition site. e.g. for AlfI (N{12}GCAN{6}TGC{12}) enter \".{12}GCA.{6}TGC.{12}\"\n";
        print "\t output\t a fasta file of those sites \n";
        print "\n"; exit;
        }
my $seqfile = $ARGV[0];
open(SEQ, $seqfile);

$patt = $ARGV[1];
@pa = split /[{}]/, $patt;
foreach $p (@pa)
	{
	if ($p =~ /\d+/) {$size+=$p; next;}
	if ($p =~ /\w/)
		{
		$p =~ s/\.//g;
		$lp = length($p);
		$size+=$lp; 
		next;
		}
	}
print "Forward sequence of recognition site: ", $patt, "\n";

my $outfile = $ARGV[2];
open(OUT, ">$outfile");

my $found = 0;
my $nseqs = 0;
$ss = "";
while(<SEQ>)
	{
	chomp;
	if ($_ =~ />/)
		{
		if ($ss ne "")
			{
			while($ss =~ /(?=$patt)/g) 
				{
				$loci = pos($ss)+$size;
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "\n";
				print OUT substr($ss, $loci-$size, $size), "\n";
				$found++;
				}
			$ss = "";
			}
		$nseqs++;
		$idi = $_; $idi =~ s/>//; $idi =~ s/\s+.*//;
		}
	else
		{
		$ss = $ss.$_;
		}
	}
if ($ss ne "")
	{
	while($ss =~ /(?=$patt)/g) 
		{
				$loci = pos($ss)+$size;
				print OUT ">", $idi, "_", $loci-$size+1, "_", $loci, "\n";
				print OUT substr($ss, $loci-$size, $size), "\n";
				$found++;
		}
	$ss = "";
	}

print $nseqs, " sequences searched\n";
print $found, " recognition sites found altogether\n";

