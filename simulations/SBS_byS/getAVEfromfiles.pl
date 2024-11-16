use strict;

my ($line, $file, $par, @par, $count, $i, $j, $count_par);

open (INFILE, "$ARGV[0]") or die "Kozyol!\n";

while (<INFILE>)
{$file = $_;
chomp($file);

open (INFILE2, "$file") or die "Kozyol!!\n";


$count=0;
<INFILE2>; <INFILE2>; <INFILE2>;
while (<INFILE2>)
{$line = $_;

$count_par=0; 
while($line)
{

($par) = ($line =~ m/(\S+)/);
($line) = ($line =~ m/\S+\s+(.*)/);

$par[$count][$count_par] += $par;
$count_par++;

}
$count++;


}

}


print "s	total_1	het2pq	homo	sample_total_1	sample_het2pq	sample_homo	total_1	het2pq	homo	sample_total_1	sample_het2pq	sample_homo\n";



for (my $i = 0; $i < 90; $i++) 
{

for (my $j = 0; $j < 13; $j++) 
{

$par[$i][$j] /= 30;

print "$par[$i][$j]\t";

}


print "\n";

}

