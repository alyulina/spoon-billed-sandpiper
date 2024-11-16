
use strict;

my ($line, $gi, $first, $outfile, $before, $after, $number, $count, $params, $average, %list, %bneck, $num, $bneck, $difference);


open (INFILE, $ARGV[0]) or die "Kozyol!\n";
open (OUTFILE, ">averages_demography.txt") or die "Kozyol!\n";
open (OUTFILE2, ">averages_bottleneck.txt") or die "Kozyol!\n";

while(<INFILE>)
{$before = $_;

chomp($before);

$bneck=$after=$before;
$after =~ s/before_/after/;
$bneck =~ s/before_/after_bneck/;


($params)=($before =~ m/before_(.*)\.txt/);
$params =~ s/_/\t/g;

print OUTFILE "$params";
print OUTFILE2 "$params";

$num=0;
open (BEFORE, $before) or die "Kozyol!!\n";
while(<BEFORE>)
{$line = $_;


if ($line =~ m/before/)
{

    if ($num != 0){$list{$num} = $average/$count;}

     $average=$count=0;
     
     $num++;
}
else
{

($line) = ($line =~ m/(\S+)/);

$average += $line;
$count++;

}

}

$list{$num} = $average/$count;

print "$before\n";


$average=$num=$count=0;
open (BNECK, $bneck) or die "Kozyol!!!\n";
while(<BNECK>)
{$line = $_;


if ($line =~ m/after/)
{

    if ($num != 0){$bneck{$num} = $average/$count;}

     $average=$count=0;
     
     $num++;
}
else
{

($line) = ($line =~ m/(\S+)/);

$average += $line;
$count++;

}

}

$bneck{$num} = $average/$count;

print "$bneck\n";



$average=$num=$count=0;
open (AFTER, $after) or die "Kozyol!!!\n";
while(<AFTER>)
{$line = $_;


if ($line =~ m/after/)
{

    if ($num != 0)
       {$average /= $count; 
       $difference = $average - $list{$num};  print OUTFILE  "\t$difference";
       $difference = 0;
       $difference = $average - $bneck{$num}; print OUTFILE2 "\t$difference";
       }

     $average=$count=0;
     
     $num++;
}
else
{

($line) = ($line =~ m/(\S+)/);

$average += $line;
$count++;
}

}

$average /= $count; 
$difference = $average - $list{$num};  print OUTFILE  "\t$difference";
$difference = $average - $bneck{$num}; print OUTFILE2 "\t$difference";

print OUTFILE  "\n";
print OUTFILE2 "\n";

print "$after\n";

}
