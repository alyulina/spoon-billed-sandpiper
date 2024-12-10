
use strict;

my ($line, $gi, $first, $outfile, $before, $after, $number, $count, $params,
$average, %list, %bneck, $num, $bneck, $difference, %ainbred, %binbred, $ainbred, $binbred);


open (INFILE, $ARGV[0]) or die "Kozyol!\n";
open (OUTFILE,  ">averages_demography.txt") or die "Kozyol!\n";
open (OUTFILE2, ">averages_bottleneck.txt") or die "Kozyol!\n";
open (OUTFILE3, ">averages_demography_inbreeding.txt") or die "Kozyol!\n";
open (OUTFILE4, ">averages_bottleneck_inbreeding.txt") or die "Kozyol!\n";

while(<INFILE>)
{$before = $_;

chomp($before);

$ainbred=$binbred=$bneck=$after=$before;
$after =~ s/before_/after/;
$bneck =~ s/before_/after_bneck/;
$ainbred =~ s/fitness_before_/inbread_after/; 
$binbred =~ s/fitness_before_/inbread_after_bneck/;


($params)=($before =~ m/before_(.*)\.txt/);
$params =~ s/_/\t/g;

print OUTFILE "$params";
print OUTFILE2 "$params";
print OUTFILE3 "$params";
print OUTFILE4 "$params";

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
open (AINBRED, $ainbred) or die "Kozyol!!!\n";
while(<AINBRED>)
{$line = $_;

if ($line =~ m/inbread/)
{
    if ($num != 0){$ainbred{$num} = $average/$count;}
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

$ainbred{$num} = $average/$count;
print "$ainbred\n";

$average=$num=$count=0;
open (BINBRED, $binbred) or die "Kozyol!!!\n";
while(<BINBRED>)
{$line = $_;

if ($line =~ m/inbread/)
{
    if ($num != 0){$binbred{$num} = $average/$count;}
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

$binbred{$num} = $average/$count;
print "$binbred\n";


$average=$num=$count=0;
open (AFTER, $after) or die "Kozyol!!!\n";
while(<AFTER>)
{$line = $_;


if ($line =~ m/after/)
{

    if ($num != 0)
       {$average /= $count; 
#       $difference = $average - $list{$num};  print OUTFILE  "\t$difference";
       $difference = $average/$list{$num} - 1;  print OUTFILE  "\t$difference";
       $difference = 0;
#       $difference = $average - $bneck{$num}; print OUTFILE2 "\t$difference";
       $difference = $average/$bneck{$num} - 1; print OUTFILE2 "\t$difference";
       $difference = 0;
       $difference = $ainbred{$num}/$list{$num} - 1;  print OUTFILE3  "\t$difference";
       $difference = 0;
       $difference = $ainbred{$num}/$binbred{$num} - 1;  print OUTFILE4  "\t$difference";
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
$difference = $ainbred{$num}/$list{$num} - 1;  print OUTFILE3  "\t$difference";
$difference = $binbred{$num}/$list{$num} - 1;  print OUTFILE4  "\t$difference";

print OUTFILE  "\n";
print OUTFILE2 "\n";
print OUTFILE3 "\n";
print OUTFILE4 "\n";

print "$after\n";

}
