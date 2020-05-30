use strict;
use warnings;

print "Warning: This test file only works on linux-like environments, with sort, diff, and rm installed."
my $sc = `SummonChimera.pl -d 500 -v 'gi|418489681|ref|NC_019488.1|' -b ../Data/Salmonella.chimeras.blast.txt > ../Data/temp.txt`;
die "SummonChimera not installed correctly\n$sc" if $?;

`sort ../Data/Salmonella.integration.txt > ../Data/Salmonella.integration.sorted.txt`;
`sort ../Data/temp.txt > ../Data/temp.sorted.txt`;
my $diff = `diff ../Data/Salmonella.integration.sorted.txt ../Data/temp.sorted.txt`;

die $diff,$/,"Output not generated properly, may need to update Perl" if $diff;

`rm ../Data/*.sorted.txt ../Data/temp.txt`;

print "No errors encountered when running SummonChimera.pl\n";

