use strict;
use warnings;
use Chimera;
use Getopt::Long;

my $pid = 90;
my $a_cut;
my $e_cut = 1e-2;
my $unique = 15;
my $distance = 500;
my @discordant_pairs;
my @chimera_blast;
my $virus_id;
my $insert = 500;
my $help;
my $cov = 0;
my $read_strands = 'fr';
my $chimeric_file;

GetOptions("s|sam=s" => \@discordant_pairs,"d|distance=i" => \$distance, 'e|evalue=s' => \$e_cut,
	'u|unique=i' => \$unique, 'v|virus=s' => \$virus_id, 'b|blast=s' => \@chimera_blast,"h|help" => \$help,
	'pid=i' => \$pid, 'ma|min_aligned=i' => \$a_cut, 'rs|read_strands=s' => \$read_strands,
	'i|insert=i' => \$insert, 'cr|chimera_report=s' => \$chimeric_file);

die Help() if $help || (!scalar(@chimera_blast) && !scalar(@discordant_pairs));

#Filter and retrieve BLASTS for chimeras
my @chimera = GetChimericReadBlasts($virus_id,$e_cut,$pid,$a_cut,\@chimera_blast);

#Retrieve SAM lines for mapped chimeras
GetChimericSAM(\@chimera,$virus_id,@discordant_pairs);

#Sort and determine chimera cases
my %chimeric_junctions = MapChimericJunctions(\@chimera,$unique,$read_strands,$insert);

#Call integrations from sorted chimeras
my %integrations = ReportIntegrations(\%chimeric_junctions,$distance);

Report(\%integrations,\%chimeric_junctions);

ReadReport(\%chimeric_junctions,$chimeric_file) if $chimeric_file;

sub Help{
	return <<HELP;
SummonChimera Copyright (C) 2014 Joshua Patrick Katz
This program comes with ABSOLUTELY NO WARRANTY; 
This is free software, and you are welcome to redistribute and/or modify 
it under the conditions outlined in the terms of the GNU General Public Licenses
as published by the Free Software Foundation either version 3 or later.
GNU General Public License version 3 is provided with this program as LICENSE.txt

SummonChimera version 1.0 by Joshua Patrick Katz joshpk105\@gmail.com
Usage:
	SummonChimera.pl [options]* -s|sam <sam> -b|blast <blast> -v|virus <virus>

	<sam>	Filename of sam file with discordant mapped reads. 
	<blast>	Filename of outfmt 6 style blast results of chimera reads.
	<virus>	Identifier of virus used in both blast and sam files all other identified treated as host.

	-s|sam, -b|blast, -v|virus can be used multiple times E.g. '-s file1.sam -s file2.sam'

Options (defaults in parentheses):
	-d|distance <integer>	
		The host nucleotide deletion size expected (DEFAULT: 500)
	-e|evalue <float>	
		The BLAST e value cut off for chimera hits (DEFAULT: 1e-2)
	-u|unique <integer>	
		Minimum number of uniquely aligned nucleotides for BLAST hits (DEFAULT: 15)
	-rs|read_strands <ff|fr|rf|rr>  
		The strands of the sequences for paired end data (DEFAULT:fr)
	-i|insert <int>         
		The insert size between paired-end reads, used for clustering (DEFAULT: 500)
	-pid <integer>		
		Minimum percent identity for BLAST chimera hits
	-ma|min_aligned <int>	
		Minimum number of total aligned nucleotides for BLAST chimera hits
	-cr|chimera_report <string>	
		Output filename for a report file of all chimeras and their reads

HELP
}
sub GetChimericSAM{
	my $chimera = shift;
	my %reads;
	my $virus = shift;
	my @discordant = @_;
	foreach my $d (@discordant){
	        open IN, $d;
	        while(<IN>){
	                next if $_ =~ /^\@/;
	                chomp;
	                my @cols = split "\t", $_;
	                $cols[0] =~ s/[\\\/]\d$//;
			
			$reads{$cols[0]}{'type'} = 'sam';
			#$cols[2] = 'chr'.$cols[2] if $cols[2] =~ /^\d+$/;
	                if($cols[2] eq $virus){
	                        $reads{$cols[0]}{'virus'} = \@cols;
	                }
	                else{
	                        $reads{$cols[0]}{'host'} = \@cols;
	                }
			if(scalar(keys(%{$reads{$cols[0]}})) == 2){
				push(@{$chimera},$reads{$cols[0]});
			}
	        }
	}
}
sub GetChimericReadBlasts{
	my ($virus,$evalue, $pid, $a_cut, $blasts) = @_;
	my %reads;
	my %amb;
	foreach my $f (@{$blasts}){
		open IN, $f;
		while(<IN>){
			chomp;
			my @cols = split "\t", $_;
			next if $cols[10] > $evalue || $cols[2] < $pid || (defined($a_cut) && $cols[3] < $a_cut);
			my $key;
			$key = 'virus' if $cols[1] eq $virus;
			$key = 'host' if $cols[1] ne $virus;
			next if !$key;
			if(!defined($reads{$cols[0]}{$key})){
        	                $reads{$cols[0]}{$key} = \@cols;
        	                $reads{$cols[0]}{'amb'} = 0;
				$reads{$cols[0]}{'type'} = 'blast';
        	        }
        	        elsif($reads{$cols[0]}{$key}[11] < $cols[11]){
        	                #print join("\t",@{$reads{$cols[0]}{$key}}),$/;
        	                $reads{$cols[0]}{$key} = \@cols;
				$reads{$cols[0]}{'amb'} = 0;
				delete $amb{$cols[0]};
        	                #print "less\n";
        	                #print join("\t",@cols),$/,$/;
        	        }
        	        elsif($reads{$cols[0]}{$key}[11] == $cols[11]){
				$reads{$cols[0]}{'amb'} += 1;
				$amb{$cols[0]} = 1;
        	        }
        	        else{
        	                $reads{$cols[0]}{$key}[13] +=1;
	                }
		}
	}
	foreach my $k (keys %amb){
		delete $reads{$k};
	}
	return map {$reads{$_}} keys %reads;
}
sub MapChimericJunctions{
	my $sites = shift;
	my $unique = shift;
	my $read_strands = shift;
	my $insert = shift;
	my %chimeras;
	my $count = 0;
	foreach my $r (@{$sites}){ #MapChimericJunctions
                next if !$r->{'virus'} || !$r->{'host'};
                my $a = $r->{'amb'};
		my $c;
		if($r->{'type'} eq 'blast'){
                	$c = new Chimera('type' => 'blast','unique' => $unique,'virus' => $r->{'virus'},
				'host' => $r->{'host'});
		}
		else{
			$c = new Chimera(type => 'sam', virus => $r->{'virus'}, insert => $insert,
				host => $r->{'host'}, read_strands => $read_strands);
		}
                $count++ if $c;
               	next if !$c; 
                my $order = $c->Get('order');
                my $ss = $c->Get('same_strand');
                my $chromo = $c->Get('chromosome');
                $chimeras{$ss}{$order}{$chromo} = [] if !$chimeras{$ss}{$order}{$chromo};
                $chimeras{$ss}{$order}{$chromo} = AddJunction($chimeras{$ss}{$order}{$chromo},$c);
        }
	return %chimeras;
}
sub PrintJunctions{
	my $chimera = shift;
	my $used = shift;
	my %chi = %{$chimera};
	foreach my $ss (keys %chi){
		foreach my $o (keys %{$chi{$ss}}){
			foreach my $c (keys %{$chi{$ss}{$o}}){
				foreach my $j (@{$chi{$ss}{$o}{$c}}){
					next if $used->{$j};
					print join("\t",'NA',$ss,$o,$c,$j->OutArray),$/;
				}	
			}
		}
	}
}
sub AddJunction{
	my $list = shift;
	my $junct = shift;
	my @toReturn;
	while(my $o = shift @{$list}){
		if($junct->CompareChimera($o)){
			$junct->CombineChimera($o);
		}
		else{
			push(@toReturn,$o);	
		}
	} 
	push(@toReturn,$junct);
	return \@toReturn;	
}
sub Report{
	my $int = shift;
	my $junct = shift;
	my %used;
	print join("\t",qw(index same_strand order host_id blast_cov sam_cov),
        	qw(amb_virus amb_host order blast_cov sam_cov amb_virus amb_host)),$/;
	for(my $i =0; $int->{$i}; $i++){
		print join("\t",$i,$int->{$i}{'HV'}->Get('same_strand'),'HV',$int->{$i}{'HV'}->Get('chromosome'),
			$int->{$i}{'HV'}->OutArray,$int->{$i}{'HV'}->Distance($int->{$i}{'VH'},'ambiguous_host'),
			'VH',$int->{$i}{'VH'}->OutArray),$/;
		$used{$int->{$i}{'HV'}} = 1;
		$used{$int->{$i}{'VH'}} = 1;
	}	
	PrintJunctions($junct,\%used);
}
sub ReadReport{
	my $chimera = shift;
	my $out= shift;
	open BLAST, ">$out.blast.txt";
	open SAM, ">$out.sam.txt";
	foreach my $ss (keys %{$chimera}){
		foreach my $o (keys %{$chimera->{$ss}}){
			foreach my $c (keys %{$chimera->{$ss}{$o}}){
				foreach my $chi (@{$chimera->{$ss}{$o}{$c}}){
					foreach my $j (@{$chi->Get('chimera')}){
						my $fh;
						if($j->Get('type') eq 'blast'){
							$fh = \*BLAST;	
						}
						else{
							$fh = \*SAM;
						}
						print $fh join("\t",@{$j->Get('virus')},$ss,$o),$/;
                                                print $fh join("\t",@{$j->Get('host')},$ss,$o),$/;
					}
				}
			}	
		}
	}
}
sub ReportIntegrations{
	my $chimeras = shift;
	my $distance = shift;
	my %used;
	my %ints;
	foreach my $ss (qw(0 1)){
		my %chrs = map {$_ => 1} keys %{$chimeras->{$ss}{'HV'}};
		my @chrs_keys = grep {$chrs{$_}} keys %{$chimeras->{$ss}{'VH'}};
		foreach my $c (@chrs_keys){
			my @donor = sort {$a->Where($b)} @{$chimeras->{$ss}{'HV'}{$c}};
                        my @acc = sort {$a->Where($b)} @{$chimeras->{$ss}{'VH'}{$c}};
                        GetIntegrations(\@donor,\@acc,\%ints,$distance);
		}
	}	
	return %ints;
}

sub GetIntegrations{
	my $donors = shift;
	my $accept = shift;
	my $ints = shift;
	my $distance = shift;
	my $i = scalar keys %{$ints};
	while(my $donor = shift @{$donors}){
                foreach my $acc (@{$accept}){
                        my $d = $donor->Distance($acc,'ambiguous_host');
                        if($d >= 0 && $d < $distance){
                                $ints->{$i}{'HV'} = $donor;
                                $ints->{$i}{'VH'} = $acc;
                                $ints->{$i}{'complete'} = 0;
                                $i++;
                        }
                        elsif($d > $distance){
                        	last;
			}
                }
        }
}
sub GetSeqs{
	my $fasta = shift;
	my %seqs;
	my $seqI = new Bio::SeqIO(-file => $fasta,-format => 'fasta');
	while(my $seq = $seqI->next_seq){
		$seqs{$seq->id} = $seq;
	}
	return %seqs;
}
