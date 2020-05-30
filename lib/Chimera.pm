package Chimera;
use strict;
use warnings;
use List::Util qw(max min);

my @ORIENTATION = qw(HV VH);

#SummonChimera Copyright (C) 2014 Joshua Patrick Katz
sub new{ 
	# 'type' => 'sam|blast', 'virus' => 'ref to result array', 
	# 'host' => 'ref to result array', 'unique' cutoff if type is blast 
	my $class = shift;
	my %params = @_;
	die "Incorrect params for Chimera" if !($params{'type'} || !$params{'virus'} || !$params{'host'});
	if($params{'type'} eq 'blast'){
		$params{'unique'} |= 15; 
	}
	elsif($params{'type'} eq 'sam'){
		die "Chimera parameter insert not set for type:sam\n" if !$params{'insert'};
		die "Chimera parameter read_strands not set/correct needs rf,fr,rr,or ff\n" 
			if $params{'read_strands'} !~ /[fr]{2}/;
	}
	else { die "Chimera parameter 'type' not set to sam or blast\n"; };

	my $self = \%params;
	$self->{'class'} = $class;
	bless $self, $class;
	return $self->newBlast if $self->{'type'} eq 'blast';
	return $self->newSAM if $self->{'type'} eq 'sam';
	return;	
}
sub newSAM{
	my $self = shift;
	my $virus = $self->{'virus'};
	my $host = $self->{'host'};
	$self->{'virus_aligned'} = [$virus->[3],$virus->[3]+length($virus->[9])];
	$self->{'host_aligned'} = [$host->[3],$host->[3]+length($host->[9])];
	$self->{'chromosome'} = $host->[2];
	my @strands = split '', $self->{'read_strands'};
	my $rs = $strands[0] eq $strands[1];
	$self->{'order'} = $self->DiscordantOrder;
	my @virus_wing = @{$self->{'virus_aligned'}};
        my @host_wing = @{$self->{'host_aligned'}};
        $self->{'wings'}{'virus'} = \@virus_wing;
        $self->{'wings'}{'host'} = \@host_wing;
	$self->{'same_strand'} = $self->DiscordantSameStrand;
	$self->{'chimera'} = [$self->Copy];
	$self->DiscordantAmbiguous;
	return $self;
}
sub DiscordantAmbiguous{ #removed range indicator
	my $self = shift;
	my $w = $self->{'wings'};
	my $i = $self->{'insert'};
	my $iv = $i-($w->{'virus'}[1]-$w->{'virus'}[0]);
	my $ih = $i-($w->{'host'}[1]-$w->{'host'}[0]);
	if($self->{'order'} eq 'HV'){
		$self->{'ambiguous_host'} = [$w->{'host'}[1],$w->{'host'}[1]];
		$self->{'ambiguous_virus'} = [$w->{'virus'}[0],$w->{'virus'}[0]] if $self->{'same_strand'};
		$self->{'ambiguous_virus'} = [$w->{'virus'}[1],$w->{'virus'}[1]] if !$self->{'same_strand'};
	}
	else{
		$self->{'ambiguous_host'} = [$w->{'host'}[0],$w->{'host'}[0]];
		$self->{'ambiguous_virus'} = [$w->{'virus'}[1],$w->{'virus'}[1]] if $self->{'same_strand'};
                $self->{'ambiguous_virus'} = [$w->{'virus'}[0],$w->{'virus'}[0]] if !$self->{'same_strand'};
	}
}
sub Copy{
	my $self = shift;
	my $copy = bless { %$self }, ref $self;
	return $copy;	
}
sub print{
	my $self = shift;
	print join("\t",'c,o,ss',$self->{'chromosome'},$self->{'order'},$self->{'same_strand'},
		$self->{'virus'}[0]),$/;
	print join("\t",'virus_wings',@{$self->{'wings'}{'virus'}}),$/;
	print join("\t",'host_wings',@{$self->{'wings'}{'host'}}),$/;
	print join("\t",'ambiguous_virus',@{$self->{'ambiguous_virus'}}),$/ if $self->{'ambiguous_virus'};
	print join("\t",'ambiguous_host',@{$self->{'ambiguous_host'}}),$/ if $self->{'ambiguous_host'};
	print join("\t",'total_chimeras', scalar(@{$self->{'chimera'}})),$/;
	print join(",",'chimera_comp:', (map{$_->Get('type')}  @{$self->{'chimera'}})),$/;
}
sub DiscordantSameStrand{ #validate this functions
	my $self = shift;
	my $vf = $self->{'virus'}->[1];
	my $hf = $self->{'host'}->[1];
	my @rs = split '', $self->{'read_strands'};
	my $ss = ((($vf & 16) xor ($hf & 16)) xor ($rs[0] eq $rs[1]));
	$ss |= 0;
	return $ss;
}
sub DiscordantOrder{ # validate this function
	my $self = shift;
	my @rs = split '', $self->{'read_strands'};
	my $hf = $self->{'host'}->[1];
	my $hs = ($hf & 128) && 1;
	#hf & 128 is true if host read is second, hf & 16 is if read was complemented to map
	return $ORIENTATION[($hf & 128) xor (($hf & 16) xor ($rs[$hs] eq 'r'))];
}
sub newBlast{
	my $self = shift;
	my $virus = $self->{'virus'};
	my $host = $self->{'host'};
	my ($RV0,$RV1,$V0,$V1) = @{$virus}[6..9];
	$self->{'virus_aligned'} = FiveToThree($V0,$V1);
	my ($RH0,$RH1,$H0,$H1) = @{$host}[6..9];
	$self->{'host_aligned'} = FiveToThree($H0,$H1);
	$self->{'chromosome'} = $host->[1];
	my ($order,$overlap) = LengthOverlap($RV0,$RV1,$RH0,$RH1);
	return if ($host->[3]-$overlap < $self->{'unique'} || $virus->[3]-$overlap < $self->{'unique'});
	
	$self->{'order'} = $ORIENTATION[!($order xor $self->{'host_aligned'}[2])];
	$self->{'overlap'} = $overlap;
	#print join("\t",@{$self->{'host_aligned'}}),$/;
	my @virus_wing = @{$self->{'virus_aligned'}};
	my @host_wing = @{$self->{'host_aligned'}};
	$self->{'wings'}{'virus'} = \@virus_wing;
	$self->{'wings'}{'host'} = \@host_wing;
	my @ambiguous = Region($V0,$V1,$H0,$H1,$order,$overlap);
	my @amb_v =  @ambiguous[0..1];
	my @amb_h = @ambiguous[2..3];
	$self->{'ambiguous_virus'} = \@amb_v;
	$self->{'ambiguous_host'} = \@amb_h;
	$self->{'same_strand'} = ($self->{'host_aligned'}[2] == $self->{'virus_aligned'}[2]);
	$self->{'same_strand'} = 0 if !$self->{'same_strand'};	
	$self->{'chimera'} = [$self->Copy];
	print join("\t",@{$self->{'virus'}}),$/ if $self->{'chromosome'} =~ 'chrUn';
	return $self;
}
sub CompareChimera{
	my $self = shift;
	my $chimera = shift;
	return 0 if $self->{'same_strand'} != $chimera->Get('same_strand') ||
		$self->{'order'} ne $chimera->Get('order');
	if($self->{'type'} eq $chimera->{'type'} && $self->{'type'} eq 'blast'){
		my @self_virus = @{$self->{'ambiguous_virus'}};
		my @self_host = @{$self->{'ambiguous_host'}};
		my @comp_virus = @{$chimera->Get('ambiguous_virus')};
		my @comp_host = @{$chimera->Get('ambiguous_host')};
		
		return (EquivCoord(@self_virus,@comp_virus) && EquivCoord(@self_host,@comp_host));
	}
	elsif(($self->{'type'} eq $chimera->{'type'} && $self->{'type'} eq 'sam') ||
		($self->{'type'} eq 'mixed' || $chimera->Get('type') eq 'mixed')){
		my @self_virus = $self->WingSpan('virus');
		my @self_host = $self->WingSpan('host');
		my $wings = $chimera->Get('wings');
		my @chimera_virus = $chimera->WingSpan('virus');
		my @chimera_host = $chimera->WingSpan('host');
		#print join("\t",join('-',@self_virus,@chimera_virus),Overlap(\@self_virus,\@chimera_virus), 
		#	join('-',@self_host,@chimera_host),Overlap(\@self_host,\@chimera_host)),$/;
		return (Overlap(\@self_virus,\@chimera_virus) && Overlap(\@self_host,\@chimera_host));
	}
	else{
		return (SAMBlastCompare($self,$chimera));
	}
}
sub WingSpan{
	my $self = shift;
	my $key = shift;
	my $wings = $self->{'wings'};
	my $i = $self->{'insert'};
	my $iw = $i-($wings->{$key}[1]-$wings->{$key}[0]);
	my @wing_span = ($wings->{$key}[0]-$iw,$wings->{$key}[1]+$iw);
	return @wing_span;
}
sub Overlap{
	my $first = shift;
	my $second = shift;
	return 1 if ($first->[0] <= $second->[0] && $first->[1] >= $second->[0]) ||
	($first->[0] <= $second->[1] && $first->[1] >= $second->[1]) ||
	($first->[0] >= $second->[0] && $first->[0] <= $second->[1]);
	return 0;
}
sub SAMBlastCompare{
	my $sam = $_[$_[1]->Get('type') eq 'sam'];
	my $blast = $_[$_[1]->Get('type') eq 'blast'];
	my @sam_virus = $sam->WingSpan('virus');
	my @sam_host = $sam->WingSpan('host');
	my $blast_virus = $blast->Get('ambiguous_virus');
	my $blast_host = $blast->Get('ambiguous_host');
	
	return 1 if Overlap(\@sam_virus,$blast_virus) && Overlap(\@sam_host,$blast_host); 
	return 0;
}
sub CombineChimera{
	my $self = shift;
	my $chimera = shift;
	push(@{$self->{'chimera'}},@{$chimera->Get('chimera')});
	$self->{'insert'} = $chimera->Get('insert') if $chimera->Get('insert');
	$self->Wings($chimera);	
	if($self->{'type'} eq $chimera->Get('type') && $self->{'type'} eq 'sam'){
		$self->DiscordantAmbiguous;		
	}
	else{ #($self->{'type'} eq $chimera->Get('type') && $self->{'type'} eq 'blast'){ 
		$self->LeastAmbiguous;
	}
	$self->{'type'} = 'mixed' if $self->{'type'} ne $chimera->Get('type');
}
sub OutArray{
	my $self = shift;
	my $blast = scalar grep {$_->Get('type') eq 'blast'} @{$self->{'chimera'}};
	my $sam = scalar grep {$_->Get('type') eq 'sam'} @{$self->{'chimera'}};
	my $amb_virus = join('-',@{$self->{'ambiguous_virus'}});
	my $amb_host = join('-',@{$self->{'ambiguous_host'}});
	my $wings_virus = join('-',@{$self->{'wings'}{'virus'}});
	my $wings_host = join('-',@{$self->{'wings'}{'host'}});
	return ($blast,$sam,$amb_virus,$amb_host,$wings_virus,$wings_host); 
}
sub LeastAmbiguous{
	my $self = shift;
	my %amb;
	my $max_chimera;
	my $max_count;
	foreach my $c (@{$self->{'chimera'}}){
		next if $c->Get('type') ne 'blast';
		my $k = join('',@{$c->Get('ambiguous_host')},@{$c->Get('ambiguous_virus')});
		if(!$max_count){
			$max_count = 1;
			$max_chimera = $c;
		}
		$amb{$k}{'total'}++;
		if($max_count < $amb{$k}{'total'}){
			$max_count = $amb{$k}{'total'};
			$max_chimera = $c;
		}
	}
	$self->{'ambiguous_virus'} = $max_chimera->Get('ambiguous_virus');
	$self->{'ambiguous_host'} = $max_chimera->Get('ambiguous_host');
}
sub Distance{
	my $self = shift;
	my $chimera = shift;
	my $key = shift;
	#print join("\t",($chimera->Get($key)->[0] - $self->{$key}[1]),$self->OutArray,$chimera->OutArray),$/;
	#return 0 if $chimera->Get($key)->[0] == $self->{$key}[0] && $chimera->Get($key)->[1] == $self->{$key}[1];
	#return $chimera->Get($key)->[0] - $self->{$key}[1]; # not sure why i originally did this comparison
	return $chimera->Get($key)->[0] - $self->{$key}[0];
}
sub Where{
	my $self = shift;
	my $chimera = shift;
	return $self->{'ambiguous_host'}[0] <=> $chimera->Get('ambiguous_host')->[0];
}
sub Wings{
	my $self = shift;
	my $chimera = shift;
	my $cwings = $chimera->Get('wings');

	$self->{'wings'}{'virus'}[0] = min($cwings->{'virus'}[0],$self->{'wings'}{'virus'}[0]);
	$self->{'wings'}{'virus'}[1] = max($cwings->{'virus'}[1],$self->{'wings'}{'virus'}[1]);
	$self->{'wings'}{'host'}[0] = min($cwings->{'host'}[0],$self->{'wings'}{'host'}[0]);
	$self->{'wings'}{'host'}[1] = max($cwings->{'host'}[1],$self->{'wings'}{'host'}[1]);
}
sub Get{
	my $self = shift;
	my $key = shift;
	return $self->{$key};
}
sub EquivCoord{
	return ($_[0] == $_[2] && $_[1] == $_[3]);
}
sub toArray{
	my $self = shift;
	return ($self->{'chromosome'},@{$self->{'ambiguous_host'}},@{$self->{'ambiguous_virus'}},
		$self->{'same_strand'});
}
sub FiveToThree{
	my @a = @_;
	return [@a,1] if $_[0] < $_[1];
	@a = reverse @a;
	return [@a,0];
}
sub LengthOverlap{
        my ($ax,$ay,$bx,$by) = @_;
        my $d1 = $ay-$bx;
        my $d0 = $by-$ax;
	my $o = $ax < $bx;
	$o = 0 if !$o;
        if($ay >= $bx && $ay <= $by){
                return ($o,$d1);
        }
        if($by >= $ax && $by <= $ay){
                return ($o,$d0);
        }
        my $dn = ($d0,$d1)[$d1 < $d0];
        return ($o,$dn);
}
sub Region{
        my ($v0,$v1,$h0,$h1,$d,$o) = @_;
        my $vs = (0,1)[$v1 > $v0];
        my $hs = (0,1)[$h1 > $h0];
        my ($v5p,$v3p,$h5p,$h3p);
        if($d && $o >= 0){ #correct
                if($vs){ ($v5p,$v3p) = ($v1-$o,$v1);}
                else{ ($v5p,$v3p) = ($v1,$v1+$o);}
                if($hs){($h5p,$h3p) = ($h0,$h0+$o);}
                else{ ($h5p,$h3p) = ($h0-$o,$h0);}
        }
        elsif($d && $o < 0){ #correct
                if($vs){ ($v5p,$v3p) = ($v1,$v1+abs($o));}
                else{ ($v5p,$v3p) = ($v1-abs($o),$v1);}
                if($hs){($h5p,$h3p) = ($h0-abs($o),$h0);}
                else{ ($h5p,$h3p) = ($h0,$h0+abs($o));}
        }
        elsif(!$d && $o >= 0){ #
                if($vs){ ($v5p,$v3p) = ($v0,$v0+$o);}
                else{ ($v5p,$v3p) = ($v0-$o,$v0);}
                if($hs){($h5p,$h3p) = ($h1-$o,$h1);}
                else{ ($h5p,$h3p) = ($h1,$h1+$o);}
        }
        elsif(!$d && $o < 0){
                if($vs){ ($v5p,$v3p) = ($v0-abs($o),$v0);}
                else{ ($v5p,$v3p) = ($v0,$v0+abs($o));}
                if($hs){($h5p,$h3p) = ($h1,$h1+abs($o));}
                else{ ($h5p,$h3p) = ($h1-abs($o),$h1);}
        }
        return ($v5p,$v3p,$h5p,$h3p);
}
1;
