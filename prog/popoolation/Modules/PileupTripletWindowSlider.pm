
{   
    package PileupTripletWindowSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use PileupTripletSlider;
    
    # annotated triplet window
    # data, chr, lower, upper, window, count_codons, $cound_valid_codons, $count_snps
    # TRIPLET definition
    # frame, pileup, chr, start, strand, valid, valid_frame, valid_coverage, valid_codon, count_snps, codon
    
    sub new
    {
        my $class=shift;
        my $pts=shift;
        my $window=shift;
        my $step=shift;
        

        
        return bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            pts=>$pts,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
    }
    
    sub nextWindow
    {
        my $self=shift;
        my $pts=$self->{pts};
        
        #get the current window, and the current chromosome
        my $curwin=$self->{curwin};
        
        #curwin contains triplets
        # TRIPLET definition
        # frame, pileup, chr, start, strand, valid, valid_frame, valid_coverage, valid_codon, count_snps, codon
        
        my $curChr="";
        $curChr=$curwin->[0]{chr} if @$curwin;
        
        my $resetchr=0;
        
        # empty unnecessary entries
        EMPTY: while(@$curwin)
        {
            my $e=shift @$curwin;
            if($e->{start}>$self->{lower})
            {
                unshift @$curwin, $e;
                last EMPTY;
            }
            
        }
        
        # fill with novel entries
        my $tri;
        FILL:while($tri=$self->_nexttriplet)
        {
            $curChr=$tri->{chr} unless $curChr;
            
            if($tri->{chr} eq $curChr && $tri->{start} <= $self->{upper})
            {
                push @$curwin,$tri;
            }
            else
            {
                $resetchr=1 if $tri->{chr} ne $curChr;
                $self->_buffertriplet($tri);
                last FILL;
            }
        }
        
        return undef unless $curChr;
        
        my $toret=_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
        if($resetchr or not defined($tri))
        {
            # we transgressed the boundaries to the next chromosome
            # reset the windows and the current buffer
            $self->{lower}=0;
            $self->{upper}=$self->{window};
            $self->{curwin}=[];
        }
        else
        {
            # next time we will still be in the same chromosome
            # increase the upper and lower boundaries by the stepsize and set the current buffer
            $self->{upper}+=$self->{step};
            $self->{lower}+=$self->{step};
            $self->{curwin}=$curwin;
        }

        return $toret;
    }
    
    sub _annotateWindow
    {
        my $curwin =shift;
        my $curChr =shift;
        my $lower  =shift;
        my $upper  =shift;
        my $window =shift;
        
        # TRIPLET definition
        # frame, pileup, chr, start, strand, valid, valid_frame, valid_coverage, valid_codon, count_snps, codon
        
        my ($count_codons, $count_valid_codons, $count_useful,$count_one_snp)=(0,0,0,0);
        for my $t (@$curwin)
        {
            $count_codons++;
            $count_one_snp++ if($t->{valid} and $t->{count_snps}==1);
            $count_valid_codons++ if $t->{valid};
        }
        
        
        # annotated triplet window
        # data, chr, lower, upper, window, count_codons, cound_valid_codons, 
        my $entry={
            data                =>$curwin,
            chr                 =>$curChr,
            lower               =>$lower,
            upper               =>$upper,
            window              =>$window,
            count_codons        =>$count_codons,
            count_valid_codons  =>$count_valid_codons,
        };
        
        return $entry;
    }
    
    
    sub _nexttriplet
    {
        my $self=shift;
        
        
        my $pts=$self->{pts};
        my $buffer=$self->{buffer};
        return shift @$buffer if @$buffer;
        return $pts->next();
    }
    
    sub _buffertriplet
    {
        my $self=shift;
        my $entry=shift;
        push @{$self->{buffer}},$entry;
    }


}
1;
