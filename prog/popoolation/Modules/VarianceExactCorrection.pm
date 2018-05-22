{
    package VarianceExactCorrection;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use VarMath;
    my $debug=0;

    sub new {
        my $class = shift;
        my $poolSize=shift;
        my $maf=shift;
        my $mincoverage=shift;
        my $maxcoverage=shift;

        
        
        # get_theta_calculator($b,$n,$snp)
        # get_pi_calculator($b,$n,$snp)

        
        my $self = bless {
                          mincoverage=>$mincoverage,
                          maxcoverage=>$maxcoverage,
                          n=>$poolSize,
                          b=>$maf,
                          pi=>get_pi_calculator(),
                          theta=>get_theta_calculator(),
                          d=>get_D_calculator()
                          }, __PACKAGE__;
        return $self;
    }
    


    sub calculate_measure
    {
        my $self=shift;
        my $measure=shift;
        my $snps=shift;
        my $covercount=shift;
        
        $measure=lc($measure);
        
        if($measure eq "pi")
        {
            return $self->_calculate_pi($snps,$covercount);
        }
        elsif($measure eq "theta")
        {
            return $self->_calculate_theta($snps,$covercount);
        }
        elsif($measure eq "d")
        {
            return $self->_calculate_d($snps,$covercount);
        }
        else
        {
            die "unknown measure to calculate $measure";
        }
    }
    
    
    sub _calculate_pi
    {
        my $self=shift;
        my $snps=shift;
        my $covercount=shift;
    
        my $measurecalculater=$self->{pi};
        my $pi_sum=$measurecalculater->($self->{b},$self->{n},$snps);

        my $toret=0;
        $toret=$pi_sum/$covercount if $covercount;
        return $toret;
        
    }
    
    
    sub _calculate_theta
    {
        my $self=shift;
        my $snps=shift;
        my $covercount=shift;
        my $measurecalculator=$self->{theta};
        my $theta=$measurecalculator->($self->{b},$self->{n},$snps);
        my $toret=0;
        $toret=$theta/$covercount if $covercount;
        return $toret;
    }
    
    
    sub _calculate_d
    {
        my $self=shift;
        my $snps=shift;
        my $covercount=shift;
        my $poolsize=$self->{n};
        my $mincount=$self->{b};
        my $mincoverage=$self->{mincoverage};
        unless($debug)
        {
            die "Corrected Tajima's D error\nMinimum count needs to be set to 2 for calculating the corrected Tajima's D;\n".
            "In case 2 is insufficient we recommend to subsample the reads to a smaller coverage" unless $mincount==2;
            die "Corrected Tajima's D error\n Poolsize >> mincoverage (as internal aproximation: 3 * minimumcoverage < poolsize)" unless 3*$mincoverage < $poolsize;
        }

        my $measurecalculator=$self->{d};
        my $d=$measurecalculator->($mincount,$poolsize,$mincoverage,$snps);
        my $toret=0;
        $toret=$d if $covercount;
        return $toret;
    }
}
1;
