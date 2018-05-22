{
    package VarianceUncorrected;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use VarMathClassical;
    
     sub new {
        my $class = shift;
        my $poolsize=shift;
        my $maf=shift;
        my $mincoverage=shift;
        my $maxcoverage=shift;
        
        die "Uncorrected population genetic estimators can only be calculated for a minor allele frequency of 1; Instead it was set to $maf" unless $maf == 1;
        # get_theta_calculator($b,$n,$snp)
        # get_pi_calculator($b,$n,$snp)
        my $self = bless {
                          pi=>get_classical_pi_calculator(),
                          theta=>get_classical_theta_calculator(),
                          d=>get_classical_D_calculator()
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
        my $pi_sum=$measurecalculater->($snps);

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
        my $theta=$measurecalculator->($snps);
        my $toret=0;
        $toret=$theta/$covercount if $covercount;
        return $toret;
    }

    sub _calculate_d
    {
        my $self=shift;
        my $snps=shift;
        my $measurecalculator=$self->{d};
        my $toret=$measurecalculator->($snps);
        return $toret;
    }   
}

1;
