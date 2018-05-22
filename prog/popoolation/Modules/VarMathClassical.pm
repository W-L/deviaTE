package VarMathClassical;
use strict;
use warnings;

use Math::BigRat; 

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT=qw(get_classical_theta_calculator get_classical_pi_calculator get_classical_D_calculator);

my $a1_buffer;
my $a2_buffer;

sub get_classical_pi_calculator
{ 
    return sub
    {
        my $snps=shift;
        my $sites =@$snps;
        return 0 unless $sites;
        my $pi_sum=0;
        foreach my $snp (@$snps)
        {
            my $pi_snp=_average_pairwise_differences($snp);
            #my $M=$snp->{eucov};
            #my $pi_snp=1;
            #$pi_snp-=($snp->{A}/$M)**2;
            #$pi_snp-=($snp->{T}/$M)**2;
            #$pi_snp-=($snp->{C}/$M)**2;
            #$pi_snp-=($snp->{G}/$M)**2;
            $pi_sum+=$pi_snp;
        }
        
        return $pi_sum;
    }
}

sub _average_pairwise_differences
{
    my $snp=shift;
    
    my $pwc =nover2($snp->{eucov});
    my $asim=nover2($snp->{A});
    my $tsim= nover2($snp->{T});
    my $csim= nover2($snp->{C});
    my $gsim= nover2($snp->{G});
    my $sim = $asim+$tsim+$csim+$gsim;
    my $pi=1-$sim/$pwc;
    return $pi;
}

sub nover2
{
    my $n=shift;
    return 0 if $n < 2;
    return ($n*($n-1))/2;
}



sub get_classical_theta_calculator
{
    my $a1c=get_a1_buffer();
    return sub
    {
        my $snps=shift;
        my $sites =@$snps;
        return 0 unless $sites;
        
        my $n=_median_coverage($snps);
        return $sites/($a1c->($n));
    }
}


sub get_classical_D_calculator
{
    my $pic=get_classical_pi_calculator();
    my $thetac=get_classical_theta_calculator();
    my $varc=_get_sqrt_variance_calculator();
    
    return sub
    {
        my $snps=shift;
        my $sites=@$snps;
        return 0 unless $sites;
        
        # pi theta and sqrt(variance)
        my $pi=$pic->($snps);
        my $theta=$thetac->($snps);
        my $var=$varc->($snps);
        
        my $smalld=$pi-$theta;        



        my $tajd=0;
        $tajd=$smalld/$var if $var;
        return $tajd;
    }
}

sub _get_sqrt_variance_calculator{
    my $a1c=get_a1_buffer();
    my $a2c=get_a2_buffer();
    return sub{
        my $snps=shift;
        my $sites=@$snps;
        my $n=_median_coverage($snps);
        
        my $a1=$a1c->($n);
        my $a2=$a2c->($n);
        my $b1=_calculate_b1($n);
        my $b2=_calculate_b2($n);
        my $c1= $b1 - (1/$a1);
        my $c2=_calculate_c2($a1,$a2,$b2,$n);
        my $e1= $c1 / $a1;
        my $e2= $c2 / (($a1**2) + $a2);
        
        my $var = $e1 * $sites + $e2 * $sites * ($sites - 1);
        $var=sqrt($var) if $var;
        return $var;
    }
}



sub _calculate_c2
{
    my $a1=shift;
    my $a2=shift;
    my $b2=shift;
    my $n=shift;
    my $t1=($n+2)/($a1 * $n);
    my $t2=($a2 / ($a1**2));
    return ($b2 - $t1 + $t2);
}


sub _calculate_b1
{
    my $n=shift;
    return ($n + 1)/(3 * ($n - 1));
}

sub _calculate_b2
{
    my $n=shift;
    my $above = 2 * (($n**2) + $n + 3);
    my $below = 9 * $n * ($n - 1);
    return $above/$below;
}


sub get_a1_buffer
{
    $a1_buffer={} unless $a1_buffer;
    return sub
    {

        my $n=shift;
        return $a1_buffer->{$n} if exists($a1_buffer->{$n});
        
        my $toret=0;
        foreach my $k (1..$n-1)
        {
            my $t1=(1/$k);
            $toret+=$t1;
        }
        #store the value in the buffer and return it 
        $a1_buffer->{$n}=$toret;
        return $toret;
    } 
}

sub get_a2_buffer
{
    $a2_buffer={} unless $a2_buffer;
    
    return sub
    {

        my $n=shift;
        return $a2_buffer->{$n} if exists($a2_buffer->{$n});
        
        my $toret=0;
        foreach my $k (1..$n-1)
        {
            my $t1=(1/($k**2));
            $toret+=$t1;
        }
        #store the value in the buffer and return it 
        $a2_buffer->{$n}=$toret;
        return $toret;
    } 
}


 sub _median_coverage
    {

        my $snps=shift;
        my $coverages= [map {$_->{eucov}} @$snps];
        $coverages=[sort {$a<=>$b} @$coverages];
        my $snpcount=@$coverages;
        die "no snps for calculating median coverage"  unless $snpcount;
        
        if($snpcount%2)
        {
            #uneven
            # 5/2 =2.5; int(2.5)=2; array starting at zero 2 is element number 3
            # int(1/2)=0 -> also true for a single element
            return $coverages->[int($snpcount/2)];
            
        }
        else
        {
            # even
            # 6/2 =3 need 2 +3
            my $i1=$snpcount/2;
            my $i2=$i1-1;
            
            my $toret=int(($coverages->[$i1]+$coverages->[$i2])/2);
            return $toret;
        }
    }

1;