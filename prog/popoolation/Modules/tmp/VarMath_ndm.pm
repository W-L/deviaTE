package VarMath;
use strict;
use warnings;

use Math::BigRat; 

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT=qw(noverk get_theta_calculator get_pi_calculator get_D_calculator);
our @EXPORT_OK = qw(get_thetadiv_buffer get_thetadiv_buffer_sqr);


# get_aMnm_calculator( M, n, m ) (M..coverage, n..poolsize,  running var ~ coverage)
# get_pidiv_calculator(b,n,M) (b..mincoverage, n..poolSize, M..coverage)
# get_thetadiv_calculator(b,n,M)
# get_theta_calculator($b,$n,$snp)
# get_pi_calculator($b,$n,$snp)


##
## BUFFER's
##
my $amnm_buffer;
my $pidiv_buffer;
my $thetadiv_buffer;
my $ddiv_buffer;
my $thetadiv_buffer_sqr;
my $amnm_buffer_sqr;
##
## BUFFER's
##


sub get_D_calculator
{
    my $thetacalc=get_theta_calculator();
    my $picalc=get_pi_calculator();
    my $dbuffer=get_ddiv_buffer();
    
    return sub
    {
        my $b=shift;
        my $n=shift;
        my $snp=shift;
	my $averagethetawindow=shift;
        
        my $pi=$picalc->($b,$n,$snp);
        my $theta=$thetacalc->($b,$n,$snp);
        
        
        my $above=$pi-$theta;
        my $below=$averagethetawindow*$dbuffer->($b,$n,$snp->{eucov});
        
        return 0 if abs($above) < 0.000000000001;
        if ($below==0)
        {
             warn "SNP: at $snp->{chr} $snp->{pos} has almost zero variance! -> Infinite Tajima's D; Details pi-theta: $above; variance: $below; SNP will be ignored\n";
             return 0;
        }
        return 0 if $below==0;
        $below=sqrt($below);
        return ($above/$below);
    }
    
}





sub get_theta_calculator
{
    
    my $thetadb=get_thetadiv_buffer();
    return sub
    {
        my $b=shift;
        my $n=shift;
        my $snp=shift;
        
        my $theta_snp=1 / ($thetadb->($b,$n,$snp->{eucov}));
        return $theta_snp;
    }
}



sub get_pi_calculator
{
    
    my $pidb=get_pidiv_buffer();
    return sub
    {
        my $b=shift;
        my $n=shift;
        my $snp=shift;
        
        my $M=$snp->{eucov};
        my $pi_snp=1;
#        print $M, "\n";
        $pi_snp-=($snp->{A}/$M)**2;
        $pi_snp-=($snp->{T}/$M)**2;
        $pi_snp-=($snp->{C}/$M)**2;
        $pi_snp-=($snp->{G}/$M)**2;
        $pi_snp*=$M/($M-1);
        
        $pi_snp/=$pidb->($b,$n,$M);
        return $pi_snp;
    }
}




sub get_ddiv_buffer
{
    #Yields the sum of the pi and theta deviations!!
    $ddiv_buffer={} unless $ddiv_buffer;
    
    my $amnm_buf=get_aMnm_buffer();
    my $pi_calc=get_pi_calculator();
    my $theta_calc=get_theta_calculator();
    
    # get_theta_calculator($b,$n,$snp)
    # get_pi_calculator($b,$n,$snp)
    
    return sub
    {
      my $b=shift;
      my $n=shift;
      my $M=shift;
      
      my $key="$b:$n:$M";
      return $ddiv_buffer->{$key} if exists($ddiv_buffer->{$key});
     
        my $div=0;
        for my $m ($b..$M-$b)
        {
            my $pibpool=$pi_calc->($b,$n,{eucov=>$M,C=>0,G=>0,A=>$m,T=>$M-$m});
            my $thetabpool=$theta_calc->($b,$n,{eucov=>$M,C=>0,G=>0,A=>$m,T=>$M-$m});
            my $term1=($pibpool-$thetabpool)**2;
            my $term2=$amnm_buf->($M,$n,$m);
            
            my $res=$term1*$term2;
            $div+=$res;
        }
        $ddiv_buffer->{$key}=$div;
        return $div;
    };
}



sub get_thetadiv_buffer
{
    $thetadiv_buffer={} unless $thetadiv_buffer;
    my $amnmcalc=get_aMnm_buffer();
    
    return sub
    {
      my $b=shift;
      my $n=shift;
      my $M=shift;
      
      my $key="$b:$n:$M";
      return $thetadiv_buffer->{$key} if exists($thetadiv_buffer->{$key});
     
        my $div=0;
        for my $m ($b..$M-$b)
        {
            my $term1=$amnmcalc->($M,$n,$m);
            $div+=$term1;
        }
      
      $thetadiv_buffer->{$key}=$div;
      return $div;
    };
}


sub get_pidiv_buffer
{
    $pidiv_buffer= {} unless $pidiv_buffer;
    my $amnmcalc=get_aMnm_buffer();
    
    return sub
    {
        my $b=shift;
        my $n=shift;
        my $M=shift;
        
        my $key="$b:$n:$M";
        return $pidiv_buffer->{$key} if exists($pidiv_buffer->{$key});
        
        # calculate the value
        my $div=0;
        for my $m ($b..$M-$b)
        {
            my $term1=(2*$m*($M-$m))/($M*($M-1));
            $term1*=$amnmcalc->($M,$n,$m);
            $div+=$term1;
        }
    
        $pidiv_buffer->{$key}=$div;
        return $div;
    }
}


sub get_aMnm_buffer
{
    $amnm_buffer={} unless $amnm_buffer;
    
    return sub
    {
        my $M=shift;
        my $n=shift;
        my $m=shift;
        
        # calculate the key
        my $key="$M:$n:$m";
        return $amnm_buffer->{$key} if exists($amnm_buffer->{$key});
        
        my $toret=0;
        foreach my $k (1..$n-1)
        {
            my $t1=binomial_term($M,$n,$m,$k);
            $t1*=(1/$k);
            $toret+=$t1;
        }
        
        #store the value in the buffer and return it 
        $amnm_buffer->{$key}=$toret;
        return $toret;
    } 
}


sub get_aMnm_buffer_sqr
# Calculates square of sum = S(M,b,r,n)^2 and stores the calculated values into buffer $amnm_buffer_sqr
# This value is used for theta^2 correction in HKA test -- in calculation of Var[S].

# Sum is defined in the following way:
#              M-b
#              ---
#              \                           
# S(M,b,r,n) = /    (1/r)*choose{M,i}*(r/n)^i*(1-r/n)^{M-i} 
#              ---
#              i=b
#
{
    $amnm_buffer_sqr={} unless $amnm_buffer_sqr;
    
    return sub
    {
        my $M=shift; # coverage -- eucov
        my $b=shift; # min-count
        my $r=shift; # summation index 1..n-1 in get_thetadiv_buffer_sqr
        my $n=shift; # pool size
        
        # calculate the key
        my $key="$M:$b:$r:$n";
        return $amnm_buffer_sqr->{$key} if exists($amnm_buffer_sqr->{$key});
        
        my $sqr=0;
        my $sumBinom=0;
        my $mM = Math::BigRat->new($M);
        foreach my $i ($b..$M-$b)
        {
        	$sumBinom+= noverk($M,$i) * ($r/$n)**($i) * (1-$r/$n)**($M-$i);
        }
        $sqr= ($sumBinom/$r)**2;
        #store the value in the buffer and return it 
        $amnm_buffer_sqr->{$key}=$sqr;
        return $sqr;
    } 
}


sub get_thetadiv_buffer_sqr
# Calculates sum R(M,b,n) and stores the calculated values into buffer $thetadiv_buffer_sqr
# This value is used for theta correction in HKA test -- in calculation of Var[S].
#
# Sum R is defined in the following way:
#            n-1
#            --- 
#            \
# R(M,b,n)=  /   S(M,b,r,n)^2
#            ---
#            r=1
{
    $thetadiv_buffer_sqr={} unless $thetadiv_buffer_sqr;
    my $amnmcalc=get_aMnm_buffer_sqr();
    
    return sub
    {
      my $b=shift;
      my $n=shift;
      my $M=shift;
      
      my $key="$M:$b:$n";
      return $thetadiv_buffer_sqr->{$key} if exists($thetadiv_buffer_sqr->{$key});
     
      my $sum=0;
      my $n_sub_1 = $n - 1;
      foreach my $r (1..$n-1)
      {
          $sum+=$amnmcalc->($M,$b,$r,$n);
      }
      
      $thetadiv_buffer_sqr->{$key}=$sum;
      return $sum;
    };
}


sub noverk
{
    my $n=shift;
    my $k=shift;
    die "n over k; n has to be larger than zero" unless $n>0;
    die "n over k; k has to be larger than zero" unless $k>0;
    die "$k mus not be larger than $n" if $k>$n;
    
    my @above=(($n-$k+1)..$n);
    my @below=(1..$k);
    
    my $val=1;
    while(@above and @below)
    {
        if($val<1)
        {
            $val*=shift @above;
        }
        else
        {
            $val/=shift @below;
        }
    }
    
    foreach(@above)
    {
        $val*=$_;
    }
    foreach(@below)
    {
        $val/=$_;
    }
    return $val
}

sub binomial_term
{
    my $M=shift; # coverage
    my $n=shift; # pool size
    my $m=shift; # running variable for $b..$M-$b
    my $k=shift; # running variable for 1..$n-1
    
    my $val=noverk($M,$m);
    die "$val is zero for M: $M and m: $m\n" unless $val;
    my $t1=($k/$n)**$m;
    my $t2=(($n-$k)/$n)**($M-$m);
    my $toret=$val*$t1*$t2;
    return $toret;
}

sub a_Mnm
{
    my $M=shift;
    my $n=shift;
    my $m=shift;
    
    my $toret=0;
    foreach my $k (1..$n-1)
    {
        my $t1=binomial_term($M,$n,$m,$k);
        $t1*=(1/$k);
        $toret+=$t1;
    }
    return $toret;
}



1;
