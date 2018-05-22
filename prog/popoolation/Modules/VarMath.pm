package VarMath;
use strict;
use warnings;

use Math::BigRat; 

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT=qw(noverk get_theta_calculator get_pi_calculator get_D_calculator);
our @EXPORT_OK = qw(get_thetadiv_buffer get_thetadiv_buffer_sqr get_nbase_buffer);


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
my $thetadiv_buffer_sqr;
my $amnm_buffer_sqr;

## Buffer for Nicolas novel Tajima's D
my $an_buffer;
my $bn_buffer;
my $astar_buffer;
my $bstar_buffer;
my $nbase_buffer;
##
## BUFFER's
##

sub get_an_buffer
{
    $an_buffer={} unless $an_buffer;
    return sub
    {
	my $n=shift;
	return $an_buffer->{$n} if(exists($an_buffer->{$n}));
	
	my $an=0;
	foreach my $i (1..($n-1))
	{
	    $an+=(1/$i);
	}
	$an_buffer->{$n}=$an;
	return $an;
    }
}

sub get_bn_buffer
{
    $bn_buffer={} unless $bn_buffer;
    return sub
    {
	my $n=shift;
	return $bn_buffer->{$n} if(exists($bn_buffer->{$n}));
	
	my $bn=0;
	foreach my $i (1..($n-1))
	{
	    $bn+=(1/($i**2));
	}
	$bn_buffer->{$n}=$bn;
	return $bn;
    }
}



sub get_D_calculator
{
    my $thetacalc=get_theta_calculator();
    my $picalc=get_pi_calculator();
    my $ddivisor=get_ddivisor();
    
    return sub
    {
        my $b=shift;
        my $n=shift;
	my $mincoverage=shift;
        my $snps=shift;
	return 0 unless @$snps;
	
	my $pi		=$picalc->($b,$n,$snps);
	my $theta	=$thetacalc->($b,$n,$snps);
	my $ddivisor	=$ddivisor->($n,$mincoverage,$snps,$theta);
	return 0 if($pi-$theta)==0;
	my $d = ($pi-$theta)/$ddivisor;
	return $d;
    }
    
}

sub get_ddivisor
{
    my $nbase_buffer;
    my $alphastarcalc=get_betastar_calculator();
    my $betastarcalc=get_betastar_calculator();
    return sub
    {
	my $n=shift;
	my $mincoverage=shift;
	my $snps=shift;
	my $theta=shift;
	$nbase_buffer=get_nbase_buffer($n) unless $nbase_buffer;
	
	
	my $snpcount=@$snps;
	
	#my $n_sum=0;
	#foreach my $snp (@$snps)
	#{
	#    my $cov=$snp->{eucov};
	#    my $nbase=$nbase_buffer->($n,$mincoverage);
	#    $n_sum+=$nbase;
	#}
	#
	#my $averagen=int($n_sum/$snpcount);
	
	my $averagen=$nbase_buffer->($n,$mincoverage);
	my $alphastar=$alphastarcalc->($averagen);
	my $betastar=$betastarcalc->($averagen);
	my $div=($alphastar/$snpcount)*$theta + $betastar* ($theta**2);
	return sqrt($div);
    }
}



sub get_nbase_buffer
{
    my $poolsize=shift;
    $nbase_buffer={} unless($nbase_buffer);
    my $pijresolver=get_nbase_matrix_resolver(3*$poolsize,$poolsize);
    return sub
    {

	my $cov=shift;

	#### shortcut
	#die"Poolsize has to be larger than coverage; maybe decreasea maximum coverage" unless $n> $cov;
	#return $cov;
	### end shortcut
	
	my $key="$poolsize:$cov";
	return $nbase_buffer->{$key} if(exists($nbase_buffer->{$key}));
	
	my $nbase=0;
	my $minj=$cov<$poolsize?$cov:$poolsize;
	
	for my $k(1..$minj)
	{
	    $nbase+=$k*$pijresolver->($cov,$k);
	}
	
	$nbase_buffer->{$key}=$nbase;
	return $nbase_buffer->{$key};
    }
}


sub _get_pij_matrix
{
    my $maxcoverage=shift;
    my $poolsize=shift;
    
    my $jboundary=$maxcoverage < $poolsize ? $maxcoverage : $poolsize;
    
    my $matrix=[];
    for my $i(1..$maxcoverage)
    {
	$matrix->[$i][0]=0;
    }
    for my $j(1..$jboundary)
    {
	$matrix->[0][$j]=0;
    }
    $matrix->[0][0]=1;
    
    for my $i(1..$maxcoverage)
    {
	for my $j (1..$jboundary)
	{
	    my $t1= ((1+$poolsize-$j)/$poolsize)*($matrix->[$i-1][$j-1]);
	    my $t2=($j/$poolsize)*($matrix->[$i-1][$j]);
	    my $pij=$t1+$t2;
	    $matrix->[$i][$j]=$pij;
	}
    }
    return $matrix;
}

sub get_nbase_matrix_resolver
{
    my $maxcoverage=shift;
    my $poolsize=shift;
    my $matrix=_get_pij_matrix($maxcoverage,$poolsize);
    
    return sub
    {
	my $C=shift;
	my $k=shift;
	unless(exists($matrix->[$C][$k]))
	{
	    $matrix=_get_pij_matrix(3*$C,$poolsize);
	}
	return $matrix->[$C][$k];
    }
}


sub get_alphastar_calculator
{
    my $anb=get_an_buffer();
    $astar_buffer= {} unless $astar_buffer;
    return sub
    {
	my $n=shift;
	die "invalid effective coverage; has to be larger than 1" unless $n>1;
	
	return $astar_buffer->{$n} if(exists($astar_buffer->{$n}));
	
	my $an=$anb->($n);
	my $fs=calculate_fstar($an,$n);
	# calculate individual terms(t) and subterms(st)
	my $t1=($fs**2)*($an-($n/($n-1)));
	my $st1=$an * ( (4*($n+1)) / (($n-1)**2) );
	my $st2=2 * (($n+3)/($n-1));
	my $t2=$fs * ($st1-$st2);
	my $t3=$an * ( (8*($n+1))/($n*(($n-1)**2)) );
	my $t4= (($n**2)+$n+60)/(3*$n*($n-1));
	my $astar= ($t1 + $t2 - $t3 + $t4);
	$astar_buffer->{$n}=$astar;
	return $astar;
    }
}

sub get_betastar_calculator
{
    my $anb=get_an_buffer();
    my $bnb=get_bn_buffer();
    $bstar_buffer={} unless $bstar_buffer;
    return sub
    {
	my $n=shift;
	die "invalid effective coverage; has to be larger than 1" unless $n>1;
	return $bstar_buffer->{$n} if(exists($bstar_buffer->{$n}));
	my $an=$anb->($n);
	my $bn=$bnb->($n);
	my $fs=calculate_fstar($an,$n);
	
	my $t1 = ($fs**2) * ($bn - ((2*($n-1)) /(($n-1)**2)));
	my $st1= $bn * (8/($n-1));
	my $st2= $an * (4/($n*($n-1)));
	my $st3= (($n**3)+12*($n**2)-35*$n+18)/($n*(($n-1)**2));
	my $t2 = $fs*($st1-$st2-$st3);
	my $t3 = $bn * (16/($n*($n-1)));
	my $t4 = $an * (8/(($n**2)*($n-1)));
	my $st4= 2*($n**4+ 110*($n**2)-255*$n+126);
	my $st5= 9*($n**2)*(($n-1)**2);
	my $t5 = $st4/$st5;
	my $bstar= ($t1 + $t2 - $t3 + $t4 + $t5);
	$bstar_buffer->{$n}=$bstar;
	return $bstar;
    }
}


sub calculate_fstar
{
    my $an=shift;
    my $n=shift;
    return (($n - 3)/($an*($n - 1) - $n));
}



sub get_theta_calculator
{
    
    my $thetadb=get_thetadiv_buffer();
    return sub
    {
        my $b=shift;
        my $n=shift;
        my $snps=shift;
	my $thetasum=0;
	foreach my $snp (@$snps)
	{
	    $thetasum+=1 / ($thetadb->($b,$n,$snp->{eucov}));
	}
        return $thetasum;
    }
}



sub get_pi_calculator
{
    
    my $pidb=get_pidiv_buffer();
    return sub
    {
        my $b=shift;
        my $n=shift;
        my $snps=shift;
	my $pisum=0;
        foreach my $snp (@$snps)
	{
	    my $M=$snp->{eucov};
	    my $pi_snp=1;
    #        print $M, "\n";
	    $pi_snp-=($snp->{A}/$M)**2;
	    $pi_snp-=($snp->{T}/$M)**2;
	    $pi_snp-=($snp->{C}/$M)**2;
	    $pi_snp-=($snp->{G}/$M)**2;
	    $pi_snp*=$M/($M-1);
	    
	    $pi_snp/=$pidb->($b,$n,$M);
	    $pisum+=$pi_snp;
	}

        return $pisum;
    }
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
