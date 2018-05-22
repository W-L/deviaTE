#!/usr/bin/perl-w
use strict;
use warnings;
use Getopt::Long; # to get options
use Pod::Usage; # to show documentations
use File::Basename; # to get the file path, file name and file extension

# Author: Ram Vinay Pandey
# Modified on 12-04-2010 to implement uniform Window sliding module.

# Define the variables
my $input;
my $output="";
my $help=0;
my $test=0;
my $verbose=1;

my $windowsize=1000;
my $step=100;
my $minCoverageFraction=0.6;


my $usage="perl $0 --input mauve-parsing-putput.txt --output dxy-output.txt --window-size 1000 --step-size 100 --min-covered-fraction 0.6\n";

GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "window-size=i"  =>\$windowsize,
    "step-size=i"   =>\$step,
    "min-covered-fraction=f"=>\$minCoverageFraction,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
Test::runTests() if $test;
#die "Tests have not yet been implemented" if $test;

   
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using window-size\t$windowsize\n";
print $pfh "Using step-size\t$step\n";
print $pfh "Using min-covered-fraction\t$minCoverageFraction\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;


      

open my $ofh, ">$output" or die "Could not open output file";


my $reader;

$reader=BpSlider->new($input,$windowsize,$step);

# population count
my $speciesCount=0;

$speciesCount=$reader->count_samples();

if($speciesCount !=3) {
	
	pod2usage(-msg=>"Only 3 species are allowed to calculate divergence\n",-verbose=>1);
}

print "\n\nCalculating Dxy started at\t". localtime() ." ....\n";

#exit;


while(my $window=$reader->nextWindow())
{
        my $chr=$window->{chr};
        my $pos=$window->{middle};
        my $win=$window->{window};
        my $above=$window->{count_covered};
        my $data=$window->{data};
	
        next unless @$data;
        
        my $coveredFrac=$above/$win;

        if ($coveredFrac>=$minCoverageFraction) {
	    
	    my($string)=Utility::calculateDxy($data,$pos,$chr,$win,$minCoverageFraction);
	    #print "$coveredFrac\t$minCoverageFraction\n";
	    print $ofh "$string\n";
            #print "$string\n";
	    
	}
        else {
		
		my $considered_pos = 0;
		foreach my $d (@$data)
		{
		    ### if condition to discard the position with gap (-) in any of 3 species
		    if (("$d->{sp1}" ne "-") and ("$d->{sp2}" ne "-") and ("$d->{sp3}" ne "-")) {
			if (($d->{sp1} =~ m/N/i) and ($d->{sp2} =~ m/N/i) and ($d->{sp3} =~ m/N/i)) {
				$considered_pos++;
			}
		    }
		}
		my $coveredFrac="0.00";
		if ($considered_pos>0) {
		    $coveredFrac=$considered_pos/$win;
		    $coveredFrac = sprintf "%.2f",$coveredFrac;
		}
		
		print $ofh "$chr\t$pos\tna\tna\tna\t$coveredFrac\n";
		#print "$chr\t$pos\tna\tna\tna\t$considered_pos\n";
        }
 
}


print "\n\nCalculating Dxy completed at\t". localtime() ." ....\n";


exit;



{
    use warnings;
    use strict;
    package BpSlider;

    sub new
    {
        my $class=shift;
        my $file=shift;
        my $window=shift;
        my $step=shift;
        
        open my $fh,"<$file" or die "Could not open file handle";
        
        return bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$fh,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
    }
    
    sub count_samples
    {
        my $self=shift;
        my $l=$self->_nextline();
        my $p=Utility::_parseLightwight($l);
        my $c=scalar(@{$p->{species}});
        $self->_bufferline($l);
        return $c;
    }
    
    sub nextWindow
    {
        my $self=shift;

	
        #get the current window, and the current chromosome
        my $curwin=$self->{curwin};
        
        my $curChr="";
        $curChr=$curwin->[0]{chr} if @$curwin;
        
        my $resetchr=0;
        
        # empty unnecessary entries
        EMPTY: while(@$curwin)
        {
            my $e=shift @$curwin;
            if($e->{pos}>$self->{lower})
            {
                unshift @$curwin, $e;
                last EMPTY;
            }
            
        }
        
        # fill with novel entries
        my $line;
        FILL:while($line=$self->_nextline)
        {
            my $e=Utility::_parseLightwight($line);
            $curChr=$e->{chr} unless $curChr;
            
            
            if($e->{chr} eq $curChr && $e->{pos} <= $self->{upper})
            {
                push @$curwin,$e;
            }
            else
            {
                $resetchr=1 if $e->{chr} ne $curChr;
                $self->_bufferline($line);
                last FILL;
            }
        }
        
        return undef unless $curChr;
        
        
        my $toret=Utility::_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
        if($resetchr or not defined($line))
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
    
    
    
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
}



{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    
      
    sub calculateDxy {
	
        my $data=shift;
	my $pos=shift;
	my $chr=shift;
        my $win=shift;
	my $minCovFraction=shift;
	
        my ($Dxy12,$Dxy13,$Dxy23) = (0,0,0);
        my ($ct12,$ct13,$ct23,$array_size) = (0,0,0,0);
        $array_size = @$data;
        
	my $considered_pos = 0;
        foreach my $d (@$data)
        {
            #print "$d->{chr}\t$d->{pos}\t$d->{astate}\t$d->{dstate}\t$d->{sp1}\t$d->{sp2}\t$d->{sp3}\n";
	    ### if condition to discard the position with gap (-) in any of 3 species
	    if (("$d->{sp1}" ne "-") and ("$d->{sp2}" ne "-") and ("$d->{sp3}" ne "-")) {
		
		if (($d->{sp1} !~ m/N/i) and ($d->{sp2} !~ m/N/i) and ($d->{sp3} !~ m/N/i)) {
		
			$considered_pos++;
			
			if ("$d->{sp1}" ne "$d->{sp2}") {
			    $ct12++;
			}
			
			if ("$d->{sp1}" ne "$d->{sp3}") {
			    $ct13++;
			}
			
			if ("$d->{sp2}" ne "$d->{sp3}") {
			    $ct23++;
			}
		}
	    }
            
        }
        
        # Calculate Dxy	    
        #$Dxy12 = $ct12/$array_size;
        #$Dxy13 = $ct13/$array_size;
        #$Dxy23 = $ct23/$array_size;
	if(($considered_pos>0) and ($ct12>0)) {
		$Dxy12 = $ct12/$considered_pos;
		$Dxy12 = sprintf "%.10f",$Dxy12;
	}
	else {
		$Dxy12 = "na";
	}
	
	if(($considered_pos>0) and ($ct13>0)) {
		$Dxy13 = $ct13/$considered_pos;
		$Dxy13 = sprintf "%.10f",$Dxy13;
	}
	else {
		$Dxy13 = "na";
	}
	
	if(($considered_pos>0) and ($ct23>0)) {
		$Dxy23 = $ct23/$considered_pos;
		$Dxy23 = sprintf "%.10f",$Dxy23;
	}
	else {
		$Dxy23 = "na";
	}
	
        
        
        
	
	my $string = "";

	my $coveredFrac="0.00";
	if ($considered_pos>0) {
	    $coveredFrac=$considered_pos/$win;
	    $coveredFrac = sprintf "%.2f",$coveredFrac;
	}
		
	if ($coveredFrac>=$minCovFraction) {
		
		$string = "$chr\t$pos\t$Dxy12\t$Dxy13\t$Dxy23\t$coveredFrac";
	}
	else {
		$string = "$chr\t$pos\tna\tna\tna\t$coveredFrac";
	}
        #print $ofh "$chr\t$middle\t$Dxy12\t$Dxy13\t$Dxy23\n";
        
        
        return $string;

    }
    
    sub _annotateWindow
    {
        my $curwin=shift;
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $window=shift;

        #my $snps=0;
        my $aboveCoverage=0;
        foreach(@$curwin)
        {
            #$snps++ if $_->{ispuresnp};
            $aboveCoverage++ if $_->{iscov};
        }

        return
        {
            chr=>$chr,
            start=>$start,
            end=>$end,
            middle=>int(($end+1+$start)/2),
            count_covered=>$aboveCoverage,
            window=>$window,
            data=>$curwin      
        };
    }
    
    
    
    sub _parseLightwight
    {
        my $line=shift;
	
        chomp $line;
        my @a=split /\s+/,$line;
        my $chr=shift @a;
        my $pos=shift @a;
	
	my $en={};

	if (($a[0] eq "-") or ($a[1] eq "-") or ($a[2] eq "-")) {

	    $en={
		chr=>$chr,
		pos=>$pos,
		iscov=>0,
		sp1=>$a[0],
		sp2=>$a[1],
		sp3=>$a[2],
		species=>\@a
	    };
        }

	else {
	    
	    $en={
		chr=>$chr,
		pos=>$pos,
		iscov=>1,
		sp1=>$a[0],
		sp2=>$a[1],
		sp3=>$a[2],
		species=>\@a
	    };
	    
	}
     
        #$en->{ispuresnp} = ($en->{issnp} and not $taintedsnp)?1:0;
        
        return $en;
    }
    
    
}





{
    package Test;
    use strict;
    use warnings;
    #use Test::More;
    #use Test::Simple;
    
    my $testcounter;
    
    sub _getBPSliderForString
    {
        my $str=shift;
        my $window=shift;
        my $step=shift;
        open my $ofh,"<",\$str or die "could not open string filehandle";
        my $cr=bless {   
	    lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$ofh,
            curwin=>[],
            buffer=>[]
	    
        },"BpSlider";
        return $cr;
    }
    
    sub runTests
    {
	$testcounter=1;
	testCalculateDxy();
	exit;
    }
    
    

    sub testCalculateDxy
    {
	
	my $teststr=
	"2L\t1\tA\t-\tA\n".
        "2L\t2\tG\tC\tG\n".
	"2L\t3\tG\t-\tG\n".
	"2L\t4\tA\tG\tA\n".
	"2L\t5\tC\tC\tG\n".
	"2L\t6\tT\t-\tT\n".
	"2L\t7\tA\tC\tA\n".
	"2L\t8\tC\tN\tC\n".
	"2L\t9\tC\tC\tC\n".
	"2L\t10\tA\tG\tA\n"; #10
	
	my $minCoverageFraction = 0.6;
	my $testfile = "test/dxytest.txt";
        my $bpsl=_getBPSliderForString($teststr,10,10);
	my $string="";
		
	#2L	249055	6	0.8333333333	na	0.8333333333
        while( my $window=$bpsl->nextWindow() ) {

		my $chr=$window->{chr};
		my $pos=$window->{middle};
		my $win=$window->{window};
		my $above=$window->{count_covered};
		my $data=$window->{data};
		
		next unless @$data;
		
		#print "TEST: $chr\t$pos\t$win\t$above\@$data\n";
		
		my $dxy = Test::calculateDxy($data,$pos,$chr,$win,$minCoverageFraction,$above);
		
		#print "$dxy->{chr}\t$dxy->{pos}\t$dxy->{considered_pos}\t$dxy->{D12}\t$dxy->{D13}\t$dxy->{D23}\n";

		is($dxy->{chr},"2L","test CalculateDxy, correct chromosome");
		is($dxy->{pos},5,"test CalculateDxy, correct window position");
		is($dxy->{considered_pos},6,"test CalculateDxy, correct number of considered base pairs to calculate Dxy");
		is($dxy->{D12},"0.6666666667","test CalculateDxy, the Dxy between species 1 and 2 is correct");
		is($dxy->{D13},"0.1666666667","test CalculateDxy, the Dxy between species 1 and 3 is correct");
		is($dxy->{D23},"0.8333333333","test CalculateDxy, the Dxy between species 2 and 3 is correct");
	
	}
  
    }
    
    
    sub calculateDxy {
	
	my ($data,$pos,$chr,$win,$minCoverageFraction,$above) = (shift,shift,shift,shift,shift,shift);
	
	my $string = "";
	my $coveredFrac=$above/$win;
	my $dxy = {};
	
        if ($coveredFrac>=$minCoverageFraction) {
	    
	    $string = Utility::calculateDxy($data,$pos,$chr,$win,$minCoverageFraction);
	    my ($chr,$pos,$D12,$D13,$D23,$considered_pos) = split("\t",$string);
	    $dxy = {
		chr=>$chr,
		pos=>$pos,
		D12=>$D12,
		D13=>$D13,
		D23=>$D23,
		considered_pos=>$considered_pos
	    };
	    
	}
        else {
		
		my $considered_pos = 0;
		foreach my $d (@$data)
		{
		    ### if condition to discard the position with gap (-) in any of 3 species
		    if (("$d->{sp1}" ne "-") and ("$d->{sp2}" ne "-") and ("$d->{sp3}" ne "-")) {
			if (($d->{sp1} =~ m/N/i) and ($d->{sp2} =~ m/N/i) and ($d->{sp3} =~ m/N/i)) {
				$considered_pos++;
			}
		    }
		}
		$string = "$chr\t$pos\tna\tna\tna\t$considered_pos";
		
		$dxy = {
		chr=>$chr,
		pos=>$pos,
		D12=>"na",
		D13=>"na",
		D23=>"na",
		considered_pos=>$considered_pos
	    };
        }
	
	#return $string;
	return $dxy;
    }

    sub is
    {

        my $a=shift;
        my $b=shift;
        my $msg=shift;
	
	
        if($a eq $b)
        {
            print "$testcounter: OK $a = $b; $msg\n";
        }
        else
        {
            print "$testcounter: FAILED $a = $b; $msg\n";
        }
        $testcounter++;
    }
    
    sub not_exists
    {
        my $a=shift;
        my $msg=shift;
        if($a)
        {
            print "$testcounter FAILED; $msg\n";
        }
        else
        {
            print "$testcounter OK: $msg\n";
        }
        $testcounter++;
    }

    sub ok
    {
        my $a=shift;
        my $msg=shift;
        if($a)
        {
            print "$testcounter OK: $msg\n";
        }
        else
        {
            print "$testcounter FAILED; $msg\n";
        }
        $testcounter++;
    }
    
    
}



    
=head1 NAME

calculate-dxy.pl - Calculates the distance between 2 species within a given window in pairwise comparison.

=head1 SYNOPSIS

 perl calculate-dxy.pl --input mauve-parser-output.txt --output Dxy.txt --window-size 1000 --step-size 100

=head1 OPTIONS

=over 4

=item B<--input>

The input file has to be mauve-parser.pl program output file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--window-size>

the size of the sliding window; default=1000

=item B<--step-size>

the size of the sliding window steps; default=100

=item B<--min-covered-fraction>

the minimum fraction of a window being between min-coverage and max-coverage in ALL populations; float; default=0.6

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

A mauve-parser.pl program output file; example:

 2L	5783	C	G	G
 2L	5784	C	C	C
 2L	5785	C	T	G

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: allelic state in species 1 (reference species)
 col 4: allelic state in species 2
 col 5: allelic state in species 3 (outgroup species)
 
=head2 Output

An output of this program looks like in the given example:

 2L	5788	0.50	0.50	0.30
 2L	5790	0.60	0.50	0.40
 2L	5791	0.44	0.10	0.33
 2L	5792	0.43	0.12	0.43

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: D12 (Distance between species1 and species2) within a given window
 col 4: D13 (Distance between species1 and species3) within a given window
 col 5: D23 (Distance between species2 and species3) within a given window
 Note: species 1 refers to reference species; species 3 refers to outgroup species in phylogenetic tree; Example: species 1 is D. melanogaster, species2 is D. simulanes and species 3 is D. yakuba

=head1 AUTHORS

Ram vinay pandey

Robert Kofler

Pablo Orozco terWengel

Christian Schloetterer

=cut
