use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use VarianceExactCorrection;
use VarianceUncorrected;
use FormatSNP;
use Pileup;
our $verbose=1;

my $pileupfile="";
my $output="";
my $snpfile="";
my $windowSize= 50000;
my $stepSize=10000;
my $minCount=2;
my $fastqtype="illumina";
my $minQual=20;
my $poolSize=0;
my $help=0;
my $test=0;
my $minCoverage=4;
my $maxCoverage=1000000;
my $minCoveredFraction=0.6;
my $region="";
my $tolerateDeletions=0;
my $uncorrected=0;


# --measure pi --help --pool-size 100 --fastq-type sanger --min-count 1 --min-coverage 4 --max-coverage 400 --min-covered-fraction 0.7 --window-size 100 --step-size 100 --input /Users/robertkofler/dev/testfiles/2R_sim_1000000.pileup --output /Users/robertkofler/dev/testfiles/output/test.pi --snp-output /Users/robertkofler/dev/testfiles/output/test.snps
# --measure theta --pool-size 100 --fastq-type sanger --min-count 1 --min-coverage 2 --max-coverage 400 --min-covered-fraction 0.5 --window-size 100 --step-size 100 --input /Volumes/Volume_4/pablo/Sakton_Merged_MQ20-resorted.pileup --output /Volumes/Volume_4/pablo/sakto-out
#--measure pi --region 2R:150-500 --pool-size 100 --fastq-type sanger --min-count 1 --min-coverage 4 --max-coverage 400 --min-covered-fraction 0.7 --window-size 100 --step-size 100 --input /Users/robertkofler/dev/testfiles/2R_sim_1000000.pileup --output /Users/robertkofler/dev/testfiles/output/test.pi 

my $measure="";


GetOptions(
    "measure=s"         =>\$measure,
    "input=s"           =>\$pileupfile,
    "output=s"          =>\$output,
    "snp-output=s"      =>\$snpfile,
    "fastq-type=s"      =>\$fastqtype,
    "window-size=i"     =>\$windowSize,
    "step-size=i"       =>\$stepSize,
    "min-count=i"       =>\$minCount,
    "min-qual=i"        =>\$minQual,
    "pool-size=i"       =>\$poolSize,
    "min-coverage=i"    =>\$minCoverage,
    "max-coverage=i"    =>\$maxCoverage,
    "min-covered-fraction=f"=>\$minCoveredFraction,
    "region=s"          =>\$region,
    "no-discard-deletions" =>\$tolerateDeletions,
    "dissable-corrections"       =>\$uncorrected,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";


pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $pileupfile;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;
pod2usage(-msg=>"Pool size not provided",-verbose=>1) unless $poolSize;
pod2usage(-msg=>"Min count not provided",-verbose=>1) unless $minCount;
pod2usage(-msg=>"Min quality not valid. Has to be between 0 and 40",-verbose=>1) if $minQual<0 || $minQual > 40;
pod2usage(-msg=>"The minimum coverage hast to be at least two times the minimum count",-verbose=>1) unless $minCoverage >= (2*$minCount);
pod2usage(-msg=>"Measure not provided",-verbose=>1) unless $measure;

# Writing the .param file
my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$pileupfile\n";
print $pfh "Using output\t$output\n";
print $pfh "Using snp-output\t$snpfile\n";
print $pfh "Using fastq-type\t$fastqtype\n";
print $pfh "Using measure\t$measure\n";
print $pfh "Using window-size\t$windowSize\n";
print $pfh "Using step-size\t$stepSize\n";
print $pfh "Using min-count\t$minCount\n";
print $pfh "Using min-qual\t$minQual\n";
print $pfh "Using pool-size\t$poolSize\n";
print $pfh "Using min-coverage\t$minCoverage\n";
print $pfh "Using max-coverage\t$maxCoverage\n";
print $pfh "Using min-covered-fraction\t$minCoveredFraction\n";
print $pfh "Using region\t$region\n";
print $pfh "Using no-discard-deletions\t$tolerateDeletions\n";
print $pfh "Dissable corrections\t$uncorrected\n";
print $pfh "Using help\t$help\n";
print $pfh "Using test\t$test\n";
close $pfh;

# qualencoding,mincount,mincov,maxcov,minqual
my $pp          = get_pileup_parser($fastqtype,$minCount,$minCoverage,$maxCoverage,$minQual,$tolerateDeletions);
my $pileslider  = PileupSlider->new($pileupfile,$windowSize,$stepSize,$pp);

# if the user provided a region a modified version of the pileupslider is used: the pileupregionslider: Decorator Pattern
if($region)
{
    $pileslider=PileupRegionSlider->new($pileslider,$region);
}

my $varianceCalculator=VarianceExactCorrection->new($poolSize,$minCount,$minCoverage,$maxCoverage);
$varianceCalculator=VarianceUncorrected->new($poolSize,$minCount,$minCoverage,$maxCoverage) if $uncorrected;

open my $ofh, ">",$output or die "Could not open output file $output\n";

my $snpwriter;
$snpwriter=get_WindowSNPFormater($snpfile) if $snpfile;


while(my $win=$pileslider->nextWindow())
{
    my $chr=$win->{chr};
    my $middle=$win->{middle};
    my $window=$win->{window};
    my $data=$win->{data};
    my $snpcount=$win->{count_snp};
    my $covercount=$win->{count_covered};
    
    print "Processing window: $chr:$win->{start}-$win->{end}\n";
    # writer snps
    $snpwriter->($win) if $snpwriter;
    
    
    my $coveredFraction=$covercount/$window;
    
    my $snps=[];
    foreach(@$data)
    {
        push @$snps,$_ if $_->{ispuresnp};
    }

    my $meas=0;
    
    $meas=$varianceCalculator->calculate_measure($measure,$snps,$covercount);
    $meas=sprintf("%.9f",$meas);
    
    $meas="na" if $coveredFraction < $minCoveredFraction;
    
    $coveredFraction=sprintf("%.3f",$coveredFraction);
    print $ofh "$chr\t$middle\t$snpcount\t$coveredFraction\t$meas\n";
    
}

close $ofh;
exit;

{
    use warnings;
    use strict;
    package PileupRegionSlider;
    
    sub new
    {
        my $class=shift;
        my $pileupslider=shift;
        my $region=shift;
        
        
        die "Region $region wrong format; has to be chr:start-end" unless $region=~m/^(\w+):(\d+)-(\d+)$/;
        my($chr,$start,$end) = ($1,$2,$3);
        
        my $item=
        {
            ps=>$pileupslider,
            chr=>$chr,
            start=>$start,
            end=>$end
        };
        
        my $self= bless $item, __PACKAGE__;
        $self->_spoolForward();
        
        return $self;
        
        
    }
    
    sub _spoolForward
    {
        my $self=shift;
        my $fh=$self->{ps}{fh};
        my $chr=$self->{chr};
        my $start=$self->{start};
        
        
        while(my $l=<$fh>)
        {
            my ($achr,$apos)=split /\t/,$l;
            
            if($achr eq $chr and $apos>=$start)
            {
                $self->{ps}->_bufferline($l);
                last;
            }
        }
        
        $self->{ps}{lower}=$start;
        $self->{ps}{upper}=$start+$self->{ps}{window};
    }
    

    
    sub nextWindow
    {
        my $self=shift;
        my $win=$self->{ps}->nextWindow();
        
        return undef unless $win;
        return undef unless $win->{chr} eq $self->{chr};
        return undef if $win->{start} >$self->{end};
        return $win;
    }
    
    
}



{
    use warnings;
    use strict;
    package PileupSlider;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use Pileup;
    sub new
    {
        my $class=shift;
        my $file=shift;
        my $window=shift;
        my $step=shift;
        my $pp=shift;
        
        open my $fh,"<$file" or die "Could not open file handle";
        
        return bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            file=>$file,
            fh=>$fh,
            pp=>$pp,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
    }
    
    
    sub nextWindow
    {
        my $self=shift;
        my $pp=$self->{pp};
        
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
            my $e=$pp->($line);
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
        
        
        my $toret=_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
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
    
    
    
    sub _annotateWindow
    {
        my $curwin=shift;
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $window=shift;
        
        my $snps=0;
        my $aboveCoverage=0;
        foreach(@$curwin)
        {
            $snps++ if $_->{ispuresnp};
            $aboveCoverage++ if $_->{iscov};
        }

        return
        {
            chr=>$chr,
            start=>$start,
            end=>$end,
            middle=>int(($end+1+$start)/2),
            count_snp=>$snps,
            count_covered=>$aboveCoverage,
            window=>$window,
            data=>$curwin      
        };
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
    package VarTest;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use VarianceExactCorrection;
    use Test::Variance;
    use Test::TClassicalVariance;
    use Test::PileupParser;
    use Pileup;
    use Test;
    
    
    sub runTests
    {
        run_PileupParserTests();
        testPileupSlider();
        run_classicalVarianceTests();
        run_VarianceTests();
        exit;
    }
    
    
    sub _getPileupSliderForString
    {
        my $str=shift;
        my $window=shift;  # win, step, mincount, mincov, minqual
        my $step=shift;
        my $mincount=shift;
        my $mincov=shift;
        my $minqual=shift;
        my $maxcov=shift;
        
        $maxcov||=1000000;
        
        # qualencoding,mincount,mincov,maxcov,minqual
        my $pp=get_pileup_parser("illumina",$mincount,$mincov,$maxcov,$minqual);
        
        open my $ofh,"<",\$str or die "could not open string filehandle";
        
        my $cr=bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$ofh,
            pp=>$pp,
            curwin=>[],
            buffer=>[]
        },"PileupSlider";
        return $cr;
    }
    

    
    sub testPileupSlider
    {
        my $str;
        my $pilsl;
        my $win;


        $str=
        "2L\t1\tA\t9\tCCCCCAAAA\tTTTTTTTTT\n".
        "2L\t2\tA\t7\tCCCCAAA\tTTTTTTT\n".
        "2L\t3\tA\t5\tGGGGT\tTTTTT\n".
        "2L\t4\tA\t3\tCCA\tTTT\n".
        "2L\t5\tA\t2\tCA\tTT\n".
        "chr1\t3\tA\t7\tCCCCCCC\tTTTTTTT\n".
        "chr2\t2\tA\t7\tCCCAAA\tTTTTTT\n".
        "chr2\t3\tA\t7\tCCCAgGAAcC\tTTTTTTTTTT\n".
        "chr4\t3\tA\t7\tCCCCAAA\tTTTTTTT\n"
        ;

        $pilsl=_getPileupSliderForString($str,3,1,2,4,20);
        $win=$pilsl->nextWindow();
        
        is($win->{chr},"2L","PileupSlider; chromosome ok");
        is($win->{window},3,"PileupSlider; window ok");
        is($win->{count_covered},3,"PileupSlider; coverage count ok");
        is($win->{count_snp},2,"PileupSlider; SNP count is ok");
        is($win->{start},0,"PileupSlider; start position is ok");
        is($win->{end},3,"PileupSlider; end position is ok");
        is($win->{data}[0]{A},4,"PileupSlider; allele count is ok");
        is($win->{data}[0]{C},5,"PileupSlider; allele count is ok");
        is($win->{data}[2]{G},4,"PileupSlider; allele count is ok");
        is($win->{data}[2]{T},0,"PileupSlider; allele count is ok");
        is($win->{data}[2]{issnp},0,"PileupSlider; issnp is ok");
        is($win->{data}[2]{iscov},1,"PileupSlider; iscov is ok");
        
        $win=$pilsl->nextWindow();
        is($win->{window},3,"PileupSlider; window ok");
        is($win->{chr},"2L","PileupSlider; chromosome ok");
        is($win->{count_covered},2,"PileupSlider; coverage count ok");
        is($win->{count_snp},1,"PileupSlider; SNP count is ok");
        is($win->{start},1,"PileupSlider; start position is ok");
        is($win->{end},4,"PileupSlider; end position is ok");
        is($win->{middle},3,"PileupSlider; end position is ok");
        
        is($win->{data}[1]{G},4,"PileupSlider; allele count is ok");
        is($win->{data}[1]{T},0,"PileupSlider; allele count is ok");
        is($win->{data}[2]{C},2,"PileupSlider; allele count is ok");
        is($win->{data}[2]{A},0,"PileupSlider; allele count is ok");
        
        $win=$pilsl->nextWindow();
        is($win->{window},3,"PileupSlider; window ok");
        is($win->{chr},"2L","PileupSlider; chromosome ok");
        is($win->{count_covered},1,"PileupSlider; coverage count ok");
        is($win->{count_snp},0,"PileupSlider; SNP count is ok");
        is($win->{start},2,"PileupSlider; start position is ok");
        is($win->{end},5,"PileupSlider; end position is ok");
        is($win->{data}[0]{G},4,"PileupSlider; allele count is ok");
        is($win->{data}[0]{T},0,"PileupSlider; allele count is ok");
        
        $win=$pilsl->nextWindow();
        is($win->{window},3,"PileupSlider; window ok");
        is($win->{chr},"chr1","PileupSlider; chromosome ok");
        is($win->{count_covered},1,"PileupSlider; coverage count ok");
        is($win->{count_snp},0,"PileupSlider; SNP count is ok");
        is($win->{start},0,"PileupSlider; start position is ok");
        is($win->{end},3,"PileupSlider; end position is ok");
        is($win->{data}[0]{C},7,"PileupSlider; allele count is ok");
        is($win->{data}[0]{T},0,"PileupSlider; allele count is ok");
        is($win->{data}[0]{issnp},0,"PileupSlider; issnp is ok");
        is($win->{data}[0]{iscov},1,"PileupSlider; iscov is ok");
        
        $win=$pilsl->nextWindow();
        is($win->{chr},"chr2","PileupSlider; chromosome ok");
        is($win->{count_covered},2,"PileupSlider; coverage count ok");
        is($win->{count_snp},2,"PileupSlider; SNP count is ok");
        is($win->{start},0,"PileupSlider; start position is ok");
        is($win->{end},3,"PileupSlider; end position is ok");
        is($win->{data}[1]{C},5,"PileupSlider; allele count is ok");
        is($win->{data}[1]{T},0,"PileupSlider; allele count is ok");
        is($win->{data}[1]{G},2,"PileupSlider; allele count is ok");
        is($win->{data}[1]{A},3,"PileupSlider; allele count is ok");
        
        $win=$pilsl->nextWindow();
        is($win->{chr},"chr4","PileupSlider; chromosome ok");
        is($win->{count_covered},1,"PileupSlider; coverage count ok");
        is($win->{count_snp},1,"PileupSlider; SNP count is ok");
        
        $win=$pilsl->nextWindow();
        not_exists($win,"PileupSlider; final window is emtpy");
        $win=$pilsl->nextWindow();
        not_exists($win,"PileupSlider; final window is still emtpy");
        
        # weird cases; * = deletion and N
        $str=
        "2L\t1\tA\t9\tCCCCCAAAA\tTTTTTTTTT\n".
        "2L\t2\tA\t8\tCCCC*AAA\tTTTTTTTT\n".
        "2L\t3\tA\t5\tGGGTTNNc\tTTTTTTTT\n".
        "2L\t4\tA\t9\tCCCCCTTTTN\tTTTTTTTTTT\n"
        ;
        $pilsl=_getPileupSliderForString($str,3,1,2,4,20);  # win, step, mincount, mincov, minqual
        $win=$pilsl->nextWindow();
        
        is($win->{chr},"2L","PileupSlider; chromosome ok");
        is($win->{window},3,"PileupSlider; window ok");
        is($win->{count_covered},3,"PileupSlider; coverage count ok");
        is($win->{count_snp},2,"PileupSlider; SNP count is ok");
        is($win->{start},0,"PileupSlider; start position is ok");
        is($win->{end},3,"PileupSlider; end position is ok");
        is($win->{data}[0]{del},0,"PileupSlider; deletion count is ok");
        is($win->{data}[1]{del},1,"PileupSlider; deletion count is ok");
        is($win->{data}[2]{del},0,"PileupSlider; deletion count is ok");
        is($win->{data}[0]{N},0,"PileupSlider; N count is ok");
        is($win->{data}[1]{N},0,"PileupSlider; N count is ok");
        is($win->{data}[2]{N},2,"PileupSlider; N count is ok");
        is($win->{data}[0]{issnp},1,"PileupSlider; snp-identification is ok");
        is($win->{data}[1]{issnp},1,"PileupSlider; snp-identificationis ok");
        is($win->{data}[2]{issnp},1,"PileupSlider; snp-identification is ok");
        is($win->{data}[0]{ispuresnp},1,"PileupSlider; pure snp-identification is ok");
        is($win->{data}[1]{ispuresnp},0,"PileupSlider; pure snp-identification is ok");
        is($win->{data}[2]{ispuresnp},1,"PileupSlider; pure snp-identification is ok");
        
        
        $win=$pilsl->nextWindow();
        is($win->{chr},"2L","PileupSlider; chromosome ok");
        is($win->{window},3,"PileupSlider; window ok");
        is($win->{count_covered},3,"PileupSlider; coverage count ok");
        is($win->{count_snp},2,"PileupSlider; SNP count is ok");
        is($win->{start},1,"PileupSlider; start position is ok");
        is($win->{end},4,"PileupSlider; end position is ok");
        is($win->{data}[0]{N},0,"PileupSlider; N count is ok");
        is($win->{data}[1]{N},2,"PileupSlider; N count is ok");
        is($win->{data}[2]{N},1,"PileupSlider; N count is ok");
        is($win->{data}[0]{eucov},7,"PileupSlider; eu-coverage is ok");
        is($win->{data}[1]{eucov},5,"PileupSlider; eu-coverage is ok");
        is($win->{data}[2]{eucov},9,"PileupSlider; eu-coverage is ok");
        is($win->{data}[0]{totcov},8,"PileupSlider; total coverage is ok");
        is($win->{data}[1]{totcov},8,"PileupSlider; total coverage is ok");
        is($win->{data}[2]{totcov},10,"PileupSlider; total coverage is ok");
        
    }
}



=head1 NAME

perl Variance-sliding.pl - A script which calculates the requested population genetics measure (pi, theta, d) along chromosomes using a sliding window approach.

=head1 SYNOPSIS

perl Variance-sliding.pl --measure pi --input input.pileup --output output.file --pool-size 500 --min-count 2 --min-coverage 4 --window-size 50000 --step-size 10000


=head1 OPTIONS

=over 4

=item B<--input>

The input file in the pileup format. A pooled population sequenced and mapped to the reference genome. Finally the mapping results have been converted to sam output format.
Using the samtools the sam format can be easily converted into the pileup format.  Mandatory.

=item B<--output>

The output file.  Mandatory.

=item B<--snp-output>

If provided, this file will contain the polymorphic sites which have been used for this analysis; default=empty

=item B<--measure>

Currently, "pi", "theta" and "D" is supported. This stands for Tajima's Pi, Watterson's Theta, Tajima's D respectively; Mandatory

=item B<--region>

Provide a subregion in which the measure should be calculated using a sliding window approach; has to be of the form chr:start-end; eg 2L:10000-20000; Example: C<--region "chr2:100-1000"> Optional

=item B<--pool-size>

The size of the pool which has been sequenced. e.g.: 500; Mandatory

=item B<--fastq-type>
The encoding of the quality characters; Must either be 'sanger' or 'illumina'; 

 Using the notation suggested by Cock et al (2009) the following applies:
 'sanger'   = fastq-sanger: phred encoding; offset of 33
 'solexa'   = fastq-solexa: -> NOT SUPPORTED
 'illumina' = fastq-illumina: phred encoding: offset of 64
 
 See also:
 Cock et al (2009) The Sanger FASTQ file format for sequecnes with quality socres,
 and the Solexa/Illumina FASTQ variants; 

default=illumina
 
=item B<--min-count>

The minimum count of the minor allele. This is important for the identification of SNPs; default=2

=item B<--min-coverage>

The minimum coverage of a site. Sites with a lower coverage will not be considered (for SNP identification and coverage estimation); default=4

=item B<--max-coverage>

the maximum coverage; used for SNP identification, the coverage in ALL populations has to be lower or equal to this threshold, otherwise no SNP will be called. default=1000000

=item B<--min-covered-fraction>

the minimum fraction of a window being between min-coverage and max-coverage in ALL populations; float; default=0.6

=item B<--min-qual>

The minimum quality; Alleles with a quality lower than this threshold will not be considered (for SNP identification and coverage estimation); default=20

=item B<--window-size>

The size of the sliding window. default=50000

=item B<--step-size>

the size of one sliding window step. If this number is equal to the --window-size the sliding window will be non overlapping (jumping window). default=10000

=item B<--no-discard-deletions>

per default sites with already a single deletion are discarded. By setting this flag, sites with deletions will be used.

=item B<--dissable-corrections>

Flag; Dissable correction factors; Calculates Pi/Theta/Tajima's D in the classical way not taking into account pooling or missing minor alleles; default=off

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

A pileup file as described here: http://samtools.sourceforge.net/pileup.shtml; example:

 2L	90131	N	11	AaAAAaaAaAA	[aUQ_a`^_\Z
 2L	90132	N	11	AaAAAaaAaAA	_bYQ_^aaT^b
 2L	90133	N	11	A$aAAAaaAaAA	_b[Xaaa__Ua
 2L	90134	N	10	tTTTttTtTT	_`aaa_a[aa
 2L	90135	N	10	a$TAAaaAaAA	aZ^a`ba`\_
 2L	90136	N	9	TTTttTtTT	`aaaaaWaa
 2L	90137	N	9	GGGggGgGG	``aaaaQaa
 2L	90138	N	9	T$TTttTtTT	[U\`a\T^_
 2L	90139	N	8	TTttTtTT	``aaaU_a
 2L	90140	N	9	CCccCcCC^FC	[aaba`aaa

=head2 Output

For every sliding window the population genetics measure will be displayed in the following format

 2L	1730000	557	0.726	0.005647899
 2L	1740000	599	0.777	0.005657701
 2L	1750000	650	0.767	0.006129462
 2L	1760000	617	0.703	0.006265200
 2L	1770000	599	0.672	0.006427032

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: number of SNPs in the sliding window; These SNPs have been used to calculate the value in col 5
 col 4: fraction of the window covered by a sufficient number of reads. Suficient means higher than min-coverage and lower than max-coverage
 col 5: population genetics estimator (pi, theta, D)

=head2 SNP Output

This output file is optional, it will only be created if the option  C<--snp-output> is provided.
For example:

 >2R:50; 2R:0-100; snps:3
 2R      13      A       17      10      0       0       7       0
 2R      30      A       34      33      1       0       0       0
 2R      37      A       40      39      1       0       0       0

 >2R:150; 2R:100-200; snps:2
 2R      105     A       76      64      0       12      0       0
 2R      118     T       75      0       72      1       2       0

 The header contains the following information
 >chr:middle chr:start-end snps:scnpcount
 chr..chromosome (contig)
 middle.. middle position of the window
 start.. start position of the window
 end.. end position of the window
 snpcount.. number of snps found in the given window

 The individual tab-delimited entries are in the following format:
 col 1: chromosome ID (contig)
 col 2: position in chromosome
 col 3: reference character
 col 4: coverage
 col 5: counts of A's
 col 6: counts of T's
 col 7: counts of C's
 col 8: counts of G's
 col 9: counts of N's


=head2 Technical detail

The pileup file may contain deletions which can be recognised by the symbol '*'. SNP sites which have a deletion are ignored.
The sum of the bases A,T,C,G is used as coverage (ignoring N's).

The script will be slow at the beginning but speed up steadily. Calculating the correction factors for Tajima's Pi, Watterson's Theta and Tajima's D is computationally intense.
However, it needs only to be done a single time. The once calculated correction factors will be stored in hashes where they can be reused when required, this speeds up computation. 
 
 
=cut
