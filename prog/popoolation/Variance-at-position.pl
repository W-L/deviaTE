use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use FormatSNP;
use VarianceExactCorrection;
use VarianceUncorrected;
use Pileup;


our $verbose=1;

my $pileupfile="";
my $gtffile="";
my $output="";
my $snpfile="";
my $minCount=2;
my $fastqtype="illumina";
my $minQual=20;
my $poolSize=0;
my $help=0;
my $test=0;
my $minCoverage=4;
my $maxCoverage=1000000;
my $makeNoise=100000;
my $tolerateDeletions=0;
my $uncorrected=0;

my $minCoveredFraction=0.6;


# --gtf /Volumes/Volume_4/analysis/genewise-fst/fb1000.gtf --pool-size 197  --measure pi --pileup /Volumes/Volume_4/analysis/genewise-fst/2Lsubsub.pileup --output /Volumes/Volume_4/analysis/genewise-fst/fem-earlygenes.out

my $measure="";


GetOptions(
    "measure=s"         =>\$measure,
    "pileup=s"          =>\$pileupfile,
    "gtf=s"             =>\$gtffile,
    "output=s"          =>\$output,
    "snp-output=s"      =>\$snpfile,
    "fastq-type=s"      =>\$fastqtype,
    "min-count=i"       =>\$minCount,
    "min-qual=i"        =>\$minQual,
    "pool-size=i"       =>\$poolSize,
    "min-coverage=i"    =>\$minCoverage,
    "max-coverage=i"    =>\$maxCoverage,
    "min-covered-fraction=f"=>\$minCoveredFraction,
    "no-discard-deletions"=>\$tolerateDeletions,
    "dissable-corrections"=>\$uncorrected,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";


pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $pileupfile;
pod2usage(-msg=>"Could not find gtf file",-verbose=>1) unless -e $gtffile;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;
pod2usage(-msg=>"Pool size not provided",-verbose=>1) unless $poolSize;
pod2usage(-msg=>"Min count not provided",-verbose=>1) unless $minCount;
pod2usage(-msg=>"Min quality not valid. Has to be between 0 and 40",-verbose=>1) if $minQual<0 || $minQual > 40;
pod2usage(-msg=>"Measure not provided",-verbose=>1) unless $measure;


my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using pileup\t$pileupfile\n";
print $pfh "Using output\t$output\n";
print $pfh "Using gtf\t$gtffile\n";
print $pfh "Using measure\t$measure\n";
print $pfh "Using snp-output\t$snpfile\n";
print $pfh "Using fastq-type\t$fastqtype\n";
print $pfh "Using min-count\t$minCount\n";
print $pfh "Using min-qual\t$minQual\n";
print $pfh "Using pool-size\t$poolSize\n";
print $pfh "Using min-coverage\t$minCoverage\n";
print $pfh "Using max-coverage\t$maxCoverage\n";
print $pfh "Using min-covered-fraction\t$minCoveredFraction\n";
print $pfh "Using no-discard-deletions\t$tolerateDeletions\n";
print $pfh "Dissable corrections\t$uncorrected\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

# get the method which should be used to calculate a feature;
#my $measureCalculator=getMeasureCalculater($measure);
my $varianceCalculator=VarianceExactCorrection->new($poolSize,$minCount,$minCoverage,$maxCoverage);
$varianceCalculator=VarianceUncorrected->new($poolSize,$minCount,$minCoverage,$maxCoverage) if $uncorrected;

# qualencoding,mincount,mincov,maxcov,minqual
my $pp=get_pileup_parser($fastqtype,$minCount,$minCoverage,$maxCoverage,$minQual,$tolerateDeletions);

my ($chrdec,$genehash)=Utility::read_gtf($gtffile);

open my $ifh, "<",$pileupfile or die "Could not open pileup file";

print "Start parsing the pileup file..\n";
open my $ofh,">",$output or die "Could not open output file";

# get a SNP writer
my $snpwriter;
$snpwriter=get_gtfSNPFormater($snpfile) if $snpfile;



my $counter=0;
while(my $line=<$ifh>)
{
    # parse the whole pileup and store every parsed pileup-entry at the corresponding
    chomp $line;
    $counter++;
    print "Processed $counter pileup entries\n" unless($counter % $makeNoise);
    
    # prefilter, without proper parsing of the pileup
    my($chr,$pos)=split/\s+/,$line;
    my $entries=$chrdec->($chr,$pos);
    next unless $entries;
    
    # proper parsing; 
    my $pu=$pp->($line);

    
    foreach my $gene (@$entries)
    {
        die "problem could not find gene" unless exists($genehash->{$gene});
        
        $genehash->{$gene}{covered}++ if $pu->{iscov};
        push @{$genehash->{$gene}{snps}}, $pu if $pu->{issnp};
        
        # in case this is the last position of the gene, handle it immediately; print it         
        if($pos==$genehash->{$gene}{last})
        {
            Utility::handle_completed_gene($ofh,$genehash,$gene,$varianceCalculator,$measure,$snpwriter);
        }
        
    }
}
print "Finished parsing the pileup file; Processing features\n";
# handle the remaining genes/features
my @keys=keys(%$genehash);
for my $key (@keys)
{
    Utility::handle_completed_gene($ofh,$genehash,$key,$varianceCalculator,$measure,$snpwriter);
}

close $ofh;

print "Done\n";
exit;

{
    package Utility;
    use strict;
    use warnings;

    
    sub handle_completed_gene
    {
        my $ofh=shift;
        my $genehash=shift;
        my $geneid=shift;
        my $varianceCalculator=shift;
        my $measure=shift;
        my $snpwriter=shift;
        
        my $temp=$genehash->{$geneid};
        
        #erase the entry from the hash -> very important memory issue
        delete($genehash->{$geneid});
        
        my $length=$temp->{length};
    
        my $snps=$temp->{snps};
        my $covered=$temp->{covered};
        
        # write the snps to the output file if the snpwriter is available
        $snpwriter->($geneid,$snps) if $snpwriter;
        
        my $snpcount=@$snps;
        my $coveredFraction=$covered/$length;
        my $meas=0;
        $meas=$varianceCalculator->calculate_measure($measure,$snps,$covered);
        $meas=sprintf("%.9f",$meas);
        
        $meas="na" if $coveredFraction < $minCoveredFraction;
        
        $coveredFraction=sprintf("%.3f",$coveredFraction);
        print $ofh "$geneid\t$snpcount\t$coveredFraction\t$meas\n";
    }
    
    
    sub _parsegtf
    {
        my $line=shift;
        my @a=split /\t/, $line;
        my $ref=$a[0];
        my $start=$a[3];
        my $end=$a[4];
        my $tfeat=$a[8];
            
        unless($ref or $start or $end or $tfeat)
        {
            die "the following line is not valid";
        }
        my $gene_id="";
        if($tfeat=~/gene_id "([^"]+)";/)
        {
            $gene_id=$1;
        }
        else
        {
            die "the following entry does not have a valid gene id: $line";
        }
        
        return
        {
            ref=>$ref,
            start=>$start,
            end=>$end,
            length=>$end-$start+1,
            gi=>$gene_id
        };
    }
    
    sub _getdefaultgenecoll
    {
        return
        {
            length=>0,
            last=>0,
            snps=>[],
            covered=>0   
        }
    }
    
    sub _getgeneint
    {
        my $geneid=shift;
        my $genemap=shift;
        my $lastcounter=shift;
        if(exists($genemap->{$geneid}))
        {
            return ($lastcounter,$genemap->{$geneid});
        }
        else
        {
            $lastcounter++;
            $genemap->{$geneid}=$lastcounter;
            return ($lastcounter,$genemap->{$geneid});
            
        }
    }
    
    sub _getDecodedGenemap
    {
        my $genemap=shift;
        
        my $decode=[];
        while(my($gene,$num)=each(%$genemap))
        {
            $decode->[$num]=$gene;
        }
        return $decode;
    }
    
    sub read_gtf
    {
        my $file=shift;
        open my $ifh,"<",$file or die "Could not open gtf-file";
        my $chrhash={};
        my $genecoll={};
        my $genemap={};
        my $lastcounter=0;

        
        print "Parsing gtf file..\n";
        while(my $line=<$ifh>)
        {
            chomp $line;
            my $ge=_parsegtf($line);
            
            my $gid=$ge->{gi};
            $genecoll->{$gid}=_getdefaultgenecoll unless exists($genecoll->{$gid});
            
            my $geneint;
            ($lastcounter,$geneint)=_getgeneint($gid,$genemap,$lastcounter);
            
            # update the chromosome hash
            my($start,$end,$ref)=($ge->{start},$ge->{end},$ge->{ref});
            
            # set the new end if it is larger than the previous one
            $genecoll->{$gid}{last} = $end if $end > $genecoll->{$gid}{last};
            
            #print "$ref $start $end\n";
            
            for(my $i=$start; $i<=$end; $i++)
            {
                if(exists($chrhash->{$ref}{$i}))
                {
                    my $ar=$chrhash->{$ref}{$i};
                    push @$ar,$geneint;
                    $ar=uniq($ar);
                    $chrhash->{$ref}{$i}=$ar;
                    
                }
                else
                {
                    $chrhash->{$ref}{$i}=[$geneint];
                }
            }
        }

        
        my $decodeGeneMap=_getDecodedGenemap($genemap);

################################################################################        
        # bless the beauty of a closure
        my $chrdecocer=sub {
            my $ref=shift;
            my $pos=shift;
            my $ta=$chrhash->{$ref}{$pos};
            return undef unless $ta;
            my $dec=[];
            for my $e (@$ta)
            {
                push @$dec,$decodeGeneMap->[$e];
            }
            return $dec;
        };
################################################################################
            
            
        #calculate the length of the features
        while(my($chr,$t)=each(%$chrhash))
        {
            while(my($pos,$genes)=each(%$t))
            {
                foreach my $g (@$genes)
                {
                    my $decg=$decodeGeneMap->[$g];
                    $genecoll->{$decg}{length}++;
                }
            }
        }
            
        
        
        
        return ($chrdecocer,$genecoll);
        #chr1 Twinscan  exon         501   650   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
    }
    
    sub uniq
    {
        my $ar=shift;
        
        my $h={};
        foreach my $a (@$ar)
        {
            $h->{$a}=1;
        }
        return [keys(%$h)];
    }
}


{
    package VarTest;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use Test::Variance;
    use Test::TClassicalVariance;
    use Test::PileupParser;
    use Test;
    
    sub runTests
    {
        run_PileupParserTests(); 
        run_classicalVarianceTests();
        run_VarianceTests();
        test_read_gtf();
        exit;
    }
    
    
    sub test_read_gtf
    {
        my $str=
        "2L\tFlyBase\texon\t2\t10\t.\t+\t.\tgene_id \"g1\"; transcript_id \"whadeva\";\n".
        "2L\tFlyBase\texon\t1\t4\t.\t+\t.\tgene_id \"g1\"; transcript_id \"whadeva\";\n".
        "2L\tFlyBase\texon\t9\t12\t.\t+\t.\tgene_id \"g2\"; transcript_id \"whadeva\";\n".
        "2R\tFlyBase\texon\t1\t5\t.\t+\t.\tgene_id \"g3\"; transcript_id \"whadeva\";\n".
        "2R\tFlyBase\texon\t11\t15\t.\t+\t.\tgene_id \"g4\"; transcript_id \"whadeva\";\n";
        
        my($dec,$gh)=Utility::read_gtf(\$str);
        is($gh->{g1}{covered},0,"test gtf_reader; covered is ok");
        is(scalar(@{$gh->{g1}{snps}}),0,"test gtf_reader; length of snps is ok");
        is($gh->{g1}{last},10,"test gtf_reader: last position is ok");
        is($gh->{g1}{length},10,"test gtf_reader: length is ok");
        
        is($gh->{g2}{last},12,"test gtf_reader: last position is ok");
        is($gh->{g2}{length},4,"test gtf_reader: length is ok");
        
        is($gh->{g3}{last},5,"test gtf_reader: last position is ok");
        is($gh->{g3}{length},5,"test gtf_reader: length is ok");

        is($gh->{g4}{last},15,"test gtf_reader: last position is ok");
        is($gh->{g4}{length},5,"test gtf_reader: length is ok");
        
        is($dec->("2L","1")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","2")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","3")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","4")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","5")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","6")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","7")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","8")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","9")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","10")[0],"g1","test gtf_reader: correct gene at given position");
        is($dec->("2L","11")[0],"g2","test gtf_reader: correct gene at given position");
        is($dec->("2L","12")[0],"g2","test gtf_reader: correct gene at given position");
        is($dec->("2L","9")[1],"g2","test gtf_reader: correct gene at given position");
        is($dec->("2L","10")[1],"g2","test gtf_reader: correct gene at given position");
        not_exists($dec->("2L","13"),"test gtf_reader: correct no gene at given position");
        not_exists($dec->("2L","14"),"test gtf_reader: correct no gene at given position");
        
        is($dec->("2R","1")[0],"g3","test gtf_reader: correct gene at given position");
        is($dec->("2R","2")[0],"g3","test gtf_reader: correct gene at given position");
        is($dec->("2R","3")[0],"g3","test gtf_reader: correct gene at given position");
        is($dec->("2R","4")[0],"g3","test gtf_reader: correct gene at given position");
        is($dec->("2R","5")[0],"g3","test gtf_reader: correct gene at given position");
        
        not_exists($dec->("2R","6"),"test gtf_reader: correct no gene at given position");
        not_exists($dec->("2R","7"),"test gtf_reader: correct no gene at given position");
        not_exists($dec->("2R","8"),"test gtf_reader: correct no gene at given position");
        not_exists($dec->("2R","9"),"test gtf_reader: correct no gene at given position");
        not_exists($dec->("2R","10"),"test gtf_reader: correct no gene at given position");

        is($dec->("2R","11")[0],"g4","test gtf_reader: correct gene at given position");
        is($dec->("2R","12")[0],"g4","test gtf_reader: correct gene at given position");
        is($dec->("2R","13")[0],"g4","test gtf_reader: correct gene at given position");
        is($dec->("2R","14")[0],"g4","test gtf_reader: correct gene at given position");
        is($dec->("2R","15")[0],"g4","test gtf_reader: correct gene at given position");
        
        
    }
}



    #"measure=s"         =>\$measure,
    #"pileup=s"          =>\$pileupfile,
    #"gtf=s"             =>\$gtffile,
    #"output=s"          =>\$output,
    #"fastq-type=s"      =>\$fastqtype,
    #"min-count=i"       =>\$minCount,
    #"min-qual=i"        =>\$minQual,
    #"pool-size=i"       =>\$poolSize,
    #"min-coverage=i"    =>\$minCoverage,
    #"max-coverage=i"    =>\$maxCoverage,
    #"min-covered-fraction=f"=>\$minCoveredFraction,
    #"test"              =>\$test,
    #"help"              =>\$help


=head1 NAME

perl Variance-at-position.pl - A script which calculates the requested population genetics measure (pi, theta, d) along chromosomes using a sliding window approach.

=head1 SYNOPSIS

perl Variance-at-position.pl --measure pi --gtf annotation.gtf --pileup input.pileup --pool-size 500 --min-count 2 --min-coverage 4

=head1 OPTIONS

=over 4

=item B<--pileup>

The input file in the pileup format. A pooled population sequenced and mapped to the reference genome. Finally the mapping results have been converted to sam output format.
Using the samtools the sam format can be easily converted into the pileup format.  Mandatory.

=item B<--gtf>

A gtf file as specified here: http://mblab.wustl.edu/GTF2.html

For more information use <--help>

=item B<--output>

The output file.  Mandatory.

=item B<--snp-output>

If provided, this file will contain the polymorphic sites which have been used for this analysis; default=empty

=item B<--measure>

Currently, "pi", "theta" and "D" is supported. This stands for Tajima's Pi, Watterson's Theta, Tajima's D respectively; Mandatory

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

=item B<--no-discard-deletions>

per default sites with already a single deletion are discarded. By setting this flag sites with deletions will be used

=item B<--dissable-corrections>

Flag; Dissable correction factors; Calculates Pi/Theta/Tajima's D in the classical way not taking into account pooling or missing minor alleles; default=off

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input pileup

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

=head2 Input gtf

A gtf-file as described here http://mblab.wustl.edu/GTF2.html

 AB000381 Twinscan  exon         501   650   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
 AB000381 Twinscan  CDS          501   650   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
 AB000381 Twinscan  exon         700   800   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
 AB000381 Twinscan  CDS          700   707   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";

The script is grouping evertying with the same C<gene_id>! The feature field is not considered, so any feature may be used. The strand is not considered;
The tag C<transcript_id> will be ignored.
Attention: Do not mix the features exon with transcript/gene since the transcripts span exons and introns!!
In case features (eg. genes) are overlapping the respective SNP is considered for every feature at a certain position once; 

=head2 Output

The output will be as in the followin example:

 FBgn0031208	21	1.000	0.001365804
 FBgn0002121	17	0.984	0.000278233
 FBgn0031209	6	1.000	0.000190078

 col 1: the id of the feature (gene id, exon id, transcription factor binding site id, etc)
 col 2: SNPs found in the feature; the SNPs of all partial features - eg. exons of a gene - are summed up
 col 3: fraction of the window covered by a sufficient number of reads. Suficient means higher than min-coverage and lower than max-coverage; all partial features are considered
 col 4: population genetics estimator (pi, theta, d); this is the weighted mean of all partial features; weighted by length of the partial feature


=head2 SNP Output

This output file is optional, it will only be created if the option  C<--snp-output> is provided. For example:

 >FBgn0248264 snps:4
 2       54273   A       62      38      0       0       24      0
 2       54274   T       60      22      38      0       0       0
 2       54339   T       58      12      46      0       0       0
 2       54382   G       72      2       1       0       69      0

 >FBgn0074145 snps:2
 2       56183   T       153     0       151     2       0       0
 2       56216   G       170     2       0       0       168     0

 The header contains the following information
 >geneid snps:scnpcount
 geneid.. the ID of the genes
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
The sum of the bases A,T,C,G is used as coverage (ignoring N)

The memory usage scales linearly with the amount of gtf entries provided. If you get an out of memory exception (eg. malloc) split your gtf file and run the script several times.
  
=cut



