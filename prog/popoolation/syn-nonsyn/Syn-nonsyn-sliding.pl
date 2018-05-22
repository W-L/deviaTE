#!/usr/bin/env perl
{
    use strict;
    use warnings;
    use Data::Dumper;
    use Getopt::Long;
    use Pod::Usage;
    use FindBin qw($RealBin);
    use lib "$RealBin/../Modules";
    use VarianceExactCorrection;
    use Pileup;
    use PileupTripletSlider;
    use PileupTripletWindowSlider;
    use SynNonSyn;
    use FormatSNP;
    our $verbose=1;
    
    # input files
    my $pileupfile="";
    my $gtf_file="";
    my $codonTableFile="";
    my $nonsynLengthTableFile="";
    my $output="";
    
    # pileup sliding
    my $windowSize= 1000;
    my $stepSize=1000;
    my $minCount=2;
    my $minCoverage=4;
    my $maxCoverage=400;
    my $minCoveredFraction=0.3;
    my $fastqtype="illumina";
    my $minQual=20;
    my $snpfile="";
    
    # general
    my $measure="";
    my $poolSize=0;
    my $maxTripletSNPs=3;
    my $suppressNa=0;
    my $uncorrected=0;
    my $help=0;
    my $test=0;
    
    # --measure pi --pool-size 500 --gtf /Users/robertkofler/dev/testfiles/2R.gff --pileup /Users/robertkofler/dev/testfiles/basepop_2R.pileup --output /Users/robertkofler/dev/testfiles/output/2R-syn-nonsyn.txt --codon-table /Users/robertkofler/dev/PopGenTools/syn-nonsyn/codon-table.txt --nonsyn-length-table /Users/robertkofler/dev/PopGenTools/syn-nonsyn/nsl_p6.txt
    # --measure D --dissable-corrections --pool-size 500 --gtf /Volumes/Volume_3/analysis/syn-nonsyn/X.gtf --pileup /Volumes/Volume_3/analysis/syn-nonsyn/x.pileup --output /Volumes/Volume_3/analysis/syn-nonsyn/testout.txt --codon-table /Users/robertkofler/dev/popoolation/syn-nonsyn/codon-table.txt --nonsyn-length-table /Users/robertkofler/dev/popoolation/syn-nonsyn/nsl_p6.txt --min-count 1   
    
    GetOptions(
        "measure=s"         =>\$measure,
        "pileup=s"          =>\$pileupfile,
        "gtf=s"             =>\$gtf_file,
        "codon-table=s"     =>\$codonTableFile,
        "nonsyn-length-table=s"=>\$nonsynLengthTableFile,
        "output=s"          =>\$output,
        "snp-output=s"      =>\$snpfile,
        "fastq-type=s"      =>\$fastqtype,
        "window-size=i"     =>\$windowSize,
        "step-size=i"       =>\$stepSize,
        "min-count=i"       =>\$minCount,
        "min-qual=i"        =>\$minQual,
        "max-triplet-snps=i"=>\$maxTripletSNPs,
        "pool-size=i"       =>\$poolSize,
        "min-coverage=i"    =>\$minCoverage,
        "max-coverage=i"    =>\$maxCoverage,
        "min-covered-fraction=f"=>\$minCoveredFraction,
        "suppress-na"       =>\$suppressNa,
        "dissable-corrections"=>\$uncorrected,
        "test"              =>\$test,
        "help"              =>\$help
    ) or die "Invalid arguments";


    pod2usage(-verbose=>2) if $help;
    SynTest::runTests() if $test;
    pod2usage(-msg=>"Could not find pileup file",-verbose=>1)                   unless -e $pileupfile;
    pod2usage(-msg=>"Could not find gtf file",-verbose=>1)                      unless -e $gtf_file;
    pod2usage(-msg=>"Could not find codon table",-verbose=>1)                   unless -e $codonTableFile;
    pod2usage(-msg=>"Could not find non-synonymous length table",-verbose=>1)   unless -e $nonsynLengthTableFile;
    pod2usage(-msg=>"Output file not provided",-verbose=>1)                     unless  $output;
    pod2usage(-msg=>"Pool size not provided",-verbose=>1)                       unless $poolSize;
    pod2usage(-msg=>"Min count not provided",-verbose=>1)                       unless $minCount;
    pod2usage(-msg=>"Min quality not valid. Has to be between 0 and 40",-verbose=>1) if $minQual<0 || $minQual > 40;
    pod2usage(-msg=>"The minimum coverage hast to be at least two times the minimum count",-verbose=>1) unless $minCoverage >= (2*$minCount);
    pod2usage(-msg=>"Measure not provided",-verbose=>1) unless $measure;
    
    
    my $paramfile=$output.".params";
    open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
    print $pfh "Using measure\t$measure\n";
    print $pfh "Using pileup\t$pileupfile\n";
    print $pfh "Using gtf\t$gtf_file\n";
    print $pfh "Using codon-talbe\t$codonTableFile\n";
    print $pfh "Using nonsyn-length-table\t$nonsynLengthTableFile\n";
    print $pfh "Using output\t$output\n";
    print $pfh "Using snp-output\t$snpfile\n";
    print $pfh "Using fastq-type\t$fastqtype\n";
    print $pfh "Using window-size\t$windowSize\n";
    print $pfh "Using step-size\t$stepSize\n";
    print $pfh "Using min-count\t$minCount\n";
    print $pfh "Using min-qual\t$minQual\n";
    print $pfh "Using max-triplet-snps\t$maxTripletSNPs\n";
    print $pfh "Using pool-size\t$poolSize\n";
    print $pfh "Using min-coverage\t$minCoverage\n";
    print $pfh "Using max-coverage\t$maxCoverage\n";
    print $pfh "Using min-covered-fraction\t$minCoveredFraction\n";
    print $pfh "Using suppress-na\t$suppressNa\n";
    print $pfh "Dissable corrections\t$uncorrected\n";
    print $pfh "Using test\t$test\n";
    print $pfh "Using help\t$help\n";
    close $pfh;

    print "Loading codon table ...\n";
    my $codonTable=load_codon_table($codonTableFile);
    print "Loading nonsynonymous length table...\n";
    my $nonsynTable=load_nonsyn_length_table($nonsynLengthTableFile);
    print "Loading gtf file...\n";
    my $chrAnotation=load_cds_gtf($gtf_file);
    print "Parsing pileup file..\n";
    
    
    # gradually building the pileup window slider
    my $pp=get_extended_parser($fastqtype,$minCount,$minCoverage,$maxCoverage,$minQual);
    my $pts=PileupTripletSlider->new($pileupfile,$chrAnotation,$pp);
    my $ptws=PileupTripletWindowSlider->new($pts,$windowSize,$stepSize);
    
    # get measure calculator
    my $meascalc = Utility::get_measure_calculator($minCount,$poolSize,$minCoverage,$maxCoverage,$measure,$nonsynTable,$uncorrected);
    
    # get snp writer
    my $snpwriter;
    $snpwriter =get_syn_nonsyn_SNPFormater($snpfile) if $snpfile;

    open my $ofh, ">", $output or die "Could not open output file";

    # annotated triplet window
    # data, chr, lower, upper, window, count_codons, cound_valid_codons, count_useful, count_onesnp_codons
    # TRIPLET definition
    # frame, pileup, chr, start, strand, valid, valid_frame, valid_coverage, valid_codon, count_snps, codon
    while(my $win=$ptws->nextWindow())
    {
        my $triplets=$win->{data};
        my $chr=$win->{chr};
        my $middle=int(($win->{lower}+$win->{upper})/2);
        my $winlength=$win->{window};
        my $count_codons=$win->{count_codons};
        print "Processing window: $chr:$win->{lower}-$win->{upper}; Codons in window: $count_codons\n";
        
        # frame is ok (123 654); valid codon (no  'N's); valid coverage of pileup (every base fully covered);
        my $valid_triplets=[grep {$_->{valid}} @$triplets];
        
        # take only triplets with an appropriate number of SNPs and incorporate the codon changes
        my $valid_snptriplets=[grep {$_->{count_snps}<=$maxTripletSNPs} @$valid_triplets];
        get_codon_changes($_,$codonTable) foreach @$valid_snptriplets;
        
        
        $snpwriter->($valid_snptriplets,$chr,$win->{lower},$win->{upper}) if $snpwriter;
        
        # calculate the measure
            # synmeasure, nonsynmeasure, synsnps, nonsynsnps, synlength, nonsynlength, validcodons, codonswithsnps
        my $measure=$meascalc->($valid_snptriplets);
        
        my $usedregions=$measure->{validcodons}*3;
        my $usedfraction=$usedregions/$winlength;
        if($usedfraction<$minCoveredFraction)
        {
            $measure->{synmeasure}="na";
            $measure->{nonsynmeasure}="na";
        }
        

        $usedfraction=sprintf("%.3f",$usedfraction);
        if($measure->{synmeasure} ne "na" or not $suppressNa)
        {
            print $ofh "$chr\t$middle\t$count_codons\t$measure->{validcodons}\t$measure->{codonswithsnps}\t$usedfraction\t$measure->{nonsynlength}\t".
            "$measure->{synlength}\t$measure->{nonsynsnps}\t$measure->{synsnps}\t$measure->{nonsynmeasure}\t$measure->{synmeasure}\n";
        }
    }
    close $ofh;
    
    
    
    exit;
}


{
    package Utility;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/../Modules";
    use VarianceExactCorrection;
    use VarianceUncorrected;
    use Test;
    use SynNonSyn;
    
    
    
    
    sub get_measure_calculator
    {
        my $mincount=shift;
        my $poolsize=shift;
        my $minCoverage=shift;
        my $maxCoverage=shift;
        my $measure=shift;
        my $nonsynTable=shift;
        my $uncorrected=shift;
        my $vec=VarianceExactCorrection->new($poolsize,$mincount,$minCoverage,$maxCoverage);
        $vec=VarianceUncorrected->new($poolsize,$mincount,$minCoverage,$maxCoverage) if $uncorrected;
        
        return sub
        {
            my $triplets=shift;
            
            my($synlength,$nonsynlength)=(0,0);
            my($synmeasure,$nonsynmeasure)=(0,0);
            my($synsnplist,$nonsynsnplist)=([],[]);
            my($count_valid,$count_codons_withsnps)=(0,0);

            foreach my $tr (@$triplets)
            {
                $count_valid++;
                
                # calculate the length
                # chr, pos, ref, strand, cstart, A, T, C, G, eucov, codon, mcodon, aa, maa, syn
                my ($sl,$nsl)=Utility::_calculate_syn_nonsynlength($tr,$nonsynTable);
                $synlength+=$sl;
                $nonsynlength+=$nsl;
                
                # only proceed with the codons that contain a SNP;
                next unless $tr->{count_snps};
                $count_codons_withsnps++;
                
                
                my $codonchanges =$tr->{cc};
                my $break;
                
                foreach my $cc (@$codonchanges)
                {
                    my $syn=$cc->{syn};
                    if($syn)
                    {
                        push @$synsnplist,$cc;
                    }
                    else
                    {
                        push @$nonsynsnplist,$cc;
                    } 
                }
            }
            my $synsnps=@$synsnplist;
            my $nonsynsnps=@$nonsynsnplist;
            $synmeasure = $vec->calculate_measure($measure,$synsnplist,$synlength);
            $nonsynmeasure=$vec->calculate_measure($measure,$nonsynsnplist,$nonsynlength);
            
            $synmeasure=sprintf("%.8f",$synmeasure);
            $nonsynmeasure=sprintf("%.8f",$nonsynmeasure);

            # synmeasure, nonsynmeasure, synsnps, nonsynsnps, synlength, nonsynlength, validcodons, codonswithsnps
            return {
                synmeasure=>$synmeasure,
                nonsynmeasure=>$nonsynmeasure,
                synsnps=>$synsnps,
                nonsynsnps=>$nonsynsnps,
                synlength=>$synlength,
                nonsynlength=>$nonsynlength,
                validcodons=>$count_valid,
                codonswithsnps=>$count_codons_withsnps
            };
        };
    }
    
    
    sub _calculate_syn_nonsynlength
    {
        my $tr=shift;
        # chr, pos, ref, strand, cstart, A, T, C, G, eucov, codon, mcodon, aa, maa, syn
        my $nonsyntable=shift;
        my $val=$nonsyntable->{$tr->{codon}};
        return ((3-$val),$val);
    }
}





{
    package SynTest;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/../Modules";
    use Test;
    use Test::TSynNonSyn;
    use Test::PileupParser;
    use Test::TPileupTripletSliding;
    use Test::Variance;
    
    
    sub runTests
    {
        run_PileupParserTests();
        run_SynNonSynTests();
        run_PileupTripletSlidingTests();
        run_VarianceTests();
        exit;
        
    }
    
    
}


    #GetOptions(
    #    "measure=s"         =>\$measure,
    #    "pileup=s"          =>\$pileupfile,
    #    "gtf=s"             =>\$gtf_file,
    #    "codon-table=s"     =>\$codonTableFile,
    #    "nonsyn-length-table=s"=>\$nonsynLengthTableFile,
    
    #    "output=s"          =>\$output,
    #    "fastq-type=s"      =>\$fastqtype,
    #    "window-size=i"     =>\$windowSize,
    #    "step-size=i"       =>\$stepSize,
    #    "min-count=i"       =>\$minCount,
    #    "min-qual=i"        =>\$minQual,
    #    "pool-size=i"       =>\$poolSize,
    #    "min-coverage=i"    =>\$minCoverage,
    #    "max-coverage=i"    =>\$maxCoverage,
    #    "min-covered-fraction=f"=>\$minCoveredFraction,
    #    "suppress-na"       =>\$suppressNa,
    #    "test"              =>\$test,
    #    "help"              =>\$help
    #) or die "Invalid arguments";
    

=head1 NAME

perl Syn-nonsyn-sliding.pl - A script which calculates the synonymousn and nonsynonymous pi/theta/d along chromosomes using a sliding window approach.

=head1 SYNOPSIS

perl Syn-nonsyn-sliding.pl --measure pi --gtf cds-annotation.gtf --pileup input.pileup --codon-table codons.txt --nonsyn-length-table nslength.txt --output output.file --pool-size 500 --min-count 2 --min-coverage 4 --window-size 50000 --step-size 10000

=head1 OPTIONS

=over 4

=item B<--pileup>

The input file in the pileup format. A pooled population sequenced and mapped to the reference genome. Finally the mapping results have been converted to sam output format.
Using the samtools the sam format can be easily converted into the pileup format.  Mandatory.

=item B<--gtf>

An annotation of the reference genome in the gtf format. Only 'CDS' entries are considered. Mandatory

=item B<--codon-table>

A file with a table containing the codons and the resulting ammino acids. A default table is provided with the software. Mandatory

=item B<--nonsyn-length-table>

A file with a table which contains the codons and the nonsynonymous length for the codon. Some default tables are provided with the software. Mandatory

=item B<--output>

The output file.  Mandatory.

=item B<--snp-output>

The SNP output; These SNPs have been used to calcualte the synonymous and the non-synonymous measure; Optional parameter

=item B<--measure>

Currently, "pi", "theta" and "d" is supported. Mandatory

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

=item B<--max-triplet-snps>

The maximum number of SNPs in a triplet to be considered. If the number of SNPs in the triplet exceeds this number it will be ignored; default=3
The algorithm will treat every SNP independently.

=item B<--window-size>

The size of the sliding window. default=50000

=item B<--step-size>

the size of one sliding window step. If this number is equal to the --window-size the sliding window will be non overlapping (jumping window). default=10000

=item B<--dissable-corrections>

Flag; Dissable correction factors; Calculates Pi/Theta/Tajima's D in the classical way not taking into account pooling or missing minor alleles; default=off

=item B<--suppress-na>

flag; If provided the output of fields containing na will be suppressed

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input Pileup

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

Example

 2R      maker   CDS     3391    3453    .       +       0       ID=maker-2R-snap-gene-0.2-mRNA-1:cds:0;Parent=maker-2R-snap-gene-0.2-mRNA-1;
 2R      maker   CDS     3507    3599    .       +       0       ID=maker-2R-snap-gene-0.2-mRNA-1:cds:1;Parent=maker-2R-snap-gene-0.2-mRNA-1;
 2R      maker   CDS     3653    4310    .       +       0       ID=maker-2R-snap-gene-0.2-mRNA-1:cds:2;Parent=maker-2R-snap-gene-0.2-mRNA-1;
 2R      maker   CDS     4369    4988    .       +       2       ID=maker-2R-snap-gene-0.2-mRNA-1:cds:3;Parent=maker-2R-snap-gene-0.2-mRNA-1;

=head2 Input codon table

Must contain the codons and the resulting amino acids. For example:

 ATT: I
 ATC: I
 ATA: I
 CTT: L 

=head2 Input nonsynonymous length

Must contain for every codon its nonsynonymous length; the synonymous length is simple calculated 3-nonsynlength; For example:

 AGT: 2.66666666666667
 TGA: 2.66666666666667
 TGT: 2.66666666666667
 CGA: 1.66666666666667


=head2 Output 

Example of output:

 2R      232500  147     147     2       0.441   315.875 125.125 1       1       0.000158349120911531    0.00231817484638718
 2R      234500  182     182     2       0.546   379.375 166.625 0       2       0       0.00113553929818699
 2R      236500  240     240     0       0.720   512.5   207.5   0       0       0       0
 2R      242500  267     267     2       0.801   559.125 241.875 1       1       0.000230523800373612    0.000861676251421832
 2R      292500  168     164     30      0.492   346.25  145.75  18      12      0.0111118663734969      0.0182644779281468

 col1: the reference contig
 col2: middle position of the sliding window
 col3: number of codons found in the window
 col4: number of codons used for the analysis
 col5: codons containing a SNP
 col6: fraction of the window covered by codons used for the analysis (useful fraction of the window)
 col7: sum of the nonsynonymous length for all codons in the window
 col8: sum of the synoymous length for all codons in the window
 col9: number of nonsynonymous SNPs in the window
 col10: number of synonymous SNPs in the window
 col11: nonsynonymous measure (pi/theta/d)
 col12: synonyomous measure (pi/theta/d)


=head2 Output SNPs

The output is provided for each window; First a definition of the window and then subequently the SNPs found in the window;
Example of output:


 >X:1448800-1449300 snps: 4
 X  1448874 G       90      0       5       0       85       non-syn -       CCA->CAA        P->Q
 X  1448909 G       75      14      0       0       61       syn     -       ACC->ACT        T->T
 X  1449241 T       105     0       96      9       0        non-syn -       AAA->GAA        K->E
 X  1449275 G       112     0       10      0       102      non-syn -       AAC->AAA        N->K

 The window position is provided in the form
 >chr:start-ent snps: snp-count
 chr:       the reference contig
 start:     the start position of the window
 end:       the end position of the window
 snp-count: the number of SNPs found in the window

 The subsequent SNP definitions are tab-delimited:
 col1: the reference contig
 col2: the position of the SNP
 col3: reference character
 col4: the coverage
 col5: the count of 'A'
 col6: the count of 'T'
 col7: the count of 'C'
 col8: the count of 'G'
 col9: SNP type, synonymous or non-synonymous
 col10: the strand of the codon
 col11: the codon and the codon change caused by the SNP in the form: from->to
 col12: the amino acid and the amino acid change caused by the SNP in the form: from->to

=cut
