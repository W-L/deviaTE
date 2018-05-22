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
    use VarianceUncorrected;
    use Pileup;
    use PileupTripletSlider;
    use SynNonSyn;
    use FormatSNP;
    our $verbose=1;
    my $shoutinterval=1_000;
    my $tripletcount=0;
    
    # input files
    my $pileupfile="";
    my $gtf_file="";
    my $codonTableFile="";
    my $nonsynLengthTableFile="";
    my $output="";
    
    # pileup sliding
    my $minCount=2;
    my $minCoverage=4;
    my $maxCoverage=400;
    my $fastqtype="illumina";
    my $minQual=20;
    my $snpfile="";
    my $usedregionfile="";
    
    # general
    my $measure="";
    my $poolSize=0;
    my $maxTripletSNPs=3;
    my $uncorrected=0;
    my $help=0;
    my $test=0;
    
    # --region-output /Users/robertkofler/dev/testfiles/output/2R-region.txt  --measure pi --pool-size 500 --gtf /Users/robertkofler/dev/testfiles/small-2R-cds.gff --pileup /Users/robertkofler/dev/testfiles/2R_sim_100000.pileup --output /Users/robertkofler/dev/testfiles/output/2R-syn-nonsyn.txt --codon-table /Users/robertkofler/dev/PopGenTools/syn-nonsyn/codon-table.txt  --snp-output /Users/robertkofler/dev/testfiles/output/2R-syn-nonsyn-snp-gene.txt --nonsyn-length-table /Users/robertkofler/dev/PopGenTools/syn-nonsyn/nsl_p6.txt

    GetOptions(
        "measure=s"         =>\$measure,
        "pileup=s"          =>\$pileupfile,
        "gtf=s"             =>\$gtf_file,
        "codon-table=s"     =>\$codonTableFile,
        "nonsyn-length-table=s"=>\$nonsynLengthTableFile,
        "output=s"          =>\$output,
        "snp-output=s"      =>\$snpfile,
        "fastq-type=s"      =>\$fastqtype,
        "min-count=i"       =>\$minCount,
        "min-qual=i"        =>\$minQual,
        "max-triplet-snps=i"=>\$maxTripletSNPs,
        "pool-size=i"       =>\$poolSize,
        "min-coverage=i"    =>\$minCoverage,
        "max-coverage=i"    =>\$maxCoverage,
        "region-output=s"   =>\$usedregionfile,
        "dissable-corrections"=>\$uncorrected,
        "test"              =>\$test,
        "help"              =>\$help
    ) or die "Invalid arguments";


    pod2usage(-verbose=>2) if $help;
    SynGeneTest::runTests() if $test;
    pod2usage(-msg => "Could not find pileup file",-verbose=>1)                   unless -e $pileupfile;
    pod2usage(-msg => "Could not find gtf file",-verbose=>1)                      unless -e $gtf_file;
    pod2usage(-msg => "Could not find codon table",-verbose=>1)                   unless -e $codonTableFile;
    pod2usage(-msg => "Could not find non-synonymous length table",-verbose=>1)   unless -e $nonsynLengthTableFile;
    pod2usage(-msg => "Output file not provided",-verbose=>1)                     unless  $output;
    pod2usage(-msg => "Pool size not provided",-verbose=>1)                       unless $poolSize;
    pod2usage(-msg => "Min count not provided",-verbose=>1)                       unless $minCount;
    pod2usage(-msg => "Min quality not valid. Has to be between 0 and 40",-verbose=>1) if $minQual<0 || $minQual > 40;
    pod2usage(-msg => "The minimum coverage hast to be at least two times the minimum count",-verbose=>1) unless $minCoverage >= (2*$minCount);
    pod2usage(-msg => "Measure not provided",-verbose=>1) unless $measure;
    
    
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
    print $pfh "Using min-count\t$minCount\n";
    print $pfh "Using min-qual\t$minQual\n";
    print $pfh "sing max-triplet-snps\t$maxTripletSNPs\n";
    print $pfh "Using pool-size\t$poolSize\n";
    print $pfh "Using min-coverage\t$minCoverage\n";
    print $pfh "Using max-coverage\t$maxCoverage\n";
    print $pfh "Using region-output\t$usedregionfile\n";
    print $pfh "Dissable corrections\t$uncorrected\n";
    print $pfh "Using test\t$test\n";
    print $pfh "Using help\t$help\n";
    close $pfh;

    print "Loading codon table ...\n";
    my $codonTable = load_codon_table($codonTableFile);
    print "Loading nonsynonymous length table...\n";
    my $nonsynTable=load_nonsyn_length_table($nonsynLengthTableFile);
    print "Encoding frame information from gtf file...\n";
    my $chrFrames = load_cds_gtf($gtf_file);
    print "Loading gene annotation from gtf file...\n";
    my $geneassigner = get_codon_to_gene_assigner($chrFrames, $gtf_file);
    
    
    # gradually building the pileup window slider
    my $pp=get_extended_parser($fastqtype,$minCount,$minCoverage,$maxCoverage,$minQual);
    my $pts=PileupTripletSlider->new($pileupfile,$chrFrames,$pp);
     
    my $vec=VarianceExactCorrection->new($poolSize,$minCount,$minCoverage,$maxCoverage);
    $vec=VarianceUncorrected->new($poolSize,$minCount,$minCoverage,$maxCoverage) if $uncorrected;
    
    # get snp writer
    my $snpwriter;
    $snpwriter =get_synnonsyngene_SNPFormater($snpfile) if $snpfile;
    
    
    
    my $genecollection={};
    open my $ofh, ">", $output or die "Could not open output file";
    print "Parsing pileup file...\n";
    while(my $tr=$pts->next())
    {
        # TRIPLET definition
        # frame, pileup, chr, start, strand, valid, valid_frame, valid_coverage, valid_codon, count_snps, codon
        my($chr,$start,$valid,$count_snps)=($tr->{chr},$tr->{start},$tr->{valid},$tr->{count_snps});
        
        $tripletcount++;
        # shout out loud
        print "Processed $tripletcount triplets\n" if($tripletcount % $shoutinterval ==0);
        
        next unless ($valid);
        next unless ($count_snps <= $maxTripletSNPs);
        
        # get genelist, nonsynlength and update the gene-info
        # for all triplets no matter if they contain snps or not
        my $genelist=$geneassigner->($chr,$start);
        my ($sl,$nsl)=Utility::calculate_syn_nonsynlength($tr,$nonsynTable);
        Utility::basicupdate_genelist($genecollection,$_,$chr,$start,$sl,$nsl) foreach (@$genelist);
        
        # from here procede only with triplets containing SNPs
        next unless ($count_snps>0);
        get_codon_changes($tr,$codonTable);
        $snpwriter->($tr,$genelist) if $snpfile;
        

        foreach my $geneid (@$genelist)
        {
            my $gene=$genecollection->{$geneid};
            my $codonchanges=$tr->{cc};
            foreach my $cc (@$codonchanges)
            {

                if($cc->{syn})
                {
                    push @{$gene->{synsnplist}},$cc;
                    $gene->{synsnps}++;
                }
                else
                {
                    push @{$gene->{nonsynsnplist}},$cc;
                    $gene->{nonsynsnps}++;                
                }
            }
        
        }
    }
    # a gene has
    # synsnplist, nonsynsnplist, synsnps,nonsynsnps,synlen,nonsynlen ,countcodon, coveredreg
    
    print "Finished parsing pielup file...\n";
    
    my $regpr = undef;
    $regpr=Utility::get_region_printer($usedregionfile) if $usedregionfile;
    
    print "Writing results to output file...\n";
    while(my($geneid,$t)=each(%$genecollection))
    {
        
        $regpr->($t) if $usedregionfile;
        #
        # synsnplist, nonsynsnplist, synsnps,nonsynsnps,synlen,nonsynlen ,countcodon, coveredreg
        my($countcodon,$synlen,$nonsynlen,$synsnps,$nonsynsnps)=($t->{countcodon},$t->{synlen},$t->{nonsynlen},$t->{synsnps},$t->{nonsynsnps});
        if(not $synlen or not $nonsynlen)
        {
            warn "No coverage for gene $t->{geneid}; Skipping gene..\n";
            next;
        }

        #$vec->calculate_measure($measure,$synsnplist,$synlength);
        my $synmeasure=$vec->calculate_measure($measure,$t->{synsnplist},$synlen);
        my $nonsynmeasure=$vec->calculate_measure($measure,$t->{nonsynsnplist},$nonsynlen);
        
        $synmeasure     =   sprintf("%.8f",$synmeasure);
        $nonsynmeasure  =   sprintf("%.8f",$nonsynmeasure);
        
        print $ofh "$t->{geneid}\t$nonsynlen\t$synlen\t$nonsynsnps\t$synsnps\t$nonsynmeasure\t$synmeasure\n";
    }
    print "DONE\n";
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
    
    sub get_region_printer
    {
        my $output=shift;
        open my $ofh, ">",$output or die "Could not open output file\n";
        return sub
        {
            my $reg=shift;
            # geneid, chr, synlen, nonsynlen, synmeasure, nonsynmeasure, synsnps, nonsynsnps, coveredreg, countcodon
            my $geneid=$reg->{geneid};
            my $covreg=$reg->{coveredreg};
            my $chr=$reg->{chr};
            print $ofh ">$geneid\n";
            foreach my $cr (@$covreg)
            {
                my $start=$cr->{start};
                my $end=$cr->{end};
                print $ofh "$chr:$start-$end\n";
                
            }
            print $ofh "\n";
        }
    }
    
    
    sub basicupdate_genelist
    {
        my $genecoll=shift;
        my $gene =shift;
        my $chr=shift;
        my $start=shift;
        my $synlength=shift;
        my $nonsynlength=shift;
        
        if(exists($genecoll->{$gene}))
        {
            # the gene already exists, only update the list of covered regions
            my $e=$genecoll->{$gene};
            my $lastcovered=@{$e->{coveredreg}}[-1];
            if($lastcovered->{end}+1==$start)
            {
                $lastcovered->{end}=$start+2;
            }
            else
            {
                push @{$e->{coveredreg}},{start=>$start,end=>$start+2};
            }
        }
        else
        {
            # gene does not exist, get a default gene entry and create one covered region
            my $defe=get_default_geneentry($gene,$chr);
            push @{$defe->{coveredreg}},{start=>$start,end=>$start+2};
            $genecoll->{$gene}=$defe;
            
        }
        $genecoll->{$gene}{nonsynlen}   += $nonsynlength;
        $genecoll->{$gene}{synlen}      +=$synlength;
        $genecoll->{$gene}{countcodon}  ++;
        # coveredreg, nonsynlen,synlen,countcodon
    }
    
    sub get_default_geneentry{
        my $geneid=shift;
        my $chr=shift;
        
        # geneid, chr, synlen, nonsynlen, synmeasure, nonsynmeasure, synsnps, nonsynsnps, coveredreg
        return
        {
          geneid        =>$geneid,
          chr           =>$chr,
          synlen        =>0,
          countcodons   =>0,
          nonsynlen     =>0,
          synmeasure    =>0,
          nonsynmeasure =>0,
          synsnps       =>0,
          nonsynsnps    =>0,
          coveredreg    =>[],
          synsnplist    =>[],
          nonsynsnplist =>[]
        };
    }
    
    sub calculate_syn_nonsynlength
    {
        my $tr=shift;
        # chr, pos, ref, strand, cstart, A, T, C, G, eucov, codon, mcodon, aa, maa, syn
        my $nonsyntable=shift;
        my $val=$nonsyntable->{$tr->{codon}};
        return ((3-$val),$val);
    }
}





{
    package SynGeneTest;
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
        run_VarianceTests();
        exit;
        
    }
    
    
}

=head1 NAME

perl Syn-nonsyn-at-position.pl - Calculates a genewise (or for any other feature) synonymousn and nonsynonymous pi/theta/d

=head1 SYNOPSIS

perl Syn-nonsyn-sliding.pl --measure pi --gtf cds-annotation.gtf --pileup input.pileup --codon-table codons.txt --nonsyn-length-table nslength.txt --output output.file --pool-size 500 --min-count 2 --min-coverage 4

=head1 OPTIONS

=over 4

=item B<--pileup>

The input file in the pileup format. A pooled population sequenced and mapped to the reference genome. The mapping results have been converted to sam output format.
Using the samtools the sam format can be easily converted into the pileup format.  Mandatory.

=item B<--gtf>

An annotation of the reference genome in the gtf format. Only 'CDS' entries are considered. The attribute gene_id "*"; is required for every CDS entry; Mandatory

=item B<--codon-table>

A file with a table containing the codons and the resulting ammino acids. A default table is provided with the software. Mandatory

=item B<--nonsyn-length-table>

A file with a table which contains the codons and the nonsynonymous length for the codon. Some default tables are provided with the software. Mandatory

=item B<--output>

The output file.  Mandatory.

=item B<--snp-output>

The SNP output; These SNPs which have been used to calcualte the synonymous and the non-synonymous measure for every gene (feature); Optional parameter

=item B<--region-output>

The regions which have been used for calculating the synonymous and the non-synonymous measure for every gene (feature). Optional parameter

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

=item B<--min-qual>

The minimum quality; Alleles with a quality lower than this threshold will not be considered (for SNP identification and coverage estimation); default=20

=item B<--max-triplet-snps>

The maximum number of SNPs in a triplet to be considered. If the number of SNPs in the triplet exceeds this number it will be ignored; default=3
The algorithm will treat every SNP independently.

=item B<--dissable-corrections>

Flag; Dissable correction factors; Calculates Pi/Theta/Tajima's D in the classical way not taking into account pooling or missing minor alleles; default=off

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

 2R      maker   CDS     3391    3453    .       +       0       gene_id "sneezeless";
 2R      maker   CDS     3507    3599    .       +       0       gene_id "sneezeless";
 2R      maker   CDS     3653    4310    .       +       0       gene_id "velourcaressasse";
 2R      maker   CDS     4369    4988    .       +       2       gene_id "promotorofRTFM";

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

 schokicraving	    683.375	1599.625    10	28	0.00631134	0.00317790 
 sleepless	    71.625	180.375	    2	9	0.00872667	0.01651562
 sneezeless	    137.5	309.5	    7	23	0.01495735	0.02248795
 senseless	    223.875	547.125	    3	13	0.00360575	0.00641439
 velourcaressasse   23.5	57.5	    3	5	0.09097559	0.05767167

 col1: gene name
 col2: non-synonymous length
 col3: synonymous length
 col4: non-synonymous SNPs
 col5: synonymous SNPs
 col6: non-synonymous measure (pi/theta/d)
 col7: synonyomous measure (pi/theta/d)

=head2 Output SNPs

The output is provided for each window; First a definition of the window and then subequently the SNPs found in the window;
Example of output:

 Senseless    X  1448874 G       90      0       5       0       85       non-syn -       CCA->CAA        P->Q
 Sleepless    X  1448909 G       75      14      0       0       61       syn     -       ACC->ACT        T->T
 Senseless    X  1449241 T       105     0       96      9       0        non-syn -       AAA->GAA        K->E
 Senseless    X  1449275 G       112     0       10      0       102      non-syn -       AAC->AAA        N->K

 The subsequent SNP definitions are tab-delimited:
 col1: name of the gene
 col2: the reference contig
 col3: the position of the SNP
 col4: reference character
 col5: the coverage
 col6: the count of 'A'
 col7: the count of 'T'
 col8: the count of 'C'
 col9: the count of 'G'
 col10: SNP type, synonymous or non-synonymous
 col11: the strand of the codon
 col12: the codon and the codon change caused by the SNP in the form: from->to
 col13: the amino acid and the amino acid change caused by the SNP in the form: from->to

=head2 Output Region

The regions which have been used to calculate the synonymous and non-synonymous measure for every gene

 >sleepless
 2R:8002-8250
 2R:8300-8302

 >sneezeless
 2R:1000-1029
 2R:1051-1137
 2R:1171-1500

 First it contains the gene id
 than for every gene the used regions in form chr:start-end

 chr: chromosome
 start: start position
 end: end position

=head1 AUTHORS

 Robert Kofler

 Christian Schloetterer

=cut