#!/usr/bin/env perl

# this script is part of PoPoolation by Robert Kofler et al.
# https://sourceforge.net/projects/popoolation/

{
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin";
use deviaTE_trim;


    

#--input /Volumes/Volume_3/analysis/clipsignal/raw/s100000.fastq --output /Users/robertkofler/tmp/filter.test.fastq

    my $processStep=100000;
    my $help=0;
    my $test=0;
    my $qualThreshold=20;
    my $minLength=40;
    my $fastqtype="illumina";
    my $discardRemainingNs=0;
    my $trimQuality=1;
    my $input1="";
    my $input2="";
    my $output="";
    my $verbose=1;
    
    GetOptions(
        "input1=s"              =>\$input1,
        "input2=s"              =>\$input2,
        "output=s"              =>\$output,
        "quality-threshold=i"   =>\$qualThreshold,
        "min-length=i"          =>\$minLength,
        "fastq-type=s"          =>\$fastqtype,
        "discard-internal-N"    =>\$discardRemainingNs,
        "no-trim-quality"       =>sub{$trimQuality=0},
        "quit"                  =>sub{$verbose=0},
        "test"                  =>\$test,
        "help"                  =>\$help
    ) or pod2usage(-msg=>"Options not valid $!",-verbose=>1);
    
    # dive into the alternative branches (if requested)
    pod2usage(-verbose=>2) if $help;
    TestTrim::runTrimTests() if $test;
    pod2usage(-msg=>"Quality must be larger than 0",-verbose=>1) if $qualThreshold<0;
    pod2usage(-msg=>"Length must be larger than 0",-verbose=>1) if $minLength<1 ;  # min length has to be 1 or larger
    pod2usage(-msg=>"At least one input file has to be provided", -verbose=>1) unless -e $input1;
    pod2usage(-msg=>"An output file has to be provided", -verbose=>1) unless $output;
    
    my $paramfile=$output.".params";
    open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
    print $pfh "Using input1\t$input1\n";
    print $pfh "Using input2\t$input2\n";
    print $pfh "Using output\t$output\n";
    print $pfh "Using quality-threshold\t$qualThreshold\n";
    print $pfh "Using min-length\t$minLength\n";
    print $pfh "Using fastq-type\t$fastqtype\n";
    print $pfh "Using discard-internal-N\t$discardRemainingNs\n";
    print $pfh "Using trim-quality (no-trim-quality)\t$trimQuality\n";
    print $pfh "Using verbose (quit)\t$verbose\n";
    print $pfh "Using test\t$test\n";
    print $pfh "Using help\t$help\n";
    close $pfh;
    
    
    my $encoder=get_quality_encoder($fastqtype);
    
    
    if(-e $input2)
    {
        print "Found an existing file for the second read; Switching to paired-read mode\n";
        
        MainProc::processPE($input1,$input2,$output,$encoder,$trimQuality,$qualThreshold,$processStep,$minLength,$discardRemainingNs,$verbose);
    }
    else
    {
        print  "Single-read mode\n";
        MainProc::processSE($input1,$output,$encoder,$trimQuality,$qualThreshold,$processStep,$minLength,$discardRemainingNs,$verbose);
        
    }
    exit;
    
}

    
{
        package MainProc;
        use strict;
        use warnings;
        
        sub processSE
        {
            my $input=shift;
            my $output=shift;
            my $encoder=shift;
            my $trimQuality=shift;
            my $qualThreshold=shift;
            my $processStep=shift;
            my $minLength=shift;
            my $discardRemainingNs =shift;
            my $verbose=shift;
            
            my $fastqr=FastqReader->new($input);
            open my $ofh, ">", $output or die "Could not open output file";
            my($count5ptrims,$count3ptrims,$countRemainingNdiscards,$countQualityTrims,$countLengthDiscard)=(0,0,0,0,0);
            my $countProcessed=0;
            my $countPasFiltering=0;
            
            my $rld=[]; #read length distribution
        
            
            while(1)
            {
                my($firstheader,$nucleotide,$secondheader,$quality)=$fastqr->nextRead();
                last unless $firstheader;
                die "Error in sequence $firstheader; nucleotide sequence and quality sequence do not have equal length $nucleotide - $quality\n" unless length($nucleotide) ==length($quality);

                
                $countProcessed++;
                print "Processed $countProcessed reads\n" if($verbose && ($countProcessed % $processStep)==0); 
                
                my($c5ptrims,$c3ptrims)=(0,0);
                ($nucleotide,$quality,$c5ptrims,$c3ptrims)=Utility::trimNs($nucleotide,$quality);
                $count5ptrims+=$c5ptrims;
                $count3ptrims+=$c3ptrims;
                
                # kick the read out if it still contains any Ns;
                 if($discardRemainingNs && $nucleotide=~m/N/i)
                 {
                    $countRemainingNdiscards++;
                    next;
                 }
                
                
                if($trimQuality)
                {
                    # do the quality trimming using the encoder chosen by the user
                    my $ctrimqual=0;
                    ($nucleotide,$quality,$ctrimqual)=Utility::trimQualityMott($nucleotide, $quality, $qualThreshold, $encoder);
                    $countQualityTrims+=$ctrimqual;
                }
                
                # kick the read out if to short
                if(length($nucleotide) < $minLength)
                {
                    $countLengthDiscard++;
                    next;
                }
                
                #update read length distribution
                $rld->[length($nucleotide)]++;
                
                $countPasFiltering++;
                Utility::printFastq($ofh,$firstheader,$nucleotide,$secondheader,$quality);
            }
            
            if($verbose)
            {
                print "\n\n";
                print "FINISHED: end statistics\n";
                print "Reads processed: $countProcessed\n";
                print "Reads passed filtering: $countPasFiltering\n";
                print "5p poly-N sequences trimmed: $count5ptrims\n";
                print "3p poly-N sequences trimmed: $count3ptrims\n";
                print "Reads discarded during 'remaining N filtering': $countRemainingNdiscards\n";
                print "Reads discarded during length filtering: $countLengthDiscard\n";
                print "Count sequences trimed during quality filtering: $countQualityTrims\n" if $trimQuality;
                print "\nRead length distribution\n";
                print "length\tcount\n";
                for(my $i=1; $i<@$rld; $i++)
                {
                    my $val=$rld->[$i];
                    next unless $val;
                    print "$i\t$val\n";
                }
            }
        }
        
        sub processPE
        {
            my $input1=shift;
            my $input2=shift;
            my $outputPrefix=shift;
            my $encoder=shift;
            my $trimQuality=shift;
            my $qualThreshold=shift;
            my $processStep=shift;
            my $minLength=shift;
            my $discardRemainingNs =shift;
            my $verbose=shift;
            
            my $output1     =$outputPrefix."_1";
            my $output2     =$outputPrefix."_2";
            my $outputse    =$outputPrefix."_SE";
            
            open my $ofh1, ">", $output1 or die "Could not open output file";
            open my $ofh2, ">", $output2 or die "Could not open output file";
            open my $ofhse, ">", $outputse or die "Could not open output file";
            
            my $fqr1=FastqReader->new($input1);
            my $fqr2=FastqReader->new($input2);
            
            
            my($count5ptr1,$count3ptr1,$countRemainingNdiscards1,$countQualityTrims1,$countLengthDiscard1) = (0,0,0,0,0);
            my($count5ptr2,$count3ptr2,$countRemainingNdiscards2,$countQualityTrims2,$countLengthDiscard2) = (0,0,0,0,0);
            my($paired,$single,$read1Passing,$read2Passing)=(0,0,0,0);

            my $countProcessed=0;
            my $rld1=[]; #read length distribution
            my $rld2=[];
        
        
            while(1)
            {
                my($fheader1,$nuc1,$sheader1,$qual1)=$fqr1->nextRead();
                my($fheader2,$nuc2,$sheader2,$qual2)=$fqr2->nextRead();
                
                if(($fheader1 and not $fheader2) or ($fheader2 and not $fheader1))
                {
                    die "paired end files do not have equal length";
                }
                last unless $fheader1;
                
                die "Error in sequence $fheader1; nucleotide sequence and quality sequence do not have equal length $nuc1 - $qual1\n" unless length($nuc1) ==length($qual1);
                die "Error in sequence $fheader2; nucleotide sequence and quality sequence do not have equal length $nuc2 - $qual2\n" unless length($nuc2) ==length($qual2);
                
                my($keep1,$keep2)=(1,1);
                
 
                
                
                $countProcessed++;
                print "Processed $countProcessed pairs\n" if($verbose && ($countProcessed % $processStep)==0); 
                
                my($c5ptr1,$c3ptr1)=(0,0);
                my($c5ptr2,$c3ptr2)=(0,0);
                ($nuc1,$qual1,$c5ptr1,$c3ptr1)=Utility::trimNs($nuc1,$qual1);
                $count5ptr1+=$c5ptr1;
                $count3ptr1+=$c3ptr1;
                ($nuc2,$qual2,$c5ptr2,$c3ptr2)=Utility::trimNs($nuc2,$qual2);
                $count5ptr2+=$c5ptr2;
                $count3ptr2+=$c3ptr2;
                
                
                # kick the read out if it still contains any Ns;
                 if($discardRemainingNs && $nuc1=~m/N/i)
                 {
                    $countRemainingNdiscards1++;
                    $keep1=0;
                 }
                if($discardRemainingNs && $nuc2=~m/N/i)
                 {
                    $countRemainingNdiscards2++;
                    $keep2=0;
                 }
                 
                
                
                if($trimQuality)
                {
                    # do the quality trimming using the encoder chosen by the user
                    my ($ctrimqual1,$ctrimqual2)=(0,0);
                    ($nuc1,$qual1,$ctrimqual1)=Utility::trimQualityMott($nuc1, $qual1, $qualThreshold, $encoder) if $keep1 and $nuc1;
                    $countQualityTrims1+=$ctrimqual1;
                    
                    ($nuc2,$qual2,$ctrimqual2)=Utility::trimQualityMott($nuc2, $qual2, $qualThreshold, $encoder) if $keep2 and $nuc2;
                    $countQualityTrims2+=$ctrimqual2;
                }
                
                # kick the read out if to short
                if(length($nuc1) < $minLength)
                {
                    $countLengthDiscard1++;
                    $keep1=0;
                }
                if(length($nuc2) < $minLength)
                {
                    $countLengthDiscard2++;
                    $keep2=0;
                }
                
                
                #update read length distribution
                $rld1->[length($nuc1)]++ if $keep1;
                $rld2->[length($nuc2)]++ if $keep2;
                $read1Passing++ if $keep1;
                $read2Passing++ if $keep2;
                
                if($keep1 and $keep2)
                {
                    $paired++;
                    Utility::printFastq($ofh1,$fheader1,$nuc1,$sheader1,$qual1);
                    Utility::printFastq($ofh2,$fheader2,$nuc2,$sheader2,$qual2);
                }
                elsif($keep1 or $keep2)
                {
                    $single++;
                    if($keep1)
                    {
                        Utility::printFastq($ofhse,$fheader1,$nuc1,$sheader1,$qual1);
                    }
                    elsif($keep2)
                    {
                        Utility::printFastq($ofhse,$fheader2,$nuc2,$sheader2,$qual2);
                    }
                    else
                    {
                        die "impossible";
                    }
                }
            }
            
            if($verbose)
            {
                print "\n\n";
                print "FINISHED: end statistics\n";
                print "Read-pairs processed: $countProcessed\n";
                print "Read-pairs trimmed in pairs: $paired\n";
                print "Read-pairs trimmed as singles: $single\n";

                #my($count5ptr1,$count3ptr1,$countRemainingNdiscards1,$countQualityTrims1,$countLengthDiscard1) = (0,0,0,0,0);
                print "\n\nFIRST READ STATISTICS\n";
                print "First reads passing: $read1Passing\n";
                print "5p poly-N sequences trimmed: $count5ptr1\n";
                print "3p poly-N sequences trimmed: $count3ptr1\n";
                print "Reads discarded during 'remaining N filtering': $countRemainingNdiscards1\n";
                print "Reads discarded during length filtering: $countLengthDiscard1\n";
                print "Count sequences trimed during quality filtering: $countQualityTrims1\n" if $trimQuality;
                print "\nRead length distribution first read\n";
                print "length\tcount\n";
                for(my $i=1; $i<@$rld1; $i++)
                {
                    my $val=$rld1->[$i];
                    next unless $val;
                    print "$i\t$val\n";
                }
                
                print "\n\nSECOND READ STATISTICS\n";
                print "Second reads passing: $read2Passing\n";
                print "5p poly-N sequences trimmed: $count5ptr2\n";
                print "3p poly-N sequences trimmed: $count3ptr2\n";
                print "Reads discarded during 'remaining N filtering': $countRemainingNdiscards2\n";
                print "Reads discarded during length filtering: $countLengthDiscard2\n";
                print "Count sequences trimed during quality filtering: $countQualityTrims2\n" if $trimQuality;
                print "\nRead length distribution second read\n";
                print "length\tcount\n";
                for(my $i=1; $i<@$rld2; $i++)
                {
                    my $val=$rld2->[$i];
                    next unless $val;
                    print "$i\t$val\n";
                }
            }
        }
}
    
    
{
    package Utility;
    use strict;
    use warnings;
        
        sub printFastq
        {
            my $ofh=shift;
            my $header1=shift;
            my $nuc=shift;
            my $header2=shift;
            my $qual=shift;
            
            print $ofh $header1."\n";
            print $ofh $nuc."\n";
            print $ofh $header2."\n";
            print $ofh $qual."\n";
            
        }
        
        sub trimNs
        {
            my $nuc=shift;
            my $qual=shift;
            
            my($t5,$t3)=(0,0);
            
            if($nuc=~m/^(N+)/i)
            {
                my $l=length($1);
                
                $nuc=~s/^.{$l}//;
                $qual=~s/^.{$l}//;
                $t5=1;
            }
            
            if($nuc=~m/(N+)$/i)
            {
                my $l=length($1);
                
                $nuc=~s/.{$l}$//;
                $qual=~s/.{$l}$//;
                $t3=1;
            }
            
            die "length of quality and nucleotide sequence must be equal" if length($nuc) != length($qual); 
            return ($nuc,$qual,$t5,$t3);
        }
        
        sub trimQuality
        {
            my $nuc=shift;
            my $qual=shift;
            my $minQual=shift;
            my $encoder=shift;
            my $trimWindow=shift;
            my $trimq=0;
            
            my $totrim = _getQualityCutoff($encoder,$qual,$minQual,$trimWindow);
            
            if($totrim>0)
            {
                $nuc=~s/.{$totrim}$//;
                $qual=~s/.{$totrim}$//;
                $trimq=1;
            }
            
            die "length of quality and nucleotide sequence must be equal $nuc vs $qual" if length($nuc) != length($qual); 
            return ($nuc,$qual,$trimq);
        }
        
        sub trimQualityMott
        {
            my $nuc=shift;
            my $qual=shift;
            die "length of quality and nucleotide sequence has to be identical $nuc vs $qual" if length($nuc)!=length($qual);
            my $minQual=shift;
            my $encoder=shift;
            
            # get the quality hsps
            my $qualhsp=_getQualityHsp($encoder,$qual,$minQual);
            return ("","",1) unless $qualhsp; #return nothing if no hsp is found
            
            my $orileng=length($nuc);
            my $newleng=$qualhsp->{end}-$qualhsp->{start}+1;

            $nuc=substr($nuc,$qualhsp->{start},$newleng);
            $qual=substr($qual,$qualhsp->{start},$newleng);
            
            my $trimq=0;
            $trimq=1 if $orileng!=$newleng;
            die "length of quality and nucleotide sequence must be equal" if length($nuc) != length($qual); 
            return ($nuc,$qual,$trimq);
        }
        
        sub _getQualityHsp
        {
            my $encoder=shift;
            my $qual=shift;
            my $minQual=shift;
            
            my @qual=split //,$qual;
            my @hsps=();
            
            
            my $highscorestart=-1;
            my $highscore=0;
            my $highscoreend=0;
            my $activescore=0;
            for(my $i=0; $i<@qual; $i++)
            {
                my $q=$encoder->($qual[$i]);
                
                my $tosub=$q-$minQual; 
                
                
                $activescore+=$tosub;
                if($activescore>0)
                {
                    if($activescore>$highscore) # has to be larger than; larger or equal leads to inconsistent behaviour at the beginning of the quality hsp.
                    {
                        $highscore=$activescore;
                        $highscoreend=$i;                        
                    }
                    
                    $highscorestart=$i if $highscorestart==-1;
                }
                else
                {
                    if($highscore>0)
                    {
                        push @hsps,{
                            start=>$highscorestart,
                            score=>$highscore,
                            end=>$highscoreend
                        };
                    }
                    $highscorestart=-1;
                    $activescore=0;
                    $highscore=0;
                    $highscoreend=0;
                }
            }
            
            if($highscore>0)
            {
                push @hsps,{
                            start=>$highscorestart,
                            score=>$highscore,
                            end=>$highscoreend
                };
            }
            
            return undef unless @hsps;
            @hsps=sort {$a->{score}<=>$b->{score}} @hsps;
            return pop @hsps;
        }
        
        
        

        
        sub _getAverageQuality
        {
            # calculate the average quality for an array of quality values
            my $quals=shift;
            
            my $toc=0;
            foreach my $q (@$quals)
            {
                die"quality for value is not defined" unless defined($q);
                $toc+=$q;
            }
            return $toc/@$quals;
        }
        
        sub _getQualityCutoff
        {
            my $encoder=shift;
            my $qual=shift;
            my $minQual=shift;
            my $window=shift;
            die"window size has to be uneven" unless($window % 2);
            
            return 0 if length($qual)<$window;
            
            my @n=split//,$qual;
            
            #start position for counting and also the length of each of the two flanking sequences of the window;
            my $start=int($window/2);
            
            my @t=();
            my $inicount=$window-1;
            while($inicount)
            {
                push @t,$encoder->(shift@n);
                $inicount--;
            }
            unshift @t,0; #also add a dummy at the qeue
            
            my $toret=-1;
            POS: for(my $i=$start; @n>0; $i++)
            {
                shift @t;
                push @t,$encoder->(shift @n);
                my $avq=_getAverageQuality(\@t);
                if($avq < $minQual)
                {
                    $toret=$i;
                    last POS;
                }

            }
            
            return 0 if $toret== -1;
            
            return length($qual)-$toret;
        }
    }

    
{
    package FastqReader;
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        open my $ofh,"<$file" or die "Could not open file handle";
        
        return bless {
            file=>$file,
            fh=>$ofh,
            buffer=>[]
        },__PACKAGE__;
    }

    
    sub nextRead
    {
        
        my $self=shift;

        my $firstheader=$self->nextLine;
        my $nucleotide=$self->nextLine;
        my $secondheader=$self->nextLine;
        my $quality=$self->nextLine;
        
        return undef unless $quality;        
        chomp $firstheader;
        chomp $nucleotide;
        chomp $secondheader;
        chomp $quality;
        
        die "first header must start with an @; line: $firstheader" unless $firstheader =~ m/^[@]/;
        die "second header must start with an +; line: $secondheader" unless $secondheader=~m/^[+]/;
        
        return ($firstheader,$nucleotide,$secondheader,$quality);
    }
    
    sub nextLine
    {
        my $self=shift;
        my $fh=$self->{fh};
        
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub bufferLine
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }


}

{
    package TestTrim;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin";
    use deviaTE_trim;
    use Test;

    

    
     
     sub runTrimTests
     {
        testNTrimming();
        #testTrimQuality();
        testTrimQualityMott();
        exit;
     }
 
     sub testNTrimming
     {
        my($t1,$t2,$c5,$c3)=Utility::trimNs("AAAN","xxxt");
        is($t1,"AAA","trim Ns: AAAN trimmed to AAA");
        is($t2,"xxx","trim Ns: quality trimmed accordingly");
        is($c5,0,"trim Ns: flag indicating 5p trim correct");
        is($c3,1,"trim Ns: flag indicating 3p trim correct");
        
        ($t1,$t2)=Utility::trimNs("AAANNNNNNNNNN","xxxtttttttttt");
        is($t1,"AAA","trim Ns: AAANNNNNNNNNN trimmed to AAA");
        is($t2,"xxx","trim Ns: quality trimmed accordingly");
        
        ($t1,$t2)=Utility::trimNs("NAAA","txxx");
        is($t1,"AAA","trim Ns: NAAA trimmed to AAA");
        is($t2,"xxx","trim Ns: quality trimmed accordingly");
        
        ($t1,$t2,$c5,$c3)=Utility::trimNs("NNNNNAAA","tttttxxx");
        is($t1,"AAA","trim Ns: NNNNNAAA trimmed to AAA");
        is($t2,"xxx","trim Ns: quality trimmed accordingly");
        is($c5,1,"trim Ns: flag indicating 5p trim correct");
        is($c3,0,"trim Ns: flag indicating 3p trim correct");
        
        ($t1,$t2,$c5,$c3)=Utility::trimNs("NNNNNAAANNNNN","tttttxxxttttt");
        is($t1,"AAA","trim Ns: NNNNNAAANNNNN trimmed to AAA");
        is($t2,"xxx","trim Ns: quality trimmed accordingly");
        is($c5,1,"trim Ns: flag indicating 5p trim correct");
        is($c3,1,"trim Ns: flag indicating 3p trim correct");
     }
     
    sub testTrimQualityMott
     {
        my $t=get_quality_encoder("illumina");
        
        my($t1,$t2,$c);
        #U=21
        #T=20
        #S=19
        
        
        ($t1,$t2,$c)=Utility::trimQualityMott("AAAA","TTTT",19,$t);
        is($t1,"AAAA","trim quality: min quality 20; AAAA with quality of TTTT not trimed");
        is($t2,"TTTT","trim quality: quality was not trimmed either");
        is($c,0,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("AAAAA","UUUTS",20,$t);
        is($t1,"AAA","trim quality: min quality 20; AAAAA with quality of UUUTS trimmed to AAA");
        is($t2,"UUU","trim quality: quality accordingly trimmed to UUU");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("AAAA","TTTT",20,$t);
        is($t1,"","trim quality: trimmed correctly");
        is($t2,"","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("TTAAAATT","TTUUUUTT",20,$t);
        is($t1,"AAAA","trim quality: sequence trimmed correctly");
        is($t2,"UUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("TTAAAATTCCCCC","SSUUUUAAUUUUU",20,$t);
        is($t1,"CCCCC","trim quality: sequence trimmed correctly");
        is($t2,"UUUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("TTAAAATTCCCCCG","SSUUUUAAUUUUUS",20,$t);
        is($t1,"CCCCC","trim quality: sequence trimmed correctly");
        is($t2,"UUUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("TTAAAAATTCCCC","SSUUUUUAAUUUU",20,$t);
        is($t1,"AAAAA","trim quality: sequence trimmed correctly");
        is($t2,"UUUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
         
         
        # can a bases with a quality of 19 be skipped if the score is high enough
        ($t1,$t2,$c)=Utility::trimQualityMott("TTAAAAATTCCCC","SSUUUTTSSUUUT",20,$t);
        is($t1,"AAAAATTCCC","trim quality: sequence trimmed correctly");
        is($t2,"UUUTTSSUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("TAAAAATTCCCCT","SUUUTTSSUUUTS",20,$t);
        is($t1,"AAAAATTCCC","trim quality: sequence trimmed correctly");
        is($t2,"UUUTTSSUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("TAAAAATTCCCCTTTTT","SUUUTTSSUUUTSSUTS",20,$t);
        is($t1,"AAAAATTCCC","trim quality: sequence trimmed correctly");
        is($t2,"UUUTTSSUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQualityMott("TTTTTAAAAATTCCCCT","UTTSSUUUTTSSUUUTS",20,$t);
        is($t1,"AAAAATTCCC","trim quality: sequence trimmed correctly");
        is($t2,"UUUTTSSUUU","trim quality: quality was also trimmed correctly");
        is($c,1,"trim quality: trimmed flag set correctly");

     }
     
     sub testTrimQuality
     {
        

        my $t=get_quality_encoder("illumina");
        
        my($t1,$t2,$c);
        
        
        ($t1,$t2,$c)=Utility::trimQuality("AAAA","TTTT",20,$t,1);
        is($t1,"AAAA","trim quality: min quality 20; AAAA with quality of TTTT not trimed");
        is($t2,"TTTT","trim quality: quality was not trimmed either");
        is($c,0,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQuality("AAAAA","TTTSS",20,$t,1);
        is($t1,"AAA","trim quality: min quality 20; AAAAA with quality of TTTSS trimmed to AAA");
        is($t2,"TTT","trim quality: quality accordingly trimmed to TTT");
        is($c,1,"trim quality: trimmed flag set correctly");

        
        ($t1,$t2)=Utility::trimQuality("AAAA","STTT",20,$t,1);
        is($t1,"","trim quality: min quality 20; AAAA with quality of STTT empty");
        is($t2,"","trim quality: quality string is empty as well");
        
        # set quality neighborhood to three
        ($t1,$t2,$c)=Utility::trimQuality("AAAAAA","TTTTSS",20,$t,3);
        is($t1,"AAA","trim quality: min quality 20; window size 3; AAAAA with quality of TTTSS trimmed to AAA");
        is($t2,"TTT","trim quality: quality accordingly trimmed to TTT");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        # set quality neighborhood to three
        ($t1,$t2,$c)=Utility::trimQuality("AAAAA","TTTSS",20,$t,3);
        is($t1,"AA","trim quality: min quality 20; window size 3; AAAAA with quality of TTTSS trimmed to AA");
        is($t2,"TT","trim quality: quality accordingly trimmed to TT");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        # set quality neighborhood to three
        ($t1,$t2,$c)=Utility::trimQuality("AAAAA","TTUSS",20,$t,3);
        is($t1,"AAA","trim quality: min quality 20; window size 3; AAAAA with quality of TTUSS trimmed to AAA");
        is($t2,"TTU","trim quality: quality accordingly trimmed to TTT");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        
        # gradient quality T=20 A=
        ($t1,$t2,$c)=Utility::trimQuality("AAAAATTTTTCCCCCGGGGGAAAAAT","ZYXWVUTSRQPONMLKJIHGFEDCBA",20,$t,1);
        is($t1,"AAAAATT","trim quality: sequence trimmed correctly; window size 1");
        is($t2,"ZYXWVUT","trim quality: quality trimmed accordingly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQuality("AAAAATTTTTCCCCCGGGGGAAAAAT","ZYXWVUTSRQPONMLKJIHGFEDCBA",20,$t,3);
        is($t1,"AAAAATT","trim quality: sequence trimmed correctly; window size 3");
        is($t2,"ZYXWVUT","trim quality: quality trimmed accordingly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQuality("AAAAATTTTTCCCCCGGGGGAAAAAT","ZYXWVUTSRQPONMLKJIHGFEDCBA",20,$t,5);
        is($t1,"AAAAATT","trim quality: sequence trimmed correctly; window size 5");
        is($t2,"ZYXWVUT","trim quality: quality trimmed accordingly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        ($t1,$t2,$c)=Utility::trimQuality("AAAAATTTTTCCCCCGGGGGAAAAAT","ZYXWVUTSRQPONMLKJIHGFEDCBA",20,$t,7);
        is($t1,"AAAAATT","trim quality: sequence trimmed correctly; window size 7");
        is($t2,"ZYXWVUT","trim quality: quality trimmed accordingly");
        is($c,1,"trim quality: trimmed flag set correctly");
        
        # test if average is calculated
        ($t1,$t2,$c)=Utility::trimQuality("AAAAATTTTT","UUSUUSSSSS",20,$t,3);
        is($t1,"AAAAA","trim quality: min quality 20; window size 3; AAAAA with quality of UUSUU trimmed to AAA");
        is($t2,"UUSUU","trim quality: quality accordingly trimmed to UUSUU");
        is($c,1,"trim quality: trimmed flag set correctly");
          
     }
}



=head1 NAME

trim-fastq.pl - Perform quality filtering of fastq files

=head1 SYNOPSIS

 # Minimum argument call; single read mode
 trim-fastq.pl --input1 input.fastq --output output.fastq
 
# Minimum argument call; single read mode
 trim-fastq.pl --input1 input_1.fastq --input2 input_2.fastq --output output_prefix

=head1 OPTIONS

=over 4

=item B<--input1>

The input file, or the input file of the first read, in fastq format. Mandatory parameter

=item B<--input2>

The input file of the second read, in fastq format. In case this file is provided the software will switch to paired read mode instead of single read mode

=item B<--output>

The output file. Will be in fastq. Mandatory parameter

=item B<--quality-threshold>

minimum average quality; A modified Mott algorithm is used for trimming; the threshold is used for calculating a score: score = quality_at_base - threshold; default=20

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

=item B<--discard-internal-N>

flag, if set reads having internal Ns will be discarded; default=off

=item B<--min-length>

The minimum length of the read after trimming; default=40


=item B<--no-trim-quality>

toggle switch: switch of trimming of quality

=item B<--quit>

suppress output to stdout (console)

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

The script removes 'N' - characters at the beginning and the end of the provided reads. If any remaining 'N' characters are found the read is discarded.
Quality removal is done using a modified Mott-algorithm;
For each base a score is calculated: score_base = quality_base - threshold

While scanning along the read a running sum of this score is calculatet; If the score drops below zero the score is set to zero;
The highest scoring region of the read is finally reported;

=head2 INPUT

Input must be a fastq file. No line breaks in the nucleotide or quality sequence are allowed. No empty lines in the file are allowed.
The fastq file produced by the Illumina pipeline can be directyl use. Example input:

    @read1
    NNNNNAAAAAAAAAATTTTTTTTTTAAAAAAAAAANNNNNNNNNN
    +read1
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    @read2
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT
    +read2
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUSSSSSSSSSS



=head2 OUTPUT

the output file will be in fastq. The quality sequence will be provided in the same format as in the input file (apart from trimming);
    
     @read1
     AAAAAAAAAATTTTTTTTTTAAAAAAAAAA
     +read1
     UUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
     @read2
     AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
     +read2
     UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

=cut