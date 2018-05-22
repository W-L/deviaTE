
{
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";


our $verbose=1;
my $maskregion=5;
my $input="";
my $output="";
my $help=0;
my $test=0;


#--gtf /Users/robertkofler/testfiles/3R-3entries.gtf --input /Users/robertkofler/testfiles/3R.pileup --output /Users/robertkofler/testfiles/output/3R-filtered.pileup

GetOptions(
    "input=s"           =>\$input,
    "mask-region=i"     =>\$maskregion,
    "output=s"          =>\$output,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";


pod2usage(-verbose=>2) if $help;
MaskTest::run_Tests() if $test;
pod2usage(-msg=>"Could not find sam file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using mask-region\t$maskregion\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;


open my $ifh, "<",$input or die "Could not open sam file";

open my $ofh,">",$output or die "Could not open output file";

my $count_total=0;
my $count_indel=0;

while(my $line=<$ifh>)
{
    chomp $line;
    if($line=~/^@/)
    {
        print $ofh "$line\n";
        next;
    }
    $count_total++;
    
    my @a=split /\t/,$line;
    unless($a[5]=~/[DI]/)
    {
        print $ofh "$line\n";
        next;
    }
    
    my $cigar=$a[5];
    my $seq=$a[9];
    
    my $maskedseq=Utility::mask_sequence($cigar,$seq,$maskregion);
    $count_indel++;
    
    $a[9]=$maskedseq;
    my $newline=join("\t",@a);
    print $ofh "$newline\n";

}
print "Processed $count_total sequences\n";
print "Masked indels in $count_indel sequences\n";
close $ofh;
exit;

}

{
    package Utility;
    use strict;
    use warnings;
    
    sub mask_sequence
    {
        my $cigar=shift;
        my $seq=shift;
        my $maskregion=shift;
        
        my $orileng=length($seq);
    
        # indelpositions
        my $indelpositions=Utility::get_indel_positions($cigar);
        
        # masking and control
        my $maskedseq=Utility::_mask_sequence_from_indels($seq,$indelpositions,$maskregion);
        my $maskleng=length($maskedseq);
        die "sequence after masking does not have its original length" unless $maskleng==$orileng;

        return $maskedseq;
    }

    sub _mask_sequence_from_indels
    {
        my $seq=shift;
        my $indels=shift;
        my $tomask=shift;
        my $last=length($seq)-1;
    
        my $masked=$seq;
        foreach my $id(@$indels)
        {
            my $pos=$id->{pos};
            my $id_len=$id->{len};

            my $start=$pos-$tomask;
            my $end=$pos+$tomask-1+$id_len;
            $start=0 if $start<0;
            $end=$last if $end>$last;
            my $leng=$end-$start+1;
            
            my $toreplace ="N"x$leng;
            
            substr($masked,$start,$leng,$toreplace);
            
        }
    return $masked;        
    }
    
    sub get_indel_positions
    {
        my $cigar=shift;
        die "Can not deal with hard clipping or padding"  if $cigar=~/[HP]/;
        my (@e)=$cigar=~/(\d+[MSDIN])/g;
        my $indels=[];
        my $runningcounter=0;
        
        foreach my $p (@e)
        {
            if($p=~/^(\d+)[SM]$/)
            {
                $runningcounter+=$1;
            }
            elsif($p=~/^(\d+)I/)
            {
             push @$indels,{
                            pos=>$runningcounter,
                            len=>$1
                           };
             $runningcounter+=$1;
            }
            elsif($p=~/^(\d+)[ND]/)
            {
                push @$indels,
                {
                    pos=>$runningcounter,
                    len=>0
                };
            }
            else
            {
                die "not allowed";
            }
        }
        return $indels;

    }
}


{
    package MaskTest;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/../Modules";
    use Test;

    sub run_Tests
    {
        test_masking_from_cigar();
        exit;
    }
    
    sub test_masking_from_cigar
    {

        my $masked;
        
        $masked=Utility::mask_sequence("3M","AAA",1);
        is($masked,"AAA","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("2M3D2M","AAAA",1);
        is($masked,"ANNA","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("2M2I2M","AAAAAA",1);
        is($masked,"ANNNNA","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("3D4M","AAAA",1);
        is($masked,"NAAA","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("4M1D","AAAA",1);
        is($masked,"AAAN","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("1I3M","AAAA",1);
        is($masked,"NNAA","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("3M1I","AAAA",1);
        is($masked,"AANN","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("5M1D5M","AAAAATTTTT",2);
        is($masked,"AAANNNNTTT","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("5M1D1M1D4M","AAAAATTTTT",2);
        is($masked,"AAANNNNNTT","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("5M1D1M1I3M","AAAAATTTTT",2);
        is($masked,"AAANNNNNNT","masking algorithm: masking ok");
        
        $masked=Utility::mask_sequence("4M1I1M1D4M","AAAAATTTTT",2);
        is($masked,"AANNNNNNTT","masking algorithm: masking ok");
        
        
        my $bla=0;
        
    }

}





=head1 NAME

perl mask-sam-indelregions.pl - Mask the indel regions in a sam file

=head1 SYNOPSIS

perl mask-sam-indelregions.pl  --input input.sam --output indelmasked.sam

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the sam format; Mandatory.

=item B<--output>

A sam file where the region around indels have been masked by 'N'-characters; Thus this regions will not be used for SNP identification from a pileup file; Mandatory

=item B<--mask-region>

The number of 5' and 3' basepairs which should be masked around the indel; The actually masked region will be 2 x C<--mask-region> (assuming the indel is not at the end of a read) ;default=5

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input sam


=head2 Output sam

  
=cut



