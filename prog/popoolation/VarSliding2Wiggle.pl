#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;
use File::Path;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use Getopt::Long;
use Pod::Usage;





my $input;
my $output;
my $help=0;
my $trackname="unknown";
my $windowsize;
my $ucscfilter="";
my $ucscprepend="";



GetOptions(
            "input=s"		=>\$input,
            "ucsc-filter=s"     =>\$ucscfilter,
            "ucsc-prepend=s"    =>\$ucscprepend,
            "output=s"		=>\$output,
            "window-size=i"      =>\$windowsize,
            "trackname=s"       =>\$trackname,
            "help"              =>\$help
        );


pod2usage(-verbose=>2) if $help;
pod2usage(-verbose=>1,-message=>"No input file specified") unless -e $input;
pod2usage(-verbose=>1, -message=>"No output file specified") unless $output;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using ucsc-filter\t$ucscfilter\n";
print $pfh "Using ucsc-prepend\t$ucscprepend\n";
print $pfh "Using window-size\t$windowsize\n";
print $pfh "Using trackname\t$trackname\n";
print $pfh "Using help\t$help\n";
close $pfh;

unless ($output=~m/\.wig/)
{
    $output.=".wig";
}

$windowsize=Utility::infer_window_size($input) unless $windowsize;
#calculate the offset; wiggle is 1-based
my $offset=int($windowsize/2)-1;

my $ucfilt;
$ucfilt= {map{ $_,1} split /\s/,$ucscfilter} if $ucscfilter;

open my $ifh,"<", $input or die "Could not open input file";
open my $ofh,">",$output or die "Could not open output file";


#print track
print $ofh "track type=wiggle_0 name=\"$trackname\" visibility=full\n";

my $activechr="";
my $chrHash={};
while(my $l=<$ifh>)
{
    chomp $l;
    my($chr,$pos,undef,undef,$measure)=split /\t/,$l;
    
    if($ucfilt)
    {
        next unless exists($ucfilt->{$chr});
    }
    if($chr ne $activechr)
    {
        my $printchr=$ucscprepend.$chr;
        print $ofh "variableStep chrom=$printchr span=$windowsize\n";
        die "Chromosome $chr occured several times; Please make sure the file is properly sorted" if exists($chrHash->{$chr});
        $chrHash->{$chr}=1;
        $activechr=$chr;
    }
    
    next if $measure eq "na";
    $pos-=$offset;
    print $ofh "$pos\t$measure\n";   
}



exit;


{
    package Utility;
    use strict;
    use warnings;
    
    sub infer_window_size
    {
        my $file=shift;
        
        print "Infering window size (~span in wiggle file)\n";
        print "Please not that wiggle only allows for non overlapping windows. This window size here will thus correspond to the step-size\n";
        open my $ifh , "<", $file or die "Could not open input file";
        
        while(1)
        {
            my $l1=<$ifh>;
            my $l2=<$ifh>;
            my ($c1,$p1)=split /\t/,$l1;
            my ($c2,$p2)=split /\t/,$l2;
            if($c1 eq $c2)
            {
                my $insertsize=$p2-$p1;
                print "Infered windowsize: $insertsize\n";
                print "If this is not correct please use the option --window-size\n";
                return $insertsize;
            }
        }
        
        die "unable to infer insert size";


    }
    
    
    
}

=head1 NAME

VarSliding2Wiggle.pl - Converts the output of the Variance-slider.pl into a wiggle file

=head1 SYNOPSIS

=head1 OPTIONS

=over 4

=item B<--input>


The input file; Mandatory parameter


=item B<--output>

The output file. Mandatory parameter

=item B<--trackname>

The name of the track. This information may for example be displayed in the IGV. default=unknown

=item B<--window-size>

The window size; Per default the script will attempt to infer the correct window-size; If this fails the user may provide the correct one; default=""

=item B<--ucsc-filter>

The UCSC genome browser only accepts certain chromosome ids. Any attempts loading a wiggle file containing not recognised chromosome ids will result in an error; this option allows to filter for certain user provided
chromosome ids; for example C<--ucsc-filter "2L 2R 3L 3R 4 X"; default=""

=item B<--ucsc-prepend>

For the UCSC genome browser, chromsomes IDs have to be formated according to the requirements of ucsc. For example chr2L has to be provided instead of 2L for D.mel.
This option allows to prepend the specified string to every chromosome id; This step is applied after the filtering using C<--ucsc-filter>; default=""

=item B<--help>

Display the help pages.

=back

=head1 DESCRIPTION

=head2 Input

Input will be the output of the variance slider. For example
 
 2L      5500    30      0.916   -0.003390146
 2L      6500    10      1.000   -0.001091447
 2L      7500    8       1.000   0.000033510

=head2 Output

A wiggle file which can be used with the Integrated Genomics Viewer (IGV).
If the names of the chromomosomes are correct this file can also be used directly with the UCSC Genome browser;
For more information see the section reformating for UCSC.

 track type=wiggle_0 name="unknown" visibility=full
 variableStep chrom=2L span=1000
 5001    -0.003390146
 6001    -0.001091447
 7001    0.000033510
 8001    -0.001735518
 9001    -0.001217919

=head2 Reformating for UCSC

The UCSC genome browser has very strict requirements for the chromosome ids.
For example UCSC requires for D. melanogaster "chom2L" instead of "2L";
Furthermore chromosomes which are unknown to UCSC will result in an error.
This may for example happen if Wolbachia is present in the D. melanogaster reference sequence.

For this reasons it is necessary to:

a.) filter for known contigs (use C<--ucsc-filter>)

b.) convert the chromosome id's into the format required by ucsc (use C<--ucsc-prepend>)





=cut


