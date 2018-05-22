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
my $windowsize;
my $filter="";



GetOptions(
            "input=s"		=>\$input,
            "filter=s"          =>\$filter,
            "output=s"		=>\$output,
            "window-size=i"      =>\$windowsize,
            "help"              =>\$help
        );


pod2usage(-verbose=>2) if $help;
pod2usage(-verbose=>1,-message=>"No input file specified") unless -e $input;
pod2usage(-verbose=>1, -message=>"No output file specified") unless $output;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using filter\t$filter\n";
print $pfh "Using window-size\t$windowsize\n";
print $pfh "Using help\t$help\n";
close $pfh;


$windowsize=Utility::infer_window_size($input) unless $windowsize;
#calculate the offset; wiggle is 1-based
my $offset=int($windowsize/2)-1;

my $filt;
$filt= {map{ $_,1} split /\s/,$filter} if $filter;

open my $ifh,"<", $input or die "Could not open input file";
open my $ofh,">",$output or die "Could not open output file";




#print track
print $ofh <<PERLSUCKS;
[general]
glyph = xyplot
graph_type=boxes
fgcolor = black
bgcolor = red
height=200
min_score=0

max_score=0.1

label=1
scale=left
key=var
PERLSUCKS



my $activechr="";
my $chrHash={};
while(my $l=<$ifh>)
{
    chomp $l;
    my($chr,$pos,undef,undef,$measure)=split /\t/,$l;
    
    if($filt)
    {
        next unless exists($filt->{$chr});
    }
    if($chr ne $activechr)
    {

        print $ofh "reference=$chr\n";
        die "Chromosome $chr occured several times; Please make sure the file is properly sorted" if exists($chrHash->{$chr});
        $chrHash->{$chr}=1;
        $activechr=$chr;
    }
    
    next if $measure eq "na";
    my $start=$pos-$offset;
    my $end=$pos+int($windowsize/2);
#2L	pi	5202..5301	score=0.004862034
#2L	pi	5302..5401	score=0.004521015
#2L	pi	5402..5501	score=0.005496933
#2L	pi	5502..5601	score=0.004151401
#2L	pi	5602..5701	score=0.003701141
#2L	pi	5702..5801	score=0.003530197
#2L	pi	5802..5901	score=0.003165151
#2L	pi	5902..6001	score=0.002726609
    print $ofh "$chr\tvar\t$start..$end\tscore=$measure\n";   
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

VarSliding2Flybase.pl - Converts the output of the Variance-slider.pl into a Flybase compatible file

=head1 SYNOPSIS

=head1 OPTIONS

=over 4

=item B<--input>

The input file; Mandatory parameter

=item B<--output>

The output file. Mandatory parameter


=item B<--window-size>

The window size; Per default the script will attempt to infer the correct window-size; If this fails the user may provide the correct one; default=""

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

A Flybase formated output file:

 [general]
 glyph = xyplot
 graph_type=boxes
 fgcolor = black
 bgcolor = red
 height=200
 min_score=na 

 max_score=na

 label=1
 scale=left
 key=pi

 reference=2L
 2L	pi	4902..5001	score=0.006027673
 2L	pi	5002..5101	score=0.005171469
 2L	pi	5102..5201	score=0.004855847
 2L	pi	5202..5301	score=0.004862034
 2L	pi	5302..5401	score=0.004521015
 2L	pi	5402..5501	score=0.005496933
 2L	pi	5502..5601	score=0.004151401
 2L	pi	5602..5701	score=0.003701141

=cut


