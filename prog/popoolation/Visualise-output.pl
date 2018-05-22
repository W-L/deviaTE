#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $tc="";
my $infile="";
my $outfile="";
my $tempr="/tmp/rinput.r";
my $lab="measure";
my $scale_equal=0;
my $outps=0;
my $ymin=undef;
my $ymax=undef;
my $help=0;

GetOptions(
    "input=s"       =>\$infile,
    "output=s"      =>\$outfile,
    "chromosomes=s" =>\$tc,
    "ylab=s"        =>\$lab,
    "ymin=f"        =>\$ymin,
    "ymax=f"        =>\$ymax,
    "scale-equal"   =>\$scale_equal,
    "ps"           =>\$outps,
    "help"          =>\$help
) or die "do not recognise option $!";

pod2usage(-verbose=>2) if $help;
pod2usage(-msg=>"Specify an input file",-verbose=>1) unless -e $infile;
pod2usage(-msg=>"Specify an output file",-verbose=>1) unless $outfile;
pod2usage(-msg=>"Specify the chromosomes",-verbose=>1) unless $tc;

my $paramfile=$outfile.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using infile\t$infile\n";
print $pfh "Using outfile\t$outfile\n";
print $pfh "Using chromosomes\t$tc\n";
print $pfh "Using ylab\t$lab\n";
print $pfh "Using ymin\t$ymin\n";
print $pfh "Using ymax\t$ymax\n";
print $pfh "Using scale-equal\t$scale_equal\n";
print $pfh "Using ps\t$outps\n";
print $pfh "Using help\t$help\n";
close $pfh;

my @chromosomes=split/\s/,$tc;

_printRCode($tempr,$infile,$outfile,\@chromosomes,$scale_equal,$lab,$ymin,$ymax,$outps);


system "R --vanilla <$tempr";
unlink $tempr;

    


exit;

sub _printRCode
{
    my $tempr=shift;
    my $input=shift;
    my $output=shift;
    my $tc=shift;
    my $scalequal=shift;
    my $ylab=shift;
    my $ymin=shift;
    my $ymax=shift;
    my $outps=shift;
    
    my $c=@$tc;
    my $columnes=3;
    my $rows= int(($c+2)/3);
    my $chromosomes=join "\",\"",@$tc;
    $chromosomes="c(\"".$chromosomes."\")";

    
open my $tfh,">$tempr" or die "Could not open output file";
print $tfh <<PERLSUCKS;
scale=$scalequal
t<-read.table("$input",sep="\t")
t\$V5[t\$V5=="na"]=NA
t\$V5=as.numeric(as.character(t\$V5))
PERLSUCKS
if($outps)
{
    print $tfh "postscript(file=\"$output\",onefile=TRUE,horizontal=FALSE)\n";
}
else
{
    print $tfh "pdf(\"$output\",width=15)\n";
}


print $tfh <<PERLSUCKS;
par(mfrow=c($rows,$columnes))
chroms<-$chromosomes
PERLSUCKS
if(defined($ymax))
{
    print $tfh "maxval=$ymax\n";
}
else
{
    print $tfh "maxval<-max(t[,5],na.rm=TRUE)\n";
}
if(defined($ymin))
{
    print $tfh "minval=$ymin\n";
}
else
{
    print $tfh "minval<-min(t[,5],na.rm=TRUE)\n";
}
print $tfh <<PERLSUCKS;
maxval
minval
xmaxtotal=max(t[,2],na.rm=TRUE)
for(chr in chroms)
{
    xmax=0
    xmaxindiv=max(t[t\$V1==chr,2])
    if(scale==1)
    {
        xmax=xmaxtotal
    }
    else{
        xmax=xmaxindiv
    }
    plot(t[t\$V1==chr,2],t[t\$V1==chr,5],type="l",xlab="position",ylab="$ylab",main=chr,xlim=c(0,xmax),ylim=c(minval,maxval))    
}
dev.off()
PERLSUCKS
close $tfh;


#dat <- read.table("$tempinput", sep="\t", header=TRUE)
#jpeg("$output");
#plot(dat)
#l<-lm(dat[,2]~dat[,1])
#abline(l,lty=1,lwd=2)
#abline(a=0,b=1, lty=2,lwd=2)  
#x<-as.character(summary(l))
#text(x[8],y=0.2,x=0.8)
#dev.off()
#PERLSUCKS

    
    
#    > t<-read.table("/Volumes/Macintosh HD/Users/robertkofler/analysis/martin/comparison_simvsmel/heterozygousity.sim-si3",sep="\t")
#> t$V3[t$V3=="na"]=NA
#> t$V3=as.numeric(as.character(t$V3))
#>par(mfrow=c(2,3))
#> plot(t[t$V1=="2L",2],t[t$V1=="2L",3],type="l",xlab="position",ylab="pi",main="2L",ylim=c(0,0.011))
#> plot(t[t$V1=="2R",2],t[t$V1=="2R",3],type="l",xlab="position",ylab="pi",main="2R",ylim=c(0,0.011))
#> summary(t)
#> plot(t[t$V1=="X",2],t[t$V1=="X",3],type="l",xlab="position",ylab="pi",main="X",ylim=c(0,0.011))
#> plot(t[t$V1=="3L",2],t[t$V1=="3L",3],type="l",xlab="position",ylab="pi",main="3L",ylim=c(0,0.011))
#> plot(t[t$V1=="3R",2],t[t$V1=="3R",3],type="l",xlab="position",ylab="pi",main="3R",ylim=c(0,0.011))
#> plot(t[t$V1=="4",2],t[t$V1=="4",3],type="l",xlab="position",ylab="pi",main="4",ylim=c(0,0.011))
#> dev.copy2eps(file="/Volumes/Macintosh HD/Users/robertkofler/tmp/sim-sim.eps")
}

    
    
=head1 NAME

perl Visualise-output.pl - A script to create a pdf displaying the extent of a population genetics measure along the selected chromosomes

=head1 SYNOPSIS

perl Visualise-output.pl --ylab tajimad --input input.file --output output.pdf --chromosomes "2 3R 3L 4"

=head1 OPTIONS


=over 4

=item B<--input>

The input file. Has to be created with the script Variance-sliding.pl; Mandatory parameter

=item B<--output>

The output file. Will be pdf. Mandatory parameter

=item B<--ylab>

The y-label used for describing the axis in the graph

=item B<--ymin>

the minimum value for the y-achsis; optional parameter

=item B<--ymax>

the maximum value for the y-achsis; optional parameter

=item B<--chromosomes>

A space separated list of chromosomes (contigs) for which the graphs should be created. Must be provided within double quotes; eg: "2L 2R X 3L 3R 4"; Mandatory parameter

=item B<--scale-equal>

Flag; Sets equal scalling of chromosomes. default=flag not set

=item B<--ps>

Flag; print a ps-file instead of a pdf

=back

=head1 Details

=head2 Input

An output file created with the script Variance-sliding.pl

 2L	1730000	557	0.726	0.005647899
 2L	1740000	599	0.777	0.005657701
 2L	1750000	650	0.767	0.006129462
 2L	1760000	617	0.703	0.006265200
 2L	1770000	599	0.672	0.006427032

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: number of SNPs in the sliding window; These SNPs have been used to calculate the value in col 5
 col 4: fraction of the window covered by a sufficient number of reads. Suficient means higher than min-coverage
 col 5: population genetics estimator (pi, theta, d)

=cut
