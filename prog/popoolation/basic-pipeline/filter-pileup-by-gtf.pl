use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use Pileup;

our $verbose=1;

my $input="";
my $gtffile="";
my $output="";
my $help=0;
my $test=0;
my $discardmode=1;

#--gtf /Users/robertkofler/testfiles/3R-3entries.gtf --input /Users/robertkofler/testfiles/3R.pileup --output /Users/robertkofler/testfiles/output/3R-filtered.pileup

GetOptions(
    "input=s"           =>\$input,
    "gtf=s"             =>\$gtffile,
    "output=s"          =>\$output,
    "keep-mode"         =>sub{$discardmode=0;},
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";

pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Could not find gtf file",-verbose=>1) unless -e $gtffile;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;


my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using gtf\t$gtffile\n";
print $pfh "Using discardmode (keep-mode)\t$discardmode\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my ($chrhash)=Utility::read_gtf($gtffile);

open my $ifh, "<",$input or die "Could not open pileup file";

print "Start parsing the pileup file..\n";
open my $ofh,">",$output or die "Could not open output file";
my $counter=0;

while(my $line=<$ifh>)
{
    # parse the whole pileup and store every parsed pileup-entry at the corresponding
    chomp $line;

    # prefilter, without proper parsing of the pileup
    my($chr,$pos)=split/\t/,$line;
    
    my $ispresent=$chrhash->{$chr}{$pos};
    
    if($discardmode)
    {
        # discard everything in the gtf
        unless($ispresent)
        {
            # if not present in the gtf
            print $ofh "$line\n";
        }
    }
    else
    {
        # keep mode
        if($ispresent)
        {
            # if you only want to keep the things in the gtf, print only when present
            print $ofh "$line\n";
        }
    }
}


close $ofh;

print "Done\n";
exit;

{
    package Utility;
    use strict;
    use warnings;
    
    
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
        
        return
        {
            ref=>$ref,
            start=>$start,
            end=>$end,
            length=>$end-$start+1
        };
    }
    
    
    sub read_gtf
    {
        my $file=shift;
        open my $ifh,"<",$file or die "Could not open gtf-file";
        my $chrhash={};
        
        print "Parsing gtf file..\n";
        while(my $line=<$ifh>)
        {
            chomp $line;
            next if $line=~m/^##/;
            my $ge=_parsegtf($line);
        
            # update the chromosome hash
            my($chr,$start,$end)=($ge->{ref},$ge->{start},$ge->{end});
            
            
            for(my $i=$start; $i<=$end; $i++)
            {
                $chrhash->{$chr}{$i}++;
            }
        }
        return $chrhash;
    }
}



=head1 NAME

perl filter-pileup-by-gtf.pl - A script which filters a pileup either discarding or keeping only the regions specified in a gtf-file

=head1 SYNOPSIS

perl filter-pileup-by-gtf.pl --gtf annotation.gtf --input input.pileup --output filtered.pileup

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the pileup format; Mandatory.

=item B<--gtf>

A gtf file as specified here: http://mblab.wustl.edu/GTF2.html; The regions specified in the gtf-file will be used for filtering of the pileup file. Mandatory

For more information use <--help>

=item B<--output>

The output file.  Mandatory.

=item B<--keep-mode>

Flag; Retain the positions provided in the gtf file instead of discarding them. Per default it is discarding the regions provided in the gtf file

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

Everything which has an entry in the gtf-file will be used for filtering the pileup file, irrespective of the source, feature, strand, gene_id etc
Therefore only the columns reference-id, start-position and end-position (1, 4, 5) will be relevant for filtering the pileup file.

=head2 Output

The output will be a filtered pileup file; Dependent on the mode (discard or keep) all regions specified in the gtf-file will be either absent or present.
Per default the discard mode is activated, i.: all regions specified in the gtf-file will be absent in the filtered pileup file.

=head2 Technical details

The memory usage scales linearly with the number of bases covered by gtf-entries. If you get an out of memory exception (eg. malloc) split your gtf file and run the script several times.
Alternatively if your gtf-file, for example, covers most of the reference genome, you could invert the gtf-file and switch the mode of this script (eg.: from discard to keep)!
  
=cut



