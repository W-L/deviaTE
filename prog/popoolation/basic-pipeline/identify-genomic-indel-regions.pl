use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
our $verbose=1;

my $input="";
my $output="";
my $indelwindow=5;
my $mincount=1;
my $help=0;
my $test=0;


#--input /Volumes/Volume_4/analysis/phix/pileup/m_q20_c100.pileup --output /Users/robertkofler/testfiles/output/indel.gtf

GetOptions(
    "input=s"           =>\$input,
    "indel-window=i"    =>\$indelwindow,
    "output=s"          =>\$output,
    "min-count=i"       =>\$mincount,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";

pod2usage(-verbose=>2) if $help;
IndelTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;
pod2usage(-msg=>"Indel window has to be larger or equal to one",-verbose=>1) unless $indelwindow>=1;


my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using indel-window\t$indelwindow\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

open my $ifh, "<",$input or die "Could not open pileup file";
open my $ofh,">",$output or die "Could not open output file";
my $counter=0;


my $indelposar=[];

my($count_entries,$count_indelpos,$count_bases_covered)=(0,0,0);

while(my $line=<$ifh>)
{
    
    chomp $line;
    # prefilter, without proper parsing of the pileup
    $count_entries++;    
    my($chr,$pos,$rc,$cov,$nucs,$qual)=split/\s+/,$line;
    
    
    next unless $nucs=~m/[-+]\d/;
    $count_indelpos++;
    my ($valid,$indelleng)=Utility::parse_pileup($nucs,$mincount);
    next unless $valid;
    
    push @$indelposar,{
            chr=>$chr,
            pos=>$pos,
            len=>$indelleng
    };
}

my $wincreator=Utility::get_region_calculator($indelposar,$indelwindow);

while(my $win=$wincreator->())
{
    print $ofh Utility::format_gtf($win);
    my $leng=$win->{end}-$win->{start}+1;
    $count_bases_covered+=$leng;
}
print "Pileup entries processed: $count_entries\n";
print "Pileup entries containing at least one indel: $count_indelpos\n";
print "How many bp of the reference are covered by indel-regions: $count_bases_covered\n";

close $ofh;


exit;

{
    package Utility;
    use strict;
    use warnings;
    
    sub parse_pileup
    {
        my $nucs=shift;
        my $mincount=shift;

        
        my $tempindels=[];
        # deletions have a length in the reference: the numbe after the '-'
        my(@dels)=$nucs=~m/[-](\d+)(??{"[ACGTNacgtn]{$1}"})/g;
        # gi|9626372|ref|NC_001422.1|     329     T       95      .$..$.$................................-1A...-1A......................................................
        # gi|9626372|ref|NC_001422.1|     330     A       94      .$...............................*..*........................................................^F.^F.
        
        foreach my $d (@dels)
        {
            push @$tempindels,$d;
        }
        
        # insertions have no length in the reference: length zero (they only have a length in the read)
        my(@ins)=$nucs=~m/[+](\d+)(??{"[ACGTNacgtn]{$1}"})/g;
        foreach my $i (@ins)
        {
            push @$tempindels,0;
        }
        
        my $tlh={}; # temporary length hash
        foreach my $t (@$tempindels)
        {
            $tlh->{$t}++
        }
        
        my $maxleng=0;
        my $valid=0;
        while(my($leng,$count)=each(%$tlh))
        {
            if($count>=$mincount)
            {
                $valid=1;
                $maxleng=$leng if $leng> $maxleng;
            }
        }
        
        return ($valid, $maxleng);
    }
    
    sub format_gtf
    {
        my $e=shift;
        # AB000381 Twinscan  exon         501   650   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
        # AB000381 Twinscan  CDS          501   650   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
        # AB000381 Twinscan  exon         700   800   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
        # AB000381 Twinscan  CDS          700   707   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
        my $str="$e->{chr}\tpileup\tindelregion\t$e->{start}\t$e->{end}\t.\t.\t.\tgene_id \".\"; transcript_id \".\";\n";
    }
    
    
    sub get_region_calculator
    {
        my $indels=shift;
        my $windowlength=shift;
        $indels=[sort {$a->{chr} cmp $b->{chr}
                          or $a->{pos} <=> $b->{pos}} @$indels];
    
        # A closure (can be used as a lightwight class!)           
        return sub
        {
          return undef unless @$indels;
          my $ini=shift @$indels;
          my $chr=$ini->{chr};
          my $start=$ini->{pos};
          
          #the end position of the current indelstretch
          my $activepos=$start+$ini->{len};
          
          while(@$indels)
          {
            my $next=shift @$indels;
            if($next->{chr} eq $chr && ($activepos+2*$windowlength)>=$next->{pos})
            {
                $activepos=$next->{pos}+$next->{len};
            }
            else
            {
                unshift @$indels,$next;
                last;
            } 
          }
          
          $start=$start-($windowlength-1);
          $start=1 if $start<1;
          $activepos+=$windowlength;
          return {
            chr=>$chr,
            start=>$start,
            end=>$activepos
          };
        };
    }

}


{
    package IndelTest;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/../Modules";
    use Test;
    
    sub runTests
    {
        test_parse_pileup();
        test_get_region_calculator();
        exit;
    }
    
    sub test_parse_pileup
    {

        my $best;
        my $valid;
        # - .. deletion in the read (they have a size in the  reference genome)
        # + .. insertion in the read (no size)
        ($valid,$best)=Utility::parse_pileup(".-1A.",1);
        is($best,1,"test parse pileup; correct size");
        is($valid,1,"test parse pileup; valid correct");
        
        ($valid,$best)=Utility::parse_pileup(".-1A...+2AA..",1);
        is($best,1,"test parse pileup; correct size",1);
        is($valid,1,"test parse pileup; valid correct");
        
        ($valid,$best)=Utility::parse_pileup(".+1A...+2AA..",1);
        is($best,0,"test parse pileup; correct size");
        is($valid,1,"test parse pileup; valid correct");
        
        ($valid,$best)=Utility::parse_pileup(".-1A...-2AA..",1);
        is($best,2,"test parse pileup; correct size");
        is($valid,1,"test parse pileup; valid correct");
        
        ($valid,$best)=Utility::parse_pileup("..-3ACG.-1A...-2AA..",1);
        is($best,3,"test parse pileup; correct size");
        is($valid,1,"test parse pileup; valid correct");
        
        ($valid,$best)=Utility::parse_pileup(".-1A...+2AA..",2);
        is($best,0,"test parse pileup; correct size");
        is($valid,0,"test parse pileup; valid correct");
        
        ($valid,$best)=Utility::parse_pileup(".+2AA...+2AA..",2);
        is($best,0,"test parse pileup; correct size");
        is($valid,1,"test parse pileup; valid correct");
        
        ($valid,$best)=Utility::parse_pileup(".-1A...-1A..",2);
        is($best,1,"test parse pileup; correct size");
        is($valid,1,"test parse pileup; valid correct");
    }
    
    sub test_get_region_calculator
    {
        my $indels;
        my $calc;
        my $r;
        
        $indels=[{chr=>2,pos=>20,len=>0}];
        $calc=Utility::get_region_calculator($indels,1);
        $r=$calc->();
        is($r->{chr},"2","region_calculator: chromosome correct");
        is($r->{start},20,"region_calculator: start correct");
        is($r->{end},21,"region_calculator: end correct");
        
        $indels=[{chr=>2,pos=>20,len=>1}];
        $calc=Utility::get_region_calculator($indels,1);
        $r=$calc->();
        is($r->{chr},"2","region_calculator: chromosome correct");
        is($r->{start},20,"region_calculator: start correct");
        is($r->{end},22,"region_calculator: end correct");
        
        $indels=[{chr=>2,pos=>20,len=>0}];
        $calc=Utility::get_region_calculator($indels,2);
        $r=$calc->();
        is($r->{chr},"2","region_calculator: chromosome correct");
        is($r->{start},19,"region_calculator: start correct");
        is($r->{end},22,"region_calculator: end correct");
        
        $indels=[{chr=>2,pos=>20,len=>2}];
        $calc=Utility::get_region_calculator($indels,3);
        $r=$calc->();
        is($r->{chr},"2","region_calculator: chromosome correct");
        is($r->{start},18,"region_calculator: start correct");
        is($r->{end},25,"region_calculator: end correct");
        $r=$calc->();
        not_exists($r,"region_calculator; correct no more entries");
        
        $indels=[{chr=>2,pos=>20,len=>2},{chr=>2,pos=>28,len=>2}];
        $calc=Utility::get_region_calculator($indels,3);
        $r=$calc->();
        is($r->{chr},"2","region_calculator: chromosome correct");
        is($r->{start},18,"region_calculator: start correct");
        is($r->{end},33,"region_calculator: end correct");
        
        $indels=[{chr=>2,pos=>20,len=>2},{chr=>2,pos=>29,len=>2}];
        $calc=Utility::get_region_calculator($indels,3);
        $r=$calc->();
        is($r->{chr},"2","region_calculator: chromosome correct");
        is($r->{start},18,"region_calculator: start correct");
        is($r->{end},25,"region_calculator: end correct");
        $r=$calc->();
        is($r->{chr},"2","region_calculator: chromosome correct");
        is($r->{start},27,"region_calculator: start correct");
        is($r->{end},34,"region_calculator: end correct");
        $r=$calc->();
        not_exists($r,"region_calculator; correct no more entries");

        
    }
}



=head1 NAME

perl identify-genomic-indel-regions.pl - A script which identifies the surrounding indels and converts these coordindates into a gtf

=head1 SYNOPSIS

perl identify-genomic-indel-regions.pl --input input.pileup --output indel-regions.gtf --indel-window 5

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the pileup format; Mandatory.

=item B<--output>

The output file.  Mandatory.

=item B<--indel-window>

Length of the window around the indel; default=5

=item B<--min-count>

Minimum count for an indel; default=1

=item B<--test>

Run the tests for the script

=item B<--help>

Display help for this script

=back

=head1 Details

This script identifies indels from a pileup file; It furthermore converts the coordinates of the positions around the indels into a gtf-formated file;
This file can for example be used to filter a pileup file (filter-pileup-by-gtf.pl)

=head2 Input pileup

A pileup file as described here: http://samtools.sourceforge.net/pileup.shtml; example:

 2L	90131	N	11	AaAAAaaAaAA	[aUQ_a`^_\Z
 2L	90132	N	11	AaAAAaaAaAA	_bYQ_^aaT^b
 2L	90133	N	11	A$aAAAaaAaAA	_b[Xaaa__Ua
 2L	90134	N	10	tTTTttTtTT	_`aaa_a[aa
 2L	90135	N	10	a$TAAaa-1TAaAA	aZ^a`ba`\_
 2L	90136	N	9	TTTttTtTT	`aaaaaWaa
 2L	90137	N	9	GGGgg+2AAGgGG	``aaaaQaa
 2L	90138	N	9	T$TTttTtTT	[U\`a\T^_
 2L	90139	N	8	TTttTtTT	``aaaU_a
 2L	90140	N	9	CCccCcCC^FC	[aaba`aaa


=head2 Output


=cut



