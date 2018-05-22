#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use SynNonSyn;


my $codontablefile;
my $outputFile;
my $transversion_penalty=6;
my $help=0;
my $test=0;

    
GetOptions(
        "codon-table=s"  => \$codontablefile,
        "output=s"       => \$outputFile,
        "transversion-penalty=i"  => \$transversion_penalty,
        "help"           => \$help,
        "test"           => \$test
    );

pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find input file",-verbose=>1) unless -e $codontablefile;
pod2usage(-msg=>"Could not find output file",-verbose=>1) unless $outputFile;

    my $paramfile=$outputFile.".params";
    open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
    print $pfh "Using codon-table\t$codontablefile\n";
    print $pfh "Using output\t$outputFile\n";
    print $pfh "Using transversion-penalyt\t$transversion_penalty\n";
    print $pfh "Using test\t$test\n";
    print $pfh "Using help\t$help\n";
    close $pfh;


my $codontable=load_codon_table($codontablefile);
my $cotrans=get_codon_translator($codontable,$transversion_penalty);
open my $ofh, ">", $outputFile or die "Could not open output file";

print $ofh "# this table contains the average non-synonymous length of a codon\n";
print $ofh "# this value is necessary to calculate the length normalised Tajima's pi (division by the length of a sequence)\n";
print $ofh "# this table has been automatically created using a transversion penalty of: $transversion_penalty\n";
print $ofh "#";
print $ofh "# Any table complying to the following format may be provided\n";
print $ofh "# codon: non-synonymous length\n";
foreach my $codon (keys(%$codontable))
{
    my $length_ns;
    my $synlength=$cotrans->($codon);
    print $ofh "$codon: $synlength\n";
}
close $ofh;


exit;


sub get_codon_translator
{
    my $codon_table=shift;
    my $trans_penalty=shift;
    
    my $score_matrix={
        AG=>$trans_penalty,
        GA=>$trans_penalty,
        CT=>$trans_penalty,
        TC=>$trans_penalty,
        GT=>1,TG=>1,AC=>1,CA=>1,AT=>1,TA=>1,GC=>1,CG=>1
    };
    
    return sub
    {
      my $codon=shift;
      my $sumscore=0;
      my $nonsynsum=0;
      my $ori_aa=$codon_table->{$codon};
      foreach my $i(0..2)
      {
        my $ori_base=substr($codon,$i,1);
        my $altbases=get_alternative_bases($ori_base);
        
        foreach my $alternat_base(@$altbases)
        {
            my $newtriplet=$codon;
            substr($newtriplet,$i,1,$alternat_base);
            my $new_aa=$codon_table->{$newtriplet};
            my $score=$score_matrix->{$ori_base.$alternat_base};
            $sumscore+=$score;
            $nonsynsum+=$score if $new_aa ne $ori_aa;
        }
        
      }
      
      my $ns_length=3*$nonsynsum/$sumscore;
      return $ns_length;
    };
}

sub get_alternative_bases
{
    my $base=shift;
    if($base eq "A")
    {
        return [qw/T C G/];
    }
    elsif($base eq "T")
    {
        return [qw/A C G/];
    }
    elsif($base eq "C")
    {
        return [qw/G A T/];
    }
    elsif($base eq "G")
    {
        return [qw/C A T/];
    }
    else
    {
        die "do not recognise base";
    }
}




