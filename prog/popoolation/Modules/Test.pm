{
package Test;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT     =qw(is not_exists ok);

my $testcounter;

sub is
    {

        my $a=shift;
        my $b=shift;
        my $msg=shift;
        $testcounter=1 unless $testcounter;
        
        
        if($a eq $b)
        {
            print "$testcounter - OK: $a = $b; $msg\n";
        }
        else
        {
            print "$testcounter - FAILED: $a = $b; $msg\n";
        }
        $testcounter++;
    }
    
    sub not_exists
    {
        my $a=shift;
        my $msg=shift;
        $testcounter=1 unless $testcounter;
        
        
        if($a)
        {
            print "$testcounter - FAILED: $msg\n";
        }
        else
        {
            print "$testcounter - OK: $msg\n";
        }
        $testcounter++;
    }

    sub ok
    {
        my $a=shift;
        my $msg=shift;
        $testcounter=1 unless $testcounter;
        
        if($a)
        {
            print "$testcounter - OK: $msg\n";
        }
        else
        {
            print "$testcounter - FAILED: $msg\n";
        }
        $testcounter++;
    }
    



}

1;