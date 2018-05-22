{
    package Test::TClassicalVariance;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use VarianceUncorrected;
    use VarMathClassical;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT=qw(run_classicalVarianceTests);
    our @EXPORT_OK = qw();
    
    
    sub run_classicalVarianceTests
    {
        test_classical_Pi();  
        test_classical_Theta();
        test_classical_D();
    }
    

    sub test_classical_Pi
    {
        my $vc;
        my $pi;
        

        $vc=get_classical_pi_calculator();
        $pi=$vc->([{eucov=>10,A=>5,T=>5,C=>0,G=>0}]);
        ok(abs($pi-0.5555555)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($pi-0.200000)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($pi-0.400000)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($pi-0.600000)<0.00001,"Testing uncorrected Pi: Pi value correct");

        $pi=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($pi-0.755556)<0.00001,"Testing uncorrected Pi: Pi value correct");
        
        $pi=$vc->([{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($pi-0.355556)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($pi-0.555556)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>10,A=>7,T=>3,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($pi-1.022222)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>10,A=>8,T=>2,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($pi-0.711111)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>10,A=>7,T=>3,C=>0,G=>0},{eucov=>10,A=>7,T=>3,C=>0,G=>0},{eucov=>10,A=>7,T=>3,C=>0,G=>0}]);
        ok(abs($pi-1.400000)<0.00001,"Testing uncorrected Pi: Pi value correct");
        
        $pi=$vc->([{eucov=>20,A=>16,T=>4,C=>0,G=>0},{eucov=>20,A=>16,T=>4,C=>0,G=>0}]);
        ok(abs($pi-0.673684)<0.00001,"Testing uncorrected Pi: Pi value correct");
        $pi=$vc->([{eucov=>30,A=>27,T=>3,C=>0,G=>0},{eucov=>30,A=>24,T=>6,C=>0,G=>0},{eucov=>30,A=>21,T=>9,C=>0,G=>0}]);
        ok(abs($pi-0.951724)<0.00001,"Testing uncorrected Pi: Pi value correct");
        
        $pi=$vc->([]);
        ok(abs($pi-0)<0.00001,"Testing uncorrected Pi: Pi value correct");
    
    }
    
    sub test_classical_Theta
    {
        my $vc;
        my $theta;
        
        $vc=get_classical_theta_calculator();
        $theta=$vc->([{eucov=>10,A=>5,T=>5,C=>0,G=>0}]);
        ok(abs($theta-0.353486)<0.00001,"Testing uncorrected Theta: Theta value correct");
        $theta=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($theta-0.353486)<0.00001,"Testing uncorrected Theta: Theta value correct");
        $theta=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($theta-0.706971)<0.00001,"Testing uncorrected Theta: Theta value correct");
        $theta=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($theta-1.060457)<0.00001,"Testing uncorrected Theta: Theta value correct");
        
        $theta=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($theta-1.060457)<0.00001,"Testing uncorrected Theta: Theta value correct");
        
        $theta=$vc->([{eucov=>30,A=>27,T=>3,C=>0,G=>0},{eucov=>30,A=>24,T=>6,C=>0,G=>0},{eucov=>30,A=>21,T=>9,C=>0,G=>0}]);
        ok(abs($theta-0.757259)<0.00001,"Testing uncorrected Theta: Theta value correct");
        
        $theta=$vc->([]);
        ok(abs($theta-0)<0.00001,"Testing uncorrected Theta: Theta value correct");
    }

    
    sub test_classical_D
    {
        my $vc;
        my $d;
        
        $vc=get_classical_D_calculator();
        $d=$vc->([{eucov=>10,A=>5,T=>5,C=>0,G=>0}]);
        ok(abs($d-1.463639)<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($d-(-1.111733))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($d-(-1.400851 ))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($d-(-1.562218))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        
        
        $d=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($d-(-1.034456))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        
        $d=$vc->([{eucov=>10,A=>2,T=>8,C=>0,G=>0}]);
        ok(abs($d-0.014992)<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>9,T=>1,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($d-(-0.690980 ))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>7,T=>3,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($d-(-0.129722))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>8,T=>2,C=>0,G=>0},{eucov=>10,A=>8,T=>2,C=>0,G=>0}]);
        ok(abs($d-(0.018891))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>7,T=>3,C=>0,G=>0},{eucov=>10,A=>7,T=>3,C=>0,G=>0},{eucov=>10,A=>7,T=>3,C=>0,G=>0}]);
        ok(abs($d-(1.151984))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>6,T=>4,C=>0,G=>0},{eucov=>10,A=>6,T=>4,C=>0,G=>0},{eucov=>10,A=>6,T=>4,C=>0,G=>0}]);
        ok(abs($d-(1.830535))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>5,T=>5,C=>0,G=>0},{eucov=>10,A=>5,T=>5,C=>0,G=>0},{eucov=>10,A=>5,T=>5,C=>0,G=>0}]);
        ok(abs($d-(2.056718))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>10,A=>5,T=>5,C=>0,G=>0},{eucov=>10,A=>9,T=>1,C=>0,G=>0}]);
        ok(abs($d-(0.221711))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        
        $d=$vc->([{eucov=>20,A=>16,T=>4,C=>0,G=>0},{eucov=>20,A=>16,T=>4,C=>0,G=>0}]);
        ok(abs($d-(0.457275))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        
        $d=$vc->([{eucov=>30,A=>27,T=>3,C=>0,G=>0},{eucov=>30,A=>24,T=>6,C=>0,G=>0},{eucov=>30,A=>21,T=>9,C=>0,G=>0}]);
        ok(abs($d-(0.604333))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        
        
        # Looking for the lower boundaries
        $d=$vc->([{eucov=>30,A=>29,T=>1,C=>0,G=>0},{eucov=>30,A=>29,T=>1,C=>0,G=>0},{eucov=>30,A=>29,T=>1,C=>0,G=>0}]);
        ok(abs($d-(-1.731783))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        $d=$vc->([{eucov=>30,A=>29,T=>1,C=>0,G=>0},{eucov=>30,A=>29,T=>1,C=>0,G=>0},{eucov=>30,A=>29,T=>1,C=>0,G=>0},{eucov=>30,A=>29,T=>1,C=>0,G=>0},{eucov=>30,A=>29,T=>1,C=>0,G=>0},{eucov=>30,A=>29,T=>1,C=>0,G=>0}]);
        ok(abs($d-(-2.099947))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
        
        
        $d=$vc->([]);
        ok(abs($d-(0))<0.00001,"Testing uncorrected Tajima's D: Tajima's D value correct");
    }
    
    

    
    


}


1;
