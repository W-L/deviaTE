{
    package Test::TPileupTripletSliding;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use SynNonSyn;
    use Test;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(run_PileupTripletSlidingTests);
    use SynNonSyn;
    use PileupTripletSlider;
    use Pileup;
    
    
    
    sub run_PileupTripletSlidingTests
    {
        test_PileupTripletSlider();
    }
    
    
    sub test_PileupTripletSlider
    {
        
        my $gtfstr;
        my $anno;
        # fastqtype, mincount, mincoverage, maxcoverage, minquality
        my $pp=get_extended_parser("illumina",2,4,1000000,0);
        my $pts;
        my $pustr;
        my $e;
        
        $gtfstr=
        "2R\tmaker\tCDS\t3\t5\t.\t+\t0\n".
        "2R\tmaker\tCDS\t8\t10\t.\t-\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pustr=
        "2R\t1\tN\t9\tAAAA\taaaa\n".
        "2R\t2\tN\t9\tTTTT\taaaa\n".
        "2R\t3\tN\t9\tTTTT\taaaa\n".
        "2R\t4\tN\t9\tCCCC\taaaa\n".
        "2R\t5\tN\t9\tCAAAA\taaaaa\n".
        "2R\t6\tN\t9\tAAAA\taaaa\n".
        "2R\t7\tN\t9\tAAAA\taaaa\n".
        "2R\t8\tN\t9\tAAAA\taaaa\n".
        "2R\t9\tN\t9\tTTTT\taaaa\n".
        "2R\t10\tN\t9\tGGGG\taaaa\n".
        "2R\t11\tN\t9\tAAAA\taaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"3","PileupTriplet slider; start ok");
        is($e->{strand},"+","PileupTriplet slider; strand ok");
        is($e->{codon},"TCA","PileupTriplet slider; codon ok");
        is($e->{valid},"1","PileupTriplet slider; Validity is ok");
        is($e->{valid_codon},"1","PileupTriplet slider; Validity is ok");
        is($e->{valid_coverage},"1","PileupTriplet slider; Validity is ok");
        is($e->{valid_consensus},"1","PileupTriplet slider; Validity is ok");
        is($e->{frame}[0],"1","PileupTriplet slider; Frame is ok");
        is($e->{frame}[1],"2","PileupTriplet slider; Frame is ok");
        is($e->{frame}[2],"3","PileupTriplet slider; Frame is ok");
        is($e->{pileup}[0]{pos},"3","PileupTriplet slider; Position is ok");
        is($e->{pileup}[1]{pos},"4","PileupTriplet slider; Position is ok");
        is($e->{pileup}[2]{pos},"5","PileupTriplet slider; Position is ok");
        is($e->{count_snps},"0","PileupTriplet slider; SNP count is ok");
        
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"8","PileupTriplet slider; start ok");
        is($e->{strand},"-","PileupTriplet slider; strand ok");
        is($e->{codon},"CAT","PileupTriplet slider; codon ok");
        is($e->{valid},"1","PileupTriplet slider; Validity is ok");
        is($e->{valid_codon},"1","PileupTriplet slider; Validity is ok");
        is($e->{valid_coverage},"1","PileupTriplet slider; Validity is ok");
        is($e->{valid_consensus},"1","PileupTriplet slider; Validity is ok");
        is($e->{frame}[0],"6","PileupTriplet slider; Frame is ok");
        is($e->{frame}[1],"5","PileupTriplet slider; Frame is ok");
        is($e->{frame}[2],"4","PileupTriplet slider; Frame is ok");
        is($e->{pileup}[0]{pos},"8","PileupTriplet slider; Position is ok");
        is($e->{pileup}[1]{pos},"9","PileupTriplet slider; Position is ok");
        is($e->{pileup}[2]{pos},"10","PileupTriplet slider; Position is ok");
        is($e->{count_snps},"0","PileupTriplet slider; SNP count is ok");
        
        $e=$pts->next();
        not_exists($e, "PileupTriplet slider; correct no more codons");
        
        # longer codon
        $gtfstr=
        "2R\tmaker\tCDS\t3\t8\t.\t+\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"3","PileupTriplet slider; start ok");
        is($e->{codon},"TCA","PileupTriplet slider; codon ok");
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"6","PileupTriplet slider; start ok");
        is($e->{codon},"AAA","PileupTriplet slider; codon ok");
        $e=$pts->next();
        not_exists($e, "PileupTriplet slider; correct no more codons");
        
        
        # invalid codon in between
        $gtfstr=
        "2R\tmaker\tCDS\t3\t5\t.\t+\t0\n".
        "2R\tmaker\tCDS\t3\t5\t.\t-\t0\n".
        "2R\tmaker\tCDS\t6\t8\t.\t+\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"6","PileupTriplet slider; start ok");
        is($e->{codon},"AAA","PileupTriplet slider; codon ok");
        $e=$pts->next();
        not_exists($e, "PileupTriplet slider; correct no more codons");
        
        $gtfstr=
        "2R\tmaker\tCDS\t3\t11\t.\t+\t0\n".
        "2R\tmaker\tCDS\t5\t7\t.\t-\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"9","PileupTriplet slider; start ok");
        is($e->{codon},"TGA","PileupTriplet slider; codon ok");
        $e=$pts->next();
        not_exists($e, "PileupTriplet slider; correct no more codons");


        # different chromosomes
        $gtfstr=
        "2R\tmaker\tCDS\t3\t5\t.\t+\t0\n".
        "3L\tmaker\tCDS\t2\t4\t.\t-\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pustr=
        "2R\t3\tN\t9\tTTTT\taaaa\n".
        "2R\t4\tN\t9\tCCCC\taaaa\n".
        "2R\t5\tN\t9\tCAAAA\taaaaa\n".
        "2R\t6\tN\t9\tAAAA\taaaa\n".
        "3L\t1\tN\t9\tAAAA\taaaa\n".
        "3L\t2\tN\t9\tAAAA\taaaa\n".
        "3L\t3\tN\t9\tTTTT\taaaa\n".
        "3L\t4\tN\t9\tGGGG\taaaa\n".
        "3L\t5\tN\t9\tAAAA\taaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"3","PileupTriplet slider; start ok");
        is($e->{codon},"TCA","PileupTriplet slider; codon ok");
        $e=$pts->next();
        is($e->{chr},"3L","PileupTriplet slider; chromosome ok");
        is($e->{start},"2","PileupTriplet slider; start ok");
        is($e->{codon},"CAT","PileupTriplet slider; codon ok");
        $e=$pts->next();
        not_exists($e, "PileupTriplet slider; correct no more codons");
        
        
        # invalid triplets
        $gtfstr=
        "2R\tmaker\tCDS\t3\t5\t.\t+\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pustr=
        "2R\t3\tN\t9\tTTTAA\taaaaa\n".
        "2R\t4\tN\t9\tCCA\taaa\n".
        "2R\t5\tN\t9\tCCCC\taaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"3","PileupTriplet slider; start ok");
        is($e->{strand},"+","PileupTriplet slider; strand ok");
        is($e->{codon},"TNC","PileupTriplet slider; codon ok");
        is($e->{valid},"0","PileupTriplet slider; Validity is ok");
        is($e->{valid_codon},"0","PileupTriplet slider; Validity is ok");
        is($e->{valid_coverage},"0","PileupTriplet slider; Validity is ok");
        is($e->{valid_consensus},"0","PileupTriplet slider; Validity is ok");
        is($e->{count_snps},"1","PileupTriplet slider; SNP count is ok");
        
        $gtfstr=
        "2R\tmaker\tCDS\t3\t5\t.\t-\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $e=$pts->next();
        is($e->{chr},"2R","PileupTriplet slider; chromosome ok");
        is($e->{start},"3","PileupTriplet slider; start ok");
        is($e->{strand},"-","PileupTriplet slider; strand ok");
        is($e->{codon},"GNA","PileupTriplet slider; codon ok");
        is($e->{valid},"0","PileupTriplet slider; Validity is ok");
        is($e->{valid_codon},"0","PileupTriplet slider; Validity is ok");
        is($e->{valid_coverage},"0","PileupTriplet slider; Validity is ok");
        is($e->{valid_consensus},"0","PileupTriplet slider; Validity is ok");
        is($e->{count_snps},"1","PileupTriplet slider; SNP count is ok");
        
    }
    
    
    
    
}

1;