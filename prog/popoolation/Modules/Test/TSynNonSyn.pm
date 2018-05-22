{
    package Test::TSynNonSyn;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use SynNonSyn;
    use Pileup;
    use PileupTripletSlider;
    use Test;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(run_SynNonSynTests);
    
    
    
    
    sub run_SynNonSynTests
    {
        test_load_cds();
        test_get_codon_changes();
    }
    
    sub test_get_codon_changes
    {
        my $gtfstr;
        my $anno;
        # fastqtype, mincount, mincoverage, maxcoverage, minquality
        my $pp=get_extended_parser("illumina",2,4,1000000,0);
        my $pts;
        my $pustr;
        my $tr;
        
        my $codonTable=_getcodontable();
        
        $gtfstr=
        "2R\tmaker\tCDS\t3\t5\t.\t+\t0\n".
        "2R\tmaker\tCDS\t8\t10\t.\t-\t0\n";
        $anno=load_cds_gtf(\$gtfstr);
        $pustr=
        "2R\t1\tN\t9\tAAAA\taaaa\n".
        "2R\t2\tN\t9\tTTTT\taaaa\n".
        "2R\t3\tN\t9\tTTTAA\taaaaa\n".
        "2R\t4\tN\t9\tCCCC\taaaa\n".
        "2R\t5\tN\t9\tCAAAA\taaaaa\n".
        "2R\t6\tN\t9\tAAAA\taaaa\n".
        "2R\t7\tN\t9\tAAAA\taaaa\n".
        "2R\t8\tN\t9\tAAAA\taaaa\n".
        "2R\t9\tN\t9\tCCCGG\taaaaa\n".
        "2R\t10\tN\t9\tGGGG\taaaa\n".
        "2R\t11\tN\t9\tAAAA\taaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $tr=$pts->next();
        get_codon_changes($tr,$codonTable);
        is($tr->{codon},"TCA","calculating amino acid change; original codon correct");
        is($tr->{chr},"2R","calculating amino acid change, chromosome is correct");
        is($tr->{count_snps},1,"calculating amino acid change, snp count is correct");
        is($tr->{cc}[0]{codon},"TCA","calculating amino acid change; original codon correct" );
        is($tr->{cc}[0]{mcodon},"ACA","calculating amino acid change; modified codon correct");
        is($tr->{cc}[0]{aa},"S","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[0]{maa},"T","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[0]{A},2,"calculating amino acid change; A count ok");
        is($tr->{cc}[0]{T},3,"calculating amino acid change; T count ok");
        is($tr->{cc}[0]{C},0,"calculating amino acid change; C count ok");
        is($tr->{cc}[0]{G},0,"calculating amino acid change; G count ok");
        is($tr->{cc}[0]{strand},"+","calculating amino acid change; strand ok");
        is($tr->{cc}[0]{syn},0,"calculating amino acid change; synonymous determination ok");
        is($tr->{cc}[0]{cstart},3,"calculating amino acid change; codon start is ok");

        $tr=$pts->next();
        get_codon_changes($tr,$codonTable);
        is($tr->{codon},"CGT","incorporating snp into codon; original codon correct");
        is($tr->{chr},"2R","calculating amino acid change, chromosome is correct");
        is($tr->{count_snps},1,"calculating amino acid change, snp count is correct");
        is($tr->{cc}[0]{codon},"CGT","calculating amino acid change; original codon correct" );
        is($tr->{cc}[0]{mcodon},"CCT","calculating amino acid change; modified codon correct");
        is($tr->{cc}[0]{aa},"R","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[0]{maa},"P","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[0]{strand},"-","calculating amino acid change; strand ok");
        is($tr->{cc}[0]{syn},0,"calculating amino acid change; synonymous determination ok");
        
        # syn change
        $pustr=
        "2R\t1\tN\t9\tAAAA\taaaa\n".
        "2R\t2\tN\t9\tTTTT\taaaa\n".
        "2R\t3\tN\t9\tAAAAA\taaaaa\n".
        "2R\t4\tN\t9\tCCCC\taaaa\n".
        "2R\t5\tN\t9\tAAAAGG\taaaaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $tr=$pts->next();
        get_codon_changes($tr,$codonTable);
        is($tr->{codon},"ACA","incorporating snp into codon; original codon correct");
        is($tr->{chr},"2R","calculating amino acid change, chromosome is correct");
        is($tr->{count_snps},1,"calculating amino acid change, snp count is correct");
        is($tr->{cc}[0]{codon},"ACA","calculating amino acid change; original codon correct" );
        is($tr->{cc}[0]{mcodon},"ACG","calculating amino acid change; modified codon correct");
        is($tr->{cc}[0]{aa},"T","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[0]{maa},"T","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[0]{strand},"+","calculating amino acid change; strand ok");
        is($tr->{cc}[0]{syn},1,"calculating amino acid change; synonymous determination ok");
        
        
        ##
        ## double SNP
        ##
        $pustr=
        "2R\t1\tN\t9\tAAAA\taaaa\n".
        "2R\t2\tN\t9\tTTTT\taaaa\n".
        "2R\t3\tN\t9\tAAAAA\taaaaa\n".
        "2R\t4\tN\t9\tCCCCAAA\taaaaaaa\n".
        "2R\t5\tN\t9\tAAAAGG\taaaaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $tr=$pts->next();
        get_codon_changes($tr,$codonTable);
        is($tr->{codon},"ACA","incorporating snp into codon; original codon correct");
        is($tr->{chr},"2R","calculating amino acid change, chromosome is correct");
        is($tr->{count_snps},2,"calculating amino acid change, snp count is correct");
        is($tr->{cc}[0]{codon},"ACA","calculating amino acid change; original codon correct" );
        is($tr->{cc}[0]{mcodon},"AAA","calculating amino acid change; modified codon correct");
        is($tr->{cc}[0]{aa},"T","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[0]{maa},"K","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[0]{strand},"+","calculating amino acid change; strand ok");
        is($tr->{cc}[0]{syn},0,"calculating amino acid change; synonymous determination ok");
        is($tr->{cc}[1]{codon},"ACA","calculating amino acid change; original codon correct" );
        is($tr->{cc}[1]{mcodon},"ACG","calculating amino acid change; modified codon correct");
        is($tr->{cc}[1]{aa},"T","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[1]{maa},"T","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[1]{strand},"+","calculating amino acid change; strand ok");
        is($tr->{cc}[1]{syn},1,"calculating amino acid change; synonymous determination ok");
        
        ##
        ## triple snp
        ## 
        $pustr=
        "2R\t1\tN\t9\tAAAA\taaaa\n".
        "2R\t2\tN\t9\tTTTT\taaaa\n".
        "2R\t3\tN\t9\tAAAAATT\taaaaaaa\n".
        "2R\t4\tN\t9\tCCCCAAA\taaaaaaa\n".
        "2R\t5\tN\t9\tAAAAGG\taaaaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $tr=$pts->next();
        get_codon_changes($tr,$codonTable);
        is($tr->{codon},"ACA","incorporating snp into codon; original codon correct");
        is($tr->{chr},"2R","calculating amino acid change, chromosome is correct");
        is($tr->{count_snps},3,"calculating amino acid change, snp count is correct");
        is($tr->{cc}[0]{codon},"ACA","calculating amino acid change; original codon correct" );
        is($tr->{cc}[0]{mcodon},"TCA","calculating amino acid change; modified codon correct");
        is($tr->{cc}[0]{aa},"T","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[0]{maa},"S","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[0]{strand},"+","calculating amino acid change; strand ok");
        is($tr->{cc}[0]{syn},0,"calculating amino acid change; synonymous determination ok");
        is($tr->{cc}[1]{codon},"ACA","calculating amino acid change; original codon correct" );
        is($tr->{cc}[1]{mcodon},"AAA","calculating amino acid change; modified codon correct");
        is($tr->{cc}[1]{aa},"T","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[1]{maa},"K","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[1]{strand},"+","calculating amino acid change; strand ok");
        is($tr->{cc}[1]{syn},0,"calculating amino acid change; synonymous determination ok");
        is($tr->{cc}[2]{codon},"ACA","calculating amino acid change; original codon correct" );
        is($tr->{cc}[2]{mcodon},"ACG","calculating amino acid change; modified codon correct");
        is($tr->{cc}[2]{aa},"T","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[2]{maa},"T","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[2]{strand},"+","calculating amino acid change; strand ok");
        is($tr->{cc}[2]{syn},1,"calculating amino acid change; synonymous determination ok");
        
        
        
        ##
        ## triple SNP minus strand
        ##
        $pustr=
        "2R\t8\tN\t9\tAAAACC\taaaaaa\n".
        "2R\t9\tN\t9\tCCCGG\taaaaa\n".
        "2R\t10\tN\t9\tGGGGCC\taaaaaa\n";
        $pts=PileupTripletSlider->new(\$pustr,$anno,$pp);
        $tr=$pts->next();
        get_codon_changes($tr,$codonTable);
        is($tr->{codon},"CGT","incorporating snp into codon; original codon correct");
        is($tr->{chr},"2R","calculating amino acid change, chromosome is correct");
        is($tr->{count_snps},3,"calculating amino acid change, snp count is correct");
        is($tr->{cc}[0]{codon},"CGT","calculating amino acid change; original codon correct" );
        is($tr->{cc}[0]{mcodon},"CGG","calculating amino acid change; modified codon correct");
        is($tr->{cc}[0]{aa},"R","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[0]{maa},"R","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[0]{strand},"-","calculating amino acid change; strand ok");
        is($tr->{cc}[0]{syn},1,"calculating amino acid change; synonymous determination ok");
        is($tr->{cc}[1]{codon},"CGT","calculating amino acid change; original codon correct" );
        is($tr->{cc}[1]{mcodon},"CCT","calculating amino acid change; modified codon correct");
        is($tr->{cc}[1]{aa},"R","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[1]{maa},"P","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[1]{strand},"-","calculating amino acid change; strand ok");
        is($tr->{cc}[1]{syn},0,"calculating amino acid change; synonymous determination ok");
        is($tr->{cc}[2]{codon},"CGT","calculating amino acid change; original codon correct" );
        is($tr->{cc}[2]{mcodon},"GGT","calculating amino acid change; modified codon correct");
        is($tr->{cc}[2]{aa},"R","calculating amino acid change; original amino acid correct" );
        is($tr->{cc}[2]{maa},"G","calculating amino acid change; modified amino acid correct");
        is($tr->{cc}[2]{strand},"-","calculating amino acid change; strand ok");
        is($tr->{cc}[2]{syn},0,"calculating amino acid change; synonymous determination ok");
    }
    
    
    sub test_load_cds
    {
        
        # illuminator:
        # 1.. fwd strand, first base of frame
        # 2.. fwd strand, second base of frame
        # 3.. fwd strand, third
        # 4.. rev strand, first base of frame
        # 5.. rev strand, second
        # 6.. rev strand, third
        # 7.. invalid frame, due to overlapping frames
        my $str;
        my $res;
        
        $str=
        "2R\tmaker\tCDS\t3\t5\t.\t+\t0\n".
        "2L\tmaker\tCDS\t3\t5\t.\t-\t0\n".
        "3R\tmaker\tCDS\t3\t6\t.\t+\t1\n".
        "3L\tmaker\tCDS\t3\t6\t.\t-\t1\n".
        "2R\tmaker\tCDS\t10\t14\t.\t+\t1\n".
        "2L\tmaker\tCDS\t10\t14\t.\t-\t1\n";
        $res=load_cds_gtf(\$str);
        is($res->{'2R'}{'3'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'4'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'5'},3,"CDS loader; frame correct");
        is($res->{'2L'}{'5'},4,"CDS loader; frame correct");
        is($res->{'2L'}{'4'},5,"CDS loader; frame correct");
        is($res->{'2L'}{'3'},6,"CDS loader; frame correct");
        
        not_exists($res->{'3R'}{'3'},"CDS loader correct, no frame at position");
        is($res->{'3R'}{'4'},1,"CDS loader; frame correct");
        is($res->{'3R'}{'5'},2,"CDS loader; frame correct");
        is($res->{'3R'}{'6'},3,"CDS loader; frame correct");
        not_exists($res->{'3R'}{'7'},"CDS loader correct, no frame at position");
        not_exists($res->{'3L'}{'6'},"CDS loader correct, no frame at position");
        is($res->{'3L'}{'5'},4,"CDS loader; frame correct");
        is($res->{'3L'}{'4'},5,"CDS loader; frame correct");
        is($res->{'3L'}{'3'},6,"CDS loader; frame correct");
        not_exists($res->{'3L'}{'2'},"CDS loader correct, no frame at position");
        
        not_exists($res->{'2R'}{'10'},"CDS loader correct, no frame at position");
        is($res->{'2R'}{'11'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'12'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'13'},3,"CDS loader; frame correct");
        not_exists($res->{'2R'}{'14'},"CDS loader correct, no frame at position");
        not_exists($res->{'2L'}{'14'},"CDS loader correct, no frame at position");
        is($res->{'2L'}{'13'},4,"CDS loader; frame correct");
        is($res->{'2L'}{'12'},5,"CDS loader; frame correct");
        is($res->{'2L'}{'11'},6,"CDS loader; frame correct");
        not_exists($res->{'2L'}{'10'},"CDS loader correct, no frame at position");
        
        #longer frames
        $str=
        "2R\tmaker\tCDS\t3\t8\t.\t+\t0\n".
        "2L\tmaker\tCDS\t3\t8\t.\t-\t0\n";
        $res=load_cds_gtf(\$str);
        not_exists($res->{'2R'}{'2'},"CDS loader correct, no frame at position");
        is($res->{'2R'}{'3'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'4'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'5'},3,"CDS loader; frame correct");
        is($res->{'2R'}{'6'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'7'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'8'},3,"CDS loader; frame correct");
        not_exists($res->{'2R'}{'9'},"CDS loader correct, no frame at position");

        not_exists($res->{'2L'}{'2'},"CDS loader correct, no frame at position");
        is($res->{'2L'}{'8'},4,"CDS loader; frame correct");
        is($res->{'2L'}{'7'},5,"CDS loader; frame correct");
        is($res->{'2L'}{'6'},6,"CDS loader; frame correct");
        is($res->{'2L'}{'5'},4,"CDS loader; frame correct");
        is($res->{'2L'}{'4'},5,"CDS loader; frame correct");
        is($res->{'2L'}{'3'},6,"CDS loader; frame correct");
        not_exists($res->{'2L'}{'9'},"CDS loader correct, no frame at position");
        
        #overlapping frames - in frame
        $str=
        "2R\tmaker\tCDS\t3\t8\t.\t+\t0\n".
        "2R\tmaker\tCDS\t6\t8\t.\t+\t0\n";
        $res=load_cds_gtf(\$str);
        not_exists($res->{'2R'}{'2'},"CDS loader correct, no frame at position");
        is($res->{'2R'}{'3'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'4'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'5'},3,"CDS loader; frame correct");
        is($res->{'2R'}{'6'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'7'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'8'},3,"CDS loader; frame correct");
        not_exists($res->{'2R'}{'9'},"CDS loader correct, no frame at position");
        
        $str=
        "2L\tmaker\tCDS\t3\t8\t.\t-\t0\n".
        "2L\tmaker\tCDS\t6\t8\t.\t-\t0\n";
        $res=load_cds_gtf(\$str);
        not_exists($res->{'2L'}{'2'},"CDS loader correct, no frame at position");
        is($res->{'2L'}{'8'},4,"CDS loader; frame correct");
        is($res->{'2L'}{'7'},5,"CDS loader; frame correct");
        is($res->{'2L'}{'6'},6,"CDS loader; frame correct");
        is($res->{'2L'}{'5'},4,"CDS loader; frame correct");
        is($res->{'2L'}{'4'},5,"CDS loader; frame correct");
        is($res->{'2L'}{'3'},6,"CDS loader; frame correct");
        not_exists($res->{'2L'}{'9'},"CDS loader correct, no frame at position");
        
        
        #overlapping frames - out of frame
        $str=
        "2R\tmaker\tCDS\t3\t11\t.\t+\t0\n".
        "2R\tmaker\tCDS\t6\t11\t.\t-\t0\n".
        "2R\tmaker\tCDS\t9\t11\t.\t+\t0\n";
        $res=load_cds_gtf(\$str);
        not_exists($res->{'2R'}{'2'},"CDS loader correct, no frame at position");
        is($res->{'2R'}{'3'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'4'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'5'},3,"CDS loader; frame correct");
        is($res->{'2R'}{'6'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'7'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'8'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'9'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'10'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'11'},7,"CDS loader; frame correct");
        not_exists($res->{'2R'}{'12'},"CDS loader correct, no frame at position");
        
        $str=
        "2R\tmaker\tCDS\t3\t11\t.\t+\t0\n".
        "2R\tmaker\tCDS\t3\t7\t.\t+\t2\n";
        $res=load_cds_gtf(\$str);
        not_exists($res->{'2R'}{'2'},"CDS loader correct, no frame at position");
        is($res->{'2R'}{'3'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'4'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'5'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'6'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'7'},7,"CDS loader; frame correct");
        is($res->{'2R'}{'8'},3,"CDS loader; frame correct");
        is($res->{'2R'}{'9'},1,"CDS loader; frame correct");
        is($res->{'2R'}{'10'},2,"CDS loader; frame correct");
        is($res->{'2R'}{'11'},3,"CDS loader; frame correct");
        not_exists($res->{'2R'}{'12'},"CDS loader correct, no frame at position");
        

    }
    
    sub _getcodontable
    {
        return {
          'GCC' => 'A',
          'AGT' => 'S',
          'TGA' => '-',
          'TGT' => 'C',
          'CGA' => 'R',
          'ATC' => 'I',
          'AAC' => 'N',
          'AGC' => 'S',
          'TAC' => 'Y',
          'TCG' => 'S',
          'ACA' => 'T',
          'CCG' => 'P',
          'CTG' => 'L',
          'GCA' => 'A',
          'AAG' => 'K',
          'GTG' => 'V',
          'CAC' => 'H',
          'GTT' => 'V',
          'AGA' => 'R',
          'ACC' => 'T',
          'CCA' => 'P',
          'TGG' => 'W',
          'CGC' => 'R',
          'CTC' => 'L',
          'TTG' => 'L',
          'TAA' => '-',
          'CAG' => 'Q',
          'ACG' => 'T',
          'AAA' => 'K',
          'ATG' => 'M',
          'GTA' => 'V',
          'CTT' => 'L',
          'TAG' => '-',
          'GGA' => 'G',
          'GTC' => 'V',
          'TGC' => 'C',
          'TCA' => 'S',
          'ATT' => 'I',
          'TAT' => 'Y',
          'AAT' => 'N',
          'ACT' => 'T',
          'CAA' => 'Q',
          'GAC' => 'D',
          'GGT' => 'G',
          'TCC' => 'S',
          'TTT' => 'F',
          'AGG' => 'R',
          'CGT' => 'R',
          'ATA' => 'I',
          'CAT' => 'H',
          'CGG' => 'R',
          'CCC' => 'P',
          'GGG' => 'G',
          'TTA' => 'L',
          'GAG' => 'E',
          'CTA' => 'L',
          'GAT' => 'D',
          'TCT' => 'S',
          'TTC' => 'F',
          'GCG' => 'A',
          'GGC' => 'G',
          'GCT' => 'A',
          'GAA' => 'E',
          'CCT' => 'P'
        };
    }
}
1;