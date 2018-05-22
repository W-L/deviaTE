{
    package FormatSNP;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT =qw(get_WindowSNPFormater get_gtfSNPFormater get_SNPwriter get_syn_nonsyn_SNPFormater get_synnonsyngene_SNPFormater);
    
    
            #pos=>$pos,
            #chr=>$chr,
            #refc=>$rc,
            #consc=>"N",
            #consc_confidence=>0,
            #totcov=>0,
            #eucov=>0,
            #A=>0,
            #T=>0,
            #C=>0,
            #G=>0,
            #del=>0,
            #N=>0,
            #iscov=>0,
            #issnp=>0,
            #ispuresnp=>0
            
    sub get_SNPwriter
    {
        my $file =shift;
        open my $ofh, ">",$file or die "Could not open SNP output file";
        return sub
        {
            my $snp=shift;
            print $ofh _formatSNP($snp);
        }
    }
            
    sub get_gtfSNPFormater
    {
        my $file =shift;
        open my $ofh, ">",$file or die "Could not open SNP output file";
        return sub
        {
            my $geneid=shift;
            my $snps=shift;
            my $snpcount=@$snps;
            print $ofh ">$geneid snps:$snpcount\n";
            foreach my $snp (@$snps)
            {
                print $ofh _formatSNP($snp);
            }
            print $ofh "\n";
            
            
        }
        
    }
    
    sub get_WindowSNPFormater
    {
        
        my $file=shift;
        open my $ofh, ">", $file or die "Could not open SNP output file";
        
        return sub
        {
            my $win=shift;
            my $chr=$win->{chr};
            my $start=$win->{start};
            my $end=$win->{end};
            my $middle=$win->{middle};
            my $window=$win->{window};
            my $data=$win->{data};
            my $snpcount=$win->{count_snp};
            
            my $snps=[];
            foreach(@$data)
            {
                push @$snps,$_ if $_->{ispuresnp};
            }
            
            print $ofh ">$chr:$middle $chr:$start-$end snps:$snpcount\n";
            foreach my $snp (@$snps)
            {
                print $ofh _formatSNP($snp);
            }
            print $ofh "\n";
        }
    }
    
    
    sub _formatSNP
    {
        my $snp=shift;
        return "$snp->{chr}\t$snp->{pos}\t$snp->{refc}\t$snp->{eucov}\t$snp->{A}\t$snp->{T}\t$snp->{C}\t$snp->{G}\t$snp->{N}\n";
    }
    
    sub get_synnonsyngene_SNPFormater
    {
        my $outfile=shift;
        open my $ofh, ">", $outfile or die "Could not open output file $outfile";
        
        # the closure
        return sub
        {
            my $tr=shift;
            my $genelist=shift;
            my $countcodons=$tr->{cc};
            next unless @$countcodons;
            
            foreach my $cc (@$countcodons)
            {
                my $formated=_format_cc($cc);
                foreach my $genename (@$genelist)
                {
                    print $ofh "$genename\t$formated\n";
                }
            }
        }
    }
    
    
    sub get_syn_nonsyn_SNPFormater
    {
        my $outfile=shift;
        open my $ofh, ">",$outfile or die "Could not open output file $outfile";
        
        return sub
        {
            my $triplets=shift;
            my $chr=shift;
            my $start=shift;
            my $end=shift;
            
            my $snpcount=0;
            foreach my $tr (@$triplets)
            {
                my $cc=$tr->{cc};
                $snpcount+=@$cc;
            }
            return unless $snpcount;
            
            print $ofh ">$chr:$start-$end snps: $snpcount\n";
            foreach my $tr (@$triplets)
            {
                my $codons=$tr->{cc};
                foreach my $cc (@$codons)
                {
                    my $formated = _format_cc($cc);
                    print $ofh $formated."\n";
                }
            }
            print $ofh "\n";
        }
    }
    
    
    sub _format_cc
    {
        my $cc=shift;
        my $codon_old=$cc->{codon};
        my $codon_novel=$cc->{mcodon};
        my $strand=$cc->{strand};
        my $aa_old=$cc->{aa};
        my $aa_new=$cc->{maa};
        my $syn = $cc->{syn};
        $syn=$syn?"syn":"non-syn";
        
        my $codon="$codon_old->$codon_novel";
        my $aa="$aa_old->$aa_new"; 
        
        
        my $toprint="$cc->{chr}\t$cc->{pos}\t$cc->{ref}\t$cc->{eucov}\t$cc->{A}\t$cc->{T}\t$cc->{C}\t$cc->{G}\t$syn\t$strand\t$codon\t$aa";
        return $toprint;
    }
}

1;