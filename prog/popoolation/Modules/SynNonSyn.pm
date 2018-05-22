{
    package SynNonSyn;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    use Data::Dumper;
    
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(load_codon_table load_nonsyn_length_table load_cds_gtf get_codon_changes get_codon_to_gene_assigner);
    
    
    sub load_codon_table
    {
        my $file=shift;
        open my $ifh,"<",$file or die "Could not open input file";
        #TTA: L
        #TTG: L
        #GTT: V
        my $count_file=0;
        
        my $codonhash={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            next unless $l;
            next if $l=~m/^#/;
            my($codon,$aa)=split /\s*:\s*/, $l;
            $codon=~s/\s+//g;
            $aa=~s/\s+//;
            die "Codons must have a length of three bases" unless length($codon)==3;
            die "Codons may only consist of ATCG" unless $codon=~m/^[ATCG]+$/;
            
            $codon=uc($codon);
            $aa=uc($aa);
            $codonhash->{$codon}=$aa;
            $count_file++;
        }
        
        die "Codon table is not complete; must contain 64 entries; provided one only contains: $count_file entries" unless $count_file==64;
        
        my $count_keys =scalar(keys(%$codonhash));
        
        die "One or more codon have been provided multiple times" unless $count_keys == 64;    
        return $codonhash;
    }
    
    sub load_nonsyn_length_table
    {
        my $file=shift;
        open my $ifh,"<",$file or die "Could not open input file";

        my $count_file=0;
        
        my $nslhash={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            next unless $l;
            next if $l=~m/^#/;
            my($codon,$nsl)=split /\s*:\s*/, $l;
            $codon=~s/\s+//g;
            $nsl=~s/\s+//;
            die "Codons must have a length of three bases" unless length($codon)==3;
            die "Codons may only consist of ATCG" unless $codon=~m/^[ATCG]+$/;
            die "Non synonymous length is not valid" if $nsl<0 or $nsl>3;
            
            $codon=uc($codon);
            $nslhash->{$codon}=$nsl;
            $count_file++;
        }
        
        die "Codon table is not complete; must contain 64 entries; provided one only contains: $count_file entries" unless $count_file==64;
        
        my $count_keys =scalar(keys(%$nslhash));
        
        die "One or more codon have been provided multiple times" unless $count_keys == 64;
        
        return $nslhash;
    }
    
    sub _parsegtf
    {
        my $line=shift;
        my @a=split /\t/, $line;
        my $ref=$a[0];
        my $feature=$a[2];
        my $start=$a[3];
        my $end=$a[4];
        my $strand=$a[6];
        my $frame=$a[7];
        my $attr=$a[8] || "";
        
            
        unless($ref or $start or $end)
        {
            die "the following line is not valid";
        }
        
        
        return
        {
            ref=>$ref,
            start=>$start,
            end=>$end,
            feature=>uc($feature),
            frame=>$frame,
            strand=>$strand,
            length=>$end-$start+1,
            attr=>$attr
        };
    }
    
    
    
    sub get_codon_to_gene_assigner
    {
        my $cdsrepresentation=shift;
        my $gtffile=shift;
        
        my $codonhash={};
        while(my($chr,$temp)=each(%$cdsrepresentation))
        {
            while(my($pos,$code)=each(%$temp))
            {
                next unless($code == 1 or $code == 6);
                $codonhash->{$chr}{$pos}={};
            }
        }
        
        
        my $geneidencoding={};
        my $currentgenenumber=1;
        open my $ifh, "<",$gtffile or die "Could not open gtf file";
        while(my $l=<$ifh>)
        {
            chomp $l;
            next if $l=~m/^##/;
            #2R	dm3_flyBaseGene	CDS	8375948	8376199	0.000000	-	0	gene_id "CG8824-RC"; transcript_id "CG8824-RC";
            my $ge=_parsegtf($l);
            next unless $ge->{feature} eq "CDS";

            # extract the values and check if appropriate
            my($chr,$start,$end,$strand,$frame,$attr)=($ge->{ref},$ge->{start},$ge->{end},$ge->{strand},$ge->{frame},$ge->{attr});
            my($geneid)=$attr=~/gene_id\s+"([^"]+)"/;
            die "gtf entry $l does not contain an attribute gene_id -> can not infer gene!" unless $geneid;
            
            unless(exists($geneidencoding->{$geneid}))
            {
                $geneidencoding->{$geneid}=$currentgenenumber;
                $currentgenenumber++;
            }
            my $genenumber=$geneidencoding->{$geneid};
            
            foreach my $i ($start..$end)
            {
                $codonhash->{$chr}{$i}{$genenumber}=1 if(exists($codonhash->{$chr}{$i}));
            }
        }
        
        # reformatcodonhash
        while(my($chr,$temp)=each(%$codonhash))
        {
            while(my($pos,$genecodes)=each(%$temp))
            {
                my $geneclist=[keys(%$genecodes)];
                $codonhash->{$chr}{$pos}=$geneclist;
            }
        }
        
        my $genenr_toid=[];
        while(my($geneid,$genenumber)=each(%$geneidencoding))
        {
            $genenr_toid->[$genenumber]=$geneid;
        }
        
        return sub
        {
            my $currentchr=shift;
            my $currentpos=shift;
            die "no gene entry found for chromosome $currentchr at position $currentpos" unless exists($codonhash->{$currentchr}{$currentpos});
            
            my $genelist=$codonhash->{$currentchr}{$currentpos};
            die "empty gene list for $currentchr at position $currentpos" unless @$genelist;
            my $toret = [map {$genenr_toid->[$_]} @$genelist];
            return $toret;
        }
    }
    
    sub load_cds_gtf
    {
        my $file=shift;
        
        open my $ifh,"<",$file or die "Could not open gtf-file";
        my $chrhash={};
        
        while(my $line=<$ifh>)
        {
            chomp $line;
            next if $line=~m/^##/;
            
            # parse gtf and filter the for CDS
            my $ge=_parsegtf($line);
            next unless $ge->{feature} eq "CDS";

            # extract the values and check if appropriate
            my($chr,$start,$end,$strand,$frame)=($ge->{ref},$ge->{start},$ge->{end},$ge->{strand},$ge->{frame});
            die "strand information not correct must be + or -; found instead: $strand" unless $strand=~m/^[-+]$/;
            die "frame information not correct must be 0,1 or 2; found instead: $frame" unless $frame=~m/^[012]$/;
            die "start position must be smaller or equal to end position" unless $start<=$end;
            
            if($strand eq "+")
            {
                my $fstart=$start+$frame;
                for(my $i=$fstart; $i<=$end-2; $i+=3)
                {
                    _set_chr_hash($chrhash,$chr,$i  ,1);
                    _set_chr_hash($chrhash,$chr,$i+1,2);
                    _set_chr_hash($chrhash,$chr,$i+2,3);
                }
            }
            elsif($strand eq "-")
            {
                my $fstart=$end-$frame;
                for(my $i=$fstart; $i>=$start+2; $i-=3)
                {
                    _set_chr_hash($chrhash,$chr,$i  ,4);
                    _set_chr_hash($chrhash,$chr,$i-1,5);
                    _set_chr_hash($chrhash,$chr,$i-2,6);                    
                }
                
            }
            else
            {
                die "Strand information is not valid";
            }
        }
        return $chrhash;
        
    }
    
    
    sub _set_chr_hash
    {
        my $chrhash=shift;
        my $chr=shift;
        my $pos=shift;
        my $value=shift;
        
        if(exists($chrhash->{$chr}{$pos}))
        {
            # so if the frame already exists and is in disagreement with the novel provided one we set it to 7 which indicates an invalid frame
            unless($chrhash->{$chr}{$pos} == $value)
            {
                $chrhash->{$chr}{$pos}=7;
            }
        }
        else
        {
            $chrhash->{$chr}{$pos}=$value;
        }
        
    }
    
    sub get_codon_changes
    {
        my $tr=shift;
        my $codonTable=shift;
        my $pileup=$tr->{pileup};
        my $codon=$tr->{codon};
        my $strand=$tr->{strand};
        my $start=$tr->{start};
        die "In order to incorporate a SNP the codon must be valid" unless $tr->{valid};
            
        my $codonchanges=[];
        foreach my $snp (@$pileup)
        {
            # pos, chr, refc, totcov, eucov, A, T, C, G, del, N, iscov, issnp, ispuresnp, alleles (list sorted by most frequent), consc, consc_confidence 
            next unless $snp->{ispuresnp};
            my $alleles=$snp->{alleles};
            die "Invalid number of alleles".Dumper($tr) unless(@$alleles ==4);
            
            my $snp_rel_pos=$snp->{pos}-$start;
            $snp_rel_pos=2-$snp_rel_pos if $strand eq "-"; # for the reverse stand the whole thing is reversed
            
            my $rc=substr($codon,$snp_rel_pos,1);
            my $alhash=_get_major_alleles($snp->{alleles},$strand);
            die "character of codon must be part of the two major alleles" unless exists($alhash->{$rc});
            
            my $notrc=(grep {$_ ne $rc} keys(%$alhash))[0];
            
            my $modcodon=$codon;
            substr($modcodon,$snp_rel_pos,1,$notrc);

            die "codon with snp is the same as codon without snps" if $modcodon eq $codon;
                    
            my $aa_ori=$codonTable->{$codon};
            my $aa_snp=$codonTable->{$modcodon};
            my $syn=$aa_ori eq $aa_snp?1:0;
            
            # chr, pos, ref, strand, cstart, A, T, C, G, eucov, codon, mcodon, aa, maa, syn
            push @$codonchanges,
            {
                chr=>       $snp->{chr},
                pos=>       $snp->{pos},
                ref=>       $snp->{refc},
                strand=>    $strand,
                cstart=>    $start,
                A=>         $snp->{A},
                T=>         $snp->{T},
                C=>         $snp->{C},
                G=>         $snp->{G},
                eucov=>     $snp->{eucov},
                codon=>     $codon,
                mcodon=>    $modcodon,
                aa=>        $aa_ori,
                maa=>       $aa_snp,
                syn=>       $syn
            };
        }
        $tr->{cc}=$codonchanges;
        return $tr;
    }
    
    #sub incorporate_snp_into_triplet
    #{
    #    my $tr=shift;
    #    my $codonTable=shift;
    #    my $pileup=$tr->{pileup};
    #    my $codon=$tr->{codon};
    #    my $strand=$tr->{strand};
    #    my $start=$tr->{start};
    #    die "In order to incorporate a SNP the codon must be valid" unless $tr->{valid};
    #        
    #    my $snp=[grep {$_->{ispuresnp}} @$pileup];
    #    die "Only one snp per codon can be incorporated" unless @$snp==1;
    #    $snp=$snp->[0];
    #    my $alleles=$snp->{alleles};
    #    
    #    die "Invalid number of alleles".Dumper($tr) unless(@$alleles ==4);
    #
    #    
    #    # if the frequency of the second and the third allele is equal; than it is not possible to estimate the snp mutated codon
    #    #return $tr if $alleles->[1]{c}== $alleles->[2]{c};
    #    
    #        
    #    # relative position in the codon
    #    my $snp_rel_pos=$snp->{pos}-$start;
    #    $snp_rel_pos=2-$snp_rel_pos if $strand eq "-"; # for the reverse stand the whole thing is reversed
    #        
    #    my $rc=substr($codon,$snp_rel_pos,1);
    #        
    #    my $alhash=_get_major_alleles($snp->{alleles},$strand);
    #    die "character of codon must be part of the two major alleles" unless exists($alhash->{$rc});
    #        
    #    my $notrc=(grep {$_ ne $rc} keys(%$alhash))[0];
    #        
    #    substr($codon,$snp_rel_pos,1,$notrc);
    #
    #    $tr->{snp_codon}=$codon;
    #    die "codon with snp is the same as codon without snps" if $tr->{snp_codon} eq $tr->{codon};
    #            
    #    my $aa_ori=$codonTable->{$tr->{codon}};
    #    my $aa_snp=$codonTable->{$tr->{snp_codon}};
    #    
    #    $tr->{syn}=($aa_ori eq $aa_snp)?1:0;
    #    $tr->{aa_old}=$aa_ori;
    #    $tr->{aa_new}=$aa_snp;
    #    $tr->{affected_snp}=$snp;
    #    return $tr;
    #    
    #}
    
    sub _get_major_alleles
    {
        my $alist=shift;
        my $strand=shift;
        
        my $alcount=@$alist;

        my $al_hash={};
        for my $i(0..1)
        {
            my $c=$alist->[$i]{a};
            $c=~tr/ATCG/TAGC/ if $strand eq "-";
            $al_hash->{$c}=1;
        }
        
        die "invalid number of keys" unless(keys(%$al_hash)==2);
        return $al_hash
    }
    
    
}
1;