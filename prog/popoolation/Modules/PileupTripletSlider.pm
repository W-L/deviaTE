{   
    package PileupTripletSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Pileup;
    use BasicUtility;
    
    my $fwd_valid={ 1=>1, 2=>1, 3=>1 };
    my $rev_valid={ 4=>1, 5=>1, 6=>1 };
    
    # disabling minimum consensus confidence
    my $consc_confidence_threshold=0.2;
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        my $annotation=shift;
        my $pp=shift;
        
        open my $fh,"<",$file or die "Could not open file handle";
        
        return bless {
            annotation=>$annotation,
            pp=>$pp,
            file=>$file,
            fh=>$fh,
            buffer=>[]
        },__PACKAGE__;
    }
    
    sub next
    {
        my $self=shift;
        
        my $codon;
        while(1)
        {
            $codon=$self->_next_codon();
            return undef unless $codon; # undef if empty
            
            my $valid=$self->_valid_frame($codon->{frame});
        
            last if $valid; # if the codon is not fucked we found a good one
        }
        
        my $toret=$self->_annotate_triplet($codon);
        return $toret;
    }
    

    
    sub _next_codon
    {
        my $self=shift;
        my $annotation=$self->{annotation};
        my $pp=$self->{pp};
        
        ##
        ## Search for the next valid entry
        ##
        my $l1;
        my($chr,$pos);
        while($l1=$self->_nextline())
        {
            chomp $l1;
            ($chr,$pos)=split/\s+/,$l1;
            next unless $annotation->{$chr}{$pos};
            my $cp=$annotation->{$chr}{$pos};
            
            # find the frame
            last if($cp==1 or $cp==6);

            #check the value
        }
        return undef unless $l1;
        # so found the first valid position in a pileup file;
        
        
        # step one get the frames;
        my $t={};
        $t->{frame} = [$annotation->{$chr}{$pos}, $annotation->{$chr}{$pos+1}, $annotation->{$chr}{$pos+2}];
        
        # step two: store the first line of the pileup and read the next two lines of the pileup

        my $p1=$pp->($l1);
        my $endpos=$pos+2;
        my $triplet=[$p1, get_empty_pileup_entry($chr,$pos+1), get_empty_pileup_entry($chr,$pos+2)];
        while(1)
        {
            my $l=$self->_nextline();
            last unless $l;
            chomp $l;
            return undef unless $l;
            my $parsed=$pp->($l);
            #check the
            if($parsed->{chr} ne $chr or $parsed->{pos}>$endpos)
            {
                $self->_bufferline($l);
                last;
            }
            my $relpos=$parsed->{pos}-$pos;
            $triplet->[$relpos]=$parsed;
        }
        $t->{pileup}=$triplet;
        
        return $t;
    }
    
    
    sub _valid_frame
    {
        my $self=shift;
        my $frame =shift;
        die "invalid number of sequences in the frame" unless @$frame==3;
        
        my $countseven=[grep {$_==7} @$frame];
        return 0 if @$countseven;
        
        # check frame
        my $countFwd=0;
        my $countRev=0;
        foreach my $f (@$frame)
        {
            if(exists($fwd_valid->{$f}))
            {
                $countFwd++;
            }
            elsif(exists($rev_valid->{$f}))
            {
                $countRev++;
            }
            else
            {
                die "invalid frame in codon $f";
            }
        }
        
        return 1 if($countFwd==3 or $countRev==3);
        return 0;
    }
    
    
    
    sub _annotate_triplet
    {
        my $self=shift;
        my $triplet=shift;
        
        my $pileups=$triplet->{pileup};
        
        
        # check chr and position
        my $start=$pileups->[0]{pos};
        my $chr=$pileups->[0]{chr};
        for(my $i=0; $i<=2; $i++)
        {
            my $act=$pileups->[$i];
            die "pileup triplet is fucked, chromsomes do not fit" unless $act->{chr} eq $chr;
            die "pileup triplet is fucked, positions are not ascending" unless($start+$i==$act->{pos});
        }
        
        # parse the pileups
        my $validCoverage=1;
        my $validConsensus=1;
        my $count_snps=0;
        my $codon="";
        foreach my $p (@$pileups)
        {
            $validCoverage=0 unless $p->{iscov};
            # consc_confidence -> is the SNP AAATTT or AAAATTT -> in the second the consensus basecall is clearly A
            $validConsensus=0 unless $p->{consc_confidence} >= $consc_confidence_threshold;
            $count_snps++ if $p->{ispuresnp};
            $codon.=$p->{consc};
        }
        
        
        my $strand;
        if($triplet->{frame}[0]==1)
        {
            $strand="+";
        }
        elsif($triplet->{frame}[0]==6)
        {
            $strand="-";
        }
        else
        {
            die "invalid frame, can not set strand";
        }
        
        # codon orientation
        $codon = reverse_complement_dna($codon) if ($strand eq "-");
        my $validCodon=$codon=~m/^[ATCG]{3}$/?1:0;
        
        # is valid
        my $valid=($validCoverage and $validCodon and $validConsensus)?1:0;

        
        $triplet->{chr}=$chr;
        $triplet->{start}=$start;
        $triplet->{strand}=$strand;
        $triplet->{valid}=$valid;
        $triplet->{valid_coverage}=$validCoverage;
        $triplet->{valid_codon}=$validCodon;
        $triplet->{valid_consensus}=$validConsensus;
        $triplet->{count_snps}=$count_snps;
        $triplet->{codon}=$codon;
        return $triplet;
    }
    
    
    # TRIPLET definition
    # frame, pileup, chr, start, strand, valid, valid_frame, valid_coverage, valid_codon, count_snps, codon
    
    
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
    
    
}
1;