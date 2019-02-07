#!/usr/bin/env perl

# this script is part of PoPoolation by Robert Kofler et al.
# https://sourceforge.net/projects/popoolation/

{
    package deviaTE_trim;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT     =qw(get_basic_parser get_pileup_parser get_extended_parser get_empty_pileup_entry get_quality_encoder get_basic_mpileupparser);

my $qualhash={
              illumina=>sub {ord(shift) - 64;},
              sanger=>sub {ord(shift) - 33;},
             };

# qualencoding,mincount,mincov,maxcov,minqual
# my $pp=get_pileup_parser();

 
 
sub get_basic_mpileupparser
{
    my $quality_encoding=shift;
    my $minqual=shift;
    my $encoder;
    $quality_encoding=lc($quality_encoding);
    die "Encoder $quality_encoding not supported" unless exists($qualhash->{$quality_encoding});
    $encoder=$qualhash->{$quality_encoding};
    
    my $subparser=_get_subnucsqualparser($encoder,$minqual);
    return sub
    {
        my $line=shift;
        die "pileup parser empty line provided $line" unless $line;
        #    chrM    1       G       76      ^F.^F.^F.^F.    hhhgf    228     ^F,^F,^F,^F,^F.  ghfhhhhh    157     F,^F,^F,^F.^F, gYheghadh
        
        my @ar=split /\s+/,$line;
        my $chr=shift @ar;
        my $pos=shift @ar;
        my $rc=shift @ar;
        $rc=uc($rc);
        my $entries=[];
        while(@ar)
        {
            my $cov=shift @ar;
            my $nucs=shift @ar;
            my $qual=shift @ar;
            die "F entry $line" unless $qual;
            my $entry=$subparser->($cov,$nucs,$qual,$rc);
            push @$entries,$entry;
        }
        
        return
        {
            chr=>$chr,
            pos=>$pos,
            refc=>$rc,
            entries=>$entries
        };
    }
}

sub get_quality_encoder
{
    my $q=shift;
    $q=lc($q);
    die "Encoder $q not supported" unless exists($qualhash->{$q});
    return $qualhash->{$q};
}

sub _get_subnucsqualparser
{
    my $encoder=shift;
    my $minqual=shift;
    
    return sub
    {
        
        my $cov=shift;
        my $nucs=shift;
        my $qual=shift;
        my $rc=shift;
        
        #first get rid of empty entries
        return {A=>0,T=>0,C=>0,G=>0, N=>0, del=>0,eucov=>0,totcov=>0} if $cov==0;
    
        
        # get rid of the crap present in the sequence line
        $nucs=~s/[-+](\d+)(??{"[ACGTNacgtn]{$1}"})//g; # I offer a beer for anyone who understand this line, I am a genius!!!!
        $nucs=~s/\^.//g;
        $nucs=~s/\$//g;
        $nucs=~s/[.]/uc($rc)/eg;
        $nucs=~s/[,]/lc($rc)/eg;
        
        die "size of sequence does not equal size of quality:  $nucs, $qual\n" unless length($nucs) == length($qual);
        my @nucs=split//,$nucs;
        my @qual=split//,$qual;
        
            # filter the pileup file by quality
        my $co=0;
        my $a=0;
        my $t=0;
        my $c=0;
        my $g=0;
        my $del=0;
        my $n=0;
        
        for(my $i=0; $i<@qual; $i++)
        {
            my $qc=$qual[$i];
            my $nc=$nucs[$i];
            
            my $quality = $encoder->($qc); 
            next if $quality <$minqual;
            $co++;
            
            if($nc=~/A/i)
            {
                $a++;
            }
            elsif($nc=~/T/i)
            {
                $t++;
            }
            elsif($nc=~/C/i)
            {
                $c++;
            }
            elsif($nc=~/G/i)
            {
                $g++;
            }
            elsif($nc=~/N/i)
            {
                $n++;
            }
            elsif($nc=~/\*/)
            {
                $del++;
            }
            elsif($nc=~/>/ or $nc=~/</)
            {
                #ignore RNA-seq shit   
            }
            else
            {
                die "Could not parse pileup; Unknown allele: $nc";
            }
        }

        # pos, chr, refc, totcov, eucov, A, T, C, G, del, N        
        my $entry={
            totcov=>$co,
            eucov=>($a+$t+$c+$g),
            A=>$a,
            T=>$t,
            C=>$c,
            G=>$g,
            del=>$del,
            N=>$n
        };
    }
}

sub get_basic_parser
{
    # qualencoding,mincount,mincov,maxcov,minqual
    my $quality_encoding=shift;
    my $minqual=shift;
    my $encoder;
    $quality_encoding=lc($quality_encoding);
    die "Encoder $quality_encoding not supported" unless exists($qualhash->{$quality_encoding});
    $encoder=$qualhash->{$quality_encoding};
    
    my $subparser=_get_subnucsqualparser($encoder,$minqual);
    
    return sub {
        my $line=shift;
        die "pileup parser empty line provided $line" unless $line;
        #4       1096    N       22      tTtttTtTtt+1gttttTttttTtt       bWZ^``_]J^_aYaa^__ZbWX
        #X        2832683 N       22      *AaAAAaaaAaaaaaaaaaaA^FA        yBaB_\aaa^a[a^____^a[^
        my($chr,$pos,$rc,$cov,$nucs,$qual)=split/\s+/,$line;
        $rc=uc($rc);
        
        my $entry=$subparser->($cov,$nucs,$qual,$rc);
        $entry->{pos}=$pos;
        $entry->{chr}=$chr;
        $entry->{refc}=$rc;
        # pos, chr, refc, totcov, eucov, A, T, C, G, del, N        
        return $entry;
    }
}




sub get_pileup_parser
{
    # qualencoding,mincount,mincov,maxcov,minqual
    my $quality_encoding=shift;
    my $mincount=shift;
    my $mincoverage=shift;
    my $maxcoverage=shift;
    my $minqual=shift;
    my $tolerateDeletions=shift || 0;
    
    
    my $pp=get_basic_parser($quality_encoding,$minqual);
    
    return sub {
        my $line=shift;
        my $pu=$pp->($line);
        
        # reset the counts to zero if not fullfilling the minimum count
        $pu->{A} =0 unless $pu->{A}>=$mincount;
        $pu->{T} =0 unless $pu->{T}>=$mincount;
        $pu->{C} =0 unless $pu->{C}>=$mincount;
        $pu->{G} =0 unless $pu->{G}>=$mincount;
        $pu->{eucov} = ($pu->{A} + $pu->{T} + $pu->{C} + $pu->{G});
        
        
        $pu->{iscov}=0;
        $pu->{issnp}=0;
        if($pu->{eucov}>=$mincoverage && $pu->{eucov}<=$maxcoverage)
        {
            $pu->{iscov}=1;
            my $alcount=0;
            $alcount++ if $pu->{A}>=$mincount;
            $alcount++ if $pu->{T}>=$mincount;
            $alcount++ if $pu->{C}>=$mincount;
            $alcount++ if $pu->{G}>=$mincount;
            $pu->{issnp}=1 if $alcount>=2;
            
            # pure snps
            $pu->{ispuresnp}=0;
            if($pu->{issnp} and ($pu->{del} < 1 or $tolerateDeletions))
            {
                $pu->{ispuresnp}=1;
            }
            
        }
        
        # pos, chr, refc, totcov, eucov, A, T, C, G, del, N, iscov, issnp, ispuresnp 
        return $pu;
    };
}



sub get_extended_parser
{
    my $quality_encoding=shift;
    my $mincount=shift;
    my $mincoverage=shift;
    my $maxcoverage=shift;
    my $minqual=shift;
    my $tolerateDeletions=shift || 0;
    
    # qualencoding,mincount,mincov,maxcov,minqual
    my $pp=get_pileup_parser($quality_encoding,$mincount,$mincoverage,$maxcoverage,$minqual,$tolerateDeletions);
    
    return sub
    {
        my $line=shift;
        my $pu=$pp->($line);
        
        
        my $ca=[{a=>"A",c=>$pu->{A}},{a=>"T",c=>$pu->{T}},{a=>"C",c=>$pu->{C}},{a=>"G",c=>$pu->{G}}];
        $ca=[sort {$b->{c}<=>$a->{c}} @$ca];
        $pu->{alleles}=$ca;

        if($pu->{iscov})
        {
            # calculate the consensus character from the pileup
            $pu->{consc}=$ca->[0]{a};
            $pu->{consc_confidence} = ($ca->[0]{c}/$pu->{eucov});
        }
        else
        {
                    $pu->{consc}="N";
                    $pu->{consc_confidence}=0;
        }
        
        # pos, chr, refc, totcov, eucov, A, T, C, G, del, N, iscov, issnp, ispuresnp, alleles (list sorted by most frequent), consc, consc_confidence 
        die "invalid number of allles for entry $line" unless scalar(@{$pu->{alleles}})==4;
        return $pu;
    }
}


    
    
    sub get_empty_pileup_entry
    {
        my $chr=shift;
        my $pos=shift;
        my $rc=shift || "N";
        
        return {
            pos=>$pos,
            chr=>$chr,
            refc=>$rc,
            consc=>"N",
            consc_confidence=>0,
            totcov=>0,
            eucov=>0,
            A=>0,
            T=>0,
            C=>0,
            G=>0,
            del=>0,
            N=>0,
            iscov=>0,
            issnp=>0,
            ispuresnp=>0
        };
    }





}

1;