use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use BasicUtility;

my $output="";
my $input="";
my $refinput="";
my $sort=1;
my $crosscheck=0;
my $help=0;
my $test=0;

GetOptions(
    "ref-input=s"	    	=>\$refinput,
    "input=s"               	=>\$input,
    "output=s"              	=>\$output,
    "no-sort"		    	=>sub{$sort=0},
    "crosscheck"		=>\$crosscheck,
    "help"                  	=>\$help,
    "test"                  	=>\$test
    ) or die "Could not parse arguments";

pod2usage(-verbose=>2) if $help;
MauveTest::runTests() if $test;
pod2usage(-msg=>"Could not open reference file",-verbose=>1) unless -e $refinput;
pod2usage(-msg=>"Could not find input file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using ref-input\t$refinput\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using sort (no-sort)\t$sort\n";
print $pfh "Using crosscheck\t$crosscheck\n";
print $pfh "Using help\t$help\n";
print $pfh "Using test\t$test\n";
close $pfh;


my $ofstrans 		= Utility::get_offset_translator($refinput);
my $mauveParser 	= Utility::get_mauve_parser($input,$ofstrans);


my $outputunsorted=$output.".unsorted";

open my $ofh,">",$outputunsorted or die "Could not open output file";

print "Parsing the mauve output\n";
while(my $al=$mauveParser->())
{
    next unless $al->{useful};
    
    my $chr=$al->{chr};
    my $pos=$al->{start};
    my $seq=$al->{seq};
    my $trseq=Utility::transpose_sequences($seq);
    
    
    foreach my $tr (@$trseq)
    {
	next if $tr->[0] eq "-";
	my $temp=join("\t",@$tr);
	print $ofh "$chr\t$pos\t$temp\n";
	$pos++;
    }
}
close $ofh;

# optional step - crosschecking the results
if($crosscheck)
{
    print "Starting to crosscheck refence character\n";
    
    print "Loading the reference sequence $refinput into hash\n";
    my $refhash=load_fasta($refinput);
    
    print "Parsing unsorted output file $outputunsorted\n";
    open my $ifh, "<", $outputunsorted or die "Could not open output file $outputunsorted\n";
    
    while(my $l=<$ifh>)
    {
	chomp $l;
	my($chr,$pos,$refc)=split /\t/,$l;
	my $refcontrol=substr($refhash->{$chr},$pos-1,1);
	die "Not fitting character states in chromosome $chr; at position $pos; $refc not equal to $refcontrol\n" unless $refc eq $refcontrol;
    }
    close $ifh;
    print "Finished crosschecking: everything is fine\n";
}

# optional step - sorting
if ($sort)
{
    print "Sorting the output file\n";
    my $command="sort -k 1,1 -k 2,2n $outputunsorted > $output";
    print "Executing command: $command\n";
    system($command);
    unlink ($outputunsorted);
}
print "Finished\n";
exit;



{
    package Utility;
    use strict;
    use warnings;
    
    sub load_reference
    {
	my $ref=shift;
	open my $ifh, "<",$ref or die "Could not open input file";
    }
    
    
    sub transpose_sequences
    {
	my $seq=shift;
	
	my $sleng=length($seq->[0]);
	my $filenr=@$seq;
	my $transposed=[];
	
	for (my $i=0; $i<$sleng; $i++)
	{
	    my $entry=[];
	    
	    for(my $k=0; $k<$filenr; $k++)
	    {
		push @$entry,substr($seq->[$k],$i,1);
	    }
	    push @$transposed,$entry;
	}
	return $transposed;	
    }

    
    sub get_mauve_parser
    {
	my $infile=shift;
	my $ofstrans=shift;
	open my $ifh,"<",$infile or die "could not open mauve reader";
	
	# read the header
	# files
	my $filecount=0;
	
	print "Reading Mauve header\n";
	while(my $l=<$ifh>)
	{
	    chomp $l;
	    die "Mauve Header must start with a #" unless $l=~m/^#/;
	    last if $l=~m/#BackboneFile/;
	    $filecount++ if $l=~/^#Sequence\d+Format/
	}
	
	return sub {
	    # > 1:85628728-85628764 + /Volumes/Volume_4/analysis/divergence/ref/dmel-short.fasta
	    # TGCGAGCATGCACCAAATACTGATCCAAGGCACTCTG
	    # > 2:102772205-102772235 + /Volumes/Volume_4/analysis/divergence/ref/dsim-short.fasta
	    # TGCGGGCATGCACCAGATACTGATCCAAGGC------
	    # > 3:0-0 + /Volumes/Volume_4/analysis/divergence/ref/dyak-short.fasta
	    # -------------------------------------
	    # =
	    # initialize with the default
	    my $entries=[];
	    
	    my $header="";
	    my $seq="";
	    while(my $l=<$ifh>)
	    {
		chomp $l;
		next if $l=~m/^#/;
		
		if($l=~m/^=/)
		{
		    my $e=_parseEntry($header,$seq);
		    $entries->[$e->{fnr}]=$e;
		    
		    # annotate the entries
		    return _annotateEntries($entries,$filecount,$ofstrans);
		}
		elsif($l=~m/^>/)
		{
		    if($header)
		    {
			my $e=_parseEntry($header,$seq);
			$entries->[$e->{fnr}]=$e;
		    }
		    $header=$l;
		    $seq="";
		}
		else
		{
		    $seq.=$l;
		}
	    }
	    return undef;
	} 
    }
    
    sub _parseEntry
    {
	my $header=shift;
	my $seq=shift;
	
	#> 2:102772205-102772235 + /Volumes/Volume_4/analysis/divergence/ref/dsim-short.fasta
	my($filenum,$start,$end,$strand)=$header=~m/^>\s*(\d+):(\d+)-(\d+)\s*([+-])/;
	
	return {
	    fnr=>	$filenum,
	    start=>	$start,
	    end=>	$end,
	    seq=>	$seq,
	    strand=>	$strand,
	    header=>	$header
	};
    }
    
    sub _annotateEntries
    {
	my $entries=shift;
	my $filecount=shift;
	my $ofstrans=shift;
	
	my $an={
	    seq=>[],
	    useful=>0
	};
	
	return $an unless ($entries->[1]);
	
	my $refe=$entries->[1];
	my $strand=$refe->{strand};
	my ($s_chr,$s_pos)=$ofstrans->($refe->{start});
	my ($e_chr,$e_pos)=$ofstrans->($refe->{end});
	
	return $an unless $s_chr eq $e_chr;
	
	# get the sequences;
	# if the reference is reverse complement, everything has to be reverse complemented
	my $refleng=length($refe->{seq});
	my $seqs=[];
	for my $nr (1..$filecount)
	{
	    my $seq;
	    if($entries->[$nr])
	    {
		$seq=$entries->[$nr]{seq};
	    }
	    else
	    {
		$seq= "-"x$refleng;
	    }
	    
	    $seq = _rc_mauve($seq) if $strand eq "-";
	    
	    die "sequence length disagreing" unless length($seq) == $refleng;
	    push @$seqs,$seq;
	}
	
	# set the sequence
	$an->{seq}=$seqs;
	$an->{useful}=1;
	$an->{chr}=$s_chr;
	$an->{start}=$s_pos;
	return $an;
    }
    
    sub _rc_mauve
    {
        my $seq=shift;
        $seq=reverse($seq);
        $seq=~tr/-ATCGNatcgn/-TAGCNtagcn/;
        return $seq;
    }
    
    
    sub get_offset_translator
    {
	my $reffile=shift;
	
	# will already be sorted properly
	# [{chr,pos},{chr,pos}]
	# offset is the first position of a new reference using a one based counting
	print "Reading reference sequence $reffile to calculate the offsets used in the Mauve multi-fasta alignment\n";
	my $ofsets=_readreference($reffile); 
	
	
	return sub
	{
	  my $position=shift;
	  my $ofs=_search($ofsets,$position);
	  
	  my $chr=$ofs->{chr};
	  
	  # ofset is one-based; first base is 1;
	  #> fasta1
	  #ATCG  1234
	  #> fasta2
	  #TGAA 5678
	  #
	  # position in mauve is also one based
	  #
	  #if I have 6 I mean second position in  > fasta2
	  # 6-5+1=2
	  my $posinchr=$position-$ofs->{pos}+1;
	  
	  return ($chr,$posinchr);
	    
	};
	
    }
    
    sub _readreference
    {
	my $reffile=shift;
	
	open my $ifh,"<", $reffile or die "Could not open reference file";
	
	# fist position in the genome is 1! so this is a one based array
	my $activecounter=1; 
	my $ofsetar=[];
	while(my $l=<$ifh>)
	{
	    chomp $l;
	    if($l=~m/^>/)
	    {
		$l=~s/^>//;
		$l=~s/\s+.*$//;
		push @$ofsetar,{
		    chr=>$l,
		    pos=>$activecounter
		};
	    }
	    else
	    {
		$l=~s/^\s+//;
		$l=~s/\s+$//;
		$activecounter+=length($l);
	    }
	}
	# ofset is the first position in the new reference, using a one based counting
	return $ofsetar;
    }
    
    sub _search
    {
        # Binary search for the correct annotation
	# [{chr,pos},{chr,pos}]
        my $ofsetar=shift;
        my $position=shift;
	
	my $i=0;
	for(;$i<@$ofsetar; $i++)
	{
	    return $ofsetar->[$i-1] if $position<$ofsetar->[$i]{pos};
	}
	
	return $ofsetar->[$i-1];
    }
}

{
    package MauveTest;
    use strict;
    use warnings;
    use Test;
    
    
    sub runTests
    {
	test_offset_translator();
	test_transpose();
	test_mauveParser();
	exit;
    }
    
    sub test_offset_translator
    {
	my $str;
	my $ot;
	my($c,$p);
	$str=
	">f1\n".
	"ATAC\n".
	">f2\n".
	"GAAT\n".
	">f3\n".
	"CACG\n";
	$ot=Utility::get_offset_translator(\$str);
	($c,$p)=$ot->(1);
	
	is($c,"f1","offset calculator, chromosome correct");
	is($p,1, "offset calculator, position is correct");
	($c,$p)=$ot->(4);
	is($c,"f1","offset calculator, chromosome correct");
	is($p,4, "offset calculator, position is correct");
	($c,$p)=$ot->(5);
	is($c,"f2","offset calculator, chromosome correct");
	is($p,1, "offset calculator, position is correct");
	($c,$p)=$ot->(8);
	is($c,"f2","offset calculator, chromosome correct");
	is($p,4, "offset calculator, position is correct");
	($c,$p)=$ot->(9);
	is($c,"f3","offset calculator, chromosome correct");
	is($p,1, "offset calculator, position is correct");
	($c,$p)=$ot->(12);
	is($c,"f3","offset calculator, chromosome correct");
	is($p,4, "offset calculator, position is correct");
    }
    
    sub test_transpose
    {
	my $str;
	my $tr;
	$str=["AATTCC","--ATCG","TT--GG"];
	$tr=Utility::transpose_sequences($str);
	is($tr->[0][0],"A","Sequence transposing correct");
	is($tr->[0][1],"-","Sequence transposing correct");
	is($tr->[0][2],"T","Sequence transposing correct");

	is($tr->[2][0],"T","Sequence transposing correct");
	is($tr->[2][1],"A","Sequence transposing correct");
	is($tr->[2][2],"-","Sequence transposing correct");

	is($tr->[3][0],"T","Sequence transposing correct");
	is($tr->[3][1],"T","Sequence transposing correct");
	is($tr->[3][2],"-","Sequence transposing correct");
	
	is($tr->[4][0],"C","Sequence transposing correct");
	is($tr->[4][1],"C","Sequence transposing correct");
	is($tr->[4][2],"G","Sequence transposing correct");

	is($tr->[5][0],"C","Sequence transposing correct");
	is($tr->[5][1],"G","Sequence transposing correct");
	is($tr->[5][2],"G","Sequence transposing correct");

    }
    
    sub test_mauveParser
    {
	my($ref,$ot,$c,$p);
	my($ms,$mp,$al);
	$ref=
	">f1\n".
	"ATAC\n".
	">f2\n".
	"GAAT\n".
	">f3\n".
	"CACG\n";
	$ot=Utility::get_offset_translator(\$ref);
	
	$ms=
	"#Sequence1Format\n".
	"#Sequence2Format\n".
	"#Sequence3Format\n".
	"#BackboneFile\n".
	"> 1:5-8 + \n".
	"GATT\n".
	"> 2:100-10000 - \n".
	"C-TT\n".
	"> 3:12-15 + \n".
	"G--T\n".
	"=\n";
	$mp= Utility::get_mauve_parser(\$ms,$ot);
	$al=$mp->();
	
	is($al->{chr},"f2","Mauve parser: chromosome correct");
	is($al->{start},1,"Mauve parser: start position is correct");
	is($al->{useful},1,"Mauve parser: is useful");
	is($al->{seq}[0],"GATT","Mauve parser: sequence is correct");
	is($al->{seq}[1],"C-TT","Mauve parser: sequence is correct");	
	is($al->{seq}[2],"G--T","Mauve parser: sequence is correct");	


	$ms=
	"#Sequence1Format\n".
	"#Sequence2Format\n".
	"#Sequence3Format\n".
	"#BackboneFile\n".
	"> 1:5-8 - \n".
	"GATT\n".
	"> 3:100-10000 - \n".
	"C-TT\n".
	"=\n";
	$mp= Utility::get_mauve_parser(\$ms,$ot);
	$al=$mp->();
	
	is($al->{chr},"f2","Mauve parser: chromosome correct");
	is($al->{start},1,"Mauve parser: start position is correct");
	is($al->{useful},1,"Mauve parser: is useful");
	is($al->{seq}[0],"AATC","Mauve parser: sequence is correct");
	is($al->{seq}[1],"----","Mauve parser: sequence is correct");	
	is($al->{seq}[2],"AA-G","Mauve parser: sequence is correct");	
    }
}



    #"ref-input=s"	    	=>\$refinput,
    #"input=s"               	=>\$input,
    #"output=s"              	=>\$output,
    #"no-sort"		    	=>sub{$sort=0},
    #"crosscheck"		=>\$crosscheck,
    #"help"                  	=>\$help,
    #"test"                  	=>\$test
=head1 NAME

perl mauve-parser.pl - Parse a multiple genome alignment as generated by Mauve.

=head1 SYNOPSIS

perl mauve-parser.pl --ref-input d.melanogaster.fa  --input d.three-way-align.mfa --output mel-guided-alignment.txt

=head1 OPTIONS

=over 4

=item B<--ref-input>

A fasta file containing the full genome of the first sequence provided to mauve! This will act as the anchor or reference of the multiple alignment; Mandatory

=item B<--input>

A multiple genome alignment created by Mauve (mfa); Mandatory

=item B<--output>

The output file; Mandatory

=item B<--no-sort>

Disable sorting of the output;

=item B<--crosscheck>

Enable crosschecking of the alignment. This step is memory intense but recommendet. It checks if every base in the C<--output> is in agreement with the corresponding base in the C<--ref-input>.

=back


=head1 Details

=head2 Output

A multiple genome alignment where the first sequence provided to Mauve acts as the reference.

 YHet    291809  T       A       T
 YHet    291810  T       T       T
 YHet    291811  T       A       A
 YHet    291812  T       C       T
 YHet    291813  A       A       A
 
 col 1: chromosome (contig); with respect to the reference genome
 col 2: position in the given chromosome (contig)
 col 3: character state in the first genome (the reference character)
 col 4: character state in the second genome
 col 5: and so on...

=cut


