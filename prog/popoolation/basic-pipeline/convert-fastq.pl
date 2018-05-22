{
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use BasicUtility;


    my $input;
    my $output;
    my $inencoding="auto";
    my $outencoding="phred";
    my $forceconversion=0;
    my $help=0;
    my $test=0;
    my $verbose=1;
    
    GetOptions(
        "input=s"               =>\$input,
        "output=s"              =>\$output,
        "in-encoding=s"         =>\$inencoding,
        "out-encoding=s"        =>\$outencoding,
        "force-conversion"      =>\$forceconversion,
        "test"                  =>\$test,
        "help"                  =>\$help
    ) or pod2usage(-msg=>"Options not valid $!",-verbose=>1);
    
    # dive into the alternative branches (if requested)
    pod2usage(-verbose=>2) if $help;
    pod2usage(-msg=>"At least one input file has to be provided", -verbose=>1) unless -e $input;
    pod2usage(-msg=>"An output file has to be provided", -verbose=>1) unless $output;
    $inencoding=lc($inencoding);
    $outencoding=lc($outencoding);
    
    
    my $paramfile=$output.".params";
    open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
    print $pfh "Using input\t$input\n";
    print $pfh "Using output\t$output\n";
    print $pfh "Using input encoding\t$inencoding\n";
    print $pfh "Using output encoding\t$outencoding\n";
    print $pfh "Using force-conversion\t$forceconversion\n";
    print $pfh "Using test\t$test\n";
    print $pfh "Using help\t$help\n";
    close $pfh;
    
    $inencoding =Utility::get_encoding($input,$inencoding);
    if($inencoding eq $outencoding)
    {
        die "Stoping Conversion\nEncoding of the input equals the encoding of the output: $inencoding = $outencoding\n" unless $forceconversion;
    }
    

    
    my $fqr=get_fastq_reader($input);
    my $fqw=get_fastq_writer($output);
    my $dec=Utility::get_qual_decoder($inencoding);
    my $enc=Utility::get_qual_encoder($outencoding);
    print "Will use input encoding $inencoding\n";
    print "Will use output encoding $outencoding\n";
    while(my $fq=$fqr->())
    {
        my $qual=$fq->{qual};
        my @ar= split  //,$qual;
        my @con=();
        
  
        foreach my $q (@ar)
        {
            my $t=$dec->($q);
            $t=$enc->($t);
            push @con,$t;
        }
        my $novqual=join("",@con);
        $fq->{qual}=$novqual;
        $fqw->($fq);
    }
    print "Finished conversion\n";
    
}

    
    
    
{
    package Utility;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/../Modules";
    use BasicUtility;
    
    
    sub get_qual_decoder
    {
        my $encoding=shift;
        
        if($encoding eq "phred")
        {
            return sub
            {
                my $c=shift;
                return ord($c)-33;
            }
        }
        elsif($encoding eq "illumina")
        {
            return sub
            {
                my $c=shift;
                return ord($c)-64;
            }
        }
        else
        {
            die "unknown encoding $encoding";
        }
    }
    
    sub  get_qual_encoder
    {
        my $encoding=shift;
        if($encoding eq "phred")
        {
            return sub
            {
                my $i=shift;
                return chr($i+33);
            }

        }
        elsif($encoding eq "illumina")
        {
            return sub
            {
                my $i=shift;
                return chr($i+64);
            }
        }
        else
        {
            die "unknown encoding $encoding";
        }
    }
    
    sub get_encoding
    {
        my $infile=shift;
        my $inencoding=shift;
        if($inencoding eq "phred" or $inencoding eq "illumina")
        {
            return $inencoding;
        }
        elsif($inencoding eq "auto")
        {
            print "Start autodetecting of encoding\n";
            my $detected=_autodetect_encoding($infile);
            die "Could not detect qualit encoding for file $infile" unless $detected;
            print "Detected encoding: $detected\n";
            return $detected;
        }
        else
        {
            die "unknow input encoding $inencoding"
        }

    }
    
    sub _autodetect_encoding
    {
        my $file=shift;
        my $fqr=get_fastq_reader($file);
        while(my $e=$fqr->())
        {
            my $qual=$e->{qual};
            my @ar=split //, $qual;
            foreach my $q (@ar)
            {
                my $phred=ord($q)-33;
                my $illumina=ord($q)-64;
                if($illumina< 0)
                {
                    return "phred";
                }
                if($phred > 40)
                {
                    return "illumina";
                }
            }
            
        }
        return undef;
    } 
    
}

    






=head1 NAME

convert-fastq.pl - Convert the quality string of fastq files

=head1 SYNOPSIS


 convert-fastq.pl --input input.fastq --output output.fastq
 


=head1 OPTIONS

=over 4

=item B<--input>

The input file in fastq format. Mandatory parameter

=item B<--output>

The output file. Will be in fastq. Mandatory parameter

=item B<--in-encoding>

The encoding of the quality in the input file; either phred/illumina/auto; default=auto

=item B<--out-encoding>

The encoding of the quality in the output file; either phred/illumina; default=phred

=item B<--force-conversion>

the script automatically stops when the input-encoding equals to the output-encoding; setting this flag an conversion may be enforced

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

=head2 INPUT

Input must be a fastq file. No line breaks in the nucleotide or quality sequence are allowed. No empty lines in the file are allowed.
The fastq file produced by the Illumina pipeline can be directyl use. Example input:

    @read1
    NNNNNAAAAAAAAAATTTTTTTTTTAAAAAAAAAANNNNNNNNNN
    +read1
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    @read2
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT
    +read2
    UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUSSSSSSSSSS



=head2 OUTPUT

the output file will also be in fastq as shown for INPUT. The quality string will be converted into the required encoding



=cut