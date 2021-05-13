#!/usr/bin/perl

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
my $transcript_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );

#my $cs_adaptor = $registry->get_adaptor( 'Human', 'Core', 'CoordSystem' );
#my $cs = $cs_adaptor->fetch_by_name('chromosome');
#printf "Coordinate system: %s %s\n", $cs->name(), $cs->version();

open(out,">>out.txt");

open(FILE,"part10.txt");
my @data=<FILE>;
close FILE;


#@data=('ENST00000354689 1 98','ENST00000603632 1 10');

$i=0;

foreach $tmp (@data)
{ 
    print("$i\n");
    $i++;		
    @stmp=();
    @stmp=split(" ",$tmp);
    my $stable_id=$stmp[0];
    my $first=(($stmp[1])*3);
    my $second=(($stmp[2])*3);
    my $transcript=();
    #$transcript = $transcript_adaptor->fetch_by_stable_id("ENSEMBLPEP",$stable_id);
    $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);

    #print"the $transcript\n";
    #Bio::EnsEMBL::Registry->set_disconnect_when_inactive();

    if($transcript)
    {
        #print "$chromosome\n";
        #printing the transcript id and corresponding protein sequence
        #print $transcript->translation()->stable_id(), "\n";
        #print">";
        #print $transcript->translation()->transcript()->stable_id(), "\n";
        #print $transcript->translate()->seq(),         "\n\n";

    #getting all the exon position for the protein sequence
    my @exon_pos=();
    foreach my $exon ( @{ $transcript->get_all_Exons() } )
    {
        my $estring = feature2string($exon);
        my @sestring=();
        @sestring=split(" ",$estring);
        push(@exon_pos,@sestring[1..2]);
        $sign=$sestring[3];
        $chromosome=$sestring[0];
        #print "\t$estring\n";
    }

    #getting 5' and 3' utr sequence and length
    #removing the utr position in the first and the last exons
    my $s5utr=0;
    if(my $n5UTR = $transcript->five_prime_utr)
    {
        #print"\n$n5UTR\n";
        my $n5UTRseq = $n5UTR->seq();
        #print" 5' UTR seq : $n5UTRseq\n";
        $s5utr=length($n5UTRseq);
        
    }
        #print"$s5utr\n";

    my $s3utr=0;
    if(my $n3UTR = $transcript->three_prime_utr)
    {
        my $n3UTRseq = $n3UTR->seq();
        #print" 3' UTR seq : $n3UTRseq\n";
        $s3utr=length($n3UTRseq);
       
    }
        #print"$s3utr\n";
    #print"@exon_pos\n";
    #print"$exon_pos[0] $one\t$exon_pos[$as-1] $two\n";
    
            #generating the chromosome position inbetween the start and end of each exon
        my @exons=();
        #print"$as\t@exon_pos\n";
        while(my ($i,$j) = splice(@exon_pos,0,2))
        {
            my @temp=();
            @temp=($i..$j);
            
            if($sign == -1)
            {
                #print"reverse successfull\n";
                @temp=reverse(@temp);
            }
            #print"@temp\n\n";
            push(@exons,@temp);
        }
    #to remove the 5utr bases in the chr positions
    if($s5utr)
    {
        my $cnt = 0;
        @exons = grep { ++$cnt > $s5utr } @exons;
    }
    
    #to remove the 3utr bases in the chr positions
    if($s3utr)
    {
        splice @exons, -$s3utr;
    }
    
        $check=$exons[$second-1];
        #$flag=0;
        if($check eq ())
        {
            $flag=1;
            #print"the check is $exons[$second-1]\n";
        }
        
        
        if($flag == 0)
        {
            print out ("$stable_id\t$chromosome:$exons[$first-3]-$exons[$second-1]\n");
        }
        elsif($exons[$second-3] eq ())
        {
            print out ("$stable_id\t$chromosome:$exons[$first-3]-$exons[$second-4]\n");
        }
        else
        {
            print out ("$stable_id\t$chromosome:$exons[$first-3]-$exons[$second-3]\n");            
        }
    }
    else
    {
        print out ("$stable_id\n");
    }
}
close(out);


sub feature2string
{
    my $feature = shift;
    my $stable_id  = $feature->stable_id();
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();
    return sprintf( "%s %d %d %+d",$seq_region, $start, $end, $strand);
}

