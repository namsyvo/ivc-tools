print "The program to extract the reads that failed to align\n";
#$fileName1=sam_file (aligned file);
#$fileName2=fastq_file (unaligned file);

$filePath=$ARGV[0];
#$fileName1=$ARGV[1];
#$fileName2=$filePath."Err/".$fileName1.".ua";
open (MYFILEO, ">>$filePath.unmap");
open (MYFILE, $filePath);
$count =0;
while (<MYFILE>) {
    chomp;
    $line =$_;
    my @values = split("\t", $line);
    if($values[1]=~/4/){
        $info = "@".$values[0];
         $len=@values;
         $read= $info."\n".$values[$len-2]."\n+\n".$values[$len-1]."\n";
         print MYFILEO $read;
         $count++;
         }
 }
close (MYFILE); 
close(MYFILEO);
print "Total unaligned reads (status 4): ".$count."\t";
print "\ndone! Check output at $fileName2\n";
