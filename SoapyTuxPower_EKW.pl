######						*SoapyTuxedo*		    			######
######	Author: H Marx @ Barker lab								######
######	Version: 16.03.2017                                	    ######
###### 	Assembly de novo transcriptome and calculate gene power #

############################ On ratel: ###############################
use File::Basename;
  
$NAME1 = "$ARGV[0]"; #sample 1 full data
$NAME2 = "$ARGV[1]"; #sample 2 full data
$NAMESP = "$ARGV[2]"; #species base name
$CPU = "$ARGV[3]";  # maximum 16 
#open NAME1 or die "No file $NAME\n";

print "\n###### SoapyTuxedo is launched ######\n";

####### Gene power with Tuxedo pipeline 

print "\n Beginning Tuxedo... \n";
my $rdlength = 150;
#my $kmer = int $rdlength/3*2;
my $kmer = int 57;

## Align (map) clean reads from each sample (time point) to reference genome (bowtie-indexed scaffold seqs) with tophat 
#system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_thout $NAMESP\.index $NAME1\_R1.clean $NAME1\_R2.clean");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_rep1_thout $NAMESP\.index $NAME1\_R1.clean.subset1 $NAME1\_R2.clean.subset1");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_rep2_thout $NAMESP\.index $NAME1\_R1.clean.subset2 $NAME1\_R2.clean.subset2");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_rep3_thout $NAMESP\.index $NAME1\_R1.clean.subset3 $NAME1\_R2.clean.subset3");

#system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_thout $NAMESP\.index $NAME2\_R1.clean $NAME2\_R2.clean");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_rep1_thout $NAMESP\.index $NAME2\_R1.clean.subset1 $NAME2\_R2.clean.subset1");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_rep2_thout $NAMESP\.index $NAME2\_R1.clean.subset2 $NAME2\_R2.clean.subset2");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; tophat2 -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_rep3_thout $NAMESP\.index $NAME2\_R1.clean.subset3 $NAME2\_R2.clean.subset3");

print "\n Finished tophat2 \n";

## Assemble genes and transcripts into transcriptome assembly for each condition (time point)
#system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_clout $NAME1\_thout/accepted_hits.bam");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_rep1_clout $NAME1\_rep1_thout/accepted_hits.bam");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_rep2_clout $NAME1\_rep2_thout/accepted_hits.bam");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME1\_rep3_clout $NAME1\_rep3_thout/accepted_hits.bam");

#system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_clout $NAME2\_thout/accepted_hits.bam");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_rep1_clout $NAME2\_rep1_thout/accepted_hits.bam");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_rep2_clout $NAME2\_rep2_thout/accepted_hits.bam");
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cufflinks -p 8 -o ~/HarvardForest/Eldridge/$NAMESP/$NAME2\_rep3_clout $NAME2\_rep3_thout/accepted_hits.bam");

print "\n Finished cufflinks \n";

## create cuffmerge file: assemblies.txt
open my $fileHandle, ">>", "./assemblies.txt" or die "Can't open './assemblies.txt'\n";
print $fileHandle "./$NAME1\_clout/transcripts.gtf\n";
print $fileHandle "./$NAME1\_rep1_clout/transcripts.gtf\n";
print $fileHandle "./$NAME1\_rep2_clout/transcripts.gtf\n";
print $fileHandle "./$NAME1\_rep3_clout/transcripts.gtf\n";
print $fileHandle "./$NAME2\_clout/transcripts.gtf\n";
print $fileHandle "./$NAME2\_rep1_clout/transcripts.gtf\n";
print $fileHandle "./$NAME2\_rep2_clout/transcripts.gtf\n";
print $fileHandle "./$NAME2\_rep3_clout/transcripts.gtf\n";
close $fileHandle;

## Merge all assemblies into single transcriptome annotation:
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cuffmerge -p 8 -s $NAMESP\-$kmer.scafSeq assemblies.txt");
print "\n Finished cuffmerge \n";

## Identify differentially expressed genes 
system ("cd ~/HarvardForest/Eldridge/$NAMESP; cuffdiff -p 8 -o ./diff_out -L $NAME1,$NAME2 -u ./merged_asm/merged.gtf \ 
./$NAME1\_thout/accepted_hits.bam,./$NAME1\_rep1_thout/accepted_hits.bam,./$NAME1\_rep2_thout/accepted_hits.bam,./$NAME1\_rep3_thout/accepted_hits.bam \
./$NAME2\_thout/accepted_hits.bam,./$NAME2\_rep1_thout/accepted_hits.bam,./$NAME2\_rep2_thout/accepted_hits.bam,./$NAME2\_rep3_thout/accepted_hits.bam");

### Remove raw reads from each sample (they are still stored in /soap_assembly_1)
#system ("cd ~/HarvardForest/power/$NAMESP; rm  $NAME1\_R*");
#system ("cd ~/HarvardForest/power/$NAMESP; rm  $NAME2\_R*");
#system ("cd ~/HarvardForest/power/$NAMESP; rm  $NAMESP\_R*"); # Remove concatenated reads

print "\n Finished $NAMESP gene power!!! \n";

