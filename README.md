# BtretaeAnalyses
Scripts and notes for population genomic analyses of B. treatae experiment

This study combines a replicated selection experiment and genotype-host associations in natural populations of Belonocnema treatae galling wasps. The data are in /uufs/chpc.utah.edu/common/home/gompert-group2/data/btreatae_gbs/.

Samples with no final character are from natural populations.

s = starting population, root gall-based sexual generation (sample from founding population) could be all female (diploid, probably HW, but maybe not) or could include haploid (males).

c = survivors on natal host (next generation).
t = survivors on non-natal host (next generation).
p2 = previous leaf gall asexual generation (probably ignore for selection experiment).

# Alignment and variant calling

*I first filtered the data for PhiX with PhiXFilterFork.pl.

````perl
#!/usr/bin/perl
#
# filter phix by aligning to the phix reference, uses fork
#

system "module load bwa\n";
system "module load samtools\n";

use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);
my $phix = "/uufs/chpc.utah.edu/common/home/u6000989/data/phixSeqIndex/NC_001422.fasta";

FILES:
foreach $fq (@ARGV){
        $pm->start and next FILES; ## fork
        $fq =~ m/^([a-zA-Z_\-0-9]+)/ or die "failed to match $fq\n";
        $base = $1;
        print "Running alignment for $fq to phix\n";
        system "bwa aln -n 5 -l 20 -k 2 -t 4 $phix $fq -f $base.sai\n";
        system "bwa samse -f $base.sam -r \'\@RG\\tID:PHIX\' $phix $base.sai $fq\n";
        system "samtools view -S -b $base.sam > $base.bam\n";
        print "Removing reads mapped to phix for $base\n";
        system "samtools view -f4 $base.bam > $base.sam\n";
        system "samtools view -S -b $base.sam > $base.bam\n";
        print "Creating filtered fastq for $base\n";
        system "samtools bam2fq $base.bam > clean_$base.fastq\n"; 
        $pm->finish;
}

$pm->wait_all_children;
````

* I then parsed and split the files using our standard perl scripts. We have data from 3072 individuals total (3072 fastq files in the Fastq subdirectory).

Next, I aligned the GBS data to a reference genome using `bwa` (version 0.7.17-r1188). This was controlled with SubBwa.sh.

````bash
#!/bin/sh
#SBATCH --time=140:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bwa
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load bwa

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/btreatae_gbs/Fastq

perl BwaFork.pl *fastq
````
````perl
#!/usr/bin/perl
#
# DNA sequence alignment
#

use Parallel::ForkManager;
my $max = 36;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group2/data/btreatae_gbs/genome/ncbi-genomes-2020-03-27/GCA_010883055.1_B_treatae_v1_genomic.fasta";

FILES:
foreach $file (@ARGV){
	$pm->start and next FILES; ## fork
        if ($file =~ m/^([A-Z0-9a-z\-_]+)/){
	        $ind = $1;
    	}
        else {
        	die "Failed to match $file\n";
	}
        system "bwa mem -t 1 -k 15 -r 1.3 -T 30 -r \'\@RG\\tID:Btreatae-"."$ind\\tPL:ILLUMINA\\tLB:Btreatae-"."$ind\\tSM:Btreatae-"."$ind"."\' $genome $file > aln"."$ind".".sam\n";

	$pm->finish;
}

$pm->wait_all_children;
````

* The sam alignments were moved to the Alignments subdirectory. I used `samtools` (version 1.10) to compress, sort and index the sam files, resulting in 3072 \*sorted.bam files. This was controlled with `SubSam2Bam.sh`.

````bash
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=samtools
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
# version 1.10

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/btreatae_gbs/Alignments

perl Sam2BamFork.pl aln*sam

#!/usr/bin/perl
#
# convert sam to bam, then sort and index 
#

use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $sam (@ARGV){
	$pm->start and next FILES; ## fork
	$sam =~ m/^aln([A-Za-z0-9_\-]+\.)sam/ or die "failed to match $sam\n";
	$base = "aln_$1";
	system "samtools view -b -O BAM -o $base"."bam $sam\n";
        system "samtools sort -O BAM -o $base"."sorted.bam $base"."bam\n";
        system "samtools index -b $base"."sorted.bam\n";
        $pm->finish;
}

$pm->wait_all_children;
````

* Variants were then called using `bcftools` (version 1.9). This was done with the script `VarCall.sh`.

````bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=BCF_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load bcftools/1.9

##BCFTOOLs mpileup version 1.9
#Usage: bcftools mpileup [options] in1.bam [in2.bam [...]]

#Input options:
#  -6, --illumina1.3+      quality is in the Illumina-1.3+ encoding
#  -A, --count-orphans     do not discard anomalous read pairs
#  -b, --bam-list FILE     list of input BAM filenames, one per line
#  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)
#  -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
# -d, --max-depth INT     max per-file depth; avoids excessive memory usage [250]
#  -E, --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs
#  -f, --fasta-ref FILE    faidx indexed reference sequence file
#      --no-reference      do not require fasta reference file
#  -G, --read-groups FILE  select or exclude read groups listed in the file
#  -q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
#  -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
#  -r, --regions REG[,...] comma separated list of regions in which pileup is generated
#  -R, --regions-file FILE restrict to regions listed in a file
#      --ignore-RG         ignore RG tags (one BAM = one sample)
#  --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset []
#  --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set
#                                            [UNMAP,SECONDARY,QCFAIL,DUP]
#  -s, --samples LIST      comma separated list of samples to include
#  -S, --samples-file FILE file of samples to include
 # -t, --targets REG[,...] similar to -r but streams rather than index-jumps
 # -T, --targets-file FILE similar to -R but streams rather than index-jumps
#  -x, --ignore-overlaps   disable read-pair overlap detection


#SNP/INDEL genotype likelihoods options:
#  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [20]
#  -F, --gap-frac FLOAT    minimum fraction of gapped reads [0.002]
#  -h, --tandem-qual INT   coefficient for homopolymer errors [100]
#  -I, --skip-indels       do not perform indel calling
#  -L, --max-idepth INT    maximum per-file depth for INDEL calling [250]
#  -m, --min-ireads INT    minimum number gapped reads for indel candidates [1]
#  -o, --open-prob INT     Phred-scaled gap open seq error probability [40]
#  -p, --per-sample-mF     apply -m and -F per-sample for increased sensitivity
#  -P, --platforms STR     comma separated list of platforms for indels [all]


#BCFTOOLs call version 1.9 options used
#v = output variant sites only
#c/m = he original calling method (conflicts with -m) or alternative model for multiallelic and rare-variant calling (conflicts with -c)
#p = variant if P(ref|D)<FLOAT with -c [0.5]
#P =  --prior <float-o mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
#O =  output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' (here it is 'v') 
#o = write output to a file [standard output]

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/btreatae_gbs/Alignments

bcftools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/gompert-group2/data/btreatae_gbs/genome/ncbi-genomes-2020-03-27/GCA_010883055.1_B_treatae_v1_genomic.fasta -q 30 -Q 20 -I -b bams -O b -o btreatae.bcf

## second run b/c of time out, can combine the two vcf files
#bcftools call -v -c -f GQ -p 0.01 -P 0.001 -r CM021346.1,CM021346.1 -O v -o btreatae_pt2.vcf btreatae.bcf 
bcftools call -v -c -f GQ -p 0.01 -P 0.001 -O v -o btreatae.vcf btreatae.bcf
````

The resulting vcf files are in the Variants subdirectory. Because of a job timeout, there are two files (the second = pt2 and has two chromosomes, with some possible overlap for one with the main file). I will combine these prior to filtering (next step). Note that I focused on the large scaffolds/chromosomes and ignored all the tiny contigs/scaffolds.

* I first combined the two vcf files into one vcf file while dropping redundant calls (from one chromosome in both files), this created /uufs/chpc.utah.edu/common/home/gompert-group2/projects/btreatae_exp/variants/btreatae_cm.vcf. Note, from here on out I am in the projects directory.

I tried two settings for variant filtering: (i) 2X coverage and 80% with data, and (ii) 1.5X coverage and 70% with data. In both cases I dropped all of the bias filters. I did this as the current bcftools versions are p-value based and with >3000 samples most p-values are tiny and not useful. Option (i) retained 8238 SNPs, whereas option (ii) retained 23,932. I think the slight decrease in stringency (especially in light of the large number of individuals and mostly allele-frequency based analyses) to triple the number of SNPs is reasonable, so I am going with (ii). The full details of the filter follow.

````perl
#!/usr/bin/perl

use warnings;
use strict;

# this program filters a vcf file based on overall sequence coverage, number of non-reference reads, number of alleles, and reverse orientation reads

# usage vcfFilter.pl infile.vcf
#
# change the marked variables below to adjust settings
#

#### stringency variables, edits as desired
my $minCoverage = 4608; # minimum number of sequences; DP, 1.5x there are 3072
#my $minCoverage = 6144; # minimum number of sequences; DP, 2x there are 3072
my $minAltRds = 10; # minimum number of sequences with the alternative allele; AC
my $notFixed = 1.0; # removes loci fixed for alt; AF
my $mq = 30;
my $miss = 922; # maximum number of individuals with no data, 70%
#my $miss = 614; # maximum number of individuals with no data
##### this set is for whole genome shotgun data
my $d;

my @line;

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z_\-]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> filtered1pt5x70pct_$1") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
	chomp;
	$flag = 1;
	if (m/^\#/){ ## header row, always write
		$flag = 1;
	}
	elsif (m/^CM/){ ## this is a sequence line, you migh need to edit this reg. expr.
		$flag = 1;
		$d = () = (m/\d\/\d:0,0,0:\d+/g); ## for bcftools call
		if ($d >= $miss){
			$flag = 0;
			##print "fail missing : ";
		}
		if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
			$flag = 0;
			#print "fail allele : ";
		}
		@line = split(/\s+/,$_);
		if(length($line[3]) > 1 or length($line[4]) > 1){
			$flag = 0;
			#print "fail INDEL : ";
		}
		m/DP=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 < $minCoverage){
			$flag = 0;
			#print "fail DP : ";
		}
		m/AC1*=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 < $minAltRds){
			$flag = 0;
			#print "fail DP : ";
		}
## bcftools call version
	
		m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error DP4 not found\n";
		if(($1 + $2) < $minAltRds){
			$flag = 0;
		}
		m/AF1*=([0-9\.e\-]+)/ or die "Syntax error, AF not found\n";
		if ($1 == $notFixed){
			$flag = 0;
		#	print "fail AF : ";
		}

		if(m/MQ=([0-9\.]+)/){
			if ($1 < $mq){
				$flag = 0;
#				print "fail MQ : ";
			}
		}
		else{
			$flag = 0;
			print "faile no MQ : ";
		}
		if ($flag == 1){
			$cnt++; ## this is a good SNV
		}
	}
	else{
		print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
		$flag = 0;
	}
	if ($flag == 1){
		print OUT "$_\n";
	}
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
````
* Next, I applied filters for excessive depth and crowded SNPs. Specifically, I removed SNPs fewer then 3 bp apart (all SNPs in such a set) and with higher coverage than the mean + 3 sd. These are both designed to deal with paralogs or repetitive regions.

````bash
perl filterSomeMore.pl filtered1pt5x70pct_btreatae_cm.vcf 
````

This filter left me with 20,446 SNPs in the file morefilter_filtered1pt5x70pct_btreatae_cm.vcf. I converted this to genotype likelihood (gl) format and split out the different populations and treatments. The 80 split files are in /uufs/chpc.utah.edu/common/home/gompert-group2/projects/btreatae_exp/popgen/glfiles.

* Covert to gl with no maf filter

````bash
perl vcf2gl.pl 0.000 morefilter_filtered1pt5x70pct_btreatae_cm.vcf
````
* split populations and treatments
 
 ```bash
 perl splitPops.pl *gl
````

# Population genomic analyses

* I estimated allele frequencies and did some initial analyses of allele frequency changes and correlations in the six experimental populations. Note that Scott pointed out we had lost one of the chromosomes. I am going back to get it and redo what I have here.

I used `popmod` to generate Bayesian allele frequency estimates. This was done from /uufs/chpc.utah.edu/common/home/gompert-group2/projects/btreatae_exp/popgen/.

````bash
#!/bin/sh
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=afest
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load hdf5

cd /uufs/chpc.utah.edu/common/home/gompert-group2/projects/btreatae_exp/popgen/glfiles

perl ../RunAFest.pl *gl

````

````perl
#!/usr/bin/perl
#
# run bayesian afreq est. 
#

use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

foreach $in (@ARGV){
	$pm->start and next; ## fork
	$out = "o_$in";
	$out =~ s/gl/hdf5/ or die "failed sub $out\n";
	system "popmod -i $in -n 10000 -b 1000 -t 3 -o $out\n";
	$pm->finish;
}

$pm->wait_all_children;
````

* I then moved the hdf5 files to the af subdirectory and generated genotype and allele frequency posterior summaries with `estpost`.

* Next I completed some initial summaries in R. The code (`btr.R`) follows, but first a few key summaries (see code for details).

* Here are a few plots to look at:

cor_dp.pdf = summarizes pairwise correlations in allele frequency change between the source (s) and control or treatment for each population. The key take home here is that patterns of change are almost entirely independent across populations, but positively correlated between treatments within populations. This means either that the two treatments evolved in parallel or the source isn't quite right. Want to think more about this (and maybe chat) but one cool possibility would be selection fluctuating between the asexual and sexual generations. Probably a long shot, but would be sweet.

cor_div.pdf = similar, but for divergence between c and t treatments. These are mostly uncorrelated across populations. Thus, we don't see much in terms of consistency in differences between treatments across experimental populations.

btr_dp.pdf = allele frequency change (c or t relative to s) for each population along the main scaffolds=chromosomes(?). Note there are two pages. The horizontal lines denote average change (mean change is low, some stuff changes a lot). You can see the similarity in change for the two treatments often.


btr_div.pdf = allele frequency difference between c and t for each population.

* I have been looking at how similar allele frequencies are for the three (s, t, c) or four (s, t, c and p2) samples associated with each experiment. First, in all cases, all three/four samples are very similar overall, it is just that some are more similar than others. Beyond that, details depend on the population.

So, where does this leave us? If we feel really confident that s is a good representative of the founders of c and t (and thus that similarities between c and t are due to parallel change), then I am inclined to fit a model of s as the ancestors of c and t. If not, I will likely fit that model anyway, as well as one with an unkown ancestor where we just look at differences between c and t. Either way, I also want to try a model for p2 to s for the cases we have it (p2 is a reasonable representative ancestor of s right). That one is maybe a bit of an aside, but could help us think about any evidence for fluctuating selection between generations especially in the context of explaining parallel change in c and t where it occurs.

````R
## get files
pf<-list.files(pattern="p_o")
L<-20446
## read allele freqs.
N<-length(pf)
P<-matrix(NA,nrow=L,ncol=N)
for(j in 1:N){
    px<-read.table(pf[j],sep=",",header=FALSE)
    P[,j]<-px[,1]
    }
    
# get ids
pids<-read.table("pids.txt",header=F)    

## pca
o<-prcomp(t(P),center=TRUE,scale=FALSE)

pdf("bt_pop_pca.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
plot(o$x[,1],o$x[,2],type='n',xlab="PC1 (51.5%)",ylab="PC2 40.6%)",cex.lab=1.5)
text(o$x[,1],o$x[,2],pids[,1],cex=.5)
dev.off()

## 3 geographic clusters
pids[o$x[,1] < -20,1]
# [1] ail    alv_c  alv_s  alv_t  bnc    chr    dia    fmn    gba    hsb   
#[11] hvf    jek    kre    kre_c  kre_p2 kre_s  kre_t  lbc    lvf    ohf   
#[21] oke_c  oke_s  oke_t  paw    pcf    pen    ppf    pry    sap    sri   


pids[o$x[,1] > 20,1]
# [1] alt bat bld bsl del egn ent gld gok hit ict inz lht lop ltx mat mcl msn ogh
#[20] pim plt poc qtz rck rut skt sll smt vdt

pids[o$x[,2] > 30,1]
# [1] abs    abs_c  abs_p2 abs_s  abs_t  ais    dps    dps_c  dps_p2 dps_s 
#[11] dps_t  ibf    llf    llf_c  llf_p2 llf_s  llf_t  och    osv    prk   

## read in groups for experiments
# s = start, c = control, t = treatment,  p2 (if present, prev. gen)
# set to columsn 1, 2, 3 and 4 resepctively
# pop order = abs, alv, dps, kre, llf, oke
grp<-read.table("groups.txt",header=FALSE)
gids<-as.data.frame(grp)
for(i in 1:6){for(j in 1:4){
    gids[i,j]<-as.character(pids[grp[i,j],1])
    }}
gids
#     V1    V2    V3
#1 abs_s abs_c abs_t abs_p2
#2 alv_s alv_c alv_t   <NA>
#3 dps_s dps_c dps_t dps_p2
#4 kre_s kre_c kre_t kre_p2
#5 llf_s llf_c llf_t llf_p2
#6 oke_s oke_c oke_t   <NA>


snps<-read.table("snps.txt",header=FALSE)
lb<-tapply(X=1:L,INDEX=snps[,1],min)
ub<-tapply(X=1:L,INDEX=snps[,1],max)

library(scales)
lgpg<-function(l=lb,u=ub,cs="gray",ss=seq(2,8,2),tt=.85,bb=-0.05){
    for(k in ss){
        polygon(c(l[k],u[k],u[k],l[k]),c(bb,bb,tt,tt),col=cs,border=NA)
    }
}

    
pop<-c("ABS","ALV","DPS","KRE","LLF","OKE")
cm<-1.4;cl<-1.4;ca<-1.2
pdf("btr_dp.pdf",width=9,height=9)
par(mfrow=c(3,2))
par(mar=c(4,5,3,1))
for(i in 1:6){
    x<-grp[i,1];y<-grp[i,2];z<-grp[i,3]
    dpy<-abs(P[,y]-P[,x])
    dpz<-abs(P[,z]-P[,x])
    plot(dpy,xlab="SNP",ylab="AF change",cex.lab=cl,cex.axis=ca,ylim=c(0,.8),type='n')
    lgpg(cs=alpha("gray",.5))
    points(dpy,pch=19,col=alpha("chocolate1",.7))
    abline(h=mean(dpy),lty=2)
    title(main=paste(pop[i]," control",sep=""),cex.main=cm)
    plot(dpz,xlab="SNP",ylab="AF change",cex.lab=cl,cex.axis=ca,ylim=c(0,.8),type='n')
    lgpg(cs=alpha("gray",.5))
    points(dpz,pch=19,col=alpha("cadetblue3",.7))
    abline(h=mean(dpz),lty=2)
    title(main=paste(pop[i]," treatment",sep=""),cex.main=cm)
    }
dev.off()

pdf("btr_div.pdf",width=9,height=9)
par(mfrow=c(3,2))
par(mar=c(4,5,3,1))
for(i in 1:6){
    x<-grp[i,1];y<-grp[i,2];z<-grp[i,3]
    div<-abs(P[,z]-P[,y])
    plot(div,xlab="SNP",ylab="AF divergence",cex.lab=cl,cex.axis=ca,ylim=c(0,.8),type='n')
    lgpg(cs=alpha("gray",.5))
    points(div,pch=19,col=alpha("firebrick1",.7))
    abline(h=mean(div),lty=2)
    title(main=paste(pop[i]," divergence",sep=""),cex.main=cm)

    }
dev.off()

## change and div cor matrix
## simple version, make a table
ctab<-matrix(NA,nrow=6,ncol=6)
myp<-matrix(c(1,2,1,3,2,3,4,1,4,2,4,3),nrow=6,ncol=2,byrow=TRUE)
for(i in 1:6){
    for(j in 1:6){
        ctab[i,j]<-cor(P[,grp[i,myp[j,1]]],P[,grp[i,myp[j,2]]])
}}

ctab
#          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#[1,] 0.9961428 0.9960964 0.9973105 0.9901029 0.9924449 0.9916113
#[2,] 0.9978572 0.9975594 0.9974115        NA        NA        NA
#[3,] 0.9981942 0.9973741 0.9979385 0.9977141 0.9979749 0.9966380
#[4,] 0.9982614 0.9964229 0.9970231 0.9977223 0.9988998 0.9965776
#[5,] 0.9968913 0.9971499 0.9961713 0.9967224 0.9983988 0.9957106
#[6,] 0.9980627 0.9973911 0.9977169        NA        NA        NA

##

div<-matrix(NA,nrow=6,ncol=L)
for(i in 1:6){
    x<-grp[i,1];y<-grp[i,2];z<-grp[i,3]
    div[i,]<-P[,z]-P[,y]
    }

dp<-matrix(NA,nrow=12,ncol=L)
    for(i in 1:6){
    x<-grp[i,1];y<-grp[i,2];z<-grp[i,3]
    dpy<-P[,y]-P[,x]
    dpz<-P[,z]-P[,x]
    j<-(i-1)*2 + 1
    k<-(i-1)*2 + 2
    dp[j,]<-dpy
    dp[k,]<-dpz
    }
    
pcc<-matrix(NA,nrow=18,ncol=L)
    for(i in 1:6){
    x<-grp[i,1];y<-grp[i,2];z<-grp[i,3]
    j<-(i-1)*3 + 1
    k<-(i-1)*3 + 2
    l<-(i-1)*3 + 3
    pcc[j,]<-P[,x]
    pcc[k,]<-P[,y]
    pcc[l,]<-P[,z]
    }    
    
pdf("cor_div.pdf",width=6,height=6)
corrgram(t(div),upper.panel=panel.cor,labels=pop)
dev.off()

pt<-paste(rep(pop,each=2),rep(c("c","t"),6),sep="_")
pdf("cor_dp.pdf",width=6,height=6)
corrgram(t(dp),upper.panel=panel.cor,labels=pt)
dev.off()    

pt<-paste(rep(pop,each=3),rep(c("s","c","t"),6),sep="_")
pdf("cor_p.pdf",width=6,height=6)
corrgram(t(pcc),upper.panel=panel.cor,labels=pt)
dev.off()   
````
