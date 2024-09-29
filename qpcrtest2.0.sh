#bin!/bin/bash
#
#PBS -N qpcr.sh
#PBS -q workq
#PBS -o fedback.log
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -V

#CDS file must be named SpeciesName.cds.fa   e.g: Bombyx_mori.cds.fa
#genome file must be named SpeciesName.fa   e.g: Bombyx_mori.fa

cd $PBS_O_WORKDIR #(All plant genomes are placed under plant folder 1, each plant genome is placed under folder 2 named after its own name, software scripts are placed under folder 1, and each run script is placed under folder 2)

#sh qprimerdb.sh Favourites class1 class2 SpeciesName.cds.fa
Favourites=Rest
class1=de 
class2=mo
cdsseq=de_mo.cds.fa
SpeciesName=`basename -s .cds.fa $cdsseq`

# Load software
module load bbmap/38.94 
module load mpprimer/1.4 
module load epcr/2.3.12 
module load primer3/2.4.0 
module load mfeprimer/2.0.0 
module load blat/37x1 
module load python/2.7.18

# The first section used for shortening gene name. You'd better don't use symbol "_" in the gene name. In the following analysis, We use the model plant ‘Arabidopsis thaliana’ as an example to demonstrate primer design pipeline.(only preserve GeneID)

perl -pe 's/\_/./g' $cdsseq | awk -F " " '/>/{$0=$1}1' > $SpeciesName.shortname.fa;

# Step 1. split the cDNA sequences into small fragment with 300-bp windows and 50-bp step.

shred.sh in=$SpeciesName.shortname.fa out=$SpeciesName.shortname.frag.fa length=300 overlap=250 minlength=200;

echo "Step 1 was done"

# Step 2. Create primer3 input file for large fasta sequences; split each transcriptome into smaller fasta files, each with 100000 splited fasta sequences

CreateMPprimerInput.py -i $SpeciesName.shortname.frag.fa -o $SpeciesName.shortname.p3in;

echo "Step 2 was done"

# Step 3. Design primers for each p3in file(the most time-consuming，#make)

primer3_core -default_version=1 < $SpeciesName.shortname.p3in > $SpeciesName.shortname.p3out;

echo "Step 3 was done"

# Step 4. Format output file to e-PCR file using perl script extr_primer3new.pl, and $SpeciesName.shortname.uniepcr will be used for e-PCR program.(#!a[$5]++ remove duplications)
perl ../scripts/extr_primer3new.pl -p $SpeciesName.shortname.p3out -l $SpeciesName.shortname.list -s $SpeciesName.shortname.blast -e $SpeciesName.shortname.epcr;
awk '{print $0"\t"$2$3}' $SpeciesName.shortname.epcr | awk  '!a[$5]++' | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $SpeciesName.shortname.uniepcr;

echo "Step 4 was done"

# Step 5. Check the specificity of all the primers, and filter out non-specific primers based on e-PCR results(#sort -u remove duplications，#make)
e-PCR -w9 -f1 -n1 -g1 -m100 $SpeciesName.shortname.uniepcr D=50-1500 $SpeciesName.shortname.fa N=1 G=1 T=3 O=$SpeciesName.shortname.STS;
awk -F "_" '{print $1"\t"$2"\t"$3"\t"$4}' $SpeciesName.shortname.STS | awk '{if($8!=0 || $9!=0) print $2"_"$3}'|sort -u > $SpeciesName.uniSTS.remove.id.tmp;
awk -F "_" '{print $1"\t"$2"\t"$3"\t"$4}' $SpeciesName.shortname.STS | awk '{if($1==$2 && $8==0 && $9==0) print $2"_"$3}' > $SpeciesName.uniSTS.id.tmp;

if [ -s $SpeciesName.uniSTS.remove.id.tmp ]
then
awk 'NR==FNR{a[$1]=1}NR!=FNR{if(a[$1]!=1) print $0}'  $SpeciesName.uniSTS.remove.id.tmp $SpeciesName.uniSTS.id.tmp > $SpeciesName.uniSTS.id;
else
cat $SpeciesName.uniSTS.id.tmp > $SpeciesName.uniSTS.id;
fi

#extract results for MFEprimer analysis
grep -F -f $SpeciesName.uniSTS.id $SpeciesName.shortname.list > $SpeciesName.uniSTS.unilist;
grep -F -f $SpeciesName.uniSTS.id $SpeciesName.shortname.blast | awk -F ">" '{print $2}' > $SpeciesName.uniSTS.uniblast.id;
#chmod +x ../faSomeRecords
../scripts/faSomeRecords $SpeciesName.shortname.blast $SpeciesName.uniSTS.uniblast.id $SpeciesName.uniSTS.uniblast;

echo "Step 5 was done"

# Step 6. Prepare primer files for each gene, and move them into individual fold.
mkdir $SpeciesName;

grep "^>" $SpeciesName.shortname.fa | awk -F '>' '{print $2}' > $SpeciesName.shortname.name;

for file1 in $(cat $SpeciesName.shortname.name); 
do
grep $file1 $SpeciesName.shortname.blast | awk -F ">" '{print$2}' > $file1.blast.in;
../scripts/faSomeRecords $SpeciesName.shortname.blast $file1.blast.in $SpeciesName.$file1.blast.fa;
rm -rf $file1.blast.in;
mv $SpeciesName.*.blast.fa ./$SpeciesName;
done

for file2 in $(cat $SpeciesName.shortname.name);
do
perl ../scripts/fastaDeal.pl -cuts 2 ./$SpeciesName/$SpeciesName.$file2.blast.fa --outdir ./$SpeciesName;
done;

rm -rf ./$SpeciesName/*.blast.fa;

echo "Step 6 was done"

# Step 7. Index original cDNA fasts file. Then, check the specificity of primers using MFPrimer

IndexDb.sh $SpeciesName.shortname.fa

echo "Step 7 was done"

# Step8. Check the specificity of primers using MFPrimer(time-consuming 700samples need 40 mins)

#cd ./$SpeciesName;
mkdir -p ./$SpeciesName/primers;

for file3 in $(cat $SpeciesName.shortname.name)
do 
    #cd ./$SpeciesName.$file3.blast.fa.cut
    echo "Processing $file3..."
    if ls ./$SpeciesName/$SpeciesName.$file3.blast.fa.cut | grep "blast.fa"
    then
        for file4 in $(ls ./$SpeciesName/$SpeciesName.$file3.blast.fa.cut)
        do
            MFEprimer.py --size_stop=1500 --size_start=50 --tab --ppc=10 \
		    -i ./$SpeciesName/$SpeciesName.$file3.blast.fa.cut/$file4 \
		    -d $SpeciesName.shortname.fa \
		    -o ./$SpeciesName/$SpeciesName.$file3.blast.fa.cut/$file4.txt

            echo "-- $file4 is done"
        done
        cat ./$SpeciesName/$SpeciesName.$file3.blast.fa.cut/*.blast.*.txt > ./$SpeciesName/primers/$file3.MFP
        #mv $file3.MFP ../primers
        rm -rf ./$SpeciesName/$SpeciesName.$file3.blast.fa.cut/*.txt;
        #cd ../
        echo "Done"
    else
        echo "The directory is empty, will ignore it"
        #cd ../
    fi
done

#cd ./primers;
cat ./$SpeciesName/primers/*.MFP > $SpeciesName.primers;
#mv $SpeciesName.primers ../../;
rm -rf ./$SpeciesName/primers/*.MFP

echo "Step 8 was done"

# Step9 format the result file to obtain gene-specific primers. l1: PPC<10; l2: FpDg <-9 && RpDG <-9; l3: FpDg <-11 && RpDG <-11
# All candidate primers
#cd ../../;
awk '$1=="AmpID" {print a} {a=$0}' $SpeciesName.primers | awk '$1 =="1" {print $0}' > $SpeciesName.primers.l1;
awk '{if($1=="AmpID" || $10 < -9) print $0}' $SpeciesName.primers | awk '{if($1=="AmpID" || $11 < -9)  print $0}' > $SpeciesName.primers.l2;
awk '{if($1=="AmpID" || $10 < -11) print $0}' $SpeciesName.primers | awk '{if($1=="AmpID" || $11 < -11)  print $0}' > $SpeciesName.primers.l3;

# obtain the unique l1 primers 
awk '{print $0"\t"$10+$11}' $SpeciesName.primers.l1 | sort -k 4,4 -k 15,15n -k 12,12nr | awk '!a[$4]++' > $SpeciesName.primers.unil1;
cut -f 4 $SpeciesName.primers.unil1 > $SpeciesName.primers.unil1.id;
# Uniq the l2 primers
awk '$1=="AmpID" {print a} {a=$0}' $SpeciesName.primers.l2 | awk '$1 =="1" {print $0}' | awk '{print $0"\t"$10+$11}' | sort -k 4,4 -k 15,15n -k 12,12nr | awk  '!a[$4]++' > $SpeciesName.primers.unil2; 
cut -f 4 $SpeciesName.primers.unil2 > $SpeciesName.primers.unil2.id;
# Uniq the l3 primers
awk '$1=="AmpID" {print a} {a=$0}' $SpeciesName.primers.l3 | awk '$1 =="1" {print $0}' | awk '{print $0"\t"$10+$11}' | sort -k 4,4 -k 15,15n -k 12,12nr | awk  '!a[$4]++' > $SpeciesName.primers.unil3; 
cut -f 4 $SpeciesName.primers.unil3 > $SpeciesName.primers.unil3.id;

echo "Step 9 was done"

# Step10 Extract uniq primers and primer parameters
# Extract l1 uniq primers and primer parameters
awk '{print $2}' $SpeciesName.primers.unil1 | sed 's/^L_//g' > $SpeciesName.primers.unil1.primer.id;
grep -F -f $SpeciesName.primers.unil1.primer.id $SpeciesName.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > $SpeciesName.primers.unil1.epcr.sorted;
sort -k 4,4 $SpeciesName.primers.unil1 | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > $SpeciesName.primers.unil1.sorted;
join -1 1 -2 1 -a 1 $SpeciesName.primers.unil1.sorted $SpeciesName.primers.unil1.epcr.sorted | perl -pe 's/ /\t/g' > $SpeciesName.primers.unil1.primers.tmp;
awk '{print $3"\tL_"$1"\t"$2"\t"$15"\t"$16"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' $SpeciesName.primers.unil1.primers.tmp > $SpeciesName.unil1.primers;
# Extract l2 uniq primer and primer parameters
awk '{print $2}' $SpeciesName.primers.unil2 | sed 's/^L_//g' > $SpeciesName.primers.unil2.primer.id;
grep -F -f $SpeciesName.primers.unil2.primer.id $SpeciesName.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > $SpeciesName.primers.unil2.epcr.sorted; 
sort -k 4,4 $SpeciesName.primers.unil2 | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > $SpeciesName.primers.unil2.sorted; 
join  -1 1 -2 1 -a 1 -o auto  $SpeciesName.primers.unil2.sorted $SpeciesName.primers.unil2.epcr.sorted | perl -pe 's/ /\t/g' > $SpeciesName.primers.unil2.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$15"\t"$16"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' $SpeciesName.primers.unil2.primers.tmp > $SpeciesName.unil2.primers;
# Extract l3 uniq primer and primer parameters
awk '{print $2}' $SpeciesName.primers.unil3 | sed 's/^L_//g' > $SpeciesName.primers.unil3.primer.id;  
grep -F -f $SpeciesName.primers.unil3.primer.id $SpeciesName.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > $SpeciesName.primers.unil3.epcr.sorted; 
sort -k 4,4 $SpeciesName.primers.unil3 | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > $SpeciesName.primers.unil3.sorted; 
join -1 1 -2 1 -a 1 -o auto  $SpeciesName.primers.unil3.sorted $SpeciesName.primers.unil3.epcr.sorted | perl -pe 's/ /\t/g' >  $SpeciesName.primers.unil3.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$15"\t"$16"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' $SpeciesName.primers.unil3.primers.tmp > $SpeciesName.unil3.primers;

# All candidate primers
awk '$1=="AmpID" {print a} {a=$0}' $SpeciesName.primers | awk '$1 =="1" {print $0}'  > $SpeciesName.primers.l1.all;
awk '$1=="AmpID" {print a} {a=$0}' $SpeciesName.primers.l2 | awk '$1 =="1" {print $0}'  > $SpeciesName.primers.l2.all; 
awk '$1=="AmpID" {print a} {a=$0}' $SpeciesName.primers.l3 | awk '$1 =="1" {print $0}' > $SpeciesName.primers.l3.all; 

#extract l1, l2 and l3 primers and primer parameters.
awk '{print $2}' $SpeciesName.primers.l1.all | sed 's/^L_//g' > $SpeciesName.primers.l1.all.primer.id;
grep -F -f $SpeciesName.primers.l1.all.primer.id $SpeciesName.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > $SpeciesName.primers.l1.all.epcr.sorted;
sort -k 4,4 $SpeciesName.primers.l1.all | cut -f 2-15 | sed 's/^L_//g' | sort -k 1,1 > $SpeciesName.primers.l1.all.sorted; 
join -1 1 -2 1 -a 1 -o auto  $SpeciesName.primers.l1.all.sorted $SpeciesName.primers.l1.all.epcr.sorted | perl -pe 's/ /\t/g' >  $SpeciesName.primers.l1.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' $SpeciesName.primers.l1.primers.tmp > $SpeciesName.l1.all.primers;
awk '{print $2}' $SpeciesName.primers.l2.all | sed 's/^L_//g' > $SpeciesName.primers.l2.all.primer.id;  
grep -F -f $SpeciesName.primers.l2.all.primer.id $SpeciesName.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > $SpeciesName.primers.l2.all.epcr.sorted; 
sort -k 4,4 $SpeciesName.primers.l2.all | cut -f 2-15 | sed 's/^L_//g'| sort -k 1,1 > $SpeciesName.primers.l2.all.sorted; 
join -1 1 -2 1 -a 1 -o auto  $SpeciesName.primers.l2.all.sorted $SpeciesName.primers.l2.all.epcr.sorted | perl -pe 's/ /\t/g' >  $SpeciesName.primers.l2.all.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' $SpeciesName.primers.l2.all.primers.tmp > $SpeciesName.l2.all.primers;
awk '{print $2}' $SpeciesName.primers.l3.all | sed 's/^L_//g' > $SpeciesName.primers.l3.all.primer.id;  
grep -F -f $SpeciesName.primers.l3.all.primer.id $SpeciesName.shortname.epcr | cut -f 1,2,3 | sort -k 1,1 > $SpeciesName.primers.l3.all.epcr.sorted; 
sort -k 4,4 $SpeciesName.primers.l3.all | cut -f 2-15 | sed 's/^L_//g'| sort -k 1,1 > $SpeciesName.primers.l3.all.sorted; 
join -1 1 -2 1 -a 1 -o auto  $SpeciesName.primers.l3.all.sorted $SpeciesName.primers.l3.all.epcr.sorted | perl -pe 's/ /\t/g' >  $SpeciesName.primers.l3.all.primers.tmp; 
awk '{print $3"\tL_"$1"\t"$2"\t"$14"\t"$15"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' $SpeciesName.primers.l3.all.primers.tmp > $SpeciesName.l3.all.primers; 

echo "Step 10 was done"

# Step 11 Obtain the best gene-specific primers, which will be put into qPrimerDB.
mkdir best.primers
mv $SpeciesName.unil1.primers ./best.primers
mv $SpeciesName.unil2.primers ./best.primers
mv $SpeciesName.unil3.primers ./best.primers
cd ./best.primers

mv $SpeciesName.unil1.primers $SpeciesName.l1;
awk '{print $1}' $SpeciesName.l1 > $SpeciesName.l1.id;
awk 'NR==FNR{a[$2]=1}NR!=FNR{if(a[$1]!=1) print $0}' $SpeciesName.l1.id $SpeciesName.unil2.primers > $SpeciesName.l2; 
awk '{print $1}' $SpeciesName.l2 > $SpeciesName.l2.id; cat $SpeciesName.l1.id $SpeciesName.l2.id > $SpeciesName.l12.id; 
awk 'NR==FNR{a[$2]=1}NR!=FNR{if(a[$1]!=1) print $0}' $SpeciesName.l12.id $SpeciesName.unil3.primers > $SpeciesName.l3; 
#rename "s/.l1/%%l1/" * 
#rename "s/.l2/%%l2/" * 
#rename "s/.l3/%%l3/" * 
rename ".l1" "%%l1" *.l1
rename ".l2" "%%l2" *.l2
rename ".l3" "%%l3" *.l3

awk '{print FILENAME"\t"$0}' $SpeciesName%%l1 > $SpeciesName.l1.tmp ; 
awk '{print FILENAME"\t"$0}' $SpeciesName%%l2 > $SpeciesName.l2.tmp; 
awk '{print FILENAME"\t"$0}' $SpeciesName%%l3 > $SpeciesName.l3.tmp;  
cat $SpeciesName.l1.tmp $SpeciesName.l2.tmp $SpeciesName.l3.tmp | awk -F "%%" '{print $1"\t"$2}' > $SpeciesName.4db.tmp;

# Obtain amplification parameters using blat against genome sequence.
cut -f 3,17 $SpeciesName.4db.tmp | awk '{print ">"$1"\n"$2}' > $SpeciesName.4db.best.tmp; 

blat ../$SpeciesName.fa $SpeciesName.4db.best.tmp -noHead $SpeciesName.4db.best.blat.tmp; 

sort -k 10,10 -k 1,1nr $SpeciesName.4db.best.blat.tmp | awk '!a[$10]++' > $SpeciesName.4db.best.blat.uniq.tmp;
sort -k 3,3 $SpeciesName.4db.tmp | awk '!a[$3]++' >  $SpeciesName.4db.best.sorted.tmp;
join -1 3 -2 10 -a 1  $SpeciesName.4db.best.sorted.tmp  $SpeciesName.4db.best.blat.uniq.tmp | perl -pe 's/ /\t/g' > $SpeciesName.4db.best.records.tmp;
awk '{print $0"\tbest"}' $SpeciesName.4db.best.records.tmp > $SpeciesName.best.4db;

awk '{if ($0~">") {sub(">", "", $0); printf "\n%s\t",$0} else {printf "%s",$0}}END{print "\n"}' ../$SpeciesName.shortname.fa | sed '1d' > ../$SpeciesName.cds

awk -v A=$Favourites -v B=$class1 -v C=$class2 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2}NR!=FNR{num = sprintf("%06d", FNR); str = $2; split(str, arr, "_"); subs1 = substr(arr[1], 1, 1); sub(arr[1], subs1, str); sub("_", ".", str); print A,B,C,str"."num"v1",$0,a[$1],length($6),length($7)}' ../$SpeciesName.cds $SpeciesName.best.4db > $SpeciesName.best.5db

awk -F "\t" 'BEGIN {print "Favourites\tClass1\tClass2\tprimerID\tGeneID\tOrganisam\tpLevel\tFpID\tRpID\tFprimer\tRprimer\tPPC\tAmpSize\tAmpGC\tFpTm\tRpTm\tFpDg\tRpDg\tBindingStart\tBindingStop\tAmpSeq\tmatches\tmisMatches\trepMatches\tNcount\tQgapcount\tQgapbases\tTgapcount\tTgapbases\tstrand\tqSize\tqStart\tqEnd\ttName\ttSize\ttStart\ttEnd\tNumExonCorss\tblockSizes\tqStarts\ttStarts\tCategory\tcdnaSeq\tFpLength\tRpLength"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36"\t"$37"\t"$38"\t"$39"\t"$40"\t"$41"\t"$42"\t"$43"\t"$44"\t"$45}'  $SpeciesName.best.5db > $SpeciesName.best.primers

echo "Step 11 was done"

#Step12 Obtain all the candidate gene-specific primers.
cd ..
mkdir all.primers
mv $SpeciesName.l1.all.primers ./all.primers
mv $SpeciesName.l2.all.primers ./all.primers
mv $SpeciesName.l3.all.primers ./all.primers

cd ./all.primers
rename ".l1.all.primers" "%%l1" *.l1.all.primers  
rename ".l2.all.primers" "%%l2" *.l2.all.primers
rename ".l3.all.primers" "%%l3" *.l3.all.primers

awk '{print FILENAME"\t"$0}' $SpeciesName%%l1 > $SpeciesName.l1.all.tmp; 
awk '{print FILENAME"\t"$0}' $SpeciesName%%l2 > $SpeciesName.l2.all.tmp; 
awk '{print FILENAME"\t"$0}' $SpeciesName%%l3 > $SpeciesName.l3.all.tmp;  
cat $SpeciesName.l1.all.tmp $SpeciesName.l2.all.tmp $SpeciesName.l3.all.tmp | awk -F "%%" '{print $1"\t"$2}' > $SpeciesName.4db.all.tmp;
sort -k 3,3 -k 2,2 $SpeciesName.4db.all.tmp | awk '!a[$4]++' > $SpeciesName.4db.all.uniq.tmp;

# Obtain amplification parameters using blat against genome sequence.
cut -f 4,17 $SpeciesName.4db.all.uniq.tmp | sed 's/^L_//g' | awk '{print ">"$1"\n"$2}' > $SpeciesName.4db.all.uniq.4blat.tmp; 
blat ../$SpeciesName.fa $SpeciesName.4db.all.uniq.4blat.tmp -noHead $SpeciesName.4db.all.blat.tmp; 
sort -k 10,10 -k 1,1nr $SpeciesName.4db.all.blat.tmp | awk '!a[$10]++' > $SpeciesName.4db.all.blat.uniq.tmp;
sort -k4,4 $SpeciesName.4db.all.uniq.tmp > $SpeciesName.4db.all.uniq.sorted.tmp;
cut -f 4 $SpeciesName.4db.all.uniq.sorted.tmp | sed 's/^L_//g' > $SpeciesName.4db.all.uniq.sorted.tmp.id;
paste  $SpeciesName.4db.all.uniq.sorted.tmp.id $SpeciesName.4db.all.uniq.sorted.tmp > $SpeciesName.4db.all.uniq.sorted.id.tmp
join -1 1 -2 10 $SpeciesName.4db.all.uniq.sorted.id.tmp  $SpeciesName.4db.all.blat.uniq.tmp | perl -pe 's/ /\t/g' > $SpeciesName.4db.all.records.tmp;
awk '{print $0"\tall"}' $SpeciesName.4db.all.records.tmp > $SpeciesName.all.4db;

awk -F "\t" 'BEGIN {print "primerID\tOrganisam\tpLevel\tGeneID\tFpID\tRpID\tFprimer\tRprimer\tPPC\tAmpSize\tAmpGC\tFpTm\tRpTm\tFpDg\tRpDg\tBindingStart\tBindingStop\tAmpSeq\tmatches\tmisMatches\trepMatches\tNcount\tQgapcount\tQgapbases\tTgapcount\tTgapbases\tstrand\tqSize\tqStart\tqEnd\ttName\ttSize\ttStart\ttEnd\tNumExonCorss\tblockSizes\tqStarts\ttStarts\tCategory"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36"\t"$37"\t"$38"\t"$39 }' $SpeciesName.all.4db > $SpeciesName.all.primers

echo "all steps were done"

# remove all tmp file
cd ..

rm -f $SpeciesName.primers*
rm -f $SpeciesName.shortname*
rm -f $SpeciesName.uniSTS*

exit 0

