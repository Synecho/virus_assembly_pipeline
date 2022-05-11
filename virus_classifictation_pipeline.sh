#!/bin/sh
#
# Defaults
rDir=0
threads=40
cDir=0
IDX=0
oDir=./
db=0
# Colors
red='\033[0;31m'
green='\033[0;32m'
nc='\033[0m'
# 
while : ; do
    case $1 in
        -r)
            shift
            rDir=$1
            shift
            ;;
        -c)
            shift
            cDir=$1
            shift
            ;;            
        -o)
            shift
            oDir=$1
            shift
            ;;
        -idx)
            shift
            idx=1
            shift
            ;;
        -t)
            shift
            threads=$1
            shift
            ;;
        -minl)
            shift
            minlength=$1
            shift
            ;;
        -co)
            shift
            compl=$1
            shift
            ;;
        -h|--help)
            echo "This script uses Bowtie2 (very-sensitive mode) and samtools to map paired Illumina reads to assembled contigs and outputs .sam, .bam and, if the -idx option is set read abuandance files."
            echo "USAGE:"
            echo "-r or -rDir [PATH]: define path to a directory containing quality quality checked reads from Illumina sequencing (mandatory)"
            echo "-c or -cDir [PATH]: define path to a directory containing corresponding assembled contigs from the reads used for mapping (mandatory)"
            echo "NOTE: file names should be: XZY_R1.fastq XZY_R2.fastq for reads and XZY_contigs.fasta"
            echo "-o [PATH]: Directory in which results will be saved. Default: ./"
            echo "-t [INT]: number of threads allocated for the process. Default: 4"
            echo "-idx: output mapped read abundance as .txt file"
            echo "-completeness [INT]: [MANDATORY] parameter for checkv completeness. Will filter >=[completeness]"
            echo "-minlength [INT]: [MANDATORY] parameter for min viral contig length. Will filter >=[minlength]"
            exit
            ;; 
        *)  
            if [ -z "$1" ]; then break; fi
            ERR=1
            shift
            ;;

    esac
done

################################## Modular Scrip #######################################
echo "Modules of the Virus Pipeline:"
echo "1: VIBRANT"
echo "2: CheckV"
echo "3: ViPTreeGen"
echo "4: Prodigal"
echo "5: Vcontact2"
echo "6: VirHostMatcher-Net"
echo "7: BBDUK repair mate-pairs in .fastq"
echo "8: Bowtie2"
echo "9: DIAMOND Blastp ViralRefSeq"


echo Run module: 
read STEP


################################ Paths to Programms ####################################
checkvdb=/home/benedikth/checkv-db-v1.0
cluster1=/bioinf/home/benedikt.heyerhoff/Resources/vcontact2/bin/cluster_one-1.0.jar
vcontact2_tax_predict=/bioinf/home/benedikt.heyerhoff/Resources/vcontact2_tax_predict.py 
VirHostMatcher=/bioinf/home/benedikt.heyerhoff/Resources/VirHostMatcher-Net/VirHostMatcher-Net.py
viralrefseq=/bioinf/home/benedikt.heyerhoff/Resources/database/viral.3.protein.faa

#hardcoded shit:
compl=90


########################################################################################

############################# Checking & editing variables #############################

if [[ "${rDir}" == 0 ]]; then
    echo "ERROR!"
    echo -e "${red}No read-directory defined. Script aborted.${nc}"
    echo -e "${red}Use -rDir /PATH/TO/READS to tell the script where they are.${nc}"
    exit
fi

if [[ "${cDir}" == 0 ]]; then
    echo "ERROR!"
    echo -e "${red}No contig directory defined. Script aborted.${nc}"
    echo -e "${red}Use -cDir /PATH/TO/READS to tell the script where they are.${nc}"
    exit
fi

if [[ ! -d $oDir ]]; then 
    echo "ERROR"
    echo -e "${red}Output directory does not exist. Set with -o or -oDIR /PATH/TO/OUTPUT/DIR${nc}"
    exit
fi

if [[ "${rDir: -1}" == "/" ]]; then
    rDir=${rDir::-1}
fi

if [[ "${oDir: -1}" == "/" ]]; then
    oDir=${oDir::-1}
fi

########################## Create list of metagenome folders ###########################
#(cd $cDir && ls *.f*) | awk 'BEGIN{FS=OFS="."}{NF--; print}'> $oDir/infiles.txt

####################################### VIBRANT ########################################
if [[ "$STEP" == 1 ]]; then
    echo -e "Next step: VIBRANT"
    for s in $(cat $oDir/infiles.txt);do
        echo -ne "\e[38;5;86mRUNNING VIBRANT WITH "${s}" CONTIGS${nc}"
        mkdir -p $oDir/vibrant/${s}
        python3 VIBRANT_run.py \
        -i $cDir/${s}.f* \
        -t $threads \
        -folder $oDir/vibrant/${s}
        echo -ne "\e[48;5;86mDONE WITH"${S}"${nc}";
    done
fi

####################################### CheckV ########################################
if [[ "$STEP" == 2 ]]; then
    echo -e "${green} Next step: CheckV ${nc}"
    echo ""
    echo ""
    for s in $(cat $oDir/infiles.txt);do
        echo -e "${green} CheckV started ${nc}"
        echo ""
        echo ""
        echo ""
        echo -e "${green}current metagenome name: "$s" ${nc}"
        #checkv end_to_end $oDir/vibrant/${s}/VIBRANT_phages_final.contigs/final.contigs.phages_combined.fna $oDir/checkv/all_metagenome/${s} -t $threads -d $checkvdb; #using full checkv pipeline
        
        #export variables for R script
        export iter=$s
        export oDir=$oDir
        export currentDir=$oDir/checkv/all_metagenome/${s}
        export compl=$compl
        export minlength=$minlength

        ######debugging messages######
        echo -e "EXPORT PARAMETERS:"
        echo -e "${green}#####"
        echo -e "1 (iteration): "$iter""
        echo -e "2 (output dir): "$oDir""
        echo -e "3 (current dir): "$currentDir""
        echo -e "4 (completeness): "$compl""
        echo -e "5 (min length): "$minlength""
        echo -e "#####${nc}"
        ######

        #Directories for filtered results from R script
        mkdir -p $oDir/checkv/q_filtered_contigs
        mkdir -p $oDir/checkv/length_filtered_contigs
        mkdir -p $oDir/checkv/high_quality_contigs
        mkdir -p $oDir/checkv/med_quality_contigs
        mkdir -p $oDir/checkv/low_quality_contigs
        mkdir -p $oDir/checkv/high_score_contigs
        mkdir -p $oDir/checkv/qc_tables
        mkdir -p $oDir/checkv/complete_phages_contigs
        mkdir -p $oDir/checkv/VirHostMatcher/${s};
        
        echo -e "${green} Quality Filtering of Phages from "${s}" Metagenome ${nc}"
        #source R script
        Rscript /bioinf/home/benedikt.heyerhoff/checkv_output_score_filtering.R
    done
fi

#################################### ViPTreeGen #######################################
if [[ "$STEP" == 3 ]]; then
    echo -e "${green} Next step: Merging length filtered contigs contigs into one file ${nc}"

    cd $oDir/checkv/length_filtered_contigs/
    mkdir -p $oDir/ViPTree
    cp *.fna $oDir/ViPTree
    cd $oDir/ViPTree

    for f in *.fna; do
        sed -i '' -e "s/^>/>${f%.fna}_/g" "${f}"; #rename fasta headers to include metagenome name
    done

    cat *.fna > $oDir/ViPTree/input/length_filtered_contigs.fna #concatenate all fastas into one 
    rm $oDir/ViPTree/*.fna
    cd $oDir

    echo -e "${green} Next step: ViPTree ${nc}"
    ViPTreeGen \
    $oDir/ViPTree/length_filtered_contigs.fna \
    $oDir/ViPTree/output \
    --notree \
    --2D $oDir/ViPTree/length_filtered_contigs.fna  \
    --ncpus $threads
fi

###################################### Prodigal ########################################
if [[ "$STEP" == 4 ]]; then
    echo -e "${green}Starting Prodigal...${nc}..."

    for s in $(cat $oDir/infiles.txt);do
        echo -e "${green} Predicting CDS of "${s}" with ${red} Prodigal ${nc}..." 
        
        mkdir -p $oDir/prodigal/
                
        prodigal -i $oDir/checkv/length_filtered_contigs/${s}.fna \
        -o $oDir/prodigal/${s}_proteins \
        -a $oDir/prodigal/${s}_proteins.faa \
        -p meta 
    done
fi
###################################### Vcontact2 ########################################
if [[ "$STEP" == 5 ]]; then
    for s in $(cat $oDir/infiles.txt);do
            echo -e "\e[44mStarting gene2genome step in ${s} Metagenome${nc}"
            mkdir -p $oDir/gene2genome/${s}
            vcontact2_gene2genome -p $oDir/prodigal/proteins/${s}_proteins.faa -o $oDir/gene2genome/${s}/${s}_viral_genomes_g2g.csv -s 'Prodigal-FAA'
            echo -e "\e[42mFinished gene2genome step ${nc}"
    done

    echo -e "${green} Next step: Vcontact2 ${nc}"
    for s in $(cat $oDir/infiles.txt);do
        mkdir -p $oDir/vcontact2/${s}
        vcontact2 --raw-proteins $oDir/prodigal/proteins/${s}_proteins.faa \
        --rel-mode "Diamond" \
        --proteins-fp $oDir/gene2genome/${s}/${s}_viral_genomes_g2g.csv \
        --db 'ProkaryoticViralRefSeq201-Merged' \
        --c1-bin /bioinf/home/benedikt.heyerhoff/Resources/vcontact2/bin/cluster_one-1.0.jar \
        -t $threads \
        --output-dir $oDir/vcontact2/${s}
    done

    echo -e "${green} Assigning taxonomy to Vcontact2 output ${nc}"
    for s in $(cat $oDir/infiles.txt);do
        mkdir -p oDir/viral_taxonomy/${s}
        $vcontact2_tax_predict \
        -i $oDir/vcontact2/${s}/genome_by_genome_overview.csv \
        -o oDir/viral_taxonomy/${s}
    done
fi

################################## VirHostMatcher-Net ###################################
if [[ "$STEP" == 6 ]]; then
    mkdir -p $oDir/VirHostMatcher/temp #tempdir for VirHostMatcher-Net
    for s in $(cat $oDir/infiles.txt);do
        python3 $VirHostMatcher \
        -q $oDir/checkv/VirHostMatcher/${s} \
        --short-contig \
        -o $oDir/VirHostMatcher/${s} \
        -n 1 \
        -t $threads \
        -l /bioinf/home/benedikt.heyerhoff/Resources/VirHostMatcher-Net/genome_list/marine_host_list.txt \
        -i $oDir/VirHostMatcher/temp \
        -d /bioinf/home/benedikt.heyerhoff/Resources/VirHostMatcher-Net/data;        
    done
fi
##################################### Fastq Repairing ##################################
if [[ "$STEP" == 7 ]]; then
    echo -e "${green} Next step: repair.sh ${nc}"
    mkdir -p $oDir/repaired/
    for s in $(cat $oDir/repair2.txt);do
        repair.sh in=$rDir/${s}_1.fastq \
        in2=$rDir/${s}_2.fastq \
        out=$oDir/repaired/${s}_1.fastq \
        out2=$oDir/repaired/${s}_2.fastq \
        outs=$oDir/repaired/${s}_singletons.fastq \
        repair;
    done
fi
####################################### Contig Mapping ##################################
if [[ "$STEP" == 8 ]]; then
    echo -e "${green} Next step: Bowtie2 ${nc}"

    contigs=$( cd $oDir/checkv/length_filtered_contigs && ls *.f* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | wc -l
    reads=$( cd $rDir && ls *.fa* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d | wc -l

    mkdir -p $oDir/bam
    mkdir -p $oDir/sam
    mkdir -p $oDir/temp
    ls $rDir

    #( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/temp/files.txt

    for s in $(cat $oDir/remap.txt);do
        echo -e "Bulding bowtie2 database from ${green}"${s}"_contigs.fasta${nc} ..."
        bowtie2-build $oDir/checkv/length_filtered_contigs/${s}.fna $oDir/temp/Bowtie2.db
        DB=$oDir/temp/Bowtie2.db

        echo
        echo -e "Mapping reads from ${green}"$s"_R1.fastq${nc} and ${green}"$s"_R2.fastq${nc} to database..." 
        bowtie2 --very-sensitive-local \
        -x $DB \
        -1 $rDir/${s}_*1.fastq \
        -2 $rDir/${s}_*2.fastq \
        -p $threads \
        -S $oDir/temp/${s}.sam

        echo -e "Reformating ${green}"$s".sam${nc} to ${green}"$s".bam${nc} ..."
        samtools view -@ $threads -b -S $oDir/temp/${s}.sam > $oDir/bam/${s}.bam
        echo -e "Deleting ${green}"$s".sam${nc}..."
        rm $oDir/temp/${s}.sam
        echo -e "Sorting ${green}"$s".bam${nc}..."
        samtools sort -@ $threads $oDir/bam/${s}.bam > $oDir/bam/${s}.sorted.bam 
        rm $oDir/bam/${s}.bam

        echo -e "Creating, sorting and indexing ${green}"$s".bam${nc} file..."
        mkdir -p $oDir/map 
        samtools index -@ $threads $oDir/bam/${s}.sorted.bam 
        samtools idxstats $oDir/bam/${s}.sorted.bam  > $oDir/map/${s}.mapped.txt
        rm $oDir/bam/${s}.sorted.bam 


        echo ""
        echo -e "Finished mapping of ${green}"$s"${nc}!";
    done 
    rm -rf $oDir/temp/
fi

####################################### ViralRefSeq Blast ##################################
if [[ "$STEP" == 9 ]]; then
    echo -e "${green} Next step: Blasting length filtered contigs against ViralRefSeq ${nc}"
    echo ""
    echo -e "${green} Creating Diamond BlastN DB ${nc}..." 
    echo ""
    mkdir -p $oDir/blastP
    diamond makedb --in $viralrefseq -d $oDir/blastP/viralrefseq_diamond_db.dmnd
    echo ""
    echo -e "${green} Diamond BlastN DB done! ${nc}..." 

    for s in $(cat $oDir/infiles.txt ); do
        echo -e "blasting ${green}"${s}"${nc} against DIAMOND database..." 
        diamond blastp --very-sensitive \
        --db $oDir/blastP/viralrefseq_diamond_db.dmnd \
        --query $oDir/prodigal/${s}/proteins/${s}_proteins.faa \
        --outfmt "6" \
        --matrix "BLOSUM45" \
        --evalue 0.0001 \
        --min-score 50 \
        --out $oDir/blastP/${s}_blastn_viralrefseq.txt;
    done
fi

####################################################################################################
echo ""
echo -e "${green} Done! ${nc}"