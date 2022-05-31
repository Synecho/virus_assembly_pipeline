#!/bin/sh
#
# Defaults
rDir=0
threads=40
cDir=0
IDX=0
oDir=./
db=0
################################################################################################################################
#                                                          Colors                                                              #
################################################################################################################################
red='\033[0;31m'
green='\033[0;32m'
blue='\e[1;34m'
yellow='\e[1;33m'
nc='\033[0m'
BgGreen='\e[42m'


################################################################################################################################
#                                                    Analysis defaults                                                         #
################################################################################################################################

#options

groups=dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae # viral groups for virsorter2
minlength=1500                                  # default minimum viral contig length in bp
id=95                                           # default minimum viral squence identity for clustering
mem=500                                         # memory for SPAdes assembly in GB
threads=80                                      # Num threads to be used                                  
################################################################################################################################
#                                                          Options                                                             #
################################################################################################################################ 
while : ; do
    case $1 in
        -r)
shift
rDir=$1
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
-id)
shift
id=1
shift
;;
-mem)
shift
mem=$1
shift
;;  
-groups)
shift
groups=$1
shift
;;                                     
-h|--help)
echo "This script uses Bowtie2 (very-sensitive mode) and samtools to map paired Illumina reads to assembled contigs and outputs .sam, .bam and, if the -idx option is set read abuandance files."
echo "USAGE:"
echo "-r or -rDir [PATH]: define path to a directory containing quality quality checked reads from Illumina sequencing (mandatory)"
echo "-c or -cDir [PATH]: define path to a directory containing corresponding assembled contigs from the reads used for mapping (mandatory)"
echo "NOTE: file names should be: XZY_R1.fastq XZY_R2.fastq for reads and XZY_contigs.fasta"
echo "-o [PATH]: Directory in which results will be saved. Default: ./"
echo "-t [INT]: number of threads allocated for the process. Default: 80"
echo "-mem [INT]: Memory for SPAdes assembly in GB (Default 500 GB)"
echo "-idx: output mapped read abundance as .txt file"
echo "-minlength [INT]: parameter for min viral contig length in base pairs. Will filter >=[minlength] Default=1500 bp"
echo "-id [INT]: Minimum sequence identity between identified phages. Default: 95 "
echo "-groups: Viral groups for virsorter2. Add comma separated after -group flag. Available: dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae. Default: All groups listed"
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
exit
;; 
*)  
if [ -z "$1" ]; then break; fi
    ERR=1
shift
;;

esac
done

################################## Modular Script #######################################
echo "Modules of the Virus Pipeline:"
echo "0: SPAdes"
echo "1: Metaviralspades"
echo "2: Virsorter2"
echo "3: VIBRANT"
echo "4: CheckV"
echo "5: FastANI"
echo "6: Bowtie2"
echo "7: VirHostMatcher-Net"
echo "8: Bacphlip"
echo "complete"
echo "10: ViPTreeGen"
echo "11: Vcontact2"



echo "12: DIAMOND Blastp ViralRefSeq"


echo Run module: 
read STEP

################################################################################################################################
#                                          Paths to third-party software and databases                                         #
################################################################################################################################
spades=/bioinf/home/benedikt.heyerhoff/Resources/SPAdes-3.15.4-Linux/bin/spades.py
metaviralspades=/bioinf/home/benedikt.heyerhoff/Resources/SPAdes-3.15.4-Linux/bin/metaviralspades.py
virsorter2=/bioinf/home/benedikt.heyerhoff/Resources/virsorter2/virsorter2.sif
checkvdb=/home/benedikth/checkv-db-v1.0
cluster1=/bioinf/home/benedikt.heyerhoff/Resources/vcontact2/bin/cluster_one-1.0.jar
vcontact2_tax_predict=/bioinf/home/benedikt.heyerhoff/Resources/vcontact2_tax_predict.py 
VirHostMatcher=/bioinf/home/benedikt.heyerhoff/Resources/VirHostMatcher-Net/VirHostMatcher-Net.py
viralrefseq=/bioinf/home/benedikt.heyerhoff/Resources/database/viral.3.protein.faa
marinehostlist=/bioinf/home/benedikt.heyerhoff/Resources/VirHostMatcher-Net/genome_list/marine_host_list.txt
virhostmatcherdata=/bioinf/home/benedikt.heyerhoff/Resources/VirHostMatcher-Net/data
vibrant=/bioinf/home/benedikt.heyerhoff/Resources/VIBRANT/VIBRANT_run.py
fastANI=/bioinf/home/benedikt.heyerhoff/Resources/FastANI/fastANI


###############################################################################################################################
#                                                      Variables check                                                        #
###############################################################################################################################

if [[ "${rDir}" == 0 ]]; then
    echo "ERROR!"
    echo -e "${red}No read-directory defined. Script aborted.${nc}"
    echo -e "${red}Use -rDir /PATH/TO/READS to tell the script where they are.${nc}"
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
(cd $rDir && ls *.f*) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq > $oDir/infiles.txt



###################################################################################################################################
# MODE: SPAdes                                                                                                                    #
###################################################################################################################################
if [[ "$STEP" == 0 || "$mode" == "complete" ]]; then
    echo -e "${green}Next step: SPAdes${nc}"
    echo ""
    echo ""
    for s in $(cat $oDir/infiles.txt);do
        echo -e "${blue}RUNNING SPAdes WITH "${s}" FASTQ files${nc}"
        mkdir -p $oDir/00_spades/
        $spades \
        -1 $rDir/${s}_*1.fastq \
        -2 $rDir/${s}_*2.fastq \
        --meta \
        -t $threads \
        -m $mem

        #move assembly from tmp to 00_spades directory
        mv $oDir/tmp/contigs.fasta $oDir/00_spades/${s}_contigs.fasta
        mv $oDir/tmp/scaffolds.fasta $oDir/00_spades/${s}_scaffolds.fasta
        mv $oDir/tmp/contigs.paths $oDir/00_spades/${s}_contigs.paths
        mv $oDir/tmp/scaffolds.paths $oDir/00_spades/${s}_scaffolds.paths
        rm -rf $oDir/tmp
    done
fi

###################################################################################################################################
# MODE: metaviralSPAdes                                                                                                           #
###################################################################################################################################

if [[ "$STEP" == 1 || "$mode" == "complete" ]]; then
    echo -e "${green}Next step: metaviralSPAdes${nc}"
    echo ""
    echo ""
    for s in $(cat $oDir/infiles.txt);do
        echo -e "${blue}RUNNING metaviralSPAdes WITH "${s}" FASTQ files${nc}"  
        mkdir -p $oDir/01_metaviralspades      
        python3 $metaviralspades \
        -1 $rDir/${s}_*1.fastq \
        -2 $rDir/${s}_*2.fastq \
        -o $oDir/tmp \
        -t $threads \
        -m $mem\

        #move assembly from tmp to 00_spades directory
        mv $oDir/tmp/contigs.fasta $oDir/01_metaviralspades/${s}_contigs.fasta
        mv $oDir/tmp/scaffolds.fasta $oDir/01_metaviralspades/${s}_scaffolds.fasta
        mv $oDir/tmp/contigs.paths $oDir/01_metaviralspades/${s}_contigs.paths
        mv $oDir/tmp/scaffolds.paths $oDir/01_metaviralspades/${s}_scaffolds.paths
        rm -rf $oDir/tmp 

        mkdir -p $oDir/03.1_all_viral_contigs
        sed 's/>.*/&|metaviralspades/' $oDir/01_metaviralspades/${s}_contigs.fasta > $oDir/03.1_all_viral_contigs/${s}/${s}.phage_contigs_metaviralspades.fna
    done
fi

###################################################################################################################################
# MODE: Virsorter2                                                                                                                #
###################################################################################################################################

if [[ "$STEP" == 2 || "$mode" == "complete" ]]; then
    echo -e "${green}Next step: Virsorter2${nc}"
    echo ""
    echo ""
    singularity run $virsorter2 config --set HMMSEARCH_THREADS=$threads #configure threads for hmmsearch
    for s in $(cat $oDir/infiles.txt);do
        echo -e "${blue}RUNNING Virsorter2 WITH ${green} "${s}" ${blue} files${nc}"
        singularity run $virsorter2 \
        run \
        -w $oDir/02_Virsorter2/tmp \
        -i $oDir/00_spades/${s}_contigs.fasta \
        --include-groups $groups \
        --min-length $minlength \
        --keep-original-seq \
        -j $threads

        mv $oDir/02_Virsorter2/tmp/*.fa* $oDir/02_Virsorter2/${s}_virsorter.fasta
        mv $oDir/02_Virsorter2/tmp/final-viral-score.tsv $oDir/02_Virsorter2/${s}_virsorter_viral_score.tsv
        #rm -rf $oDir/02_Virsorter2/tmp

        mkdir -p $oDir/03.1_all_viral_contigs    
        sed 's/>.*/&|virsorter2/' $oDir/02_Virsorter2/${s}_virsorter.fasta > $oDir/03.1_all_viral_contigs/${s}/${s}.phage_contigs_virsorter.fna
    done
fi

###################################################################################################################################
# MODE: VIBRANT                                                                                                                   #
###################################################################################################################################

if [[ "$STEP" == 3 || "$mode" == "complete" ]]; then
    echo -e "${green}Next step: VIBRANT${nc}"
    echo ""
    echo ""
    for s in $(cat $oDir/infiles.txt);do
        echo -ne "${blue}RUNNING VIBRANT WITH ${green} "${s}" ${blue}CONTIGS${nc}"
        mkdir -p $oDir/03_vibrant/${s}
        mkdir -p $oDir/03_vibrant/tmp
        python3 $vibrant \
        -i $oDir/00_spades/${s}_contigs.fasta \
        -t $threads \
        -folder $oDir/03_vibrant/tmp \
        -no_plot \
        -l $minlength

        #keep important files only
        mv $oDir/03_vibrant/tmp/*.phages_combined.fna $oDir/03_vibrant/${s}/${s}.phages_combined_vibrant.fna
        mv $oDir/03_vibrant/tmp/*.phages_combined.faa $oDir/03_vibrant/${s}/${s}.phages_combined_vibrant.faa
        mv $oDir/03_vibrant/tmp/VIBRANT_HMM_tables_parsed_${s}_contigs $oDir/03_vibrant/${s}/
        mv $oDir/03_vibrant/tmp/VIBRANT_genome_quality_${s}_contigs.tsv $oDir/03_vibrant/${s}/VIBRANT_genome_quality_${s}_contigs.tsv
        mv $oDir/03_vibrant/tmp/VIBRANT_annotations_${s}_contigs.tsv $oDir/03_vibrant/${s}/VIBRANT_annotations_${s}_contigs.tsv
        mv $oDir/03_vibrant/tmp/VIBRANT_genbank_table_${s}_contigs.tsv $oDir/03_vibrant/${s}/VIBRANT_genbank_table_${s}_contigs.tsv
        mv ${s}_contigs.phages_combined.txt $oDir/03_vibrant/${s}/${s}_contigs.phages_combined.txt
        rm -rf $oDir/03_vibrant/tmp

        mkdir -p $oDir/03.1_all_viral_contigs/
        sed 's/>.*/&|vibrant/' $oDir/03_vibrant/${s}/${s}.phages_combined.fna > $oDir/03.1_all_viral_contigs/${s}/${s}.phages_contigs_vibrant.fna
    done
fi

###################################################################################################################################
# MODE: CheckV                                                                                                                    #
###################################################################################################################################

if [[ "$STEP" == 4 || "$mode" == "complete" ]]; then
    echo -e "${green} Next step: CheckV ${nc}"
    echo ""
    echo ""

    (cd $oDir/03.1_all_viral_contigs && ls -d */ | cut -f1 -d'/' > $oDir/03.1_all_viral_contigs/infiles.txt) #list all folders in target directory
    
    for i in $(cat $oDir/03.1_all_viral_contigs/infiles.txt); do #copy all potential virus contigs of iteration into one directory
        cat $oDir/03.1_all_viral_contigs/${i}/*.f* > $oDir/03.1_all_viral_contigs/${i}_combined_contigs.fna #combine all contigs into one multifasta per metagenome    
    done
    #rm -rf $oDir/03.1_all_viral_contigs/${i}
    
    for i in $(cat $oDir/03.1_all_viral_contigs/infiles.txt); do
        echo -e "${blue}checking quality of viruses in${green} "${i}" ${blue}metagenome${nc}"
        checkv end_to_end \
        $oDir/03.1_all_viral_contigs/${i}_combined_contigs.fna \
        $oDir/04_checkv/${i} \
        -t $threads \
        -d $checkvdb
        rm -rf $oDir/04_checkv/${i}/tmp

        echo -e "${blue}CheckV finished${nc}"

        ############################################################
        #Quality filtering and binning
        ############################################################
        
        #export variables for R script
        export iter=$i
        export oDir=$oDir
        export minlength=$minlength
        ######debugging messages######
        echo -e "${green}R EXPORT PARAMETERS:"
        echo -e "${green}#####"
        echo -e "1 (iteration): "$iter""
        echo -e "2 (output dir): "$oDir""
        echo -e "5 (min length): "$minlength""
        echo -e "#####${nc}"
        ######

        
        echo -e "${green} Quality Filtering of Phages from "${i}" Metagenome ${nc}"
        #source R script
        Rscript /bioinf/home/benedikt.heyerhoff/checkv_output_score_filtering.R
    done
fi

###################################################################################################################################
# MODE: FastANI                                                                                                                   #
###################################################################################################################################

if [[ "$STEP" == 5 ]]; then
    echo -e "${green} Next step: FastANI ${nc}"
    echo ""
    echo ""

    mkdir -p $oDir/05_fastANI #make output directory
    (cd $oDir/04_checkv && ls -d */ | cut -f1 -d'/' > $oDir/05_fastANI/infiles.txt) #list all folders in target directory

    for s in $(cat $oDir/05_fastANI/infiles.txt); do
        cd $oDir/04_checkv/${s}/phages && ls | xargs readlink -f | uniq > $oDir/05_fastANI/fastANI_phage_files.txt
    done
       
    echo -e "${blue}Comparing Phages in "${i}" with FastANI ${nc}"
    for i in $(cat $oDir/05_fastANI/infiles.txt); do
        mkdir -p $oDir/05_fastANI/${i}
        $fastANI --ql $oDir/05_fastANI/fastANI_phage_files.txt \
        --rl $oDir/05_fastANI/fastANI_phage_files.txt \
        -o $oDir/05_fastANI/${i}/${i}.fastANI_out.txt \
        -t $threads

        #export variables for R script
        export iter=$i
        export oDir=$oDir
        export id=$id
        ######debugging messages######
        echo -e "${green}R EXPORT PARAMETERS:"
        echo -e "${green}#####"
        echo -e "1 (iteration): "$iter""
        echo -e "2 (output dir): "$oDir""
        echo -e "#####${nc}"
        ######
    done
fi

####################################### Contig Mapping ##################################
if [[ "$STEP" == 6 ]]; then
    echo -e "${green} Next step: Bowtie2 ${nc}"
    echo ""
    echo ""

    contigs=$( cd $oDir/04_checkv/length_filtered_contigs && ls *.f* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | wc -l
    reads=$( cd $rDir && ls *.fa* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d | wc -l

    mkdir -p $oDir/bam
    mkdir -p $oDir/sam
    mkdir -p $oDir/temp
    ls $rDir

    #( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/temp/files.txt

    for s in $(cat $oDir/remap.txt);do
        echo -e "Bulding bowtie2 database from ${green}"${s}"_contigs.fasta${nc} ..."
        bowtie2-build $oDir/04_checkv/length_filtered_contigs/${s}.fna $oDir/temp/Bowtie2.db
        DB=$oDir/temp/Bowtie2.db

        echo
        echo -e ${blue}"Mapping reads from ${green}"$s"_R1.fastq${nc} and ${green}"$s"_R2.fastq${nc} to database..." 
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
        samtools idxstats $oDir/bam/${s}.sorted.bam > $oDir/06_map/${s}.mapped.txt
        rm $oDir/bam/${s}.sorted.bam 


        echo ""
        echo -e "Finished mapping of ${green}"$s"${nc}!";
    done 
    rm -rf $oDir/temp/
fi

################################## VirHostMatcher-Net ###################################
if [[ "$STEP" == 7 ]]; then
    echo -e "${green} Next step: VirHostMatcher-net ${nc}"
    echo ""
    echo ""
    mkdir -p $oDir/07_VirHostMatcher/temp #tempdir for VirHostMatcher-Net
    for s in $(cat $oDir/infiles.txt);do
        python3 $VirHostMatcher \
        -q $oDir/04_checkv/VirHostMatcher/${s} \
        --short-contig \
        -o $oDir/07_VirHostMatcher/${s} \
        -n 1 \
        -t $threads \
        -l $marinehostlist \
        -i $oDir/07_VirHostMatcher/temp \
        -d $virhostmatcherdata;    
    done
fi

##################################### BACPHLIP ########################################
if [[ "$STEP" == 8 ]]; then
    echo -e "${green} Next step: BACPHLIP ${nc}"

    contigs=$( cd $oDir/04_checkv/length_filtered_contigs && ls *.f* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | wc -l
    for s in $(cat $oDir/infiles.txt);do
        bacphlip -i /valid/path/to/a/multigenome.fasta \
        --multi_fasta
    done



fi



###################################### Vcontact2 ########################################
if [[ "$STEP" == 9 ]]; then
    echo -e "${green} Next step: Vcontact2 ${nc}"
    echo ""
    echo ""
    echo -e "${blue} predicting ORFs with Prodigal${nc}"
    echo ""
    echo ""
    for s in $(cat $oDir/infiles.txt);do        
        mkdir -p $oDir/prodigal/
        prodigal -i $oDir/checkv/length_filtered_contigs/${s}.fna \
        -o $oDir/prodigal/${s}_proteins \
        -a $oDir/prodigal/${s}_proteins.faa \
        -p meta 
    done
    
    echo -e "${blue}Starting gene2genome step in ${s} Metagenome${nc}"
    echo ""
    echo ""

    for s in $(cat $oDir/infiles.txt);do
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
        --c1-bin $cluster1 \
        -t $threads \
        --output-dir $oDir/09_vcontact2/${s}
    done

    echo -e "${green} Assigning taxonomy to Vcontact2 output ${nc}"
    for s in $(cat $oDir/infiles.txt);do
        mkdir -p oDir/viral_taxonomy/${s}
        $vcontact2_tax_predict \
        -i $oDir/vcontact2/${s}/genome_by_genome_overview.csv \
        -o oDir/viral_taxonomy/${s}
    done
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

#################################### ViPTreeGen #######################################
if [[ "$STEP" == 8 ]]; then
    echo -e "${green} Next step: Merging length filtered contigs contigs into one file ${nc}"
    echo ""
    echo ""

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
    $oDir/08_ViPTree/length_filtered_contigs.fna \
    $oDir/08_ViPTree/ \
    --notree \
    --2D $oDir/ViPTree/length_filtered_contigs.fna  \
    --ncpus $threads
fi


echo ""
echo ""
echo -e "############################################"
echo -e "###################${BgGreen}Done!${nc}####################"
echo -e "############################################"

