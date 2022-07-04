# virus_assembly_pipeline
++Important: Work in progress!++
This pipeline will be integrated into a larger pipeline soon (https://github.com/LeonDlugosch/MetaSeq-Toolkit). For now this pipeline does no QC and quality assessment prior to assembly as it will be added later!
This pipeline assembles metagenomic reads, detects viruses via metaviralspades, Virsorter2 and VIBRANT, assesses quality of viral assemblies via checkv, calculates ANI and eliminates detected duplicate viruses, determines abundance of viral contigs via Bowtie2, determines viral lifestyle of contigs, clusters by protein similarity to infer taxonomy via vcontact2, and blasts viral contigs to viralrefseq to infer taxonomy.
