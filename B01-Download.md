# download

```
mkdir -p {1_fq/log,2_cellranger/log,3_velocyto/log}
```



```
cat filereport_read_run_PRJNA671310_tsv.txt | cut -f1 | sed 1d | tr ";" "\n" | sed 1,2d | awk '{print "ascp -k 1 -QT -l 300m -P33001 -i ~/miniconda3/envs/py3_velocity/etc/asperaweb_id_dsa.openssh era-fasp@"$1" . "}' > download.sh
nohup bash download.sh &>log/download.log &
```



```
cat filereport_read_run_PRJNA545094_tsv.txt | grep -v bonemarrow | cut -f1 | sed 1d | tr ";" "\n" | awk '{print "ascp -k 1 -QT -l 300m -P33001 -i ~/miniconda3/envs/py3_velocity/etc/asperaweb_id_dsa.openssh era-fasp@"$1" . "}' > download.sh
nohup bash download.sh &>log/download.log &
```



# Test md5

```
cat filereport_read_run_PRJNA545094_tsv.txt | grep -v bonemarrow | cut -f5,1 | sed 1d > raw_md5.txt
cat filereport_read_run_PRJNA671310_tsv.txt | cut -f5,1 | sed 1d > raw_md5.txt

bash ena2md5.sh raw_md5.txt > md5.txt
```

```
cat > ena2md5.sh
paste <(cat $1 | \
 cut -f1 | awk -F ';' '{print$1,$2}' | \
 tr ' ' '\n') <(cat $1  | \
 cut -f2 | while read id; do   echo ${id}_1.fastq.gz;   echo ${id}_2.fastq.gz; done) | \
 awk -F '\t' '{print$1,$2}' | sed 's/[ ]/  /'
```

```
md5sum -c md5.txt
```



# Move file

```
cat sample.txt | while read sample
do
mkdir $sample
ls | grep $sample | while read id
do
mv $id ${sample}/${id}
done
done
```





# cellranger

```
cat > submit.sh
cat $1 | while read id
do
    if ((i%$2==$3))
    then
        echo $id | bash
    fi
i=$((i+1))
done

```

```
cat ../fq/filereport_read_run_PRJNA671310_tsv.txt | cut -f3 | sed 1d | sort | uniq | awk 'BEGIN{p=4} {print "( cellranger count --id="$1" --transcriptome=/home/pengyz/3-data/10x/refdata-gex-GRCh38-2020-A --fastqs=../fq/"$1" --sample="$1" --nosecondary --localcores "p" ) &>log/"$1".log"}' > cellranger.sh
```

```
for i in {0..2}; do nohup bash submit.sh cellranger2.sh 3 $i 2>&1 & done
```



```
nohup cellranger count --id=UPN2 --transcriptome=/home/pengyz/3-data/10x/refdata-gex-GRCh38-2020-A --fastqs=../1_fq/UPN2 --sample=UPN2 --nosecondary --localcores 16 &>log/UPN2.log &

nohup cellranger count --id=UPN3 --transcriptome=/home/pengyz/3-data/10x/refdata-gex-GRCh38-2020-A --fastqs=../1_fq/UPN3 --sample=UPN3 --nosecondary --localcores 16 &>log/UPN2.log &

nohup cellranger count --id=UPN5 --transcriptome=/home/pengyz/3-data/10x/refdata-gex-GRCh38-2020-A --fastqs=../1_fq/UPN5 --sample=UPN5 --nosecondary --localcores 16 &>log/UPN5.log &
```



```
nohup samtools sort -@ 16 -t CB -O BAM -o UPN2/outs/cellsorted_possorted_genome_bam.bam UPN2/outs/possorted_genome_bam.bam &
nohup samtools sort -@ 16 -t CB -O BAM -o UPN3/outs/cellsorted_possorted_genome_bam.bam UPN3/outs/possorted_genome_bam.bam &
nohup samtools sort -@ 16 -t CB -O BAM -o UPN5/outs/cellsorted_possorted_genome_bam.bam UPN5/outs/possorted_genome_bam.bam &
```





# velocity

```
nohup velocyto run10x -m /home/pengyz/3-data/10x/hg38_rmsk.gtf ../2_cellranger/UPN2 /home/pengyz/3-data/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf &>log/UPN2.log &
nohup velocyto run10x -m /home/pengyz/3-data/10x/hg38_rmsk.gtf ../2_cellranger/UPN3 /home/pengyz/3-data/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf &>log/UPN3.log &
nohup velocyto run10x -m /home/pengyz/3-data/10x/hg38_rmsk.gtf ../2_cellranger/UPN5 /home/pengyz/3-data/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf &>log/UPN5.log &
```

```
mv 2_cellranger/UPN2/velocyto/UPN2.loom 3_velocyto
mv 2_cellranger/UPN3/velocyto/UPN3.loom 3_velocyto
mv 2_cellranger/UPN5/velocyto/UPN5.loom 3_velocyto

```





# Extract 10x mtx results into Aliyunpan

```
ls 2_cellranger | grep UPN | while read id
do
aliyunpan upload 2_cellranger/${id}/outs/web_summary.html GastricAtlas/E-MTAB-9221/${id}/
aliyunpan upload 2_cellranger/${id}/outs/filtered_feature_bc_matrix/* GastricAtlas/E-MTAB-9221/${id}/
aliyunpan upload 3_velocyto/${id}.loom GastricAtlas/E-MTAB-9221/${id}/
done
```



```
nohup unrar x Raw.rar &
```



- tcr & bcr HRA000704 multi

```
ls | grep -v csv | grep GC | while read id
do
aliyunpan upload ${id}/outs/per_sample_outs/${id}/web_summary.html GastricAtlas/HRA000704_vdj/${id}/
aliyunpan upload ${id}/outs/per_sample_outs/${id}/vdj_t/filtered_contig_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_t
aliyunpan upload ${id}/outs/per_sample_outs/${id}/vdj_t/consensus_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_t
aliyunpan upload ${id}/outs/per_sample_outs/${id}/vdj_t/clonotypes.csv GastricAtlas/HRA000704_vdj/${id}/vdj_t
aliyunpan upload ${id}/outs/per_sample_outs/${id}/vdj_b/filtered_contig_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_b
aliyunpan upload ${id}/outs/per_sample_outs/${id}/vdj_b/consensus_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_b
aliyunpan upload ${id}/outs/per_sample_outs/${id}/vdj_b/clonotypes.csv GastricAtlas/HRA000704_vdj/${id}/vdj_b
done
```



- vdj 

```
ls | grep -v csv | grep GC | while read idd
do
id=${idd:0:5}
aliyunpan upload ${id}_vdj_b/outs/web_summary.html GastricAtlas/HRA000704_vdj/${id}/vdj_b 
aliyunpan upload ${id}_vdj_b/outs/filtered_contig_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_b
aliyunpan upload ${id}_vdj_b/outs/consensus_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_b
aliyunpan upload ${id}_vdj_b/outs/clonotypes.csv GastricAtlas/HRA000704_vdj/${id}/vdj_b
aliyunpan upload ${id}_vdj_t/outs/web_summary.html GastricAtlas/HRA000704_vdj/${id}/vdj_t
aliyunpan upload ${id}_vdj_t/outs/filtered_contig_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_t
aliyunpan upload ${id}_vdj_t/outs/consensus_annotations.csv GastricAtlas/HRA000704_vdj/${id}/vdj_t
aliyunpan upload ${id}_vdj_t/outs/clonotypes.csv GastricAtlas/HRA000704_vdj/${id}/vdj_t
done
```



