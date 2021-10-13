#Author: Ziyan Lin
#Description: this is the qiime2 script
module purge
module load miniconda3/4.6.14
conda activate qiime2-2020.2

##################################
#qiime2
##################################
#import sequence (paired-end)
qiime tools import \
--type SampleData[PairedEndSequencesWithQuality] \
--input-path /gpfs/home/zl1887/SaxenaLab/2021-03-09_IECOH/2021-03-09/rawdata/fastq \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path demux_seqs.qza

#demux visualization
qiime demux summarize \
  --i-data ./demux_seqs.qza \
  --o-visualization ./demux_seqs.qzv

#dada2 denoising
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ./demux_seqs.qza \
--p-trim-left-f 17 \
--p-trim-left-r 21 \
--p-trunc-len-f 300 \
--p-trunc-len-r 196 \
--o-table ./dada2_table.qza \
--o-representative-sequences ./dada2_rep_seq.qza \
--o-denoising-stats ./dada2_denoise_stats.qza \
--p-n-threads 0 #all available core would be used

#visualization, denosie stats
qiime metadata tabulate \
  --m-input-file ./dada2_denoise_stats.qza  \
  --o-visualization ./dada2_denoise_stats.qzv

#visualize metadata
qiime metadata tabulate \
  --m-input-file metadata.tsv \
  --o-visualization metadata.qzv

#Feature table summary
qiime feature-table summarize \
  --i-table ./dada2_table.qza \
  --m-sample-metadata-file ./metadata.tsv \
  --o-visualization ./dada2_table.qzv
  
#visualize rep-seq
qiime feature-table tabulate-seqs \
  --i-data dada2_rep_seq.qza \
  --o-visualization dada2_rep_seq.qzv

####################
#rarefied
qiime feature-table rarefy \
  --i-table dada2_table.qza \
  --p-sampling-depth 1001 \
  --o-rarefied-table rarefied_table.qza
  
qiime feature-table summarize \
  --i-table ./rarefied_table.qza \
  --m-sample-metadata-file ./metadata.tsv \
  --o-visualization ./rarefied_table.qzv

####################
#Phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences dada2_rep_seq.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

####################
#import taxonomic assignment (HOMD database)
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ./REF/HOMD_16S_rRNA_RefSeq_V15.2.fasta \
  --output-path homd_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ./REF/HOMD_16S_rRNA_RefSeq_V15.2.qiime.taxonomy \
  --output-path ref-taxonomy.qza
  
#extract reference reads
qiime feature-classifier extract-reads \
  --i-sequences homd_otus.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --o-reads ref-seqs.qza

#train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

#test the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads dada2_rep_seq.qza \
  --o-classification taxonomy.qza

#visualize taxonomy
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#taxo barplot
qiime taxa barplot \
  --i-table rarefied_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

#Creating a TSV BIOM table
qiime tools export \
  --input-path ./rarefied_table.qza \
  --output-path ./rarefied-biom_table

biom convert -i rarefied-biom_table/feature-table.biom -o rarefied-biom_table/feature-table.tsv --to-tsv

#Creating taxonomy table
qiime tools export \
  --input-path taxonomy.qza \
  --output-path ./rarefied-biom_table
  
mv taxonomy.tsv biom-taxonomy.tsv

#vim biom-taxonomy.tsv
#change first line: #OTUID    taxonomy    confidence

#exporting tree
qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path exported-tree/unrooted
  
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported-tree/rooted


