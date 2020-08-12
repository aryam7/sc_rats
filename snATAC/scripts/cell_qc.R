#!/usr/bin/env Rscript

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(AnnotationHub))

args = commandArgs(trailingOnly = TRUE)
sample = args[1]
out = args[2]

# Files from CellRanger:
# 1) peak/cell matrix - counts of the Tn5 cut sites within each peak
raw = paste0(sample, "/filtered_peak_bc_matrix.h5")
# 2) some metadata
metadata = paste0(sample, "/singlecell.csv")
# 3) a list of unique fragments across all single cells
fragments = paste0(sample, "/fragments.tsv.gz")

EnsDb.Rnorvegicus.v98 = query(AnnotationHub(), pattern = c("Rattus Norvegicus", "EnsDb", 98))[[1]]


counts = Read10X_h5(filename = raw)
metadata = read.csv(
  file = metadata,
  header = TRUE,
  row.names = 1
)

amygdala = CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

fragment.path = fragments
amygdala = SetFragments(
  object = amygdala,
  file = fragment.path
)

# note we currently only do this for chr1 because it can take a while
amygdala <- NucleosomeSignal(object = amygdala, region = "1-1-282763074")

amygdala$pct_reads_in_peaks <- amygdala$peak_region_fragments / amygdala$passed_filters * 100

p1 <- VlnPlot(amygdala, c('pct_reads_in_peaks'), pt.size = 0.1) + NoLegend()
p2 <- VlnPlot(amygdala, c('peak_region_fragments'), pt.size = 0.1) + NoLegend() + ylim(0, 7000)
p3 <- VlnPlot(amygdala, c('nucleosome_signal'), pt.size = 0.1, y.max=20) + ylim(0,5)

p1 | p2 | p3
ggsave(paste0(out, "/vln.pdf"), width=11.15, height=6.69, units='in')

amygdala$nucleosome_group <- ifelse(5 < amygdala$nucleosome_signal | amygdala$nucleosome_signal < 0.2, '5 < NS or NS < 0.2', '0.2 < NS < 5')
FragmentHistogram(object = amygdala, group.by = 'nucleosome_group', region='1-1-282763074')
ggsave(paste0(out, "/fragment_hist.pdf"), width=11.15, height=6.69, units='in')

# create granges object with TSS positions
# not sure which reference genome is the right one to use here, but I chose Rnorvegicus from EnsDB
# we grab the TSS positions from the reference genome
gene.ranges <- genes(EnsDb.Rnorvegicus.v98)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding',]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse', species='Rattus_norvegicus')

tss.ranges <- resize(gene.ranges, width = 1, fix = "start")
# drop the mt genome because it is not in our tabix index
tss.ranges <- dropSeqlevels(tss.ranges, "MT", pruning.mode = 'coarse')
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse', species='Rattus_norvegicus')
amygdala <- TSSEnrichment(object = amygdala, tss.positions = tss.ranges)

amygdala$high.tss <- ifelse(amygdala$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(amygdala, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
ggsave(paste0(out, "/tss.pdf"), width=11.15, height=6.69, units='in')

amygdala <- subset(
  x = amygdala,
  subset = peak_region_fragments > 600 &
    peak_region_fragments < 6000 &
    pct_reads_in_peaks > 10 &
    nucleosome_signal < 5 &
    nucleosome_signal > 0.2 &
    TSS.enrichment > 2
)

save(amygdala, file=paste0(out, "/cells.rda"))
