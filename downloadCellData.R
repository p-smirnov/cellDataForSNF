library(PharmacoGx)

options(timeout=600)

outputDir <- "./output"

if(!file.exists(outputDir)) dir.create(outputDir)

GDSC <- downloadPSet("GDSC_2020(v2-8.2)", saveDir="./dwl/")

cell.info <- cellInfo(GDSC)[, c("cellid", "tissueid", "Cellosaurus.Disease.Type", "CellLine.Type")]
tissue.disease.map <- unique(cell.info[,c("tissueid", 'Cellosaurus.Disease.Type')])

tissue.disease.map <- tissue.disease.map[order(tissue.disease.map[[1]], tissue.disease.map[[2]]),]


rna.data <- summarizeMolecularProfiles(GDSC, "rna")

rna.data <- t(SummarizedExperiment::assay(rna.data))


cnv.data <- summarizeMolecularProfiles(GDSC, "cnv")

cnv.data <- t(SummarizedExperiment::assay(cnv.data))

mut.data <- summarizeMolecularProfiles(GDSC, "mutation_exome", summary.stat="or")

mut.data <- t(SummarizedExperiment::assay(mut.data))


write.csv(cell.info, file=file.path(outputDir, "GDSC_cell_info.csv"))
write.csv(tissue.disease.map, file=file.path(outputDir, "GDSC_tissue_table.csv"))
write.csv(rna.data, file=file.path(outputDir, "GDSC_rna_data.csv"))
write.csv(cnv.data, file=file.path(outputDir, "GDSC_cnv_data.csv"))
write.csv(mut.data, file=file.path(outputDir, "GDSC_mut_data.csv"))


