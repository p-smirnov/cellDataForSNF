library(PharmacoGx)

options(timeout=600)

outputDir <- "./output"

if(!file.exists(outputDir)) dir.create(outputDir)

GDSC <- downloadPSet("GDSC_2020(v2-8.2)", saveDir="./dwl/")

GDSC.cell.info <- cellInfo(GDSC)[, c("cellid", "tissueid", "Cellosaurus.Disease.Type", "CellLine.Type")]
GDSC.tissue.disease.map <- unique(GDSC.cell.info[,c("tissueid", 'Cellosaurus.Disease.Type')])

GDSC.tissue.disease.map <- GDSC.tissue.disease.map[order(GDSC.tissue.disease.map[[1]], GDSC.tissue.disease.map[[2]]),]


GDSC.rna.data <- summarizeMolecularProfiles(GDSC, "Kallisto_0.46.1.rnaseq")

GDSC.rna.data <- t(SummarizedExperiment::assay(GDSC.rna.data))


GDSC.cnv.data <- summarizeMolecularProfiles(GDSC, "cnv")

GDSC.cnv.data <- t(SummarizedExperiment::assay(GDSC.cnv.data))

GDSC.mut.data <- summarizeMolecularProfiles(GDSC, "mutation_exome", summary.stat="or")

GDSC.mut.data <- t(SummarizedExperiment::assay(GDSC.mut.data))

GDSC.aac.data <- t(summarizeSensitivityProfiles(GDSC, "aac_recomputed"))

write.csv(GDSC.cell.info, file=file.path(outputDir, "GDSC_cell_info.csv"))
write.csv(GDSC.tissue.disease.map, file=file.path(outputDir, "GDSC_tissue_table.csv"))
write.csv(GDSC.rna.data, file=file.path(outputDir, "GDSC_rna_data.csv"))
write.csv(GDSC.cnv.data, file=file.path(outputDir, "GDSC_cnv_data.csv"))
write.csv(GDSC.mut.data, file=file.path(outputDir, "GDSC_mut_data.csv"))
write.csv(GDSC.aac.data, file=file.path(outputDir, "GDSC_aac_data.csv"))



CCLE <- downloadPSet("CCLE_2015", saveDir="./dwl/")
CTRPv2 <- downloadPSet("CTRPv2_2015", saveDir="./dwl/")


CCLE.cell.info <- cellInfo(CCLE)[, c("cellid", "tissueid", "Cellosaurus.Disease.Type", "CellLine.Type")]
CCLE.tissue.disease.map <- unique(CCLE.cell.info[,c("tissueid", 'Cellosaurus.Disease.Type')])

CCLE.tissue.disease.map <- CCLE.tissue.disease.map[order(CCLE.tissue.disease.map[[1]], CCLE.tissue.disease.map[[2]]),]


CCLE.rna.data <- summarizeMolecularProfiles(CCLE, "Kallisto_0.46.1.rnaseq")

CCLE.rna.data <- t(SummarizedExperiment::assay(CCLE.rna.data))


CCLE.cnv.data <- summarizeMolecularProfiles(CCLE, "cnv")

CCLE.cnv.data <- t(SummarizedExperiment::assay(CCLE.cnv.data))

CCLE.mut.data <- summarizeMolecularProfiles(CCLE, "mutation", summary.stat="or")

CCLE.mut.data <- t(SummarizedExperiment::assay(CCLE.mut.data))

CTRPv2.aac.data <- t(summarizeSensitivityProfiles(CTRPv2, "aac_recomputed"))


write.csv(CCLE.cell.info, file=file.path(outputDir, "CCLE_cell_info.csv"))
write.csv(CCLE.tissue.disease.map, file=file.path(outputDir, "CCLE_tissue_table.csv"))
write.csv(CCLE.rna.data, file=file.path(outputDir, "CCLE_rna_data.csv"))
write.csv(CCLE.cnv.data, file=file.path(outputDir, "CCLE_cnv_data.csv"))
write.csv(CCLE.mut.data, file=file.path(outputDir, "CCLE_mut_data.csv"))
write.csv(CTRPv2.aac.data, file=file.path(outputDir, "CTRPv2_aac_data.csv"))



gCSI <- downloadPSet("gCSI_2017", saveDir="./dwl/")
gCSI.cell.info <- cellInfo(gCSI)[, c("cellid", "tissueid", "Cellosaurus.Disease.Type", "CellLine.Type")]
gCSI.tissue.disease.map <- unique(gCSI.cell.info[,c("tissueid", 'Cellosaurus.Disease.Type')])

gCSI.tissue.disease.map <- gCSI.tissue.disease.map[order(gCSI.tissue.disease.map[[1]], gCSI.tissue.disease.map[[2]]),]


gCSI.rna.data <- summarizeMolecularProfiles(gCSI, "Kallisto_0.46.1.rnaseq")

gCSI.rna.data <- t(SummarizedExperiment::assay(gCSI.rna.data))


gCSI.cnv.data <- summarizeMolecularProfiles(gCSI, "cnv")

gCSI.cnv.data <- t(SummarizedExperiment::assay(gCSI.cnv.data))

gCSI.mut.data <- summarizeMolecularProfiles(gCSI, "mutation", summary.stat="or")

gCSI.mut.data <- t(SummarizedExperiment::assay(gCSI.mut.data))

gCSI.aac.data <- t(summarizeSensitivityProfiles(gCSI, "aac_recomputed"))


write.csv(gCSI.cell.info, file=file.path(outputDir, "gCSI_cell_info.csv"))
write.csv(gCSI.tissue.disease.map, file=file.path(outputDir, "gCSI_tissue_table.csv"))
write.csv(gCSI.rna.data, file=file.path(outputDir, "gCSI_rna_data.csv"))
write.csv(gCSI.cnv.data, file=file.path(outputDir, "gCSI_cnv_data.csv"))
write.csv(gCSI.mut.data, file=file.path(outputDir, "gCSI_mut_data.csv"))
write.csv(gCSI.aac.data, file=file.path(outputDir, "gCSI_aac_data.csv"))


