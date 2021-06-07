library(PharmacoGx)

options(timeout=600)

outputDir <- "./output"

if(!file.exists(outputDir)) dir.create(outputDir)

GDSC2 <- downloadPSet("GDSC_2020(v2-8.2)", saveDir="./dwl/")

GDSC2.cell.info <- cellInfo(GDSC2)[, c("cellid", "tissueid", "Cellosaurus.Disease.Type", "CellLine.Type")]
GDSC2.tissue.disease.map <- unique(GDSC2.cell.info[,c("tissueid", 'Cellosaurus.Disease.Type')])

GDSC2.tissue.disease.map <- GDSC2.tissue.disease.map[order(GDSC2.tissue.disease.map[[1]], GDSC2.tissue.disease.map[[2]]),]


GDSC2.rna.data <- summarizeMolecularProfiles(GDSC2, "Kallisto_0.46.1.rnaseq")

GDSC2.rna.data <- t(SummarizedExperiment::assay(GDSC2.rna.data))


GDSC2.cnv.data <- summarizeMolecularProfiles(GDSC2, "cnv")

GDSC2.cnv.data <- t(SummarizedExperiment::assay(GDSC2.cnv.data))

GDSC2.mut.data <- summarizeMolecularProfiles(GDSC2, "mutation_exome", summary.stat="or")

GDSC2.mut.data <- t(SummarizedExperiment::assay(GDSC2.mut.data))

GDSC2.aac.data <- t(summarizeSensitivityProfiles(GDSC2, "aac_recomputed"))

write.csv(GDSC2.cell.info, file=file.path(outputDir, "GDSC2_cell_info.csv"))
write.csv(GDSC2.tissue.disease.map, file=file.path(outputDir, "GDSC2_tissue_table.csv"))
write.csv(GDSC2.rna.data, file=file.path(outputDir, "GDSC2_rna_data.csv"))
write.csv(GDSC2.cnv.data, file=file.path(outputDir, "GDSC2_cnv_data.csv"))
write.csv(GDSC2.mut.data, file=file.path(outputDir, "GDSC2_mut_data.csv"))
write.csv(GDSC2.aac.data, file=file.path(outputDir, "GDSC2_aac_data.csv"))


GDSC1 <- downloadPSet("GDSC_2020(v1-8.2)", saveDir="./dwl/")

GDSC1.cell.info <- cellInfo(GDSC1)[, c("cellid", "tissueid", "Cellosaurus.Disease.Type", "CellLine.Type")]
GDSC1.tissue.disease.map <- unique(GDSC1.cell.info[,c("tissueid", 'Cellosaurus.Disease.Type')])

GDSC1.tissue.disease.map <- GDSC1.tissue.disease.map[order(GDSC1.tissue.disease.map[[1]], GDSC1.tissue.disease.map[[2]]),]


GDSC1.rna.data <- summarizeMolecularProfiles(GDSC1, "Kallisto_0.46.1.rnaseq")

GDSC1.rna.data <- t(SummarizedExperiment::assay(GDSC1.rna.data))


GDSC1.cnv.data <- summarizeMolecularProfiles(GDSC1, "cnv")

GDSC1.cnv.data <- t(SummarizedExperiment::assay(GDSC1.cnv.data))

GDSC1.mut.data <- summarizeMolecularProfiles(GDSC1, "mutation_exome", summary.stat="or")

GDSC1.mut.data <- t(SummarizedExperiment::assay(GDSC1.mut.data))

GDSC1.aac.data <- t(summarizeSensitivityProfiles(GDSC1, "aac_recomputed"))

write.csv(GDSC1.cell.info, file=file.path(outputDir, "GDSC1_cell_info.csv"))
write.csv(GDSC1.tissue.disease.map, file=file.path(outputDir, "GDSC1_tissue_table.csv"))
write.csv(GDSC1.rna.data, file=file.path(outputDir, "GDSC1_rna_data.csv"))
write.csv(GDSC1.cnv.data, file=file.path(outputDir, "GDSC1_cnv_data.csv"))
write.csv(GDSC1.mut.data, file=file.path(outputDir, "GDSC1_mut_data.csv"))
write.csv(GDSC1.aac.data, file=file.path(outputDir, "GDSC1_aac_data.csv"))


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


