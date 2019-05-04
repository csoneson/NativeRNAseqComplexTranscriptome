suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

## Define remapping of sample names
remap <- c(Srpk_1 = "ONT-NSK007-HAP_srpk_1", 
           Srpk_2 = "ONT-NSK007-HAP_srpk_2", 
           wt1 = "ONT-NSK007-HAP_wt_1", 
           wt2 = "ONT-NSK007-HAP_wt_2",
           wt_4_FGCZ = "ONT-PCS108-HAP_wt_4", 
           `20171207_1645_p2557_4017_2_ALLREADS.pass` = "ONT-PCS108-HAP_wt_4",
           srpk_1_DCS108 = "ONT-DCS108-HAP_srpk_1", 
           srpk_2_DCS108 = "ONT-DCS108-HAP_srpk_2",
           wt_1_DCS108 = "ONT-DCS108-HAP_wt_1", 
           wt_2_DCS108 = "ONT-DCS108-HAP_wt_2",
           srpk_1_RNA001 = "ONT-RNA001-HAP_srpk_1", 
           srpk_2_RNA001 = "ONT-RNA001-HAP_srpk_2",
           srpk_3_RNA001 = "ONT-RNA001-HAP_srpk_3", 
           srpk_4_RNA001 = "ONT-RNA001-HAP_srpk_4",
           srpk_5_RNA001 = "ONT-RNA001-HAP_srpk_5", 
           srpk_6_RNA001 = "ONT-RNA001-HAP_srpk_6",
           wt_1_RNA001 = "ONT-RNA001-HAP_wt_1", 
           wt_2_RNA001 = "ONT-RNA001-HAP_wt_2",
           wt_3_RNA001 = "ONT-RNA001-HAP_wt_3", 
           wt_4_RNA001 = "ONT-RNA001-HAP_wt_4",
           wt_5_RNA001 = "ONT-RNA001-HAP_wt_5", 
           wt_6_RNA001 = "ONT-RNA001-HAP_wt_6",
           SS2_srpk_1 = "TempSwitch-HAP_srpk_1", 
           SS2_srpk_2 = "TempSwitch-HAP_srpk_2",
           SS2_wt_1 = "TempSwitch-HAP_wt_1", 
           SS2_wt_2 = "TempSwitch-HAP_wt_2",
           `20170918.A-dSprk1_1` = "Illumina_srpk_1",
           `20170918.A-dSprk1_2` = "Illumina_srpk_2",
           `20170918.A-dSprk1_3` = "Illumina_srpk_3",
           `20170918.A-dSprk1_4` = "Illumina_srpk_4",
           `20170918.A-WT_1` = "Illumina_wt_1",
           `20170918.A-WT_2` = "Illumina_wt_2",
           `20170918.A-WT_3` = "Illumina_wt_3", 
           `20170918.A-WT_4` = "Illumina_wt_4",
           HEK_1 = "ONT-RNA001-HEK_wt_1", 
           HEK_2 = "ONT-RNA001-HEK_wt_2", 
           HEK_3 = "ONT-RNA001-HEK_wt_3",
           HEK_4 = "ONT-RNA001-HEK_wt_4", 
           HEK_5 = "ONT-RNA001-HEK_wt_5",
           `20181026_0911_p2557_4959_1.pass` = "ONT-PCS109-HAP_wt_4",
           wt4 = "ONT-SQK-PCS109-HAP_wt_4",
           p2557_5265_1 = "ONT-SQK-PCS109-HAP_wt_1",
           p2557_5265_2 = "ONT-SQK-PCS109-HAP_wt_2",
           p2557_5265_5 = "ONT-SQK-PCS109-HAP_srpk_1",
           p2557_5265_6 = "ONT-SQK-PCS109-HAP_srpk_2")

## Create a sample annotation table
sample_annotation <- data.frame(
  sample_orig = names(remap),
  sample_remap = remap,
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(sample_remap, into = c("dataset", "condition", "sample_nbr"),
                  sep = "_", remove = FALSE)

## Define remapping of data set names
remapds <- c(pilot = "ONT-NSK007-HAP",
             NSK007 = "TempSwitch-HAP",
             DCS108 = "ONT-DCS108-HAP",
             RNA001 = "ONT-RNA001-HAP",
             HEK293RNA = "ONT-RNA001-HEK",
             FGCZ = "ONT-PCS108-HAP", 
             Illumina = "Illumina",
             FGCZ_PCS109 = "ONT-PCS109-HAP",
             FGCZ_SQK_PCS109 = "ONT-SQK-PCS109-HAP",
             FGCZ_PCS109_GridION = "ONT-SQK-PCS109-HAP")

## Define color schemes
ds_colors <- c(`ONT-DCS108-HAP`= "#66C2A5", `ONT-NSK007-HAP` = "#FC8D62", 
               `ONT-RNA001-HEK` = "#8DA0CB", `TempSwitch-HAP` = "darkgrey", 
               `ONT-PCS108-HAP` = "#A6D854", `ONT-RNA001-HAP` = "#FFD92F", 
               Illumina = "#B3B3B3", `ONT-PCS109-HAP` = "#E78AC3",
               `ONT-SQK-PCS109-HAP` = "darkblue")

ds_order <- c("Illumina", "ONT-NSK007-HAP", "ONT-DCS108-HAP", "ONT-RNA001-HAP", 
              "ONT-RNA001-HEK")

## Help function to remove the dataset name from the sample name
removeDatasetFromSample <- function(sample, dataset) {
  sapply(seq_along(sample), function(i) gsub(paste0(dataset[i], "_"), "", sample[i]))
}
