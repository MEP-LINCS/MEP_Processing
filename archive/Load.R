#Loader for shiny image QA app
library(data.table)
# dt <- fread("/graylab/share/dane/MEP-LINCS/MEP_LINCS/AnnotatedData/MCF10A_SS2_Level1.txt")
# cDT <- dt[,.(Barcode,Well,Spot,PrintSpot,ArrayRow,ArrayColumn,ImageID,ECMp,Ligand,Nuclei_CP_AreaShape_Area,Nuclei_CP_AreaShape_Center_X,Nuclei_CP_AreaShape_Center_Y,Spot_PA_SpotCellCount,Nuclei_PA_Gated_EdUPositive,Nuclei_PA_Cycle_State,Nuclei_PA_AreaShape_Neighbors)]
# save(cDT,file = "/graylab/share/dane/MEMAImageQA/MCF10A_SS2.RData")
load(file = "/graylab/share/dane/MEMAImageQA/MCF10A_SS2.RData")
