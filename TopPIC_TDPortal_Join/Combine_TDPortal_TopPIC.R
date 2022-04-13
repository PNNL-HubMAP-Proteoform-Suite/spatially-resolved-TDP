library(tidyverse)

TDPortal <- read.csv("All_Hits_Rat_Brain.csv", header = TRUE,check.names = FALSE)
 

load("x53m.RData")


### Change headers to match TDPortal
### Remove redundant or extra columns 
w_cluster <- w_cluster %>%
  mutate(`Scan#` = `Scan(s)`) %>%
  select(-`Scan(s)`) %>%
  select(-NormRecalMass, -NormRTalign, -cluster, -cluster_new, -Obs, -ProtLen, -AccMap, -AnnType,
         -`First residue`, -`Last residue`, -`Retention time`, -Charge, -`Protein accession`,
         -`#variable PTMs`, -Mass, -RecalMass,-RTest, -RTalign)


############  Need to fix the multiple identifications per spectrum issue with TDPortal
TDPortal <- TDPortal %>%
  mutate(Dataset = `File Name`) %>%
  mutate(Dataset = substr(Dataset,15,25)) %>%
  mutate(CV = case_when(grepl("cv50", `File Name`) ~ "-50",
                        grepl("cv40", `File Name`) ~ "-40",
                        grepl("cv30", `File Name`) ~ "-30")) %>%
  separate_rows(`Fragment Scans`) %>%
  mutate(`Scan#` = as.integer(`Fragment Scans`)) %>%
  group_by(Dataset, CV, `Scan#`) %>%
  slice_min(order_by = `Global Q-value`, n = 1) %>%
  ungroup() %>%
  select(-PFR, -`Uniprot Id`,
          -`% Cleavages`, -`P-score`, -`C-score`)

names(w_cluster) <- paste0("TopPIC.", names(w_cluster))

names(TDPortal) <- paste0("TDPortal.", names(TDPortal))

combined <- full_join(w_cluster, TDPortal, by = c("TopPIC.Dataset" = "TDPortal.Dataset",
                                                  "TopPIC.CV" = "TDPortal.CV",
                                                  "TopPIC.Scan#" = "TDPortal.Scan#"))


combined <- combined %>%
  mutate(SeqConflict = case_when(TopPIC.cleanSeq != TDPortal.Sequence ~ "Yes",
                                 TopPIC.cleanSeq == TDPortal.Sequence ~ "No",
                                 TRUE ~ "NA"))

combined <- combined %>%
  mutate(Software = case_when(is.na(`TDPortal.Proteoform Level`) ~ "TopPIC Only",
                              is.na(`TopPIC.Proteoform Level`) ~ "TDPortal Only",
                              SeqConflict != "NA" ~ "Both"))


write.csv(combined, file = "TDportal_TopPIC_overlappingIDs.csv")
