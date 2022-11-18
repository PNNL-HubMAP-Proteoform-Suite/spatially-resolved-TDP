## HubMAP
## David Degnan, Pacific Northwest National Laboratory
## Last Edited: 10/06/2022

# Load libraries
library(dplyr)
library(purrr)
library(tibble)

setwd("~/Git_Repos/spatially-resolved-TDP/")

######################################################
## STEP 1: Collapse down TopPIC and TDPortal Output ##
######################################################

# REASONING:
#-------------------------------------------------------------------------------
# There are duplicate proteoform "IDs" in both the TopPIC and TDPortal output.
# To collapse down TopPIC output, a "GCC" column was created to identify unique
# groups and the lowest e-score selected per GCC. In TDPortal, the lowest e-score
# was selected per TDPortal Accession and Monoisotopic Mass. 

# Load proteoform data
Proteoform <- read.csv("ProMexAlign_Proteoforms/Inputs/Proteoforms/TDportal_TopPIC_overlappingIDs.csv")

## TopPIC-----------------------------------------------------------------------

# Pull Modification Names
ModName <- lapply(Proteoform$TopPIC.Proteoform, function(Name) {
  
  Name %>% 
    gsub(pattern = "[", replacement = "^[", fixed = T) %>% 
    gsub(pattern= "]", replacement = "]^", fixed = T) %>% 
    strsplit("^", fixed = T) %>% 
    unlist() %>%
    lapply(function(x) {if (grepl("[", x, fixed = T)) {x}}) %>%
    unlist() %>%
    paste0(collapse = " & ")
  
}) %>% unlist()

# Pull TopPIC-Relevant Data. "GCC" is a unique ID based on a soon-to-be published method.
# Pull the lowest score per GCC. 
TopPIC <- Proteoform %>%
  select(c(TopPIC.gcc, TopPIC.UniProtAcc, TopPIC.Precursor.mass, 
           TopPIC.RTmin, TopPIC.First_AA, TopPIC.Last_AA,
           TopPIC.E.value)) %>%
  mutate(
    "Unique.ID" = TopPIC.gcc,
    "Mass" = TopPIC.Precursor.mass,
    "Retention.Time" = TopPIC.RTmin, 
    "Proteoform.Details" = paste0(TopPIC.UniProtAcc, " ", " (", 
      TopPIC.First_AA, "-", TopPIC.Last_AA, ") ", ModName),
    "E.Value" = TopPIC.E.value
  ) %>%
  select(Unique.ID, Mass, Retention.Time, Proteoform.Details, E.Value) %>%
  group_by(Unique.ID) %>%
  slice(which.min(E.Value)) %>%
  arrange(Mass)

## TDPortal---------------------------------------------------------------------

# Pull TDPortal relevant data, group by accession number, and pull out the 
# lowest e-score. Group by Accession and Mass for TDPortal
TDPortal <- Proteoform %>%
  select(c(TDPortal.Accession, TDPortal.Monoisotopic.Mass, TDPortal.Average.Retention.Time,
           TDPortal.Start.Index, TDPortal.End.Index, TDPortal.Modifications,
           TDPortal.Modification.Codes, TDPortal.E.value)) %>%
  filter(!is.na(TDPortal.Accession)) %>%
  mutate(
    "Unique.ID" = paste0(TDPortal.Accession, "_", TDPortal.Monoisotopic.Mass),
    "Mass" = TDPortal.Monoisotopic.Mass,
    "Retention.Time" =  TDPortal.Average.Retention.Time,
    "Proteoform.Details" = paste0(TDPortal.Accession, " ", " (", 
      TDPortal.Start.Index, "-", TDPortal.End.Index, ") ", TDPortal.Modifications,
      " ", TDPortal.Modification.Codes),
    "E.Value" = TDPortal.E.value
  ) %>%
  select(Unique.ID, Mass, Retention.Time, Proteoform.Details, E.Value) %>%
  group_by(Unique.ID) %>%
  slice(which.min(E.Value)) %>%
  arrange(Mass)

################################################################
## STEP 2: Zero-filling new intensities in ProMexAlign output ##
################################################################

# REASONING:
#-------------------------------------------------------------------------------
# ProMexAlign will try to identify features that were not present in the original
# MS1 Feature Tables (MS1FTs). These intensities may be noise. In this approach, 
# we decided to remove them. 

# Remove intensities identified by ProMexAlign
cleanData <- function(filepath, sample) {
  
  # Read file 
  File <- read.table(filepath, header = T)
  
  # Rename columns
  colnames(File) <- colnames(File) %>% gsub(pattern = paste0(sample, "_"), replacement = "")
  
  # Create an identification matrix. Any non-zero is converted to 1.
  Identified <- File[,22:30]
  Identified[Identified > 0] <- 1
  
  # Multiply identification matrix to the abundance matrix
  Abundance <- File[,4:12] * Identified
  
  # Remove abundance from the column names
  colnames(Abundance) <- colnames(Abundance) %>% gsub(pattern = "_Abundance", replacement = "")
  
  # Return dataframe
  return(data.frame(File[,1:3], Abundance))
  
}

# Load the data 
RatBrain <- rbind(
  cleanData("~/Git_Repos/spatially-resolved-TDP/ProMexAlign_Proteoforms/Inputs/ProMexAlign/Rat_Brain_CV30_crosstab.tsv", "CV30"),
  cleanData("~/Git_Repos/spatially-resolved-TDP/ProMexAlign_Proteoforms/Inputs/ProMexAlign/Rat_Brain_CV40_crosstab.tsv", "CV40"),
  cleanData("~/Git_Repos/spatially-resolved-TDP/ProMexAlign_Proteoforms/Inputs/ProMexAlign/Rat_Brain_CV50_crosstab.tsv", "CV50") 
) %>%
  arrange(MonoMass) 

# Add average RT
RatBrain$AverageRT <- lapply(1:nrow(RatBrain), function(row) {
  mean(c(RatBrain$MinElutionTime[row], RatBrain$MaxElutionTime[row]))
}) %>% unlist()
  

#######################################################################
## STEP 3:  Merge with TopPIC and TDPortal Proteoform Annotations ##
#######################################################################

# REASONING:--------------------------------------------------------------------
# Proteoforms should be within +/-1 Da and +/-4 minutes RT. If there's more than
# one match, return the closest in M/Z. 

# Build a function to iterate through each proteoform table and identify
# matches to the RatBrain table

#' @param Proteoforms The cleaned TDPortal or TopPIC dataframes with redundancies removed.
#' @param ProMexTable The combined ProMexAlign tables. 
#' @param PPMWindow The ppm window for collapsing proteoforms. Default is 15 ppm.
#' @param RTWindow The retention time for collapsing proteoforms. Default is 2 min. 
FindMatches <- function(Proteoforms, ProMexTable, PPMWindow = 15, RTWindow = 2) {

  NoMatches <- c()
  
  Tolerance <- ProMexTable$MonoMass * (PPMWindow / 1e6)
  Upper <- ProMexTable$MonoMass + Tolerance
  Lower <- ProMexTable$MonoMass - Tolerance
  
  Matches <- lapply(1:nrow(Proteoforms), function(row) {
    
    # Filter down to options within the PPM and RT Window
    Matches <- which(Upper >= Proteoforms$Mass[row] - (PPMWindow/1e6) &
                     Lower <= Proteoforms$Mass[row] + (PPMWindow/1e6) &
                     abs(ProMexTable$AverageRT - Proteoforms$Retention.Time[row]) <= RTWindow)
    
    # If no matches, return message
    if (length(Matches) == 0) {
      NoMatches <<- c(NoMatches, Proteoforms$Unique.ID[row])
      return(NA)
    } else if (length(Matches) > 1) {
      Matches <- Matches[which.min(abs(ProMexTable[Matches, "MonoMass"] - Proteoforms$Mass[row]))]
    } else if (length(Matches) > 1) {
      warning("There are still multiple matches that need to be addressed.")
    } 
    
    return(Matches)
    
  }) %>% unlist()
  
  return(list(Matches, NoMatches))
    
}

# Get matches
TopPIC_Matches <- FindMatches(TopPIC, RatBrain)
TDPortal_Matches <- FindMatches(TDPortal, RatBrain)

# Add rows to ProMexAlign table
TopPIC_Cols <- paste0("TopPIC.", colnames(TopPIC))
TDPortal_Cols <- paste0("TDPortal.", colnames(TDPortal))
for (col in c(TopPIC_Cols, TDPortal_Cols)) {
  RatBrain[[col]] <- NA
}

# Add TopPIC data
Out <- lapply(1:nrow(TopPIC), function(row) {
  
  # Add to the proteoform entry so that no information replaces
  Match <- TopPIC_Matches[[1]][row]
  
  if (!is.na(Match)) {
    Current <- RatBrain[Match, TopPIC_Cols]
    Current <- ifelse(is.na(Current), "", Current)
    NewEntry <- ifelse(Current == "", TopPIC[row,], paste(unlist(Current), TopPIC[row,], sep = " + "))
    RatBrain[Match, TopPIC_Cols] <<- NewEntry
  }
  
})

# Add TDPortal Data
Out <- lapply(1:nrow(TDPortal), function(row) {
  
  # Add to the proteoform entry so that no information replaces
  Match <- TDPortal_Matches[[1]][row]
  
  if (!is.na(Match)) {
    Current <- RatBrain[Match, TDPortal_Cols]
    Current <- ifelse(is.na(Current), "", Current)
    NewEntry <- ifelse(Current == "", TDPortal[row,], paste(unlist(Current), TDPortal[row,], sep = " + "))
    RatBrain[Match, TDPortal_Cols] <<- NewEntry
  }
  
})
rm(Out)

# Rearrange rat brain rows 
RatBrain <- RatBrain[,c(1:3, 14:23, 4:12, 13)]

######################################
## STEP 4: Group redundant features ##
######################################

# REASONING:
#-------------------------------------------------------------------------------
# There are redundant features within MS1FT files and across the CVs that need
# to be collapsed down to +/- 1 Da, +/- 4 Mean Retention Time groups. This is 
# not a simple binning approach, rather a 3 step grouping approach. 
# 1. Arrange by lowest to highest mass. Determine if each subsequent mass falls 
# within +/- 1 Da of the previous mass. Assign an ID for each "mass group"
# 2. Then, per "mass group", sort by mean retention time and determine if each 
# subsequent mean retention time is within +/- 4 minutes of the previous mean retention time. 
# Assign an ID for each of these "RT groups" within the "mass groups."
# 3. Combine mass groups and RT groups to generate a unique group ID. 
# Note that this only applies to cases where there are no proteoform annotations. 

# Split RatBrain data into proteoform identified, and not 
RatBrain_NA <- RatBrain[is.na(RatBrain$TopPIC.Unique.ID) & is.na(RatBrain$TDPortal.Unique.ID),]
RatBrain_Proteoform <- RatBrain[!is.na(RatBrain$TopPIC.Unique.ID) | !is.na(RatBrain$TDPortal.Unique.ID),]

# Define a collapse function for each subset of data 
large_collapse <- function(x, proteoform_mode, PPMWindow = 15, RTWindow = 2) {

  # Determine group numbers based on whether each subsequent mass falls within a
  # +/- 1 Da range of the previous mass.
  MassGroups <- x %>%
    mutate(
      Tolerance = MonoMass * (PPMWindow / 1e6),
      LowMass = MonoMass - Tolerance,
      HighMass = MonoMass + Tolerance,
      MassFlag = MonoMass >= lag(LowMass, default = first(LowMass)) &
        MonoMass <= lag(HighMass, default = first(HighMass))
    )
  
  # The first mass will always be TRUE so set to FALSE
  MassGroups$MassFlag[1] <- FALSE
  
  # Create new IDs by for each "FALSE" mass flag
  ID <- 0
  MassGroups$MassID <- lapply(1:nrow(MassGroups), function(row) {
    if (MassGroups$MassFlag[row] == FALSE) {ID <<- ID + 1}
    return(ID)
  }) %>% unlist()
  
  # Remove Low and High Mass
  MassGroups <- MassGroups %>% select(-c(LowMass, HighMass))
  
  # Group by the mass ID, sort retention times within a group, 
  # and assign a retention time within a -/+4 retention time window. 
  ID <- 0
  RTGroups <- do.call(rbind, lapply(1:max(MassGroups$MassID), function(ID) {
    
    sub <- MassGroups[MassGroups$MassID == ID,]
    
    if (nrow(sub) == 1) {
      
      ID <<- ID + 1
      sub$RTFlag <- FALSE
      sub$RTID <- ID
      return(sub)
      
    } else {
      
      RTGroups <- sub %>% 
        arrange(AverageRT) %>%
        mutate(
          LowRT = AverageRT - RTWindow,
          HighRT = AverageRT + RTWindow,
          RTFlag = AverageRT >= lag(LowRT, default = first(LowRT)) &
            AverageRT <= lag(HighRT, default = first(HighRT))
        )
      
      RTGroups$RTFlag[1] <- FALSE
      
      RTGroups$RTID <- lapply(1:nrow(RTGroups), function(row) {
        if (RTGroups$RTFlag[row] == FALSE) {ID <<- ID + 1}
        return(ID)
      }) %>% unlist()
      
      return(RTGroups %>% select(-c(LowRT, HighRT)))
      
    }
    
  })) %>% select(-c(MassFlag, RTFlag))
  
  # Generate a new ID and return maximum mass, mean min and max RT, and summed 
  # intensity
  median2 <- function(x) {
    x <- x[!is.na(x)]
    median(x, na.rm = T)
  }
  
  # Generate function to collapse non-NA values
  collapse_fun <- function(x, y) {
    vals <- unlist(x[[y]])
    vals[!is.na(vals)] %>% paste(collapse = " & ")
  }
  
  if (proteoform_mode == FALSE) {
  
    # Get final groupings for NA data 
    Final <- RTGroups %>%
      mutate(
        ID = factor(paste0(MassID, "_", RTID), 
                    levels = unique(paste0(MassID, "_", RTID))) %>% as.numeric()
      ) %>% 
      relocate(ID) %>%
      select(-c(MassID, RTID)) %>% 
      group_by(ID) %>%
      tidyr::nest() %>%
      mutate(
        NumFeatures = purrr::map(data, function(x) {nrow(x)}) %>% unlist(),
        MonoMass = purrr::map(data, function(x) {x$MonoMass %>% paste0(collapse = ", ")}) %>% unlist(),
        MinElutionTime = purrr::map(data, function(x) {
          x$MinElutionTime %>% min() %>% round(3)}) %>% unlist(),
        MaxElutionTime = purrr::map(data, function(x) {
          x$MaxElutionTime %>% max() %>% round(3)}) %>% unlist(),
        H3_C1 = purrr::map(data, function(x) {sum(x$H3_C1 %>% unlist())}) %>% unlist(),
        H3_C2 = purrr::map(data, function(x) {sum(x$H3_C2 %>% unlist())}) %>% unlist(),
        H3_C3 = purrr::map(data, function(x) {sum(x$H3_C3 %>% unlist())}) %>% unlist(),
        H3_C4 = purrr::map(data, function(x) {sum(x$H3_C4 %>% unlist())}) %>% unlist(),
        H3_C5 = purrr::map(data, function(x) {sum(x$H3_C5 %>% unlist())}) %>% unlist(),
        H4_A1 = purrr::map(data, function(x) {sum(x$H4_A1 %>% unlist())}) %>% unlist(),
        H4_A2 = purrr::map(data, function(x) {sum(x$H4_A2 %>% unlist())}) %>% unlist(),
        H4_A3 = purrr::map(data, function(x) {sum(x$H4_A3 %>% unlist())}) %>% unlist(),
        H4_A4 = purrr::map(data, function(x) {sum(x$H4_A4 %>% unlist())}) %>% unlist()
      ) %>%
      select(-data) %>%
      arrange(MonoMass)
    
  } else {
    
    # Get final groupings for NA data 
    Final <- RTGroups %>%
      mutate(
        ID = factor(paste0(MassID, "_", RTID), 
                    levels = unique(paste0(MassID, "_", RTID))) %>% as.numeric()
      ) %>% 
      relocate(ID) %>%
      select(-c(MassID, RTID)) %>% 
      group_by(ID) %>%
      tidyr::nest() %>%
      mutate(
        NumFeatures = purrr::map(data, function(x) {nrow(x)}) %>% unlist(),
        MonoMass = purrr::map(data, function(x) {
          if (nrow(x) == 1) {x$MonoMass} else {
            choose <- x %>% mutate(theSum = H3_C1 + H3_C2 + H3_C3 + H3_C4 + H3_C5 + H4_A1 + H4_A2 + H4_A3 + H4_A4) %>% 
              dplyr::select(theSum) %>% unlist() %>% which.max()
            return(x$MonoMass[choose])
          }
        }) %>% unlist(),
        MinElutionTime = purrr::map(data, function(x) {
          x$MinElutionTime %>% min() %>% round(3)}) %>% unlist(),
        MaxElutionTime = purrr::map(data, function(x) {
          x$MaxElutionTime %>% max() %>% round(3)}) %>% unlist(),
        TopPIC.Unique.ID = purrr::map2(data, "TopPIC.Unique.ID", collapse_fun) %>% unlist(),
        TopPIC.Mass = purrr::map2(data, "TopPIC.Mass", collapse_fun) %>% unlist(),                 
        TopPIC.Retention.Time = purrr::map2(data, "TopPIC.Retention.Time", collapse_fun) %>% unlist(),      
        TopPIC.Proteoform.Details = purrr::map2(data, "TopPIC.Proteoform.Details", collapse_fun) %>% unlist(), 
        TopPIC.E.Value = purrr::map2(data, "TopPIC.E.Value", collapse_fun) %>% unlist(),              
        TDPortal.Unique.ID = purrr::map2(data, "TDPortal.Unique.ID", collapse_fun) %>% unlist(),         
        TDPortal.Mass = purrr::map2(data, "TDPortal.Mass", collapse_fun) %>% unlist(),
        TDPortal.Retention.Time = purrr::map2(data, "TDPortal.Retention.Time", collapse_fun) %>% unlist(),     
        TDPortal.Proteoform.Details = purrr::map2(data, "TDPortal.Proteoform.Details", collapse_fun) %>% unlist(),
        TDPortal.E.Value = purrr::map2(data, "TDPortal.E.Value", collapse_fun) %>% unlist(),
        H3_C1 = purrr::map(data, function(x) {sum(x$H3_C1 %>% unlist())}) %>% unlist(),
        H3_C2 = purrr::map(data, function(x) {sum(x$H3_C2 %>% unlist())}) %>% unlist(),
        H3_C3 = purrr::map(data, function(x) {sum(x$H3_C3 %>% unlist())}) %>% unlist(),
        H3_C4 = purrr::map(data, function(x) {sum(x$H3_C4 %>% unlist())}) %>% unlist(),
        H3_C5 = purrr::map(data, function(x) {sum(x$H3_C5 %>% unlist())}) %>% unlist(),
        H4_A1 = purrr::map(data, function(x) {sum(x$H4_A1 %>% unlist())}) %>% unlist(),
        H4_A2 = purrr::map(data, function(x) {sum(x$H4_A2 %>% unlist())}) %>% unlist(),
        H4_A3 = purrr::map(data, function(x) {sum(x$H4_A3 %>% unlist())}) %>% unlist(),
        H4_A4 = purrr::map(data, function(x) {sum(x$H4_A4 %>% unlist())}) %>% unlist()
      ) %>%
      select(-data) %>%
      arrange(MonoMass)

  }
  
  return(Final)
    
}

Final_NA_Groups <- large_collapse(RatBrain_NA, FALSE)
Final_Proteoform_Groups <- large_collapse(RatBrain_Proteoform, TRUE)

############################################
## STEP 5: Identify and collapse isotopes ##
############################################

# REASONING:
#-------------------------------------------------------------------------------
# Isotopes within 1 or 2 Da (with a range established by the PPM window) should
# be collapsed down and not counted as separate protein/proteoform identifications.
# They should also be within the retention time window.

PPMWindow <- 15

Proteoform_Isotopes <- Proteoform_Groups %>%
  arrange(MonoMass) %>%
  ungroup() %>%
  mutate(
    Tolerance = MonoMass * (PPMWindow / 1e6),
    FirstIsoLow = MonoMass - Tolerance + 1,
    FirstIsoHigh = MonoMass + Tolerance + 1,
    SecondIsoLow = MonoMass - Tolerance + 2,
    SecondIsoHigh = MonoMass + Tolerance + 2,
    FirstIsoLow = lag(FirstIsoLow, default = first(FirstIsoLow)),
    FirstIsoHigh = lag(FirstIsoHigh, default = first(FirstIsoHigh)),
    SecondIsoLow = lag(SecondIsoLow, default = first(SecondIsoLow)),
    SecondIsoHigh = lag(SecondIsoHigh, default = first(SecondIsoHigh)),
    MinRT = lag(MinElutionTime, default = first(MinElutionTime)),
    MaxRT = lag(MaxElutionTime, default = first(MaxElutionTime)),
    MeanRT = purrr::map2(MinElutionTime, MaxElutionTime, function(x,y) {mean(x,y)}) %>% unlist(),
    IsotopeFlag1 = MonoMass >= FirstIsoLow & MonoMass <= FirstIsoHigh & MeanRT >= MinRT & MeanRT <= MaxRT,
    IsotopeFlag2 = MonoMass >= SecondIsoLow & MonoMass <= SecondIsoHigh & MeanRT >= MinRT & MeanRT <= MaxRT,
    ID = ifelse(IsotopeFlag1 | IsotopeFlag2, ID-1, ID)
  ) %>%
  dplyr::select(-c(FirstIsoLow, FirstIsoHigh, SecondIsoLow, SecondIsoHigh, MeanRT, MinRT, MaxRT, Tolerance)) %>%
  dplyr::relocate(IsotopeFlag1, IsotopeFlag2)


# Fill out NAs for unmatched features 
for (col in c(TopPIC_Cols, TDPortal_Cols)) {
  Final_NA_Groups[[col]] <- NA
}

# Merge final tables 
FinalTable <- rbind(
  Final_NA_Groups[,c(2:5, 15:24, 6:14)],
  Final_Proteoform_Groups[,c(2:24)]
) %>% 
  data.frame() %>%
  arrange(MonoMass)

write.table(FinalTable, "~/Git_Repos/spatially-resolved-TDP/ProMexAlign_Proteoforms/Outputs_and_Code/Proteoforms_and_Collapsed_Features.txt", 
            quote = F, row.names = F, sep = "\t", na = "")
  




  
