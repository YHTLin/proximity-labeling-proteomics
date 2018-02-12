### THIS SUITE OF CODE CLEANS AND REDUCES MAXQUANT OUTPUT ON PROTEIN GROUPS FOR DENIZ

# Step 1 - Set working directory and increase memory size for xlsx file writing
options(java.parameters = "- Xmx8000m")
setwd("C:/Users/Tony Lin/Desktop/Wiita Lab/Projects/Proteomics Project/BioID2-APEX-AGO2 Deniz/R script")


# Step 2 - Read in protein groups CSV data file
raw.data = read.csv("20170616 AGO2 BioID2 APEX Maxquant proteinGroups.csv", header = TRUE, colClasses = "character")


# Step 3 - Data clean-up
## Remove reverse proteins and contaminations
data = raw.data
data = data[data$Only.identified.by.site != "+", ]
data = data[data$Reverse != "+", ]
data = data[data$Potential.contaminant != "+", ]


# Step 4 - Data-specific manipulations
## Cast LFQ intensity columns as numeric type
LFQ.names = grep("^LFQ", colnames(data), value = TRUE)
data[, LFQ.names] = lapply(data[, LFQ.names], function(x) as.numeric(x))

## Calculate log2(LFQ intensity)
LOG.names = gsub("^LFQ.intensity", "LOG2", LFQ.names)
data[, LOG.names] = lapply(data[, LFQ.names], function(x) log2(x))

## Filter columns by names
filter_cols = function(df, exp_grp) {
  exp_grp.names = grep(exp_grp, colnames(df), value = TRUE)
  LFQ.names = grep("^LFQ", exp_grp.names, value = TRUE)
  IDtype.names = grep("^Identification.type", exp_grp.names, value = TRUE)
  LOG2.names = grep("^LOG2", exp_grp.names, value = TRUE)
  keepCols = c("Protein.IDs", "Protein.names", "Gene.names", "Number.of.proteins", 
             "Peptides", LOG2.names, LFQ.names, IDtype.names)
  return(df[, keepCols])
}
APEX.data = filter_cols(data, "APEX")
BirA.data = filter_cols(data, "BirA")

## Group filtering to remove -Inf values
NA_filter = function(df, groupList = list("AGO2", "ctrl"), 
                     min_number = 1, impute = TRUE) {
  # df = data frame containing log2 intensity data
  # groupList = a list of vectors containing replicate groups
  # min_number = filter out rows that do not meet the minimum number of valid values 
  #       in ALL groups
  
  #Count up all valid values in each group
  validSum = lapply(groupList,
                    function(group) {
                      group.names = grep(paste0("LOG2.*", group), 
                                         colnames(df), value = TRUE)
                      return(apply(data.frame(df[, group.names]), 1, 
                                          function(x) sum(x != -Inf)))
                    })
  #Mark rows that do not pass the min_number filter
  filter = apply(data.frame(validSum), 1, function(x) all(x >= min_number))
  
  return(df[filter, ])
}
APEX.data = NA_filter(APEX.data)
BirA.data = NA_filter(BirA.data)

## Impute missing data
impute_data = function(df, width = 0.3, downshift = 1.8) {
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  LOG2.names = grep("^LOG2", colnames(df), value = TRUE)
  impute.names = gsub("^LOG2", "impute", LOG2.names)
  
  # Create new columns indicating whether the values are imputed
  df[, impute.names] = lapply(LOG2.names, function(x) is.infinite(df[, x]))
  
  # Imputation
  set.seed(1)
  df[, LOG2.names] = lapply(LOG2.names, 
                            function(x) {
                              temp = df[, x]
                              temp[is.infinite(temp)] = NA
                              temp.sd = width * sd(temp, na.rm = TRUE)
                              temp.mean = mean(temp, na.rm = TRUE) - downshift * sd(temp, na.rm = TRUE)
                              n.missing = sum(is.na(temp))
                              temp[is.na(temp)] = sample(rnorm(n.missing, mean = temp.mean, sd = temp.sd),
                                                         n.missing, replace = FALSE)
                              return(temp)
                            })
  return(df)
}
APEX.data = impute_data(APEX.data)
BirA.data = impute_data(BirA.data)

## Calculate LOG2 difference (ONLY WORKS WHEN LENGTH OF groupList EQUALS 2)
LOG2_diff = function(df, groupList = list("AGO2", "ctrl"), identifier = "") {
  LOG2.names = grep("^LOG2", colnames(df), value = TRUE)
  mean.names = paste0("mean.LOG2.", groupList, identifier)
  df[, mean.names] = lapply(groupList, 
                            function(x) {
                              group.names = grep(x, LOG2.names, value = TRUE)
                              return(rowMeans(data.frame(df[, group.names])))
                            })
  df[, paste0("diff.mean.LOG2", identifier)] = df[, mean.names[1]] - df[, mean.names[2]]
  return(df)
}
APEX.data = LOG2_diff(APEX.data, identifier = ".APEX")
BirA.data = LOG2_diff(BirA.data, identifier = ".BirA")

## Merge APEX and BirA data to check for 
merge.data = merge(APEX.data, BirA.data)
mean.names = grep("mean", names(merge.data), value = TRUE)
diff.names = grep("diff.mean", mean.names, value = TRUE)
merge.data = cbind(merge.data[, 1:5], merge.data[, mean.names])
merge.data[, "bothHigh"] = apply(merge.data[, diff.names], 1, function(x) all(x > 1))
merge.data[, "bothLow"] = apply(merge.data[, diff.names], 1, function(x) all(x < -1))
merge.data[, "bothUp"] = apply(merge.data[, diff.names], 1, function(x) all(x > 0))
merge.data[, "bothDown"] = apply(merge.data[, diff.names], 1, function(x) all(x < 0))

## List of AGO2 interactors from UniProt
AGO2_interactors = c("DICER1", "TARBP2", "AGO1", "DDB1", "DDX5", "DDX6", "DDX20", "DHX30", 
                     "DHX36", "DDX47", "DHX9", "ELAVL", "FXR1", "GEMIN4", "HNRNPF", "IGF2BP1", 
                     "ILF3", "IMP8", "MATR3", "PABPC1", "PRMT5", "P4HA1", "P4HB", "RBM4", 
                     "SART3", "TNRC6A", "TNRC6B", "UPF1", "YBX1", "DCP1A", "XRN1", "LIMD1", 
                     "WTIP", "AJUBA", "TRIM71", "APOBEC3G", "APOBEC3A", "APOBEC3C", "APOBEC3F",
                     "APOBEC3H", "EIF6", "MOV10", "RPL7A", "FMR1", "ZFP36")
is_interactor = function(df) {
  AGO2.regex = paste0(rep("\\b", length(AGO2_interactors)), AGO2_interactors, 
                      rep("\\b", length(AGO2_interactors))) #search string must be preceded and followed by nonwords
  AGO2.index = sapply(AGO2.regex, grep, df$Gene.names)
  AGO2.index = unique(unlist(AGO2.index))
  df[, "is.interactor"] = 1:nrow(df) %in% AGO2.index
  return(df)
}
BirA.data = is_interactor(BirA.data)
APEX.data = is_interactor(APEX.data)
merge.data = is_interactor(merge.data)


# Step 5 - Output Excel file with processed data
if(FALSE) {
require(xlsx)
wb <- createWorkbook()

addDataFrame(raw.data, row.names = FALSE, 
             createSheet(wb, sheetName="raw_data"))
addDataFrame(data, row.names = FALSE, 
             createSheet(wb, sheetName="cleaned_data"))
addDataFrame(APEX.data, row.names = FALSE, 
             createSheet(wb, sheetName="APEX_data"))
addDataFrame(BirA.data, row.names = FALSE, 
             createSheet(wb, sheetName="BirA_data"))
addDataFrame(merge.data, row.names = FALSE, 
             createSheet(wb, sheetName="merged_data"))
addDataFrame(AGO2_interactors, row.names = FALSE, 
             createSheet(wb, sheetName="AGO2_interactors"))

saveWorkbook(wb, "AGO2_BioID2_APEX_MaxQuant_proteinGroups processed.xlsx")
}