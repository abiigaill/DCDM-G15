sop <- read.csv("IMPC_SOP_2026.csv", header=T)

head(sop, 10) 

allFields <- sop[,1]
allMins <- sop[,3]
allMaxs <- sop[,4]

analysis_id <- allFields[1]
anaid_Min <- 12   
anaid_Max <- 12

gene_accession_id <- allFields[2]
geneacc_Min <- allMins[2]
geneacc_Max <- allMaxs[2]

gene_symbol <- allFields[3]
genesym_Min <- allMins[3]
genesym_Max <- allMaxs[3]

strain_name <- allFields[4]
strname_Min <- allMins[4]
strname_Max <- 26 

mouse_life_stage <- allFields[5]
mls_Min <- allMins[5]
mls_Max <- allMaxs[5]

parameter_id <- "parameter_stable_id"
paraid_Min <- allMins[6]
paraid_Max <- allMaxs[6]

parameter_name <- allFields[7]
paraname_Min <- allMins[7]
paraname_Max <- allMaxs[7]

pvalue <- "p_value"
pval_Min <- allMins[8]
pval_Max <- allMaxs[8]

allele_symbol <- allFields[9]
allsym_Min <- allMins[9]
allsym_Max <- allMaxs[9]

lacz <- allFields[10]
lacz_Min <- allMins[10]
lacz_Max <- 20 

library(tidyverse)

files <- list.files("C:/Users/Utente/Desktop/DCDM Coursework/Group15 data/group15_data/group15_data", pattern="\\.csv$", full.names=T)
files <- files[!startsWith(basename(files), "._")]

simple <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

read_file <- function(i) {
  lines <- readLines(i, warn=F)
  lines <- lines[lines != ""]
  lines <- gsub('"', "", lines)
  
  field <- sub(",.*", "", lines)
  value <- sub("^[^,]*,?", "", lines)
  value[value %in% c("", "NA")] <- NA
  
  data.frame(field, value)
}

template <- read_file(files[1])
correct_fields <- template[,1]

get_values <- function(i) {
  theFile <- read_file(i)
  pos <- match(simple(correct_fields), simple(theFile$field))
  theFile$value[pos]
}

flagged_files <- c()

for (i in files) {
  theFile = read_file(i)
  theField = theFile[,1]
  
  if (!identical(simple(theField), simple(correct_fields))) {
    flagged_files = c(flagged_files, basename(i))
  }
}

length(flagged_files)
flagged_files

dirty_DF <- map_dfr(files, function(i) {
  allValues <- get_values(i)
  data.frame(file=basename(i), setNames(as.list(allValues), correct_fields), check.names=F)
})

dirty_DF$flagged <- dirty_DF$file %in% flagged_files
write_csv(dirty_DF, "dirty_dataDF.csv")

strainNames <- c("B6Brd;B6Dnk;B6N-Tyr<c-Brd>", "C57BL/6Brd-Tyr<c-Brd>", "C57BL/6Dnk", "C57BL/6JTyr;C57BL/6N", "C57BL/6N", "C57BL/6NCrl", "C57BL/6NJ", "C57BL/6NTac", "C57BL/6NTacDen")
lifeStages <- c("E12.5", "E15.5", "E18.5", "E9.5", "Early adult", "Late adult", "Middle aged adult")
laczValues <- c("No data", "ambiguous", "expression", "imageOnly", "no expression", "tissue not available")

clean_data <- "C:/Users/Utente/Desktop/DCDM Coursework/Group15 data/clean_data"
dir.create(clean_data, showWarnings=F, recursive=T)

for (i in files) {
  theField = correct_fields
  allValues = trimws(get_values(i))
  
  pos = which(theField == strain_name)
  
  for (p in pos) {
    match_value = match(simple(allValues[p]), simple(strainNames))
    if (!is.na(match_value)) allValues[p] = strainNames[match_value]
  }
  
  pos = which(theField == mouse_life_stage)
  
  for (p in pos) {
    match_value = match(simple(allValues[p]), simple(lifeStages))
    if (!is.na(match_value)) allValues[p] = lifeStages[match_value]
  }
  
  pos = which(grepl("^lacz_", theField))
  
  for (p in pos) {
    match_value = match(simple(allValues[p]), simple(laczValues))
    if (!is.na(match_value)) allValues[p] = laczValues[match_value]
  }
  
  pos = which(theField == gene_symbol)
  
  for (p in pos) {
    allValues[p] = gsub("[^A-Za-z0-9]", "", allValues[p])
    allValues[p] = str_to_title(allValues[p])
  }
  
  pos = which(theField == parameter_id)
  
  for (p in pos) {
    allValues[p] = gsub("[!@#$%^&*+?=().[:space:]]", "", allValues[p])
  }
  
  gene_pos = which(theField == gene_symbol)
  pos = which(theField == allele_symbol)
  
  for (p in pos) {
    allele = allValues[p]
    
    if (!is.na(allele) && grepl("<", allele) && !is.na(allValues[gene_pos])) {
      allele = gsub("[!@#$%^&*+?=_[:space:]-]", "", allele)
      inner = sub("^[^<]*<", "", allele)
      inner = gsub("[<>]", "", inner)
      inner = gsub("[()]", "", inner)
      
      for (lab in c("KOMP", "EUCOMM", "NCOM")) {
        inner = sub(lab, paste0("(", lab, ")"), inner, ignore.case=T)
      }
      
      allValues[p] = paste0(allValues[gene_pos], "<", inner, ">")
    }
  }
  
  theFile = data.frame(theField, allValues)
  write.table(theFile, file=file.path(clean_data, basename(i)), sep=",", row.names=F, col.names=F, quote=T)
}

valid_brackets <- function(x) {
  if (is.na(x) || x == "") return(FALSE)
  
  old <- ""
  while (x != old) {
    old <- x
    x <- gsub("\\([^()<>]*\\)|<[^()<>]*>", "", x)
  }
  
  if (grepl("[<>()]", x)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

QC_passed <- c()
QC_failed <- c()

clean_files <- file.path(clean_data, basename(files))

for (i in clean_files) {
  
  theFile = read.csv(i, header=F)
  theField = theFile[,1]
  allValues = as.character(theFile[,2])
  
  QC_checks = 0
  
  pos = which(theField == analysis_id)
  
  if (length(pos) == 1 &&
      !is.na(allValues[pos]) &&
      nchar(allValues[pos]) >= anaid_Min &&
      nchar(allValues[pos]) <= anaid_Max) {
    QC_checks = QC_checks + 1
  }
  
  pos = which(theField == gene_accession_id)
  
  if (length(pos) == 1 &&
      !is.na(allValues[pos]) &&
      nchar(allValues[pos]) >= geneacc_Min &&
      nchar(allValues[pos]) <= geneacc_Max) {
    QC_checks = QC_checks + 1
  }
  
  pos = which(theField == gene_symbol)
  
  if (length(pos) == 1 &&
      !is.na(allValues[pos]) &&
      nchar(allValues[pos]) >= genesym_Min &&
      nchar(allValues[pos]) <= genesym_Max) {
    QC_checks = QC_checks + 1
  }
  
  pos = which(theField == strain_name)
  
  if (length(pos) == 1 &&
      !is.na(allValues[pos]) &&
      nchar(allValues[pos]) >= strname_Min &&
      nchar(allValues[pos]) <= strname_Max &&
      allValues[pos] %in% strainNames) {
    QC_checks = QC_checks + 1
  }
  
  pos = which(theField == mouse_life_stage)
  
  if (length(pos) == 1 &&
      !is.na(allValues[pos]) &&
      nchar(allValues[pos]) >= mls_Min &&
      nchar(allValues[pos]) <= mls_Max &&
      allValues[pos] %in% lifeStages) {
    QC_checks = QC_checks + 1
  }
  
  pos = which(theField == parameter_id)
  
  if (length(pos) == 1 &&
      !is.na(allValues[pos]) &&
      nchar(allValues[pos]) >= paraid_Min &&
      nchar(allValues[pos]) <= paraid_Max) {
    QC_checks = QC_checks + 1
  }
  
  pos = which(theField == parameter_name)
  
  if (length(pos) == 1 &&
      !is.na(allValues[pos]) &&
      nchar(allValues[pos]) >= paraname_Min &&
      nchar(allValues[pos]) <= paraname_Max) {
    QC_checks = QC_checks + 1
  }
  
  pos = which(theField == pvalue)
  
  if (length(pos) == 1) {
    value = suppressWarnings(as.numeric(allValues[pos]))
    
    if (!is.na(value) &&
        value >= pval_Min &&
        value <= pval_Max) {
      QC_checks = QC_checks + 1
    }
  }
  
  pos = which(theField == allele_symbol)
  gene_pos = which(theField == gene_symbol)
  
  if (length(pos) == 1 && length(gene_pos) == 1) {
    allele = allValues[pos]
    gene = allValues[gene_pos]
    
    if (!is.na(allele) &&
        !is.na(gene) &&
        nchar(allele) >= allsym_Min &&
        nchar(allele) <= allsym_Max &&
        startsWith(allele, gene) &&
        valid_brackets(allele)) {
      QC_checks = QC_checks + 1
    }
  }
  
  pos = which(grepl("^lacz_", theField))
  lacz_valid = 0
  
  for (p in pos) {
    
    if (!is.na(allValues[p]) &&
        nchar(allValues[p]) >= lacz_Min &&
        nchar(allValues[p]) <= lacz_Max &&
        allValues[p] %in% laczValues) {
      lacz_valid = lacz_valid + 1
    }
  }
  
  if (length(pos) > 0 && lacz_valid == length(pos)) {
    QC_checks = QC_checks + 1
  }
  
  if (QC_checks == 10) {
    QC_passed = c(QC_passed, i)
  } else {
    QC_failed = c(QC_failed, i)
  }
}

QC_results <- data.frame(
  file = basename(c(QC_passed, QC_failed)),
  status = c(rep("PASS", length(QC_passed)),
             rep("FAIL", length(QC_failed)))
)

QC_results$status[QC_results$file %in% flagged_files] <- "FLAGGED"


clean_DF <- map_dfr(clean_files, ~ read_csv(.x, col_names=c("field", "value"), col_types="cc") %>%
                      mutate(file=basename(.x)))
clean_DF <- clean_DF %>% 
  pivot_wider(names_from=field, values_from=value) %>% 
  relocate(file)

clean_DF$status <- QC_results$status[match(clean_DF$file, QC_results$file)]

clean_DF <- clean_DF[order(rowSums(is.na(clean_DF))), ]
clean_DF <- clean_DF[is.na(clean_DF$analysis_id) | !duplicated(clean_DF$analysis_id), ]

write_csv(clean_DF, "clean_dataDF.csv")
write_csv(QC_results, "QC_results.csv")

bespoke_task <- read.table("bespokeTask.txt", sep="\t", header=F, fill=T)
bespoke_task <- bespoke_task[,-6]
colnames(bespoke_task) <- c("Human Marker Symbol", "Human Entrez Gene ID", "Mouse Marker Symbol", "MGI Marker Accession ID", "High-level Mammalian Phenotype ID")

chrom1 <- read.table("chrom1.txt", sep="\t", header=T)

disease_info <- read.table("Disease_information.txt", sep="\t", header=F, fill=T, comment.char="#", quote="")
colnames(disease_info) <- c("Human Gene Symbol", "Human Gene Name", "HGNC ID for human gene", "DO term name associated with human gene", "DO term ID associated with human gene", "Mouse genotype IDs associated with DO term", "Mouse Gene", "MGI Mouse Gene ID")


IMPC_parameter <- read.table("IMPC_parameter_description.txt", sep= ",", header=T, fill=T)
IMPC_procedure <- read.table("IMPC_procedure.txt", sep=",", header=T, fill=T)

library(tidyr)

disease_info <- separate_rows(disease_info,
                              `Mouse genotype IDs associated with DO term`, 
                              sep="\\|")

bespoke_task <- separate_rows(bespoke_task,
                              `High-level Mammalian Phenotype ID`,
                              sep=",")

bespoke_task[] <- lapply(bespoke_task, trimws)
chrom1[] <- lapply(chrom1, trimws)
disease_info[] <- lapply(disease_info, trimws)
IMPC_parameter[] <- lapply(IMPC_parameter, trimws)
IMPC_procedure[] <- lapply(IMPC_procedure, trimws)

bespoke_task[bespoke_task == "" | bespoke_task == "NA"] <- NA
chrom1[chrom1 == "" | chrom1 == "NA"] <- NA
disease_info[disease_info == "" | disease_info == "NA"] <- NA
IMPC_parameter[IMPC_parameter == "" | IMPC_parameter == "NA"] <- NA
IMPC_procedure[IMPC_procedure == "" | IMPC_procedure == "NA"] <- NA

bespoke_task <- unique(bespoke_task)
chrom1 <- unique(chrom1)
disease_info <- unique(disease_info)
IMPC_parameter <- unique(IMPC_parameter)
IMPC_procedure <- unique(IMPC_procedure)

write.csv(bespoke_task, "bespoke_task_clean.csv", row.names=F)
write.csv(chrom1, "chrom1_clean.csv", row.names=F)
write.csv(disease_info, "disease_info_clean.csv", row.names=F)
write.csv(IMPC_parameter, "IMPC_parameter_description_clean.csv", row.names=F)
write.csv(IMPC_procedure, "IMPC_procedure_clean.csv", row.names=F)

