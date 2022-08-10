#library(reticulate) # so that python can be run from within R
library(stringr)
library(tidyr)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please supply directory of mpwt output", call.=FALSE)
}

# check if args[1] is a directory??

## first get the directories:
#python_ex <- system("which python3", intern = T)
#use_python(python_ex)


arg <- "5_mpwt_output/"

dirs <- list.dirs(path = arg, full.names = TRUE, recursive = FALSE)


## create dataframe to rbind to at the end:
all_genomes_dataframe <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(all_genomes_dataframe) <- c("genome","REACTIONS.PRESENT","Pathway.Frame.id","Pathway.Name","Pathway.Class.Name","Pathway.Class.Frame.id",
                                     "Pathway.Score","Pathway.Frequency.Score","Pathway.Abundance","Reason.to.Keep","Pathway.URL",
                                     "KEEP.","gene")


for (dir in dirs) {
  print(dir)
  ### find the reports file:
  reports_path <- paste0(dir,"/1.0/reports",sep="")
  files <- list.files(path = reports_path, full.names = TRUE)

  ##file in pathways-report file:
  for (file in files) {
    if (grepl("pathways-report", file, fixed = TRUE)) {
      pathways_file <- file
    }
  }


  ## now we have the file - we can read it in:

  pathways <- read.table(pathways_file,
                         comment.char = "#",
                         sep="|",
                         header = T,
                         quote = "")



  #### other data files:


  #pwy-inference-description.data
  #has all the reactions present for the pathway - the reactions that are missing - and a list of the core reactions that are needed for that pathway
  #importantly - has instances where the pathways weren't kept as "present" because key reactions are

  #parse this using python
  #py_run_file("~/github/pathwaytools_gff2gbk/pathway_inference_parser.py --input pwy-inference-description.data --output pwy_inference_description.csv")
  input_pwf_inf <- paste0(dir,"/1.0/reports/pwy-inference-description.data",sep="")
  output_pwf_inf <- paste0(dir,"/1.0/reports/pwy_inference_description.csv",sep="")
  system(paste("python3 ~/github/pathwaytools_gff2gbk/pathway_inference_parser.py","--input",input_pwf_inf,"--output",output_pwf_inf,sep=" "))
  #now read the output csv file in:
  pathway_inf <- read.csv(output_pwf_inf,
                          header = T,
                          quote = "")

  #pwy-inference-report_2020-05-04.txt
  #has the pathwasy imported, those not imported initially, those imported initially, those pruned, and those kept
  #"determine-pathways-with-cf" - pruned according to confidence value??

  #name-matching-report_2020-05-04.txt
  #has the individual reactions linked to GO erms, EC-numbers, and the gene name that corresponds with that reaction

  ###again, parse this usng python:
  #first find the name file:
  ##file in pathways-report file:
  for (file in files) {
    if (grepl("name-matching-report", file, fixed = TRUE)) {
      names_file <- file
    }
  }

  #set output name
  output_gene_names <- paste0(dir,"/1.0/reports/gene_names_reactions.csv",sep="")
  #run script
  system(paste("python3 ~/github/pathwaytools_gff2gbk/name_matching_parser.py","--input",names_file,"--output",output_gene_names,sep=" "))

  ##read this in
  gene_names <- read.csv(output_gene_names,
                          header = T,
                          quote = "")

  #pwy-evidence-list.dat
  #list of all the inferred pathways and super-pathways, and the list of reactions for that pathway where there is evidence for that reaction (i.e. the enzyme exists)

  # so - step 1 was to read in the pathways from pathways_file (file containing "pathways-report")

  ##now:
  #step 2- get the number of reactions in each pathway and number not present into the df:

  # so take all of the pathways from the pathway_inf, and along with it, take the "REACTIONS-PRESENT":

  pathways_w_reactions_present <- subset.data.frame(pathway_inf, select=c("PATHWAY","REACTIONS.PRESENT","KEEP."))


  ###merge with pathways df for pathways in this pgdb:
  merged_pathways <- merge(pathways,pathways_w_reactions_present,by.x = "Pathway.Frame.id", by.y = "PATHWAY")

  ## go through each pathway - take out the pathways/ inflate and match against the different gene names

  ## oR - possibly inflate the df in the first place for the diff. reactions - and merge with all possible genes with those reactions?
  ## first make the reactions as strings:
  #merged_pathways$REACTIONS.PRESENT <- toString(merged_pathways$REACTIONS.PRESENT)

  ## then split and use unnest?:
  unnested_pathways <- merged_pathways %>%
    mutate(REACTIONS.PRESENT = strsplit(as.character(REACTIONS.PRESENT), " ")) %>%
    unnest(REACTIONS.PRESENT)

  ## now merge in the gene names for those reagions:
  pathways_with_genes <- merge(unnested_pathways, gene_names, by.x = "REACTIONS.PRESENT", by.y = "reaction")
  pathways_with_genes$genome <- str_split_fixed(dir, '/',2)[2]

  ##now rbind:
  all_genomes_dataframe <- rbind(all_genomes_dataframe,pathways_with_genes)

}

####write out the csv for all of the information:

write.table(x = all_genomes_dataframe,
          file = "metabolic_pathways_rxns_with_genes_only_gene_evidence.tsv",
          sep="\t",
          quote = F,
          row.names = F)
