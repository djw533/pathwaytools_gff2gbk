library(reticulate) # so that python can be run from within R
library(stringr)
library(tidyr)
library(dplyr)


## first get the directories:
#python_ex <- system("which python3", intern = T)
#use_python(python_ex)

dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)


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
          file = "metabolic_pathways_rxns_with_genes_all_serratia_only_gene_evidence.tsv",
          sep="\t",
          quote = F,
          row.names = F)

### look at how common the pathways are (the ones that are true):

pathways_and_genomes <- unique(subset.data.frame(all_genomes_dataframe, select = c("Pathway.Frame.id","genome")))


#plot histogram
library(ggplot2)
library(ggtree)
library(phytools)

pathways_and_genomes <- within(pathways_and_genomes, Pathway.Frame.id <- factor(Pathway.Frame.id, levels=names(sort(table(Pathway.Frame.id),  decreasing=TRUE))))

ggplot(pathways_and_genomes, aes(Pathway.Frame.id)) + geom_bar() +
  scale_x_discrete()

#### plot it with the tree:


# input_tree <- args[1]

input_tree <- "~/Documents/202002-Feb/20200214-BAPS_clustering_panaroo_serratia_alignment/panaroo_serratia.newick"

strains_tree <- read.tree(input_tree)

###need to change the tip labels to reflect changes:
strains_tree$tip.label <- gsub("#", "_",strains_tree$tip.label)
strains_tree$tip.label <- gsub("[.]", "_",strains_tree$tip.label)

# if (is.null(strains_tree)) {
#   stop("Please provide a treefile for plotting", call.=FALSE)
# }
# 
# if (length(args)== 2) {
#   rooted_strains_tree <- root(strains_tree, outgroup = gsub("#","_",args[2]))
# } else if (length(args) == 1) {
#   rooted_strains_tree <- midpoint.root(strains_tree)
# }

rooted_strains_tree <- midpoint.root(strains_tree)
g1 <- ggtree(rooted_strains_tree, ladderize = T, right = T)


####addd serratia colouring:

serratia_metadata <- read.csv("~/Documents/thesis/chapter_1/5_All_serratia_tr/tree_metadata.csv", header = T,
                              quote = "",
                              stringsAsFactors = F,
                              comment.char = "")

###set groups
genus_groups <- list()
i <- 0
###################### Loop to set species groups in tree
for (species_cluster in as.data.frame(table(unlist(serratia_metadata$Cluster.y)))$Var1){
  i <- i + 1
  new_list <- c(subset.data.frame(serratia_metadata, Cluster.y == species_cluster)$File_prefix)
  print(species_cluster)
  genus_groups[[species_cluster]] <- new_list
}


genus_groups_tree <- groupOTU(rooted_strains_tree, genus_groups) #, overlap='abandon',connect = T)
tree_cols <- c("TssA" = "#3cb44b",
               "TssB" = "#ffe119",
               "TssC" = "#e6194b",
               "TssD" = "#4363d8",
               "TssE" = "#ff1493",
               "TssF" = "#911eb4",
               "TssG" = "#46f0f0",
               "TssH" = "#f032e6",
               "TssI" = "#bcf60c",
               "TssJ" = "#fabebe",
               "TssK" = "#008080",
               "TssL" = "#e6beff",
               "TssM" = "#9a6324",
               "Non-model" = "#ffffff",
               "DUF4150" = "#808000",
               "FHA" = "red",
               "ImpE" = "blue",
               "PAAR_motif" = "pink",
               "Pkinase" = "green",
               "PP2C" = "purple",
               "TagF_N" = "yellow",
               "black",
               "0" = "black",
               "i1" = "#3cb44b",
               "i2" = "#ffe119",
               "i3" = "#e6194b",
               "Other_i3" = "#76b7b2",
               "i3_A" = "#ff9da7",
               "i3_B" = "#f28e2b",
               "i4a" = "#4363d8",
               "i4b" = "#ff1493",
               "i5" = "#911eb4",
               "black",
               "14" = "black",
               "0" = "black",
               "ref" = "#4C4C4C")


### plot tree with the groups ####

genus_t1 <- ggtree(genus_groups_tree, ladderize = T, right = T, size=0.5, aes(color=group)) +
  scale_color_manual(values = c(tree_cols))  +
  geom_treescale(x=0,y=0,width = 0.05, offset = 4) +
  geom_rootedge(rootedge = 0.04) + 
  theme(legend.position = "right")

### separate ut so can use gheatmap with ggtree:

df2 <- pathways_and_genomes %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Pathway.Frame.id, values_from = present)

df3 <- as.data.frame(df2)

##set colnames:
rownames(df3) <- df3$genome
drops = "genome"
df3 <- df3[ , !(names(df3) %in% drops)]

### now sort colnames so that they are in descending order of prevalence:
orders <- subset.data.frame(as.data.frame(table(unlist(pathways_and_genomes$Pathway.Frame.id))), Freq > 0)

orders_to_sort <- as.character(unlist(orders$Var1))

df4<- df3[, orders_to_sort]

#gheatmap:
genus_t1 %>% gheatmap(df4, color = NULL,
                colnames_position = "top",
                colnames_offset_y = 20,
                width = 2,
                colnames_angle = 0,
                offset = 0,
                hjust = 0.5)  +
  scale_color_manual(values = c("blue","yellow","orange","red","black","grey","pink","purple","cyan","green","brown"))

## facet plot:
#facet_plot(genus_t1, panel="Barplot", data=pathways_and_genomes, geom = geom_point,
#           mapping = aes(x=Pathway.Frame.id, group=genome), shape = '.')


