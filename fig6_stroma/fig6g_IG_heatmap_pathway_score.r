library(qusage,rlist,Category)
library(Category)
library(rlist)
library(pipeR)
library(stringr)


{
  #get full genelist from all go terms
  table<- read.csv("/Users/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/downs_top_10_pathways_summary/temp_merged_13_summary_summary_node.csv",stringsAsFactors = F)
  try =  as.character(table$EnrichmentMap..Formatted_name)
  try = gsub("[c(]|[)]","",try)
  try = gsub("[\\\n]","",try)
  #try = str_replace_all(try,"[:punct:]","")
  names(try) <- table$name
  try = (strsplit(try,","))
  names(try)
  library(biomaRt)
  for (z in names(try)){
    print(z)
    go_terms<-unlist(try[z])
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
    #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
    gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                       filters = 'go', values = go_terms, mart = ensembl)
    gene_data<-unique(gene.data$hgnc_symbol)
    try[z]<-list(gene_data)
  }
}

##Or just use the rncihment map genes
#get genes from enrichment map
table<- read.csv("/Users/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/downs_top_10_pathways_summary/temp_merged_13_summary_summary_node.csv",stringsAsFactors = F)
try =  as.character(table$EnrichmentMap..Genes)
try = gsub("[c(]|[)]","",try)
try = gsub("[\\\n]","",try)
try = str_replace_all(try,"[:punct:]","")
try = (strsplit(try," "))
names(try) <- table$name
names(try)



###group the list
#response = try[c("Response to type I interferon signalling","Response to interleukin 6,7,12","tumor necrosis factor mediated pathway","interleukin-1 response","response to interferon gamma")]
#production = try[c("tumor necrosis factor production","leukocyte chemotaxis","interferon gamma production","interleukin 1,2 secretion","monocyte and neutrophil chemotaxis")]
#response_score_results<- data.frame(row.names = names(response))

production_score_results<- data.frame(row.names = names(try))
production<-try
degs_dir<-"/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/DEGs"
#Get all DE files
for (z in list.files(degs_dir)){
  print(z)
  f_name<-paste0(degs_dir,"/",z)
  de = read.csv(f_name)
  de_name = gsub("DEGs_|_DEGs.csv","",z)
  colnames(de)<-c("X","p","fc","genes")

  {#start of loop
    production_score_results[de_name] <- NA
    for (i in names(production)){
      print(i)
      temp <- de[(de$genes%in%unlist(production[i])),]
      temp_med <- mean(as.numeric(temp$fc))
      production_score_results[de_name][(rownames(production_score_results)==i),]<-temp_med
    }
  }#End of loop
}
write.csv(production_score_results,"/Users/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/top_10_pathway_production_scores.csv")


##Plot heatmap
library(ggplot2)
library(hrbrthemes)

#Remove natural log scale
production_score_results<-exp(production_score_results)

#Scale the scores for each column
library(tidyverse)
production_score_results<-as.data.frame(scale(production_score_results))
#trial<-production_score_results %>% mutate(zscore = (B_lineage - mean(B_lineage))/sd(B_lineage))

#plot production
library(tidyverse)
production_map <- production_score_results %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(production_map)
colnames(production_map)<-c("pathway","cell_state","logfc")
production_map$logfc[production_map$logfc=="NaN"]<-0
production_map$logfc[is.na(production_map$logfc)]<-0

#Reorder factors:
production_map$pathway<-factor(production_map$pathway)
production_map$pathway<- ordered(production_map$pathway, levels = c("monocyte and lymphocyte chemotaxis", 
                                                                    "Interleukin 12 production", 
                                                                    "type 1 interferon signalling pathway",
                                                                    "wnt signalling pathway",
                                                                    "CD4 tcell differentiation",
                                                                    "response to oxygen radical",
                                                                    "G1 DNA damage checkpoint",
                                                                    "interleukin 1 secretion",
                                                                    "tumor necrosis factor production and response",
                                                                    "telomere lengthening via telomerase"))

ggplot(production_map, aes(pathway, cell_state, fill= logfc)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = median(production_map$logfc)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55) 




###############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu95av2.db")

library(qusage)
library(rlist,Category)
library(Category)
library(rlist)
library(pipeR)
library(biomaRt)

################################

##Get genesets for GO
database <- read.gmt("/Users/Downloads/gprofiler_hsapiens.name/hsapiens.GO:BP.name.gmt")
database

#Make it so that the temp database only contains GO terms for each pathway
#get full genelist from all go terms
table<- read.csv("/Users/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/downs_top_10_pathways_summary/temp_merged_13_summary_summary_node.csv",stringsAsFactors = F)
try =  as.character(table$EnrichmentMap..Formatted_name)
try = gsub("[c(]|[)]","",try)
try = gsub("[\\\n]","",try)
#try = str_replace_all(try,"[:punct:]","")
names(try) <- table$name
try = (strsplit(try,","))
names(try)

db = list()
for(p in names(try)){
  print(p)
  temp_db <- database[names(database)%in%unlist(try[p])]
  temp_db<-as.character(unlist(temp_db,recursive = TRUE))
  temp_db<-list(temp_db)
  names(temp_db)<-p
  db<-append(db,temp_db)
}
database <- db

##Prepare qeury from DEGS
production_score_results<- list()
production <- db
degs_dir<-"/Users/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/DEGs"
#Get all DE files
for (z in list.files(degs_dir)){
  print(z)
  f_name<-paste0(degs_dir,"/",z)
  de = read.csv(f_name)
  de_name = gsub("DEGs_|_DEGs.csv","",z)
  colnames(de)<-c("X","p","fc","genes")
  temp <- de[(de$genes%in%unlist(production)),]
  temp_med <- list(as.character(temp$genes))
  names(temp_med)<-z
  production_score_results<-append(production_score_results,temp_med)
}


##Main loop to assign enrichment score to each cellstate

#Initialise empty frame
results_concat <- data.frame(module=character(),tf=character(),assayed=character(), significant=character(), universe=character(), representation=character(), p=character(), odds=character(), expected=character(),stringsAsFactors=FALSE)
for (p in names(production_score_results)){
  print(p)
  sigsets <- Map(function(x, y) x[x %in% y], db, MoreArgs=list(y=as.character(unlist(production_score_results[p]))))
  universe <- unlist(db, use.names=FALSE)
  result <- hyperg(db, sigsets, universe)
  head(result)
  result$module <- p
  result$tf <- rownames(result)
  #result<-result[(result$p<0.05)&(result$p!=0),]
  row.names(result) <- NULL
  results_concat<- rbind(result,results_concat)
}

#results_export <- results_concat[, c("module","tf","assayed", "p", "significant", "universe", "representation", "odds", "expected")]
#results_export <- as.data.frame(t(apply(results_export, 1, unlist)))

#plot heatmap
production_map<-results_concat[,c("tf","odds","module")]
colnames(production_map)<-c("pathway","odds","cell_state")
production_map$odds<-unlist(production_map$odds)

#Reorder factors:
production_map$pathway<-factor(production_map$pathway)
production_map$pathway<- ordered(production_map$pathway, levels = c("monocyte and lymphocyte chemotaxis", 
                                                                    "Interleukin 12 production", 
                                                                    "type 1 interferon signalling pathway",
                                                                    "wnt signalling pathway",
                                                                    "CD4 tcell differentiation",
                                                                    "response to oxygen radical",
                                                                    "G1 DNA damage checkpoint",
                                                                    "interleukin 1 secretion",
                                                                    "tumor necrosis factor production and response",
                                                                    "telomere lengthening via telomerase"))

ggplot(production_map, aes(pathway, cell_state, fill= odds)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = median(production_map$odds)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55) 
