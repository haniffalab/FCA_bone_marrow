#Get all genes from a list og GO terms
library(qusage,rlist,Category)
library(Category)
library(rlist)
library(pipeR)
library(stringr)
library(biomaRt)
{
  #get full genelist from all go terms
  #table<- read.csv("/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/downs_top_10_pathways_summary/temp_merged_13_summary_summary_node.csv",stringsAsFactors = F)
  #try =  as.character(table$EnrichmentMap..Formatted_name)
  #try = gsub("[c(]|[)]","",try)
  #try = gsub("[\\\n]","",try)
  #try = str_replace_all(try,"[:punct:]","")
  #names(try) <- table$name
  #try = (strsplit(try,","))
  #names(try)
  #library(biomaRt)
  #for (z in names(try)){
  #  print(z)
  #  go_terms<-unlist(try[z])
  #  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
  #  #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
  #  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
  #                     filters = 'go', values = go_terms, mart = ensembl)
  #  gene_data<-unique(gene.data$hgnc_symbol)
  #  try[z]<-list(gene_data)
  #}
#}

#TNF response genes
go_term <- c("GO:0034612","GO:0033209","GO:0034612","GO:0071356","GO:0071847","GO:1903265")

#Hypoxia response genes
go_term<-c("GO:0071456","GO:0061418","GO:1900039","GO:0097411","GO:1990144","GO:1900037")

#Erythropoeitin response genes
go_term<-c("GO:0038162","GO:0004900","GO:0036018","GO:0038162","")
#go_term<-c("GO:0004900")

go_terms<-unlist(go_term)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go', values = go_terms, mart = ensembl)
gene_data<-unique(gene.data$hgnc_symbol)
gene.data
#Get DEG and gene_list intersect
DEG<-read.csv("/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/all_DEGS_concat.csv")

#DEG<-read.csv("/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/downs_all_pathways/DEGs/DEGs_HSC_DEGs.csv",stringsAsFactors = F)
intersect<-toupper(DEG$gene)[toupper(DEG$gene)%in%toupper(gene_data)]

DEG_intersect<-DEG[(DEG$gene%in%intersect),]
DEG_intersect<-DEG_intersect[!duplicated(DEG_intersect$gene), ]
#DEG_intersect$cell<-"HSC_downs"

DEG_intersect
##Plot heatmap
library(ggplot2)

##plot response
#library(tidyverse)
#response_map <- response_score_results %>%
#  rownames_to_column() %>%
#head(response_map)
##  gather(colname, value, -rowname)

#DEG_intersect$global<-"downs"

#Plot global
ggplot(DEG_intersect, aes(celltype, gene, fill= logfc)) + ggtitle("global erythropoietin response genes by DE") +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "purple", midpoint = median(DEG_intersect$logfc),limits = c((0), max(as.vector(DEG_intersect$logfc)))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55)

#write.csv(DEG_intersect$gene,"~/Downloads/hypoxia_de.csv")
DEG_intersect<-DEG_intersect[(DEG_intersect$celltype=="HSC"),]
ggplot(DEG_intersect, aes(celltype, gene, fill= logfc)) + ggtitle("HSC TNF response genes by DE") +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "purple", midpoint = median(DEG_intersect$logfc),limits = c((0), max(as.vector(DEG_intersect$logfc)))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55)

DEG_intersect<-DEG_intersect[(DEG_intersect$celltype=="HSC"),]
ggplot(DEG_intersect, aes(celltype, gene, fill= logfc)) + ggtitle("HSC TNF response genes by DE") +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white",high="red",limits = c((0), max(as.vector(DEG_intersect$logfc)))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55)


library(scales)
ggplot(DEG_intersect, aes((celltype %in% c("HSC")), gene, fill= logfc)) + ggtitle("Global TNF response genes by DE") +
  geom_tile(color = "white")+
  scale_fill_gradient2(high = "red",low="white",limits = c((0), max(as.vector(DEG_intersect$logfc))))+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55) 

midpoint = mean(DEG_intersect$logfc[DEG_intersect$celltype=="HSC"])

ggplot(DEG_intersect, aes(global, gene, fill= logfc)) + ggtitle("Global Hypoxia response genes by DE") +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "purple", midpoint = median(DEG_intersect$logfc)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55) 

#plot production
library(tidyverse)
production_map <- production_score_results %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(production_map)
colnames(production_map)<-c("pathway","cell_state","fc")

ggplot(production_map, aes(pathway, cell_state, fill= fc)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = median(response_map$value)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip() + theme(aspect.ratio = 0.55) 
