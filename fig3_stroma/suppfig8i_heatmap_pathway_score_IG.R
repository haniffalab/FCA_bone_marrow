library(qusage,rlist,Category)
library(Category)
library(rlist)
library(pipeR)
library(stringr)

{#get full genelist from all go terms
    #This table contains enrichment terms from the cytoscape app "enricr" that predicts enriched genesets. Use the "network predcition" script to get this input
  table<- read.csv("/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/Stroma_bm_downs/downs_Fbm_stroma_comp_merged_3_summary/downs_Fbm_stroma_comp_merged_3_summary_summary_node.csv",stringsAsFactors = F)
  try =  as.character(table$EnrichmentMap_Formatted_name)
  try = gsub("[c(]|[)]","",try)
  try = gsub("[\\\n]","",try)
  #try = str_replace_all(try,"[:punct:]","")
  try = (strsplit(try,","))
  names(try) <- table$name
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
table<- read.csv("/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/Stroma_bm_downs/downs_Fbm_stroma_comp_merged_3_summary/downs_Fbm_stroma_comp_merged_3_summary_summary_node.csv",stringsAsFactors = F)
try =  as.character(table$EnrichmentMap_Genes)
try = gsub("[c(]|[)]","",try)
try = gsub("[\\\n]","",try)
try = str_replace_all(try,"[:punct:]","")
try = (strsplit(try," "))
names(try) <- table$name
names(try)



##group the list
response = try[c("Response to type I interferon signalling","Response to interleukin 6,7,12","tumor necrosis factor mediated pathway","interleukin-1 response","response to interferon gamma")]
production = try[c("tumor necrosis factor production","leukocyte chemotaxis","interferon gamma production","interleukin 1,2 secretion","monocyte and neutrophil chemotaxis")]
response_score_results<- data.frame(row.names = names(response))
production_score_results<- data.frame(row.names = names(production))


#Get DE file 
de = read.csv("~/Downloads/DEGS_strom_endo_downs.csv")
de_name = "endo"
de = read.csv("~/Downloads/DEGS_strom_osteo_downs.csv")
de_name = "osteo"
de = read.csv("~/Downloads/DEGS_strom_mac_downs.csv",stringsAsFactors = F)
de_name = "mac"

{
de <- de[5:7]
colnames(de)<-c("genes","p","fc")


{#start of loop
production_score_results[de_name] <- NA
for (i in names(production)){
  print(i)
  temp <- de[(de$genes%in%unlist(production[i])),]
  temp_med <- median(as.numeric(temp$fc))
  production_score_results[de_name][(rownames(production_score_results)==i),]<-temp_med
}

response_score_results[de_name] <- NA
for (i in names(response)){
  print(i)
  temp <- de[(de$genes%in%unlist(response[i])),]
  temp_med <- median(as.numeric(temp$fc))
  response_score_results[de_name][(rownames(response_score_results)==i),]<-temp_med
}
}#End of loop
}
write.csv(production_score_results,"/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/Stroma_bm_downs/production_score_results.csv")
write.csv(response_score_results,"/Users/issac/Documents/Projects/Fetal Bone Marrow/network_gsea/Stroma_bm_downs/response_score_results.csv")

rbind(production_score_results,response_score_results)



##Plot heatmap
library(ggplot2)
library(hrbrthemes)

#plot response
library(tidyverse)
response_map <- response_score_results %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
head(response_map)

ggplot(response_map, aes(rowname, colname, fill= value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = median(response_map$value)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
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
