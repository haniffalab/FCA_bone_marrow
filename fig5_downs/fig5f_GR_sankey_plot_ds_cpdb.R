library(networkD3)
library(dplyr)
setwd("~/Desktop/bm_data")
###build ggplot heatmap used to derive colours
y<- paste0("var", seq(1:length(value)))
col_plot<-data.frame(value = value, place = as.character(y))
ggplot(col_plot, aes(place, value, fill= value)) + 
   geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.003)

###build object
links <- data.frame(
source=c("Myeloid_DC_GRN", "Neut_GRN",
         "Mono_TNF", "NK_TNF", "Neut_TNF", "Myeloid_DC_TNF",
         "B_LTA",
         "Neut_TNFSF13", "Myeloid_DC_TNFSF13", "Mono_TNFSF13",
         "B_LTA",
         "Neut_TNFSF13", "Myeloid_DC_TNFSF13", "Mono_TNFSF13",
         "Stroma_DLK1", "Myeloid_pro_DLK1"
         ),
target=c("TNFRSF1A", "TNFRSF1A",
         "TNFRSF1A", "TNFRSF1A","TNFRSF1A", "TNFRSF1A",
         "TNFRSF1A",
         "TNFRSF1A", "TNFRSF1A","TNFRSF1A",
         "TNFRSF14",
         "TNFRSF14", "TNFRSF14", "TNFRSF14",
         "NOTCH1", "NOTCH1"
         ), 
value=c(46, 53,
        43,78,51,78,
        30,
        48,40,33,
        30,
        48,40,33,
        133,33))
links$IDsource<-c(0,1,
                  2,3,4,5,
                  6,
                  7,8,9,
                  6,
                  7,8,9,
                  10,11)
links$IDtarget<-c(12,12,
                  12,12,12,12,
                  12,
                  12,12,12,
                  13,
                  13,13,13,
                  14,14)
nodes <- data.frame(
  name=c(as.character(links$source), 
  as.character(links$target)) %>% unique()
)
links$IDsource<-as.numeric(links$IDsource)
links$IDtarget<-as.numeric(links$IDtarget)
my_color <- 'd3.scaleOrdinal() 
.domain(["Myeloid_DC_GRN", "Neut_GRN",
"Mono_TNF", "NK_TNF", "Neut_TNF", "Myeloid_DC_TNF", "B_LTA",
"Neut_TNFSF13", "Myeloid_DC_TNFSF13", "Mono_TNFSF13", "Stroma_DLK1",  "Myeloid_pro_DLK1",
"TNFRSF1A", "TNFRSF14", "NOTCH1",
"type_a", "type_b", "type_c"]) 
.range(["#f9c2b0", "#f7b29e", 
"#f9c2b0", "#f08268", "#f7b29e", "#f08268",
"#fad1c4",
"#f7b29e", "#f9c2b0", "#fad1c4",
"#e52021", "#fad1c4",
"454c9b", "c5b6da", "786aae",
"#75a55e", "#7fccbe", "#cccc7a"])'
links$group <- as.factor(c("type_a","type_a","type_a",
                           "type_a","type_a","type_a",
                           "type_a","type_a","type_a","type_a",
                           "type_b","type_b","type_b", "type_b",
                           "type_c", "type_c"))
nodes$group <- as.factor(c("my_unique_group"))
p<-sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
             sinksRight=FALSE, colourScale=my_color, LinkGroup="group"
             )
p