# Rationale: using supervised learning, we would like to distinguish between 30 cell types using the 198 features (antibodies) in CITE-seq data.
# Inspired by implementation of decision tree in March 2021 bioarchive paper: https://doi.org/10.1101/2021.03.18.435922

#----------
# 1. BUILD DECISION TREE (DT) 
#----------

#Install required libraries and load data
library(rpart) # version 4.1-15
library(rattle) # version 5.4.0
library(rpart.plot) # version 3.0.9
library(RColorBrewer) # version 1.1-2
library(caret) # version 6.0-86
library(rlang) # version 0.4.10
library(reshape2) # version 1.4.3 
library(scales) # version 1.1.0 

# Load the Ab matrix for the FBM MNC CITE-seq data in as 'train'. This is DSB normalised data, with each cell type subset to k=10 (to account for class imbalance). 
train <- read.csv('/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs_mk3/resources_for_pipelines/figs3a_fbm_mnc_citeseq_dsb_counts_train_20210421.csv', 
                  header=TRUE, sep=",", row.names="X")
train[1:10]

# Build a DT using rpart and name the obj 'mytree'.
mytree <- rpart(
  cell.labels ~ ., 
  data = train, 
  method = "class", 
  cp = -1, #-ve complexity ensures tree growth isn't terminated in favour of reduced complexity. can lead to overfitting, but can be pruned later
  minsplit = 2, # minumim n(obs) in a node before split attempted (this is pre-pruning). Set to =2 to ensure that tree growth is unrestricted
  xval = 10 #  This runs 10-fold cross-validation on my training data (as default)
)

# Have a look at the DT obj
mytree

# Have a look at a table of the DT model (shows xerror/xstd by complexity)
df <- as.data.frame(mytree$cptable)
colnames(df) <- c("CP", "nsplit", "relError", "xerror", "xstd")
df["CP_plus_xerror"] <- df["CP"] + df["xerror"]
df
#             CP nsplit  relError    xerror       xstd
#1   0.033437827      0 1.0000000 1.0344828 0.00000000
#2   0.032392894      3 0.8996865 0.9749216 0.01326506
#3   0.031347962      6 0.8025078 0.9529781 0.01534178
#4   0.028213166      7 0.7711599 0.8871473 0.01990190
#5   0.025078370     13 0.5987461 0.7993730 0.02386454
#6   0.021943574     14 0.5736677 0.7492163 0.02544907
#7   0.018808777     19 0.4639498 0.7178683 0.02624405
#8   0.015673981     21 0.4263323 0.6739812 0.02713442
#9   0.012539185     22 0.4106583 0.6645768 0.02729363
#10  0.009404389     24 0.3855799 0.6426332 0.02762385
#11  0.006269592     32 0.3103448 0.6206897 0.02789793
#12  0.005224660     48 0.2100313 0.6144201 0.02796616
#13  0.003134796     51 0.1943574 0.6175549 0.02793259
#14 -1.000000000    113 0.0000000 0.6175549 0.02793259

# Have a look at graph of the DT model (shows xerror/xstd by size of tree). Save this graph as pdf for context of later pruning
plotcp(mytree)

#----------
# 2. AUTOMATICALLY PRUNE DT AND VISUALISE FINAL TREE
#----------
# Step one: Automatically adjust pruning level based on: greatest complexity for which [xerror < min(xerror + xstd)] - step one
pruneLevel <- df["CP"][df["xerror"] < min(df["xerror"] + df["xstd"])]  # Define pruning equation                                                            # View nsplits which pass automatic pruning criteria
df <- df[df$CP %in% pruneLevel,]                                       # Subset df for complexities which pass pruning step one
df                                                                     # View complexities which pass pruning step one
# Step two: Select for top 3 pruning levels which have smallest(cp+xerror)
df <- df[order(df["CP_plus_xerror"]),][1:4,]
df
#              CP nsplit  relError    xerror       xstd CP_plus_xerror
# 14 -1.000000000    113 0.0000000 0.6332288 0.02774809     -0.3667712
# 11  0.006269592     32 0.3103448 0.6175549 0.02793259      0.6238245
# 13  0.003134796     51 0.1943574 0.6332288 0.02774809      0.6363636
# 12  0.005224660     48 0.2100313 0.6332288 0.02774809      0.6384535
# Step three: Select the pruning level with lowest nsplit with all classes present as terminal node in decision tree 
# Built decision trees for all +ve complexities in df above: nsplit=32 omits 5 classes. nsplit=48 contains all classes > auto_cp must be set at this threshold)
df <- as.data.frame(mytree$cptable)
colnames(df) <- c("CP", "nsplit", "relError", "xerror", "xstd")
auto_cp <- df["CP"][df["nsplit"]==48]  # Define pruning equation   

# Prune and visualise final tree
mytree <- prune(mytree, cp = auto_cp) 
printcp(mytree)
df <- as.data.frame(mytree$cptable)
colnames(df) <- c("CP", "nsplit", "relError", "xerror", "xstd")
df
#             CP nsplit  relError    xerror       xstd
# 1  0.033437827      0 1.0000000 1.0344828 0.00000000
# 2  0.032392894      3 0.8996865 0.9780564 0.01293201
# 3  0.031347962      6 0.8025078 0.9592476 0.01478831
# 4  0.028213166      7 0.7711599 0.8620690 0.02122267
# 5  0.025078370     13 0.5987461 0.7868339 0.02429980
# 6  0.021943574     14 0.5736677 0.7711599 0.02480616
# 7  0.018808777     19 0.4639498 0.7115987 0.02638663
# 8  0.015673981     21 0.4263323 0.6739812 0.02713442
# 9  0.012539185     22 0.4106583 0.6677116 0.02724175
# 10 0.009404389     24 0.3855799 0.6489028 0.02753532
# 11 0.006269592     32 0.3103448 0.6394984 0.02766640
# 12 0.005224660     48 0.2100313 0.6426332 0.02762385

# Save the automatically pruned DT plot as A4 landscape pdf
rpart.plot(mytree, type = 3, clip.right.labs = FALSE, 
           branch = .3, under = FALSE, extra=FALSE, box.palette=0)

#----------
# 4. ASSESS DT FEATURE IMPORTANCE
#----------

# Visually inspect the feature (Ab) importance in the DT model with final complexity we've selected (cp = 0.0096774)
# Make dataframe of Abs and their importance
df <- as.data.frame(mytree$variable.importance)
df["ab"] = rownames(df)

# Plot the Abs (ordered by importance) as a barplot
ggplot(df, aes(x=reorder(ab, mytree$variable.importance), y=mytree$variable.importance)) + # 
  geom_bar(stat = "identity", width=0.5)  +
  theme(axis.text.y = element_text(hjust=1, size=3)) +
  coord_flip()

#----------
# 4. TEST DT MODEL
#----------

# Time to evaluate the DT model for how generalisable it is (incl. purity and recall) using a test dataset

# Load the test data in
test_data <- read.csv('/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs_mk3/resources_for_pipelines/figs3a_fbm_mnc_citeseq_dsb_counts_test_20210421.csv', 
                  header=TRUE, sep=",", row.names="X")
test_data[1:5]

# Run class predictions for the test data
preds <- factor(predict(mytree, newdata=test_data, type="class"), levels=levels(test_data$cell.labels))

# reorder the levels for preds (so that the confusion matrix looks nice)
celltype_list <- c( 'CD38- pro.', 'CD38+ pro.',
                   'eosinophil', 'basophil', 'mast cell', 
                   'early MK', 'MK', 'early erythroid', 'mid erythroid', 'late erythroid', 
                   'CD4 T cell',  'CD56 bright NK',  
                   'pre pro B progenitor', 'pro B progenitor', 'pre B progenitor', 'immature B cell',  'naive B cell',
                   'pDC', 'DC1', 'DC2', 'DC3', 
                   'MOP', 'promonocyte', 'CD14 monocyte', 'promyelocyte', 'neutrophil',
                   'stromal macrophage', 'osteoclast', 'sinusoidal EC',  'tip EC')

levels(preds) <- celltype_list
levels(test_data$cell.labels) <- celltype_list
preds

#----------
# 5. SAVE DT STATS + CONFUSION MATRIX OUTPUT 
#----------

# To evaluate both true positives (sensitivity) and true negatives (specificity), I generate confusion matrix
c_mat <- confusionMatrix(data=preds, reference=test_data$cell.labels) # default beta=1 
confusion_matrix_df <- melt(data.frame(c_mat$table))
head(confusion_matrix_df)

# Plot heatmap of confusion matrix (save as pdf in R)
ggplot(confusion_matrix_df, aes(x=Prediction, y=Reference)) + 
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  scale_fill_gradient2(low = "white", high = muted("midnightblue"), midpoint = 0) + # colour
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1, size=8, face="bold"),
        plot.title = element_text(size=20, face="bold"),
        axis.text.y = element_text(size=8, face="bold")) + 
  ggtitle("Confusion matrix for CITE-seq decision tree (DT)") + 
  theme(legend.title=element_text(face="bold", size=8)) + 
  xlab("Test data predicted label based on DT") + ylab("Test data actual label") +
  labs(fill="")

# Save the overall DT statistics for use as a supplementary table
write.csv(c_mat$byClass, "/Users/b8058304/Documents/PhD_work/Coding/manuscript_figs_mk3/resources_for_pipelines/confusion_matrix_stats_20210421.csv")

