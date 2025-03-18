
#BIOL 358 BIOINFORMATICS LAB
#SPRING 2024. DATE: APRIL 18, 2024
#AUTHOR: MARTHA I. NATUKUNDA
#DATASET USED FOR ANALYSIS: READ COUNT DATA FROM NATUKUNDA ET AL., 2022 (https://link.springer.com/article/10.1186/s12864-021-08251-4) (ONLY ONE SORGHUM GENOTYPE - BTx623)


#*********************************************************************************************************************#
#START OF CODE FOR SECTION 1

#STEP 1: INSTALLING THE R PACKAGES NEEDED FOR DIFFERENTIAL GENE EXPRESSION ANALYSIS AND DRAWING A HEATMAP

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("tidyverse")
BiocManager::install("pheatmap")

#STEP 2: LOAD THE INSTALLLED PACKAGES
library(edgeR) #loading required library edgeR
library(tidyverse) #loading required library tidyverse
library(pheatmap) #loading required library pheatmap

#STEP 3: SET THE WORKING DIRECTORY

#Note: Change the directory to match the one on your computer before you run the next line
setwd("C:/Users/mibore/Desktop/BIOL 358 Spring 2024/4. Files for the Bioinformatics Lab")

#STEP 4: READ DATA INTO R AND CLEAN IT UP A BIT TO REMOVE GENES WITH NO DATA
data <- read.csv("BIOL358_Count_Table_BTx623_Control.csv", header = TRUE) #reading count data in to R
data2 <- data[apply(data[, 2:8], 1, sum) > 0, ] #tells R to only keep genes that have counts >0 in at least one sample in a file called data2
dim(data2) #checking dimensions for data2
head(data2) #checking the top rows and columns of data2

group <- factor(c(1, 1, 2, 2, 3, 3, 3)) #defining groups for each column in data2
group

y <- DGEList(counts = data2[,2:8], group = group)
dim(y)
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~ group) #contrast. Creating your design matrix. This design uses group 1 as a reference or intercept.
design
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)


#STEP 5: FITTING THE GENERALIZED LINEAR MODEL_USED BY edgeR######
fit <- glmFit(y, design)


#STEP 6:STARTING THE CANOPY LAYER COMPARISONS####

###COMPARISON 1: GROUP2 (middle layer) VERSUS GROUP1 (lower layer) COMPARISON##### GROUP2 is the middle layer (leaf 8) AND GROUP1 is the lower layer (leaf 5). This is the MvsL comparison.
lrt.2vs1 <- glmLRT(fit, coef = 2)
topTags(lrt.2vs1)

# find out the number of DE genes that have higher or lower expression when canopy layers were compared. 
summary(dt <- decideTestsDGE(lrt.2vs1))

##### Plot DE genes for GROUP2 VERSUS GROUP1 IN AN MA PLOT#####  
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]
plotSmear(lrt.2vs1, de.tags=DEnames)
abline(h=c(-1,1), col="blue")


###COMPARISON 2: GROUP3 VERSUS GROUP2 COMPARISON##### GROUP3 IS CPFL AND GROUP2 IS C8
###COMPARISON 2: GROUP3 (upper layer) VERSUS GROUP2 (middle layer) COMPARISON##### GROUP3 is the upper layer (pre-flag leaf at the top of the plant) AND GROUP2 is the middle layer (leaf 8). This is the UvsM comparison.
lrt.3vs2 <- glmLRT(fit, contrast = c(0, -1, 1))#contrast checks if the group3 mean is equal to that for group2
topTags(lrt.3vs2)

# find out the number of DE genes that have higher or lower expression when canopy layers were compared. 
summary(dt <- decideTestsDGE(lrt.3vs2))

##### Plot DE genes for GROUP3 VERSUS GROUP2 IN AN MA PLOT#####  
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]
plotSmear(lrt.3vs2, de.tags=DEnames)
abline(h=c(-1,1), col="blue")


###COMPARISON 3: GROUP3 VERSUS GROUP1 COMPARISON##### GROUP3 IS CPFL AND GROUP1 IS C5
###COMPARISON 3: GROUP3 (upper layer) VERSUS GROUP1 (lower layer) COMPARISON##### GROUP3 is the upper layer (pre-flag leaf at the top of the plant) AND GROUP1 is the lower layer (leaf 5). This is the UvsL comparison.
lrt.3vs1 <- glmLRT(fit, coef = 3)
topTags(lrt.3vs1)

# find out the number of DE genes that have higher or lower expression when canopy layers were compared. 
summary(dt <- decideTestsDGE(lrt.3vs1))

##### Plot DE genes for GROUP3 VERSUS GROUP1 IN AN MA PLOT#####  
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]
plotSmear(lrt.3vs1, de.tags=DEnames)
abline(h=c(-1,1), col="blue")



######*************END OF CODE FOR SECTION 1 i.e R CODE FOR DIFFERENTIAL GENE EXPRESSION ANALYSIS**********************


###START OF SECTION 2 CODE (CODE FOR DRAWING THE HEATMAP)

#NOTE: WE ARE USING A DATASET THAT IS PROVIDED IN THE FOLDER TO DRAW A HEATMAP
#Drawing a heatmap for genes common among all three comparisons for BTx623 (N = 392 genes)

#STEP 1: load library, read in data and look at top 10 values in dataset
library(pheatmap) #loading library
df = read.csv("BTx623_Common genes_For Heatmap.csv", header= TRUE) #reading in data
head(df) #looking at the top values in the dataset called "df" just read in

#STEP 2: Transfer only the numerical part of dataset (columns 2 to 4) into a matrix
df_num = as.matrix(df[,2:4])
df_num

#STEP 3: Create the heatmap. First define the values for the heatmap colors to ensure that they are well distributed and the 0 is in the middle of the color bar, then draw the heatmap.
brk1 = seq(from = -8.5, to = 0, by = 0.1)
brk2 = seq(0.1, 8.5, by = 0.1)
brk = c(brk1, brk2)
my_palette <- c(colorRampPalette(colors = c("blue", "white"))(n = length(brk1)-1),
                c(colorRampPalette(colors = c("white", "red"))(n = length(brk2)-1)))
pheatmap(df_num, color = my_palette, breaks = brk, scale = "none", 
         cluster_rows = T, cluster_cols = F, display_numbers = F, margin = c(5,5))


#NOTE: THE HEATMAP PLOT WILL SHOW UP ON THE LOWER RIGHT. SAVE IT AS A PDF USING THE EXPORT FUNCTION ON THE LOWER RIGHT. WHEN SAVING THE PLOT, YOU CAN USE THE DEFAULT DIMENSIONS PROVIDED BY R (5.13 x 3.65 inches) OR MODIFY AS YOU WISH.


#####END OF CODE FOR SECTION 2 i.e. DRAWING HEATMAP#####***************************************************************




