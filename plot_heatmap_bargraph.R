# script to normalize count data from FeatureCounts and then plot a heatmap
library(dplyr)
setwd("C:/Users/")

### normalize count data

##import count data

input_filename <- "gene_counts_amoA.txt"
counts <- read.table(input_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
counts
#colnames(counts)[1] <- "Sample" #change column name of Chr to Sample
#change sample names to get rid of X infront
colnames(counts)[3] <- "1FNE"
colnames(counts)[4] <- "8FNE"
colnames(counts)[5] <- "1JNE"
colnames(counts)[6] <- "8JNE"
colnames(counts)[7] <- "1SNE"
colnames(counts)[8] <- "8SNE"
colnames(counts)

order_groups <- c( "RBC_Group_I", "RBC_Group_E", "RBC_Group_C", "RBC_Group_D", "RBC_Group_F", "RBC_Group_B", "RBC_Group_M", "RBC_Group_L", "RBC_Group_A", "RBC_Group_G", "RBC_Group_H", "RBC_Group_K", "RBC_Group_J")
#change row order in counts so it matches that in tree see https://stackoverflow.com/questions/26548495/reorder-rows-using-custom-order
counts <- counts %>% 
  slice(match(order_groups, Chr))
counts<- as.data.frame(counts)
counts

##normalize by length of amoA gene mapped to

samples <- counts[,c(3:16)] #get sample names
gene_length <- counts [,c(2)] #get gene/amoA length, save as a variable
norm_length <- counts %>%
  mutate_at(vars(3:16), funs(./gene_length)) #divide by length of gene mapped to 
norm_length <- norm_length [,c(1,3:16), drop=FALSE] 
# see https://stackoverflow.com/questions/28123098/multiply-many-columns-by-a-specific-other-column-in-r-with-data-table
norm_length

##import raw reads data
reads_filename <- "all_read_counts.tsv"
reads_file <- read.table(reads_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

reads <- reads_file %>%
  filter(Step=="QC") %>%
  select("Sample", "Step", "Total_Reads") #select and filter so only rows and columns of interest to normalize by

reads
order <- c( "1FNE", "8FNE", "1JNE", "8JNE", "1SNE", "8SNE", "NE1", "NE8", "NW1", "NW8", "SE2", "SE8", "SW1", "SW8")
#change row order in reads so it matches that in other table and other figures see https://stackoverflow.com/questions/26548495/reorder-rows-using-custom-order
reads <- reads %>% 
  slice(match(order, Sample))
reads <- as.data.frame(reads)

norm_length2 <- data.frame(t(norm_length[-1]))
colnames(norm_length2) <- norm_length[,1]
# flip columns and rows, see https://stackoverflow.com/questions/33643181/how-do-i-flip-rows-and-columns-in-r
norm_length2


row.names(norm_length2) <- reads$Sample #add sample column and lable row by sample name
norm_length2["Sample"] <- reads$Sample
norm_length2

join <- full_join(reads, norm_length2, by='Sample') #join the two tables together
join

##normalize by number of reads per sample

Tot_reads <- join [,c(3)] #get total read numbers (total reads column), save as variable
norm_by_reads <- join %>%
  mutate_at(vars(4:16), funs(./Tot_reads)) #divide by total number of reads per sample
norm_by_reads
norm_counts <- norm_by_reads %>%
  mutate_at(vars(4:16), funs(.*100000000000)) %>% #multiply by 100 billion to get non-decimal numbers
  select(c(1,4:16)) #select only columns needed for heatmap
norm_counts
norm_counts2 <- data.frame(t(norm_counts[-1])) #transpose rows and columns (switch them)
colnames(norm_counts2) <- norm_counts[,1]
norm_counts2



## heatmap part
#source("http://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
#install.packages("RColorBrewer")
library("RColorBrewer")
Heatmap(counts_matrix, cluster_rows = FALSE, cluster_columns = FALSE,  name = "counts", heatmap_legend_param = list(title_position = "topcenter", color_bar = "continuous"), col =brewer.pal(9,"Blues"), column_title = "RBC Samples", column_title_side = "bottom", row_title = "RBC Group", row_title_side = "right")

###make bargraph to summarize total mapped reads to Nitrospira amoA per sample
library(ggplot2)
## get total reads mapped to sample (add up all normalized counts)
tot_counts <- colSums(norm_counts2)
tot_counts
tot_counts <- data.frame(tot_counts) #convert to a data frame
tot_counts
row.names(tot_counts) <- reads$Sample #add sample column and lable row by sample name
tot_counts["Sample"] <- reads$Sample
tot_counts

## make bar graph
tot_counts$Sample <- factor(tot_counts$Sample, levels = tot_counts$Sample) #lock in order of samples, see https://stackoverflow.com/questions/38131596/ggplot2-geom-bar-how-to-keep-order-of-data-frame

ggplot(tot_counts, aes(x=Sample, y=tot_counts)) + geom_col(fill="#08306B", color="black") + ylab("total normalized mapped counts") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) #plot graph

#see http://felixfan.github.io/ggplot2-remove-grid-background-margin/, 

#export tables
write.table(tot_counts, "normalized_counts.tsv", sep="\t")
write.table(counts_matrix, "counts_all_samples.tsv", sep="\t")
