# R-scripts-for-read-mapping
R scripts  to make heatmaps and bargraphs summarizing mapped reads.

This script was used in a Master's thesis project, and is intended to be used specifically for the data generated in that project. The code is here to allow for reproducibility. This script is not intended to be used for generic sequencing data. You can use it as a template to plot your own data, but it might not work exactly 'as is' with other data, depending on how it is formatted.

This script is intended to be used to visualize data obtained from read-mapping-code-amoA script. 

Needed: R

R packages: ComplexHeatmap from BiocLite, circlize, ggplot2, RcolorBrewer, dplyr

Normalization of count data: First divide each count by its respective gene length. Then divide each count by its respective sample total read count. This gives you Reads per kilobase per 100 million unassembled reads.

