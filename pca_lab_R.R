
# Install and load packages

install.packages(c("factoextra", "FactoMineR", "RColorBrewer"))
devtools::install_github("kassambara/ggpubr")

library("factoextra")
library("FactoMineR")
library("ggpubr")
library("RColorBrewer")


# read in data
data <- read.csv("pca_data.csv", row.names = 1)

# run PCA and visualize
pca.data <- PCA(data[, -1], scale.unit = TRUE, graph = FALSE)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))


# Create a correlation plot 
fviz_pca_var(pca.data, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)


# Run PCA using cell types. Use t() to flip the table so cell types are rows
pca.data <- PCA(t(data[,-1]), scale.unit = TRUE, graph = FALSE)


#fviz_pca_ind is used to visualize the data in a PCA plot.
fviz_pca_ind(pca.data, col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

# Store PCA plot in "a" and use ggpar to add labels
a <- fviz_pca_ind(pca.data, col.ind = "cos2", 
                  gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
                  repel = TRUE)

ggpar(a, title = "Principal Component Analysis",
      xlab = "PC1",
      ylab = "PC2",
      legne.title = "Cos2",
      legend.position = "top",
      ggtheme = theme_minimal())


# Run pca using genes instead of cell types
pca.data <- PCA(data[,-1], scale.unit = TRUE, ncp = 2, graph = FALSE)

# Convert lineages to factor data type
data$Lineage <- as.factor(data$Lineage)


#Create a 3 color palette (for each cell type)
nb.cols <- 3
mycolors <- colorRampPalette(brewer.pal(3, "Set1"))(nb.cols)
```


# Create PCA plot
fviz_pca_ind(pca.data, col.ind = data$Lineage,
             palette = mycolors, addEllipses = TRUE)

# Store PCA plot as "a" and use ggpar to add labels
a <- fviz_pca_ind(pca.data, col.ind = data$Lineage,
                  palette = mycolors, addEllipses = TRUE)

ggpar(a,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cell type", legend.position = "top",
      ggtheme = theme_minimal())
