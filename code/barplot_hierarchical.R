library(microbiome)
library(dendextend)
library(phyloseq)
library(patchwork)
library(purrr) 
library(ggdendro)

# Create tax_table for other's phyloseq
tax_table_other <- function(colnames){
  df <- data.frame(matrix('other',nrow=1,ncol = length(colnames)))
  colnames(df) <- colnames
  rownames(df) <- 'other'
  return(df)
}

# Preprocessing: sort by detection and prevalence 
preprocessing <- function(ps,taxrank,detection,prevalence){
  taxas <- core_members(ps_obj, detection = detection, prevalence = prevalence)
  ps_obj <- prune_taxa(taxas, ps_obj)
  ps_obj <- speedyseq::tax_glom(ps_obj, taxrank=taxrank, NArm = FALSE)
  
  return(ps_obj)
}

barplot_hierarchical <- function(ps,taxrank,detection,prevalence,transformation = 'compositional',feature,top ,colors_barplot,colors_line,dist,size,
                                 decreasing,tax_table_other){
  ps_obj <- ps
  ps_obj <- microbiome::transform(ps_obj, transformation)
  
  #Create sample data with the desired characteristic
  metadata <- ps_obj@sam_data[,feature]

  
  taxa_matrix <- as(phyloseq::tax_table(ps_obj), 'data.frame')
  taxa_matrix <- cbind(ASV=rownames(taxa_matrix), taxa_matrix) 
  
  # Choose top taxa

  topx <- top_taxa(ps_obj, n = top )
  # Choose others taxa
  not_topx <- !(taxa_matrix$ASV %in% topx)

  # Create two phyloseq objects
  ps_obj_top <- prune_taxa(topx, ps_obj)
  ps_obj_top@sam_data <- metadata
  ps_obj_other <-   prune_taxa(not_topx, ps_obj)

  # Let's summarize the columns and create a separate phyloseq object containing only information about the others 
  otu_table_other <- as.matrix(rowSums(otu_table(ps_obj_other)))
  colnames(otu_table_other) <- "other"
  
  tax_tab_oth <- tax_table_oth
  tax_tab_oth <- as.matrix(tax_tab_oth)
  TAX = phyloseq::tax_table(tax_tab_oth)
  
  OTU = otu_table(otu_table_other, taxa_are_rows = FALSE)
  
  samples = sample_data(metadata)
  
  ps_other <- phyloseq(OTU,TAX, samples)
  # Create wide plot
  data_plot <- psmelt(ps_obj_top)
  data_plot_other <- psmelt(ps_other)

  colnames_to_keep <- colnames(data_plot)
  data_plot_other <- data_plot_other[,colnames_to_keep]
  df <- rbind(data_plot,data_plot_other)
  
  # Create dendogram
  ps_obj <- ps
  
  if (dist == 'bray'){
    ps_obj <- microbiome::transform(ps_obj, 'compositional')
    dend <- as.dendrogram(hclust(vegdist(ps_obj@otu_table, method="bray")))
  }
  if (dist == 'aitchinson'){
    ps_obj <- microbiome::transform(ps_obj, 'clr')
    dend <- as.dendrogram(hclust(dist(ps_obj@otu_table)))
  }
  
  hc <- dendro_data(as.dendrogram(dend))
  
  df[,taxrank] <- factor(df[,taxrank], levels = unique(df[,taxrank][order(df[,taxrank] == 'other',decreasing = decreasing)]))
  df$Sample <- factor(df$Sample, levels = labels(dend))
  
  sample_data <- sample_data(ps_obj)
  
  title <- paste('detection:',detection,'prevalence:',prevalence,'tax_rank',taxrank,sep='__')
  #df$taxrank <- str_wrap(df$taxrank, width = 2)
  # Stacked Percent
  p1 <- ggplot(data = segment(hc)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    scale_x_discrete(labels=label(hc)$label) +
    ylab(dist) + theme_bw()+ ggtitle(title)+
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank())
  #p1 <- ggplot(dend, horiz = F) 
  p2 <- ggplot(df, aes(fill=!!sym(taxrank), y=Abundance, x=Sample)) + 
    geom_bar( stat="identity", position="fill",width = 1,color="black",size=0.1) +scale_fill_manual("legend", values = colors_barplot, name = taxrank)+
    theme(axis.text.x = element_blank(),
          axis.line.x = element_line(size = 1, colour = "black"),
          axis.title.x = element_text(size = 10, hjust = 1),
          legend.text =  element_text(size = 12),
          legend.title = element_text(size = 10,face = "bold"),
          legend.position = 'right',
          legend.box = "vertical")+
    guides(fill = guide_legend(ncol = 2))
  # Phenotype's metadata
  p3<-ggplot(df,aes(Sample,y=1,fill=!!sym(feature)))+geom_tile(color='black',size=0.1)+
    theme(axis.title=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          legend.text =  element_text(size = 10),
          legend.title = element_text(size = 10,face = "bold"))+
    scale_fill_manual(values=colors_line)
  
  p4 <- p1 /p3 / p2 + plot_layout(heights = size) 
  
  return(p4)
}

