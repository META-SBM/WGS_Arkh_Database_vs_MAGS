library(phyloseq)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)

library(phyloseq)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lemon)
library(ggsignif)
library(rstatix)
library(vegan)
library(ggrepel)

# Define the function to create the ordination plots
create_ordination_plots <- function(ps_obj, method = "PCoA", distance_method = "wunifrac", group, size = 10,palette,level_factors,taxa_plot = TRUE,
                                    len_axis =10) {
  
  #ps_obj <- ps_obj_transform
  # Calculate distance matrix
  dist <- phyloseq::distance(ps_obj, method = distance_method)
  #dist <- vegdist(ps_obj, method = distance_method)
  
  # Perform ordination
  ordination <- ordinate(ps_obj, method = method, distance = dist)
  P <- cbind(as.data.frame(ordination$vectors), as.data.frame(as.matrix(sample_data(ps_obj))))
  
  # Check if the group column exists
  if (!group %in% colnames(P)) {
    stop(paste("Column", group, "does not exist in the sample data."))
  }
  
  # Calculate percentage of variance explained by the first two axes
  variance_explained <- 100 * ordination$values$Relative_eig[1:2]
  variance_explained <- 100 * ordination$values$Relative_eig[1:2]
  # Calculate the means for each unique value of the group
  means <- P %>%
    group_by(!!rlang::sym(group)) %>%
    summarise(mean_Axis1 = mean(Axis.1), mean_Axis2 = mean(Axis.2))
  
  # Merge means back into the original data
  P <- P %>%
    dplyr::left_join(means, by = group)
  P[,group] <- factor(P[,group],level = level_factors)
  anosim_res <- anosim(dist, P[[group]],parallel = 30)
  anosim_text <- paste("ANOSIM R:", round(anosim_res$statistic, 3), "p-value:", round(anosim_res$signif, 3))
  
  if(taxa_plot != F){
    ps_ph <- speedyseq::tax_glom(ps_obj,taxrank = 'Genus', NArm=T )
    env_data <- as.data.frame(otu_table(ps_ph))  # Replace with your environmental data
    # Fit environmental vectors using envfit
    env_fit <- envfit(as.data.frame(ordination$vectors), env_data, permutations = 999, na.rm = TRUE)
    # Extract significant vectors (e.g., p-value < 0.05)
    significant_vectors <- as.data.frame(env_fit$vectors$arrows)
    
    # Add p-values to the data frame
    significant_vectors$pvals <- env_fit$vectors$pvals
    # Add correlation coefficients (r values) to the data frame
    significant_vectors$r <- env_fit$vectors$r
    
    # Filter for significant vectors (p-value < 0.05)
    significant_vectors <- significant_vectors %>%
      filter(pvals < 0.05)
    
    # Select top 10 vectors with the largest r values
    top_significant_vectors <- significant_vectors %>%
      arrange(desc(r)) %>%
      head(10)
    
    tax_df <- as.data.frame(tax_table(ps_obj))
    tax_df <- tax_df[,c('Phylum','Genus'),drop = F]
    top_significant_vectors <- merge(top_significant_vectors,tax_df,by = 0,all.x=T)
    rownames(top_significant_vectors) <- top_significant_vectors$Row.names
    top_significant_vectors$taxa <- paste0(top_significant_vectors$Phylum,'_',top_significant_vectors$Genus)
  }
  
  
  # Create scatter plot with ellipses and means
  pl <-  ggplot(P, aes_string(x = "Axis.1", y = "Axis.2", color = group)) +
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse() +
    geom_segment(aes(x = Axis.1, y = Axis.2, xend = mean_Axis1, yend = mean_Axis2), alpha = 0.3) +
    geom_label(data = means, aes(x = mean_Axis1, y = mean_Axis2, label = !!rlang::sym(group)), fill = "white", color = "black", fontface = "bold", size = 3) +
    geom_text(x = Inf, y = Inf, label = anosim_text, hjust = 1.1, vjust = 1.1, size = size * 0.3, color = "black") +
    scale_color_manual(values=palette) +
    xlab(label = paste0('Axis.1 [', round(variance_explained[1], 1), '%]')) +
    ylab(label = paste0('Axis.2 [', round(variance_explained[2], 1), '%]')) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size, color = 'black'),
      legend.position = "none",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      text = element_text(size = size, color = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    )
  #pl
  if(taxa_plot != F){
    pl <-pl + geom_segment(data = top_significant_vectors, 
                    aes(x = 0, y = 0, 
                        xend = Axis.1 * len_axis, 
                        yend = Axis.2 * len_axis),  # Scale vectors for better visualization
                    arrow = arrow(length = unit(0.2, "cm")), 
                    color = "black", alpha = 0.8) +
      geom_text_repel(data = top_significant_vectors, 
                      aes(x = Axis.1 * len_axis, 
                          y = Axis.2 * len_axis, 
                          label = taxa),  # Add taxonomic labels
                      color = "black",
                      alpha = 0.9,
                      size = 4, 
                      vjust = 0.5, 
                      hjust = 1,
                      box.padding = 0.5,   # Optional: Adjust padding around text
                      point.padding = 0.5,
                      max.overlaps = Inf)
  }
  # Perform Wilcoxon test for Axis.1
  # Perform pairwise comparisons for Axis.1
  # Define pairwise comparisons
  unique_vals <- unique(P[[group]])
  my_comparisons <- lapply(combn(unique_vals, 2, simplify = FALSE), as.vector)
  anno_df <- ggpubr::compare_means(as.formula(paste("Axis.1 ~", group)), data = P, method = "wilcox.test",p.adjust.method = "BH")%>%
    add_significance("p.adj") %>%
    add_x_position()%>%
    add_y_position( data=P,formula = as.formula(paste("Axis.1 ~", group)),step.increase = 0.2)
  
  for (i in 1:nrow(anno_df)) {
    # Проверяем group1
    if (anno_df$group1[i] %in% level_factors) {
      anno_df$xmin[i] <- which(level_factors == anno_df$group1[i]) # Получаем порядковый номер для xmin
    }
    
    # Проверяем group2
    if (anno_df$group2[i] %in% level_factors) {
      anno_df$xmax[i] <- which(level_factors == anno_df$group2[i]) # Получаем порядковый номер для xmax
    }
  }
  
  # Create violin plot for Axis.1
  x_dens <- ggplot(P, aes_string(x = group, y = "Axis.1",color=group)) +
    geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75, position = position_dodge(0.9)) +
    geom_jitter(size = 1.5, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9)) +
    scale_color_manual(values=palette) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size, color = 'black'),
      legend.position = "none",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      panel.background = element_rect(fill = "white", color = "black"),
      text = element_text(size = size, color = 'black'),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    )+
    stat_pvalue_manual(
      anno_df,  label = "p.adj.signif", tip.length = 0.02,
      step.increase = 0.05,coord.flip = TRUE
    )+
    coord_flip()
  
  # Perform Wilcoxon test for Axis.2
  test_result_axis2 <- ggpubr::compare_means(as.formula(paste("Axis.2 ~", group)), data = P, method = "wilcox.test")
  anno_df <- ggpubr::compare_means(as.formula(paste("Axis.2 ~", group)), data = P, method = "wilcox.test",p.adjust.method = "BH")%>%
    add_significance("p.adj") %>%
    add_x_position()%>%
    add_y_position( data=P,formula = as.formula(paste("Axis.2 ~", group)),step.increase = 0.2)
  
  for (i in 1:nrow(anno_df)) {
    # Проверяем group1
    if (anno_df$group1[i] %in% level_factors) {
      anno_df$xmin[i] <- which(level_factors == anno_df$group1[i]) # Получаем порядковый номер для xmin
    }
    
    # Проверяем group2
    if (anno_df$group2[i] %in% level_factors) {
      anno_df$xmax[i] <- which(level_factors == anno_df$group2[i]) # Получаем порядковый номер для xmax
    }
  } 
  # Create violin plot for Axis.2
  y_dens <- ggplot(P, aes_string(x = group, y = "Axis.2", color = group)) +
    geom_violin(trim = FALSE, alpha = 0.1) +
    geom_boxplot(width = 0.5, alpha = 0.75, position = position_dodge(0.9)) +
    geom_jitter(size = 1.5, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9)) +
    scale_color_manual(values=palette) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 45, hjust = 1, size = size, color = 'black'),
      legend.position = "none",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      plot.title = element_text(size = 25),
      panel.background = element_rect(fill = "white", color = "black"),
      text = element_text(size = size, color = 'black'),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.1, linetype = 'solid', color = "gray")
    ) +
    stat_pvalue_manual(
      anno_df,  label = "p.adj.signif", tip.length = 0.02,
      step.increase = 0.05,coord.flip = FALSE
    )
  
  l <- ggplot(P, aes_string(x = "Axis.1", y = "Axis.2", color = group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse() +
    geom_segment(aes(xend = mean(Axis.1), yend = mean(Axis.2)), alpha = 0.3) +
    geom_label(aes(x = mean(Axis.1), y = mean(Axis.2), label = !!rlang::sym(group)), fill = "white", color = "black", fontface = "bold", size = 3) +
    scale_color_manual(values=palette) +
    xlab(label = paste0('Axis.1 [', round(variance_explained[1], 1), '%]')) +
    ylab(label = paste0('Axis.2 [', round(variance_explained[2], 1), '%]')) +
    theme(
      axis.text.y = element_text(color = "black", size = size),
      axis.text.x = element_text(angle = 90, hjust = 1, size = size, color = 'black'),
      legend.position = "right",
      axis.title.y = element_text(color = "black", size = size, angle = 90),
      axis.title.x = element_text(color = "black", size = size),
      text = element_text(size = size, color = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.major = element_line(size = 0.01, linetype = 'solid', color = "gray"),
      panel.grid.minor = element_line(size = 0.01, linetype = 'solid', color = "gray")
    ) 
  # Create a blank plot
  legend <- g_legend(l)
  
  # Make grobs
  gA <- ggplotGrob(x_dens)
  gB <- ggplotGrob(pl)
  gD <- ggplotGrob(y_dens)
  gL <-  legend
  # Get width
  xWidth = unit.pmax( gA$widths[2:4], gB$widths[2:4])
  yHeight = unit.pmax(gB$heights[4:5], gD$heights[4:5])
  
  # Set the widths
  gA$widths[2:3] <- xWidth
  gB$widths[2:3] <- xWidth
  # Set the heights
  gB$heights[4:5] <- yHeight
  gD$heights[3:5] <- yHeight
  
  p =grid.arrange(gD,gB,gL ,gA,ncol=2, nrow=2, widths=c(2, 5), heights=c(5, 2))
  
}

plot_alpha_div <- function(
    ps_obj,
    group ,
    color ,
    measure 
) {
  # Generate boxplot of alpha diversity without transformation
  p <- plot_richness(ps_obj, x = group, color = color, measures = measure)
  
  return(p)
}
create_alpha_plots <- function(ps_obj,col,measure = 'Shannon',method = 'wilcox.test',color,my_comparisons,size=10,level_factors){
  
  p <- plot_alpha_div(ps_obj, group = col, color = col, measure = measure) +
  geom_violin(trim=F, alpha=0.1) +
  geom_boxplot(width=0.5, alpha=0.75, position=position_dodge(0.9)) +
  geom_jitter(size=1.5, alpha=0.5, position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) + scale_color_manual(values=color)+
  theme_bw(base_size=20)+
  theme(
    plot.title = element_text(size = 25),
    axis.text.y = element_text(color = "black", size = size),
    axis.text.x = element_text(angle=45, hjust=1,size=15,color = 'black'),
    legend.position = "none",
    axis.title.y  = element_text(color = "black", size = size,angle=90),
    axis.title.x  = element_text(color = "black", size = size),
    #legend.key.size = unit(0.5, 'cm'),
    text = element_text(size = size,colour ='black'))+
    scale_color_manual(values=color)
  
  group <- col
  df <- as.data.frame(p$data)
  anno_df <- ggpubr::compare_means(as.formula(paste("value ~", group)), data = df, method = "wilcox.test",p.adjust.method = "BH")%>%
    add_significance("p.adj") %>%
    add_x_position()%>%
    add_y_position( data=df,formula = as.formula(paste("value ~", group)),step.increase = 0.2)
  
  for (i in 1:nrow(anno_df)) {
    # Проверяем group1
    if (anno_df$group1[i] %in% level_factors) {
      anno_df$xmin[i] <- which(level_factors == anno_df$group1[i]) # Получаем порядковый номер для xmin
    }
    
    # Проверяем group2
    if (anno_df$group2[i] %in% level_factors) {
      anno_df$xmax[i] <- which(level_factors == anno_df$group2[i]) # Получаем порядковый номер для xmax
    }
  }
  
  p <- p +
    stat_pvalue_manual(
    anno_df,  label = "p.adj.signif", tip.length = 0.02,
    step.increase = 0.05,coord.flip = FALSE
  )
  return(p)
}