library(tibble)
library(ggtext)
library(glue)
library(microViz)
library(PNWColors)
library(showtext)
font_add_google("Montserrat", "montserrat")

plot_permanova <- function(ps_obj, formula, transformation = "compositional", method = "bray", show_plot = TRUE,level,det,prev,taxa_num,year,lab.size=lab.size,size=15,cols_meta,
                           palette,threshold,threshold_loc_y) {
  
  # Convert sample data to a data frame
  metadf <- as.data.frame(as.matrix(ps_obj@sam_data))
  
  if(method == 'euclidean'){
    distance <- 'aitchison'
  }else{
    distance <- method 
  }
  bd1 <- ps_obj %>%
    dist_calc(distance) %>%
    dist_bdisp(variables = cols_meta$new_name) %>%
    bdisp_get()
  
  if (method == 'euclidean'){
    transformation <- 'clr'
    ps_obj <- microbiome::transform(ps_obj,transformation)
  }
  if ((method == 'unifrac') | (method =='wunifrac')){
    transformation <- "None"
  }
  
  if ((transformation != "None" )& (method != 'euclidean') &(method != 'unifrac')) {
    ps_obj <- microbiome::transform(ps_obj, transformation)
  }
  
  # Calculate distance matrix based on the specified method
  dist <- phyloseq::distance(ps_obj, method =method)
  
  # Create formula from the string
  formula1 <- as.formula(paste("dist ~", formula))
  
  # Perform adonis analysis with the specified formula
  permanova_result <- adonis2(formula1, data = metadf,by='term',parallel = 30)
  
  # Process adonis results for plotting
  res <- as.data.frame(permanova_result)
  res$meta <- row.names(res)
  res <- subset(res, res$meta != 'Total')
  res$meta <- gsub('Residuals', 'other', res$meta)
  res[which(res$`Pr(>F)` > 0.05), 'signif'] <- ' '
  res[which(res$`Pr(>F)` <= 0.05), 'signif'] <- '*'
  res[which(res$`Pr(>F)` <= 0.01), 'signif'] <- '**'
  res[which(res$`Pr(>F)` <= 0.001), 'signif'] <- '***'
  res$Feature <- as.character(res$meta)
  res <- res %>%
    mutate(  Feature = ifelse((signif == '*') | (signif == '**') | (signif == '***') ,
                              glue("**{Feature}**"),  # Use glue to format as bold
                              Feature))
  res$Feature <- factor(res$Feature, levels = res$Feature[order(res$R2, decreasing = FALSE)])
  #res$Feature <- as.character(res$Feature)
  #res$Feature <- ifelse(nchar(res$Feature) > 10, 
  #                      str_replace_all(res$Feature, "_", "\n"), 
  #                      res$Feature)
  res$Feature <- str_wrap(res$Feature, width = 1)
  res <- tibble::rownames_to_column(res, var = "new_name")
  res <- res %>%
    left_join(cols_meta, by = 'new_name')
  res <- res %>%
    filter(new_name != "Residual")
  permdisp_df <- data.frame(new_name = character(), p_val = numeric(), stringsAsFactors = FALSE)
  
  # Проходим по всем элементам в cols_meta$new_name
  for (elem in cols_meta$new_name) {
    # Извлекаем p-value из bd1 для текущего элемента
    p_value <- bd1[[elem]]$anova$`Pr(>F)`[1]  # Используем двойные квадратные скобки для доступа к элементу
    
    # Добавляем новую строку в results
    permdisp_df <- rbind(permdisp_df, data.frame(new_name = elem, p_val = p_value))
  }
  permdisp_df[which(permdisp_df$p_val > 0.05), 'signif_disp'] <- ' '
  permdisp_df[which(permdisp_df$p_val <= 0.05), 'signif_disp'] <- '*'
  permdisp_df[which(permdisp_df$p_val <= 0.01), 'signif_disp'] <- '**'
  permdisp_df[which(permdisp_df$p_val <= 0.001), 'signif_disp'] <- '***'
  
  #print(permdisp_df[1,])
  #print(res[1,])
  res <- merge(res,permdisp_df,by.x='meta',by.y='new_name',all.x=T)
  #res <- res %>%
  #  dplyr::left_join(permdisp_df,by='new_name')
  #print(res[1,])
  # Create the ggplot object for plotting
  res <- res %>%
    group_by(category) %>%          # Group by 'category'
    arrange(desc(R2), .by_group = TRUE) %>%  # Sort within each group by 'R2' in descending order
    ungroup() 
  res$Feature <- factor(res$Feature, levels = res$Feature[order(-res$R2)])
  #print(res)  
  main_plot <- ggplot(res, aes(x = Feature, y = R2, fill = category,order =R2)) +
    geom_col() +  # Use geom_col() to plot actual values
    labs(title = "",
         x = "",
         y = "R2") +
    theme_classic() +  # Optional: use a minimal theme
    coord_flip()+
    geom_label(data = res %>% filter(signif != ' '),
               aes(label = signif,color='PERMANOVA \n significance'),  # Set fontface based on significance
               size = 3, hjust = -0.2,fill='white')+
    geom_label(data = res %>% filter(signif_disp != ' '),
               aes(label = signif_disp,color = 'PERMDISP \n significance'),  # Set fontface based on significance
               size = 3,fill='white', hjust = 1)+
    theme(
      panel.background = element_rect(fill="white"),
      legend.position = 'top',               # Position legend at top left corner
      legend.title = element_blank(), 
      legend.text = element_text(size=14,family = "montserrat",color='black'),
      legend.direction = 'horizontal',
      legend.justification = "left",
      axis.title.y = element_blank(),
      axis.line.x = element_line(color="black"),
      axis.line.y = element_blank(),
      axis.text.y = element_markdown(size=14, hjust=1,family = "montserrat",color='black',linewidth = 1.2),
      axis.text.x = element_text(size=14, hjust=1,family = "montserrat",color='black'),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.text.x = element_blank(),  # Remove x strip labels
      strip.text.y = element_blank()
    )+ 
    facet_grid(category ~ ., scales = "free", space = "free")+
    guides(fill = guide_legend(nrow = 1))+
    
    # Use pnw_colors for filling
    scale_fill_manual(values = palette)+
    scale_color_manual(values=c("PERMANOVA \n significance"="red", "PERMDISP \n significance"="black"))+
    #scale_fill_manual(values = c("Significance"="white", "Dispersion"="black"))+
    geom_hline(yintercept = threshold, col = "grey30", lty = "dashed")
  #ggtitle(paste('method:',distance,'\nformula',formula))
  #annotate("text", x = 0.6, y = 0.6, label = "major contribution",
  #           family = "Fira Sans", size = 3, hjust = 0) 
  plot(main_plot)
  
  total_effect_size <- round(sum(res$R2)*100,3)
  unexplained_variance <- 100 - total_effect_size
  
  
  # Prepare data for circular barplot
  circular_data <- data.frame(
    category = c("Explained Variance", "Unexplained Variance"),
    value = c(total_effect_size, unexplained_variance)
  )
  
  
  circular_data$fraction <- circular_data$value / sum(circular_data$value)
  
  # Compute the cumulative percentages (top of each rectangle)
  circular_data$ymax <- cumsum(circular_data$fraction)
  
  # Compute the bottom of each rectangle
  circular_data$ymin <- c(0, head(circular_data$ymax, n=-1))
  
  # Compute label position
  circular_data$labelPosition <- (circular_data$ymax + circular_data$ymin) / 2
  
  # Compute a good label
  circular_data$label <- paste0(circular_data$category, "\n value: ", circular_data$value)
  
  # Make the plot
  circular_plot <- ggplot(circular_data, aes(ymax=ymax, ymin=ymin, xmax=3.5, xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=4, aes(y=labelPosition, label=label), size=5) +
    scale_fill_manual(values = pnw_palette("Shuksan2", n=2))+
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none")
  #plot(circular_plot)
  combined_plot <- main_plot + circular_plot + plot_layout(ncol=2,widths = c(2,1)) 
  
  # Check if plot should be generated
  if (show_plot) {
    print(combined_plot)
  }
  # Return both adonis results and the ggplot object
  return(list(permanova_result, combined_plot,res))
}