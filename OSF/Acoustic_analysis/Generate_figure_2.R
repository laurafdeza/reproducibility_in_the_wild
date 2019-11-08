
## Description: Extraction of pitch information to plot Figure 2 of Roettger & Cole
## Authors: Timo Roettger
## Date: June, 20th, 2018

#######################
#### Preprocessing ####
#######################

## libraries
library(tidyverse)
library(ggbeeswarm)

###########
## legend #
###########

## specify path
setwd(getwd())
ax <- read.csv("./OSF/Acoustic_analysis/all_acoustics.csv")


###########################
#### acoustic analysis ####
###########################

## aggregate accented vowel duration of subject and object
ax_agg <- ax %>%
  spread(position, max_f0) %>%
  group_by(focus, target) %>%
  summarise(f0_max_subj = mean(subject, na.rm = T),
            f0_max_obj = mean(object, na.rm = T),
  ) %>%
  mutate(rel_peak = f0_max_subj - f0_max_obj)
  

## F0max on subject
f0max_subj <- 
  ggplot(ax_agg, aes(color = focus, x = focus, y = f0_max_subj)) +
  geom_violin(fill = "grey95", color = "grey95", trim = F) +
  geom_quasirandom(size = 3, width = 0.2) + 
  scale_color_manual(labels = c("b" = "broad", "n" = "narrow", "c" = "contrastive", "g" = "given"),
                     guide = guide_legend(title = "Intended focus type"),
                     values = c("#CC6666", "#999999", "#66CC99", "#9999CC")) +
  scale_y_continuous(expand = c(0, 0), breaks = (c(150,200,250,300)), limits = c(150,300)) +
  ylab("Maximum f0 \n of sentence subject\n") +
  xlab("\n Focus type") +
  theme_classic() + 
  theme(legend.position = "none",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"))

## F0 max on object
f0max_obj <- 
  ggplot(ax_agg, aes(color = focus, x = focus, y = f0_max_obj)) +
  geom_violin(fill = "grey95", color = "grey95", trim = F) +
  geom_quasirandom(size = 3, width = 0.2) + 
  scale_color_manual(labels = c("b" = "broad", "n" = "narrow", "c" = "contrastive", "g" = "given"),
                     guide = guide_legend(title = "Intended focus type"),
                     values = c("#CC6666", "#999999", "#66CC99", "#9999CC")) +
  scale_y_continuous(expand = c(0, 0), breaks = (c(110,150,190,230)), limits = c(110,230)) +
  ylab("Maximum f0 \n of sentence object\n") +
  xlab("\n Focus type") +
  theme_classic() + 
  theme(legend.position = "none",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.line.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"))  


## function for multiple plots
multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(f0max_subj, f0max_obj)  
