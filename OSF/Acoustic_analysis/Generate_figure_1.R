
## Description: Extraction of pitch contours to plot Figure 1 of Roettger & Cole
## Authors: Timo Roettger
## Date: June, 20th, 2018


#######################
#### Preprocessing ####
#######################

library(tidyverse)

# load in data
setwd(getwd())
data <-  read.csv2("./OSF/Acoustic_analysis/Output_pitch_tracks.csv", 
                   header = F, row.names=NULL)
data_overall <-  read.csv2("./OSF/Acoustic_analysis/Output_pitch_tracks_old.csv", 
                   header = F, row.names=NULL)

# rename colnames
colnames(data) <- c("label","word","position","start","end","duration","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100")
colnames(data_overall) <- c("label","start","end","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100")

# create tidy data frames
cols = c(7:107) 
data[,cols] = apply(data[,cols], 2, function(x) as.numeric(as.character(x)))

cols2 = c(4:103) 
data_overall[,cols2] = apply(data_overall[,cols2], 2, function(x) as.numeric(as.character(x)))

data_tidy <- data %>%
  gather(key = "pos", value = "f0", c(-label, -word, -position, -start, -end, -duration), na.rm = T)

data_tidy2 <- data_overall %>%
  gather(key = "pos", value = "f0", c(-label, -start, -end), na.rm = T)

data_tidy$pos <-  as.numeric(data_tidy$pos)
data_tidy2$pos <-  as.numeric(data_tidy2$pos)

data_tidy$label <- as.character(data_tidy$label)
data_tidy$item <- as.factor(sapply(strsplit(data_tidy$label, split = '_'), function(x){(x[1])}))
data_tidy$focus <- as.factor(sapply(strsplit(data_tidy$label, split = '_'), function(x){(x[2])}))

data_tidy2$label <- as.character(data_tidy2$label)
data_tidy2$item <- as.factor(sapply(strsplit(data_tidy2$label, split = '_'), function(x){(x[1])}))
data_tidy2$focus <- as.factor(sapply(strsplit(data_tidy2$label, split = '_'), function(x){(x[2])}))

data_tidy$duration <- as.numeric(as.character(data_tidy$duration))

data_tidy$start <- as.numeric(as.character(data_tidy$start))
data_tidy$end <- as.numeric(as.character(data_tidy$end))

data_tidy2$start <- as.numeric(as.character(data_tidy2$start))
data_tidy2$end <- as.numeric(as.character(data_tidy2$end))

data_tidy <- data_tidy %>%
  mutate(duration = end - start)

data_tidy2 <- data_tidy2 %>%
  mutate(duration = end - start)


xagg <-
  data_tidy2 %>%
  group_by(focus, pos, item) %>%
  summarise(mean_f0 = mean(f0, na.rm = T))

yagg <-
  data_tidy2 %>%
  group_by(focus, pos) %>%
  summarise(mean_f0 = mean(f0, na.rm = T),
            mean_dur_o = mean(duration, na.rm = T))

dur.agg <-
  data_tidy %>%
  group_by(focus, position) %>%
  summarise(mean_dur = mean(duration, na.rm = T))

dur.agg <- full_join(dur.agg,yagg) 

dur.agg2 <-
  dur.agg %>%
  mutate(rel_dur = mean_dur/mean_dur_o) %>%
  group_by(focus, position) %>%
  summarise(mean_rel = 100 * mean(rel_dur))



##################
#### Plotting ####
##################

# define theme
theme.ED <- 
  theme_classic() + 
  theme(legend.justification = c(0.05,0.05), 
        legend.position = "right",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        strip.text.x = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

# plot indidual plots
broad <- 
  ggplot(yagg[yagg$focus == "broad",], aes(x = pos, y = mean_f0, group = focus)) +
  geom_line(data = xagg[xagg$focus == "broad",], aes(x = pos, y = mean_f0, group = item), lwd = 1, colour = "grey", alpha = 0.5) +
  geom_line(lwd = 3, colour = "#CC6666") +
  annotate("text", x = 15, y = 260, label = "L+H*", size = 10) + 
  annotate("text", x = 73, y = 220, label = "H*", size = 10) +
  annotate("text", x = dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 1,]$mean_rel / 2
             , y = 90, label = "SUBJECT", size = 4, hjust = 0.5, vjust = 0.5) + 
  annotate("text", x = dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 2,]$mean_rel / 2 + 
             dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 1,]$mean_rel
             , y = 90, label = "VERB", size = 4, hjust = 0.5, vjust = 0.5) +
  annotate("text", x = dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 4,]$mean_rel / 2 +
             dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 1,]$mean_rel + 
             dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 2,]$mean_rel
             , y = 90, label = "OBJECT", size = 4, hjust = 0.5, vjust = 0.5) +
  geom_segment(x = dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 1,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 1,]$mean_rel, yend = 0, lty = "dotted", size = 1, colour = "grey") +
  geom_segment(x = dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 2,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "broad" & dur.agg2$position == 2,]$mean_rel,
               yend = 0, lty = "dotted", size = 1, colour = "grey") +
  #facet_wrap(  ~ focus, ncol = 2) +
  labs(x = "\nnormalized time", y = "fundamental\nfrequency in Hz\n") + 
  scale_y_continuous(expand = c(0, 0), breaks = (c(75,300)), limits = c(75,300)) + 
  scale_x_continuous(expand = c(0, 0), breaks = (c(0,50,100)), limits = c(0,100)) + 
  ggtitle("Broad\n") +
  theme.ED +
  theme(legend.position = "none")


given <- 
  ggplot(yagg[yagg$focus == "given",], aes(x = pos, y = mean_f0, group = focus)) +
  geom_line(data = xagg[xagg$focus == "given",], aes(x = pos, y = mean_f0, group = item), lwd = 1, colour = "grey", alpha = 0.5) +
  geom_line(lwd = 3, color = "#66CC99") +
  annotate("text", x = 15, y = 220, label = "H*", size = 10) + 
  annotate("text", x = dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 1,]$mean_rel / 2
           , y = 90, label = "SUBJECT", size = 4, hjust = 0.5, vjust = 0.5) + 
  annotate("text", x = dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 2,]$mean_rel / 2 + 
             dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 1,]$mean_rel
           , y = 90, label = "VERB", size = 4, hjust = 0.5, vjust = 0.5) +
  annotate("text", x = dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 4,]$mean_rel / 2 +
             dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 1,]$mean_rel + 
             dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 2,]$mean_rel
           , y = 90, label = "OBJECT", size = 4, hjust = 0.5, vjust = 0.5) +
  geom_segment(x = dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 1,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 1,]$mean_rel, yend = 0, lty = "dotted", size = 1, colour = "grey") +
  geom_segment(x = dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 2,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "given" & dur.agg2$position == 2,]$mean_rel,
               yend = 0, lty = "dotted", size = 1, colour = "grey") +
  #facet_wrap(  ~ focus, ncol = 2) +
  labs(x = "\nnormalized time", y = "fundamental\nfrequency in Hz\n") + 
  scale_y_continuous(expand = c(0, 0), breaks = (c(75,300)), limits = c(75,300)) + 
  scale_x_continuous(expand = c(0, 0), breaks = (c(0,50,100)), limits = c(0,100)) + 
  ggtitle("Given\n") +
  theme.ED +
  theme(legend.position = "none")


contrastive <- 
  ggplot(yagg[yagg$focus == "contrastive",], aes(x = pos, y = mean_f0, group = focus)) +
  geom_line(data = xagg[xagg$focus == "contrastive",], aes(x = pos, y = mean_f0, group = item), lwd = 1, colour = "grey", alpha = 0.5) +
  geom_line(lwd = 3, color = "#999999") +
  annotate("text", x = 20, y = 280, label = "L+H*", size = 10) + 
  annotate("text", x = 75, y = 170, label = "L*", size = 10) + 
  annotate("text", x = dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 1,]$mean_rel / 2
           , y = 90, label = "SUBJECT", size = 4, hjust = 0.5, vjust = 0.5) + 
  annotate("text", x = dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 2,]$mean_rel / 2 + 
             dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 1,]$mean_rel
           , y = 90, label = "VERB", size = 4, hjust = 0.5, vjust = 0.5) +
  annotate("text", x = dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 4,]$mean_rel / 2 +
             dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 1,]$mean_rel + 
             dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 2,]$mean_rel
           , y = 90, label = "OBJECT", size = 4, hjust = 0.5, vjust = 0.5) +
  geom_segment(x = dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 1,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 1,]$mean_rel, yend = 0, lty = "dotted", size = 1, colour = "grey") +
  geom_segment(x = dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 2,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "contrastive" & dur.agg2$position == 2,]$mean_rel,
               yend = 0, lty = "dotted", size = 1, colour = "grey") +
  #facet_wrap(  ~ focus, ncol = 2) +
  labs(x = "\nnormalized time", y = "fundamental\nfrequency in Hz\n") + 
  scale_y_continuous(expand = c(0, 0), breaks = (c(75,300)), limits = c(75,300)) + 
  scale_x_continuous(expand = c(0, 0), breaks = (c(0,50,100)), limits = c(0,100)) + 
  ggtitle("Contrastive\n") +
  theme.ED +
  theme(legend.position = "none")


narrow <- 
  ggplot(yagg[yagg$focus == "narrow",], aes(x = pos, y = mean_f0, group = focus)) +
  geom_line(data = xagg[xagg$focus == "narrow",], aes(x = pos, y = mean_f0, group = item), lwd = 1, colour = "grey", alpha = 0.5) +
  geom_line(lwd = 3, color = "#9999CC") +
  annotate("text", x = dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 1,]$mean_rel / 2
           , y = 90, label = "SUBJECT", size = 4, hjust = 0.5, vjust = 0.5) + 
  annotate("text", x = dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 2,]$mean_rel / 2 + 
             dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 1,]$mean_rel
           , y = 90, label = "VERB", size = 4, hjust = 0.5, vjust = 0.5) +
  annotate("text", x = dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 4,]$mean_rel / 2 +
             dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 1,]$mean_rel + 
             dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 2,]$mean_rel
           , y = 90, label = "OBJECT", size = 4, hjust = 0.5, vjust = 0.5) +
  geom_segment(x = dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 1,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 1,]$mean_rel, yend = 0, lty = "dotted", size = 1, colour = "grey") +
  geom_segment(x = dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 2,]$mean_rel, 
               y = +Inf, xend = dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 1,]$mean_rel +
                 dur.agg2[dur.agg2$focus == "narrow" & dur.agg2$position == 2,]$mean_rel,
               yend = 0, lty = "dotted", size = 1, colour = "grey") +
  #facet_wrap(  ~ focus, ncol = 2) +
  labs(x = "\nnormalized time", y = "fundamental\nfrequency in Hz\n") + 
  scale_y_continuous(expand = c(0, 0), breaks = (c(75,300)), limits = c(75,300)) + 
  scale_x_continuous(expand = c(0, 0), breaks = (c(0,50,100)), limits = c(0,100)) + 
  ggtitle("Narrow\n") +
  theme.ED +
  theme(legend.position = "none")



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

# plot four plots together
multiplot(given, broad, narrow, contrastive)
