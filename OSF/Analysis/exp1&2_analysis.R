
## Description: Bayesian analysis and plotting posterior distributions of Experiment 1 and 2 of Roettger & Cole
## authors: Timo Roettger

#######################
#### Preprocessing ####
#######################

## libraries
library(tidyverse)
library(brms)
library(rstan)
library(ggbeeswarm)
library(rstudioapi)
options(mc.cores=parallel::detectCores()) # Run on multiple cores for Bayesian regressions

###########
## legend #
###########

# stimulus_id:                  the unique id given to a stimulus pair
# participant_id:               the unique id given to each participant
# task_id:                      always 'ident' for experiment1
# pronoun:                      unique id of experimental sentences
# correct_ans:                  the correct answer for this stimulus_id
# incorrect_ans:                the incorrect answer for this stimulus_id
# stim_pair_id:                 category pair information
# correct_response:             1 character form of /correct_ans/ column       
# actual_response:              user's response in 1 character form (c, b, g, n)
# binned_response:              binary correct (1) or incorrect (0) response label.  1 if correct_response = actual_response and 0 otherwise.

## specify path

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

setwd("../Data/")

## load in data
x2 <- read.csv("experiment1&2.csv")
x2$binned_response <- as.numeric(x2$binned_response)


###############
#### Model ####
###############

# set priors for intercept and coefficiants of fixed effects
prior_x2 <- c(prior(normal(0,5), class = b), prior(normal(0,5), class = Intercept))

# model Exp1 
xmdl_x2_discr <- brm(binned_response ~  correct_ans * incorrect_ans 
                 + (1 + correct_ans * incorrect_ans | pronoun)
                 + (1 | participant_id)
                 , 
                 prior = prior_x2, 
                 family = binomial("logit"),
                 control = list(adapt_delta = 0.8),
                 data = x2[x2$task_id == "discr",])

# model Exp2
xmdl_x2_ident <- brm(binned_response ~  correct_ans * incorrect_ans 
                     + (1 + correct_ans * incorrect_ans | pronoun)
                     + (1 | participant_id)
                     , 
                     prior = prior_x2, 
                     family = binomial("logit"),
                     control = list(adapt_delta = 0.8),
                     data = x2[x2$task_id == "ident",])

## save models for later use
setwd("../Analysis/")
save(xmdl_x2_discr, xmdl_x2_ident, file = "CMR_Bayesian_model_x1_and_x2.RData")

# load bayesian model
load(file = "CMR_Bayesian_model_x1_and_x2.RData")


####################################
#### extract descriptive values ####
####################################

# aggregate raw values for experimental stimuli
xfit2_discr_agg_item <- x2[x2$task_id == "discr",] %>%
  group_by(incorrect_ans, correct_ans, pronoun) %>%
  summarise(est = mean(binned_response))

levels(xfit2_discr_agg_item$incorrect_ans) <- c("Broad", "Contrastive", "Given", "Narrow")
levels(xfit2_discr_agg_item$correct_ans) <- c("Broad", "Contrastive", "Given", "Narrow")

xfit2_ident_agg_item <- x2[x2$task_id == "ident",] %>%
  group_by(incorrect_ans, correct_ans, pronoun) %>%
  summarise(est = mean(binned_response))

levels(xfit2_ident_agg_item$incorrect_ans) <- c("Broad", "Contrastive", "Given", "Narrow")
levels(xfit2_ident_agg_item$correct_ans) <- c("Broad", "Contrastive", "Given", "Narrow")


############################################
#### extract posterior values for discr ####
############################################

# posterior distribution against chance in discr 
psamples_d = posterior_samples(xmdl_x2_discr, pars = rownames(brms:::fixef.brmsfit(xmdl_x2_discr))) %>% 
  mutate(bc = b_Intercept + b_incorrect_anscontrastive,
         bg = b_Intercept + b_incorrect_ansgiven,
         bn = b_Intercept + b_incorrect_ansnarrow,
         cb = b_Intercept + b_correct_anscontrastive,
         cg = b_Intercept + b_correct_anscontrastive + b_incorrect_ansgiven + `b_correct_anscontrastive:incorrect_ansgiven`,
         cn = b_Intercept + b_correct_anscontrastive + b_incorrect_ansnarrow + `b_correct_anscontrastive:incorrect_ansnarrow`,
         gb = b_Intercept + b_correct_ansgiven,
         gc = b_Intercept + b_correct_ansgiven + b_incorrect_anscontrastive + `b_correct_ansgiven:incorrect_anscontrastive`,
         gn = b_Intercept + b_correct_ansgiven + b_incorrect_ansnarrow + `b_correct_ansgiven:incorrect_ansnarrow`,
         nb = b_Intercept + b_correct_ansnarrow,
         nc = b_Intercept + b_correct_ansnarrow + b_incorrect_anscontrastive + `b_correct_ansnarrow:incorrect_anscontrastive`,
         ng = b_Intercept + b_correct_ansnarrow + b_incorrect_ansgiven + `b_correct_ansnarrow:incorrect_ansgiven`
  )


# 95% HDI of difference
correct_ans <- c(rep("Broad",3),rep("Contrastive",3),rep("Given",3),rep("Narrow",3))
incorrect_ans <- c("Contrastive", "Given", "Narrow","Broad", "Given", "Narrow", "Broad","Contrastive", "Narrow","Broad","Contrastive","Given")

bc_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$bc))[1])
bg_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$bg))[1])
bn_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$bn))[1])
cb_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$cb))[1])
cg_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$cg))[1])
cn_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$cn))[1])
gb_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$gb))[1])
gc_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$gc))[1])
gn_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$gn))[1])
nb_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$nb))[1])
nc_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$nc))[1])
ng_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$ng))[1])

lci <- c(bc_lci,bg_lci,bn_lci,cb_lci,cg_lci,cn_lci,gb_lci,gc_lci,gn_lci,nb_lci,nc_lci,ng_lci)

bc_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$bc))[2])
bg_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$bg))[2])
bn_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$bn))[2])
cb_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$cb))[2])
cg_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$cg))[2])
cn_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$cn))[2])
gb_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$gb))[2])
gc_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$gc))[2])
gn_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$gn))[2])
nb_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$nb))[2])
nc_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$nc))[2])
ng_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_d$ng))[2])

uci <- c(bc_uci,bg_uci,bn_uci,cb_uci,cg_uci,cn_uci,gb_uci,gc_uci,gn_uci,nb_uci,nc_uci,ng_uci)

bc_m <- plogis(mean(psamples_d$bc))
bg_m <- plogis(mean(psamples_d$bg))
bn_m <- plogis(mean(psamples_d$bn))
cb_m <- plogis(mean(psamples_d$cb))
cg_m <- plogis(mean(psamples_d$cg))
cn_m <- plogis(mean(psamples_d$cn))
gb_m <- plogis(mean(psamples_d$gb))
gc_m <- plogis(mean(psamples_d$gc))
gn_m <- plogis(mean(psamples_d$gn))
nb_m <- plogis(mean(psamples_d$nb))
nc_m <- plogis(mean(psamples_d$nc))
ng_m <- plogis(mean(psamples_d$ng))

m <- c(bc_m,bg_m,bn_m,cb_m,cg_m,cn_m,gb_m,gc_m,gn_m,nb_m,nc_m,ng_m)

# probability that coeff. is bigger than 0 (i.e. performance is above chance)
bc_above_chance <- length(which(psamples_d$bc > 0)) / length(psamples_d$bc)
bg_above_chance <- length(which(psamples_d$bg > 0)) / length(psamples_d$bg)
bn_above_chance <- length(which(psamples_d$bn > 0)) / length(psamples_d$bn)
cb_above_chance <- length(which(psamples_d$cb > 0)) / length(psamples_d$cb)
cg_above_chance <- length(which(psamples_d$cg > 0)) / length(psamples_d$cg)
cn_above_chance <- length(which(psamples_d$cn > 0)) / length(psamples_d$cn)
gb_above_chance <- length(which(psamples_d$gb > 0)) / length(psamples_d$gb)
gc_above_chance <- length(which(psamples_d$gc > 0)) / length(psamples_d$gc)
gn_above_chance <- length(which(psamples_d$gn > 0)) / length(psamples_d$gn)
nb_above_chance <- length(which(psamples_d$nb > 0)) / length(psamples_d$nb)
nc_above_chance <- length(which(psamples_d$nc > 0)) / length(psamples_d$nc)
ng_above_chance <- length(which(psamples_d$ng > 0)) / length(psamples_d$ng)

above_chance <- c(bc_above_chance,bg_above_chance,bn_above_chance,
                  cb_above_chance,cg_above_chance,cn_above_chance,
                  gb_above_chance,gc_above_chance,gn_above_chance,
                  nb_above_chance,nc_above_chance,ng_above_chance)

post_data_discr <- data.frame(correct_ans,incorrect_ans,lci,uci,m, above_chance)
#write.table(post_data_discr, "exp2_posterior.csv")


#############################################
#### plot posteriors for Exp1 - Figure 3 ####
#############################################

## plot discrimination
x2_discr <- 
  ggplot(post_data_discr, aes(x = incorrect_ans, y = m, colour = incorrect_ans, fill = incorrect_ans, shape = incorrect_ans)) +
  geom_segment(x = 0, y = 0.5, xend = 4, yend = 0.5, lty = "dotted", size = 1, colour = "black") +  
  geom_point(size = 3, stroke = 1) +
  geom_quasirandom(data = xfit2_discr_agg_item, aes(x = incorrect_ans, y = est), alpha = 0.2, size = 3, width = 0.2, dodge.width = 0.2) +
  #geom_line(data = xfit2_discr_agg_item, aes(x = incorrect_ans, y = est, group = pronoun), alpha = 0.3, colour = "gray") +  
  geom_errorbar(aes(ymin = lci, ymax = uci), 
                colour = "black", width = .2, position = position_dodge(0.2), size = 1) +
  geom_point(size = 3,  stroke = 1) +
  scale_fill_manual("Focus competitor", labels = c("B" = "Broad", "N" = "Narrow", "C" = "Contrastive", "G" = "Given"),
                    guide = guide_legend(title = "Focus competitor"),
                    values = c("#CC6666", "#999999", "#66CC99", "#9999CC")) +
  scale_colour_manual("Focus competitor", labels = c("B" = "Broad", "N" = "Narrow", "C" = "Contrastive", "G" = "Given"),
                      guide = guide_legend(title = "Focus competitor"),
                    values = c("#CC6666", "#999999", "#66CC99", "#9999CC")) +
  scale_shape_manual(values = c(21,22,24,25)) +
  scale_y_continuous(expand = c(0, 0), breaks = (c(0,0.25,0.5,0.75,1)), limits = c(0,1.0)) +
  facet_wrap(~ correct_ans, ncol = 2) +
  labs(title = "Results Experiment 1: one context - two prosodic forms",
       subtitle = "posterior means and 95% CIs; transparent points are item averages\n") +
  ylab("Predicted probability of \n correct response\n") +
  xlab("\n Focus competitor") +
  theme_classic() + 
  theme(legend.position = "none",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold", hjust = 0),
        title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "plain"),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"))


setwd("../Figures/")
ggsave(filename = "x2_discr.jpeg",
       plot = x2_discr,
       device = "jpeg",
       width = 200, 
       height = 200,
       units = "mm", 
       bg = "transparent",
       dpi = 1200)


############################################
#### extract posterior values for ident ####
############################################

# posterior distribution against chance in ident 
psamples_i = posterior_samples(xmdl_x2_ident, pars = rownames(brms:::fixef.brmsfit(xmdl_x2_ident))) %>% 
  mutate(bc = b_Intercept + b_incorrect_anscontrastive,
         bg = b_Intercept + b_incorrect_ansgiven,
         bn = b_Intercept + b_incorrect_ansnarrow,
         cb = b_Intercept + b_correct_anscontrastive,
         cg = b_Intercept + b_correct_anscontrastive + b_incorrect_ansgiven + `b_correct_anscontrastive:incorrect_ansgiven`,
         cn = b_Intercept + b_correct_anscontrastive + b_incorrect_ansnarrow + `b_correct_anscontrastive:incorrect_ansnarrow`,
         gb = b_Intercept + b_correct_ansgiven,
         gc = b_Intercept + b_correct_ansgiven + b_incorrect_anscontrastive + `b_correct_ansgiven:incorrect_anscontrastive`,
         gn = b_Intercept + b_correct_ansgiven + b_incorrect_ansnarrow + `b_correct_ansgiven:incorrect_ansnarrow`,
         nb = b_Intercept + b_correct_ansnarrow,
         nc = b_Intercept + b_correct_ansnarrow + b_incorrect_anscontrastive + `b_correct_ansnarrow:incorrect_anscontrastive`,
         ng = b_Intercept + b_correct_ansnarrow + b_incorrect_ansgiven + `b_correct_ansnarrow:incorrect_ansgiven`
  )


# 95% HDI of difference

correct_ans <- c(rep("Broad",3),rep("Contrastive",3),rep("Given",3),rep("Narrow",3))
incorrect_ans <- c("Contrastive", "Given", "Narrow","Broad", "Given", "Narrow", "Broad","Contrastive", "Narrow","Broad","Contrastive","Given")

bc_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$bc))[1])
bg_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$bg))[1])
bn_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$bn))[1])
cb_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$cb))[1])
cg_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$cg))[1])
cn_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$cn))[1])
gb_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$gb))[1])
gc_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$gc))[1])
gn_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$gn))[1])
nb_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$nb))[1])
nc_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$nc))[1])
ng_lci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$ng))[1])

lci <- c(bc_lci,bg_lci,bn_lci,cb_lci,cg_lci,cn_lci,gb_lci,gc_lci,gn_lci,nb_lci,nc_lci,ng_lci)

bc_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$bc))[2])
bg_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$bg))[2])
bn_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$bn))[2])
cb_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$cb))[2])
cg_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$cg))[2])
cn_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$cn))[2])
gb_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$gb))[2])
gc_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$gc))[2])
gn_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$gn))[2])
nb_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$nb))[2])
nc_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$nc))[2])
ng_uci <- plogis(coda::HPDinterval(as.mcmc(psamples_i$ng))[2])

uci <- c(bc_uci,bg_uci,bn_uci,cb_uci,cg_uci,cn_uci,gb_uci,gc_uci,gn_uci,nb_uci,nc_uci,ng_uci)

bc_m <- plogis(mean(psamples_i$bc))
bg_m <- plogis(mean(psamples_i$bg))
bn_m <- plogis(mean(psamples_i$bn))
cb_m <- plogis(mean(psamples_i$cb))
cg_m <- plogis(mean(psamples_i$cg))
cn_m <- plogis(mean(psamples_i$cn))
gb_m <- plogis(mean(psamples_i$gb))
gc_m <- plogis(mean(psamples_i$gc))
gn_m <- plogis(mean(psamples_i$gn))
nb_m <- plogis(mean(psamples_i$nb))
nc_m <- plogis(mean(psamples_i$nc))
ng_m <- plogis(mean(psamples_i$ng))

m <- c(bc_m,bg_m,bn_m,cb_m,cg_m,cn_m,gb_m,gc_m,gn_m,nb_m,nc_m,ng_m)

# probability that coeff. is bigger than 0 (i.e. performance is above chance)
bc_above_chance <- length(which(psamples_i$bc > 0)) / length(psamples_i$bc)
bg_above_chance <- length(which(psamples_i$bg > 0)) / length(psamples_i$bg)
bn_above_chance <- length(which(psamples_i$bn > 0)) / length(psamples_i$bn)
cb_above_chance <- length(which(psamples_i$cb > 0)) / length(psamples_i$cb)
cg_above_chance <- length(which(psamples_i$cg > 0)) / length(psamples_i$cg)
cn_above_chance <- length(which(psamples_i$cn > 0)) / length(psamples_i$cn)
gb_above_chance <- length(which(psamples_i$gb > 0)) / length(psamples_i$gb)
gc_above_chance <- length(which(psamples_i$gc > 0)) / length(psamples_i$gc)
gn_above_chance <- length(which(psamples_i$gn > 0)) / length(psamples_i$gn)
nb_above_chance <- length(which(psamples_i$nb > 0)) / length(psamples_i$nb)
nc_above_chance <- length(which(psamples_i$nc > 0)) / length(psamples_i$nc)
ng_above_chance <- length(which(psamples_i$ng > 0)) / length(psamples_i$ng)

above_chance <- c(bc_above_chance,bg_above_chance,bn_above_chance,
                  cb_above_chance,cg_above_chance,cn_above_chance,
                  gb_above_chance,gc_above_chance,gn_above_chance,
                  nb_above_chance,nc_above_chance,ng_above_chance)

post_data_ident <- data.frame(correct_ans,incorrect_ans,lci,uci,m, above_chance)
#write.table(post_data_ident, "exp2_posterior_ident.csv")


## save stuff for later use
setwd("../Data/")
save(xmdl_x2_discr, xmdl_x2_ident, post_data_discr, post_data_ident, file = "CMR_Bayesian_model_x2.RData")
#write.table(post_data_discr, "exp2_posterior_discr.csv")
#write.table(post_data_ident, "exp2_posterior_ident.csv")


#############################################
#### plot posteriors for Exp2 - Figure 4 ####
#############################################

x2_ident <- 
  ggplot(post_data_ident, aes(x = incorrect_ans, y = m, fill = incorrect_ans, colour = incorrect_ans, shape = incorrect_ans)) +
  geom_segment(x = 0, y = 0.5, xend = 4, yend = 0.5, lty = "dotted", size = 1, colour = "black") +  
  geom_point(size = 3, stroke = 1) +
  geom_quasirandom(data = xfit2_ident_agg_item, aes(x = incorrect_ans, y = est), alpha = 0.2, size = 3, width = 0.2, dodge.width = 0.1) +
  #geom_line(data = xfit2_discr_agg_item, aes(x = incorrect_ans, y = est, group = pronoun), alpha = 0.3, colour = "gray") +  
  geom_errorbar(aes(ymin = lci, ymax = uci), 
                colour = "black", width = .2, position = position_dodge(0.2), size = 1) +
  geom_point(size = 3,  stroke = 1) +
  scale_fill_manual("Focus competitor", labels = c("B" = "Broad", "N" = "Narrow", "C" = "Contrastive", "G" = "Given"),
                    guide = guide_legend(title = "Focus competitor"),
                    values = c("#CC6666", "#999999", "#66CC99", "#9999CC")) +
  scale_colour_manual("Focus competitor", labels = c("B" = "Broad", "N" = "Narrow", "C" = "Contrastive", "G" = "Given"),
                    guide = guide_legend(title = "Focus competitor"),
                    values = c("#CC6666", "#999999", "#66CC99", "#9999CC")) +
  scale_shape_manual(values = c(21,22,24,25)) +
  scale_y_continuous(expand = c(0, 0), breaks = (c(0,0.25,0.5,0.75,1)), limits = c(0,1.0)) +
  facet_wrap(~ correct_ans, ncol = 2) +
  labs(title = "Results Experiment 2: two contexts - one prosodic form",
       subtitle = "posterior means and 95% CIs; transparent points are item averages\n") +
  ylab("Predicted probability of \n correct response\n") +
  xlab("\n Focus competitor") +
  theme_classic() + 
  theme(legend.position = "none",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold", hjust = 0),
        title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "plain"),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = unit(c(1,1,1.5,1.2),"cm"))

setwd("../Figures/")

ggsave(filename = "x2_ident.jpeg",
       plot = x2_ident,
       device = "jpeg",
       width = 200, 
       height = 200,
       units = "mm", 
       bg = "transparent",
       dpi = 1200)
