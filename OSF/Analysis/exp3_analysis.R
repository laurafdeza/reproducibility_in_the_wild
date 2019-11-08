
## Description: Bayesian analysis and plotting posterior distributions of Experiment 3 of Roettger & Cole
## authors: Timo Roettger

#######################
#### Preprocessing ####
#######################

## libraries
library(tidyverse)
library(brms)
library(rstan)
library(stringr)
library(rstudioapi)
options(mc.cores = parallel::detectCores()) # Run on multiple cores for Bayesian regressions


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

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

## specify path
setwd("../Data/")

## load in data
x3 <- read.csv("experiment3.csv")
x3 <- x3[,1:24]
x3$participant_id <- as.character(x3$participant_id)
dummy <- as.factor(sapply(strsplit(x3$participant_id, split = '_'), function(x){paste0(x[1],x[2],x[3])}))
x3$participant_id_n <- dummy

# Reshaping function
mapFunc <- function(pair, correctResponse, userResponse) {
  pair <- toString(pair)
  correctResponse <- toString(correctResponse)
  userResponse <- toString(userResponse)
  #print(toString(userResponse));
  valueA = substr(pair, 1, 1);
  if(correctResponse == valueA) {
    valueB = substr(pair, 3, 3);
  }
  else {
    valueB = valueA;
    valueA = substr(pair, 3, 3);
  }
  
  if (substr(userResponse, 1, 1) == valueA) {
    if (nchar(userResponse) == 1)
    {
      retVal = 5;
    }
    else
    {
      retVal = 4;
    }
  }
  else
  {
    if (substr(userResponse, 1, 1) == valueB) {
      if (nchar(userResponse) == 1)
      {
        retVal = 1;
      }
      else
      {
        retVal = 2;
      }
    }
    else
    {
      retVal = 3;
    }
  }
  
  return(retVal);
}


## transform data to tidy 5-point responses
x3$five_point_response <- mapply(mapFunc, x3$stim_pair_id, x3$correct_response, x3$mapped_response_full)
x3$five_point_response <- as.ordered(x3$five_point_response)


###############
#### Model ####
###############

# define priors
prior_x3 <- c(prior(student_t(5,0,1), class = b), prior(normal(0,1), class = Intercept))

# model that shit
xmdl_x3_discr <- brm(five_point_response ~  correct_ans * incorrect_ans 
                     + (1 + correct_ans * incorrect_ans | pronoun)
                     #+ (1  | pronoun)
                     + (1 | participant_id_n)
                     , 
                     prior = prior_x3, 
                     family = cumulative("logit"),
                     control = list(adapt_delta = 0.9),
                     data = x3[x3$task_id == "discr",])

xmdl_x3_ident <- brm(five_point_response ~  correct_ans * incorrect_ans 
                     + (1 + correct_ans * incorrect_ans | pronoun)
                     + (1 | participant_id_n)
                     , 
                     prior = prior_x3, 
                     family = cumulative("logit"),
                     control = list(adapt_delta = 0.9),
                     data = x3[x3$task_id == "ident",])


## save models for later use
# save(xmdl_x3_discr, 
#      xmdl_x3_ident,
#      file = "CMR_Bayesian_model_x3.RData")

############################################
#### extract posterior values for discr ####
############################################

# load bayesian model
setwd("../Analysis/")
load(file = "CMR_Bayesian_model_x3.RData")

psamples_3_discr = posterior_samples(xmdl_x3_discr) %>% 
  mutate(bc0  = plogis(`b_Intercept[1]` - b_incorrect_anscontrastive), 
         bc_01 = plogis(`b_Intercept[2]` - b_incorrect_anscontrastive),
         bc_02 = plogis(`b_Intercept[3]` - b_incorrect_anscontrastive),
         bc_03 = plogis(`b_Intercept[4]` - b_incorrect_anscontrastive),
         bc1  = bc_01 - bc0,
         bc2  =  bc_02 - bc_01,
         bc3  =  bc_03 - bc_02,
         bc4  =  1 - bc_03,
         bg0  = plogis(`b_Intercept[1]` - b_incorrect_ansgiven), 
         bg_01 = plogis(`b_Intercept[2]` - b_incorrect_ansgiven),
         bg_02 = plogis(`b_Intercept[3]` - b_incorrect_ansgiven),
         bg_03 = plogis(`b_Intercept[4]` - b_incorrect_ansgiven),
         bg1  = bg_01 - bg0,
         bg2  =  bg_02 - bg_01,
         bg3  =  bg_03 - bg_02,
         bg4  =  1 - bg_03,
         bn0  = plogis(`b_Intercept[1]` - b_incorrect_ansnarrow), 
         bn_01 = plogis(`b_Intercept[2]` - b_incorrect_ansnarrow),
         bn_02 = plogis(`b_Intercept[3]` - b_incorrect_ansnarrow),
         bn_03 = plogis(`b_Intercept[4]` - b_incorrect_ansnarrow),
         bn1  = bn_01 - bn0,
         bn2  =  bn_02 - bn_01,
         bn3  =  bn_03 - bn_02,
         bn4  =  1 - bn_03,
         cb0  = plogis(`b_Intercept[1]` - b_correct_anscontrastive), 
         cb_01 = plogis(`b_Intercept[2]` - b_correct_anscontrastive),
         cb_02 = plogis(`b_Intercept[3]` - b_correct_anscontrastive),
         cb_03 = plogis(`b_Intercept[4]` - b_correct_anscontrastive),
         cb1  = cb_01 - cb0,
         cb2  =  cb_02 - cb_01,
         cb3  =  cb_03 - cb_02,
         cb4  =  1 - cb_03,
         cg0  = plogis(`b_Intercept[1]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`), 
         cg_01 = plogis(`b_Intercept[2]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`),
         cg_02 = plogis(`b_Intercept[3]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`),
         cg_03 = plogis(`b_Intercept[4]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`),
         cg1  = cg_01 - cg0,
         cg2  =  cg_02 - cg_01,
         cg3  =  cg_03 - cg_02,
         cg4  =  1 - cg_03,
         cn0  = plogis(`b_Intercept[1]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`), 
         cn_01 = plogis(`b_Intercept[2]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`),
         cn_02 = plogis(`b_Intercept[3]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`),
         cn_03 = plogis(`b_Intercept[4]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`),
         cn1  = cn_01 - cn0,
         cn2  =  cn_02 - cn_01,
         cn3  =  cn_03 - cn_02,
         cn4  =  1 - cn_03,
         gb0  = plogis(`b_Intercept[1]` - b_correct_ansgiven), 
         gb_01 = plogis(`b_Intercept[2]` - b_correct_ansgiven),
         gb_02 = plogis(`b_Intercept[3]` - b_correct_ansgiven),
         gb_03 = plogis(`b_Intercept[4]` - b_correct_ansgiven),
         gb1  = gb_01 - gb0,
         gb2  =  gb_02 - gb_01,
         gb3  =  gb_03 - gb_02,
         gb4  =  1 - gb_03,
         gc0  = plogis(`b_Intercept[1]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`), 
         gc_01 = plogis(`b_Intercept[2]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`),
         gc_02 = plogis(`b_Intercept[3]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`),
         gc_03 = plogis(`b_Intercept[4]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`),
         gc1  = gc_01 - gc0,
         gc2  =  gc_02 - gc_01,
         gc3  =  gc_03 - gc_02,
         gc4  =  1 - gc_03,
         gn0  = plogis(`b_Intercept[1]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`), 
         gn_01 = plogis(`b_Intercept[2]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`),
         gn_02 = plogis(`b_Intercept[3]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`),
         gn_03 = plogis(`b_Intercept[4]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`),
         gn1  = gn_01 - gn0,
         gn2  =  gn_02 - gn_01,
         gn3  =  gn_03 - gn_02,
         gn4  =  1 - gn_03,
         nb0  = plogis(`b_Intercept[1]` - b_correct_ansnarrow), 
         nb_01 = plogis(`b_Intercept[2]` - b_correct_ansnarrow),
         nb_02 = plogis(`b_Intercept[3]` - b_correct_ansnarrow),
         nb_03 = plogis(`b_Intercept[4]` - b_correct_ansnarrow),
         nb1  = nb_01 - nb0,
         nb2  =  nb_02 - nb_01,
         nb3  =  nb_03 - nb_02,
         nb4  =  1 - nb_03,
         nc0  = plogis(`b_Intercept[1]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`), 
         nc_01 = plogis(`b_Intercept[2]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`),
         nc_02 = plogis(`b_Intercept[3]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`),
         nc_03 = plogis(`b_Intercept[4]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`),
         nc1  = nc_01 - nc0,
         nc2  =  nc_02 - nc_01,
         nc3  =  nc_03 - nc_02,
         nc4  =  1 - nc_03,
         ng0  = plogis(`b_Intercept[1]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`), 
         ng_01 = plogis(`b_Intercept[2]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`),
         ng_02 = plogis(`b_Intercept[3]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`),
         ng_03 = plogis(`b_Intercept[4]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`),
         ng1  = ng_01 - ng0,
         ng2  =  ng_02 - ng_01,
         ng3  =  ng_03 - ng_02,
         ng4  =  1 - ng_03)


# 95% HDI of difference

correct_ans <- c(rep("Broad",15),rep("Contrastive",15),rep("Given",15),rep("Narrow",15))

incorrect_ans <- c(rep("Contrastive",5), rep("Given",5), rep("Narrow",5), 
                   rep("Broad",5), rep("Given",5), rep("Narrow",5),
                   rep("Broad",5), rep("Contrastive",5), rep("Narrow",5),
                   rep("Broad",5), rep("Contrastive",5), rep("Given",5))


fpc <- rep(c("1: never",
             "2: dispreferred",
             "3: equal",
             "4: preferred",
             "5: always"), 12)

fpc_n <- rep(c(1:5), 12)

bc0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc0))[1]
bc1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc1))[1]
bc2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc2))[1]
bc3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc3))[1]
bc4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc4))[1]
bg0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg0))[1]
bg1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg1))[1]
bg2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg2))[1]
bg3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg3))[1]
bg4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg4))[1]
bn0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn0))[1]
bn1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn1))[1]
bn2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn2))[1]
bn3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn3))[1]
bn4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn4))[1]
cb0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb0))[1]
cb1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb1))[1]
cb2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb2))[1]
cb3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb3))[1]
cb4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb4))[1]
cg0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg0))[1]
cg1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg1))[1]
cg2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg2))[1]
cg3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg3))[1]
cg4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg4))[1]
cn0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn0))[1]
cn1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn1))[1]
cn2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn2))[1]
cn3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn3))[1]
cn4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn4))[1]
gb0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb0))[1]
gb1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb1))[1]
gb2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb2))[1]
gb3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb3))[1]
gb4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb4))[1]
gc0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc0))[1]
gc1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc1))[1]
gc2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc2))[1]
gc3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc3))[1]
gc4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc4))[1]
gn0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn0))[1]
gn1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn1))[1]
gn2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn2))[1]
gn3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn3))[1]
gn4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn4))[1]
nb0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb0))[1]
nb1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb1))[1]
nb2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb2))[1]
nb3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb3))[1]
nb4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb4))[1]
nc0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc0))[1]
nc1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc1))[1]
nc2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc2))[1]
nc3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc3))[1]
nc4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc4))[1]
ng0_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng0))[1]
ng1_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng1))[1]
ng2_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng2))[1]
ng3_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng3))[1]
ng4_lci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng4))[1]

lci <- c(bc0_lci,bc1_lci,bc2_lci,bc3_lci,bc4_lci,
         bg0_lci,bg1_lci,bg2_lci,bg3_lci,bg4_lci,
         bn0_lci,bn1_lci,bn2_lci,bn3_lci,bn4_lci,
         cb0_lci,cb1_lci,cb2_lci,cb3_lci,cb4_lci,
         cg0_lci,cg1_lci,cg2_lci,cg3_lci,cg4_lci,
         cn0_lci,cn1_lci,cn2_lci,cn3_lci,cn4_lci,
         gb0_lci,gb1_lci,gb2_lci,gb3_lci,gb4_lci,
         gc0_lci,gc1_lci,gc2_lci,gc3_lci,gc4_lci,
         gn0_lci,gn1_lci,gn2_lci,gn3_lci,gn4_lci,
         nb0_lci,nb1_lci,nb2_lci,nb3_lci,nb4_lci,
         nc0_lci,nc1_lci,nc2_lci,nc3_lci,nc4_lci,
         ng0_lci,ng1_lci,ng2_lci,ng3_lci,ng4_lci
         )

bc0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc0))[2]
bc1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc1))[2]
bc2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc2))[2]
bc3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc3))[2]
bc4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bc4))[2]
bg0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg0))[2]
bg1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg1))[2]
bg2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg2))[2]
bg3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg3))[2]
bg4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bg4))[2]
bn0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn0))[2]
bn1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn1))[2]
bn2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn2))[2]
bn3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn3))[2]
bn4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$bn4))[2]
cb0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb0))[2]
cb1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb1))[2]
cb2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb2))[2]
cb3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb3))[2]
cb4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cb4))[2]
cg0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg0))[2]
cg1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg1))[2]
cg2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg2))[2]
cg3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg3))[2]
cg4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cg4))[2]
cn0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn0))[2]
cn1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn1))[2]
cn2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn2))[2]
cn3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn3))[2]
cn4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$cn4))[2]
gb0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb0))[2]
gb1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb1))[2]
gb2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb2))[2]
gb3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb3))[2]
gb4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gb4))[2]
gc0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc0))[2]
gc1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc1))[2]
gc2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc2))[2]
gc3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc3))[2]
gc4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gc4))[2]
gn0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn0))[2]
gn1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn1))[2]
gn2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn2))[2]
gn3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn3))[2]
gn4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$gn4))[2]
nb0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb0))[2]
nb1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb1))[2]
nb2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb2))[2]
nb3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb3))[2]
nb4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nb4))[2]
nc0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc0))[2]
nc1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc1))[2]
nc2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc2))[2]
nc3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc3))[2]
nc4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$nc4))[2]
ng0_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng0))[2]
ng1_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng1))[2]
ng2_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng2))[2]
ng3_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng3))[2]
ng4_uci <- coda::HPDinterval(as.mcmc(psamples_3_discr$ng4))[2]

uci <- c(bc0_uci,bc1_uci,bc2_uci,bc3_uci,bc4_uci,
         bg0_uci,bg1_uci,bg2_uci,bg3_uci,bg4_uci,
         bn0_uci,bn1_uci,bn2_uci,bn3_uci,bn4_uci,
         cb0_uci,cb1_uci,cb2_uci,cb3_uci,cb4_uci,
         cg0_uci,cg1_uci,cg2_uci,cg3_uci,cg4_uci,
         cn0_uci,cn1_uci,cn2_uci,cn3_uci,cn4_uci,
         gb0_uci,gb1_uci,gb2_uci,gb3_uci,gb4_uci,
         gc0_uci,gc1_uci,gc2_uci,gc3_uci,gc4_uci,
         gn0_uci,gn1_uci,gn2_uci,gn3_uci,gn4_uci,
         nb0_uci,nb1_uci,nb2_uci,nb3_uci,nb4_uci,
         nc0_uci,nc1_uci,nc2_uci,nc3_uci,nc4_uci,
         ng0_uci,ng1_uci,ng2_uci,ng3_uci,ng4_uci
)

bc0_m <- mean(psamples_3_discr$bc0)
bc1_m <- mean(psamples_3_discr$bc1)
bc2_m <- mean(psamples_3_discr$bc2)
bc3_m <- mean(psamples_3_discr$bc3)
bc4_m <- mean(psamples_3_discr$bc4)
bg0_m <- mean(psamples_3_discr$bg0)
bg1_m <- mean(psamples_3_discr$bg1)
bg2_m <- mean(psamples_3_discr$bg2)
bg3_m <- mean(psamples_3_discr$bg3)
bg4_m <- mean(psamples_3_discr$bg4)
bn0_m <- mean(psamples_3_discr$bn0)
bn1_m <- mean(psamples_3_discr$bn1)
bn2_m <- mean(psamples_3_discr$bn2)
bn3_m <- mean(psamples_3_discr$bn3)
bn4_m <- mean(psamples_3_discr$bn4)

cb0_m <- mean(psamples_3_discr$cb0)
cb1_m <- mean(psamples_3_discr$cb1)
cb2_m <- mean(psamples_3_discr$cb2)
cb3_m <- mean(psamples_3_discr$cb3)
cb4_m <- mean(psamples_3_discr$cb4)
cg0_m <- mean(psamples_3_discr$cg0)
cg1_m <- mean(psamples_3_discr$cg1)
cg2_m <- mean(psamples_3_discr$cg2)
cg3_m <- mean(psamples_3_discr$cg3)
cg4_m <- mean(psamples_3_discr$cg4)
cn0_m <- mean(psamples_3_discr$cn0)
cn1_m <- mean(psamples_3_discr$cn1)
cn2_m <- mean(psamples_3_discr$cn2)
cn3_m <- mean(psamples_3_discr$cn3)
cn4_m <- mean(psamples_3_discr$cn4)
gb0_m <- mean(psamples_3_discr$gb0)
gb1_m <- mean(psamples_3_discr$gb1)
gb2_m <- mean(psamples_3_discr$gb2)
gb3_m <- mean(psamples_3_discr$gb3)
gb4_m <- mean(psamples_3_discr$gb4)
gc0_m <- mean(psamples_3_discr$gc0)
gc1_m <- mean(psamples_3_discr$gc1)
gc2_m <- mean(psamples_3_discr$gc2)
gc3_m <- mean(psamples_3_discr$gc3)
gc4_m <- mean(psamples_3_discr$gc4)
gn0_m <- mean(psamples_3_discr$gn0)
gn1_m <- mean(psamples_3_discr$gn1)
gn2_m <- mean(psamples_3_discr$gn2)
gn3_m <- mean(psamples_3_discr$gn3)
gn4_m <- mean(psamples_3_discr$gn4)
nb0_m <- mean(psamples_3_discr$nb0)
nb1_m <- mean(psamples_3_discr$nb1)
nb2_m <- mean(psamples_3_discr$nb2)
nb3_m <- mean(psamples_3_discr$nb3)
nb4_m <- mean(psamples_3_discr$nb4)
nc0_m <- mean(psamples_3_discr$nc0)
nc1_m <- mean(psamples_3_discr$nc1)
nc2_m <- mean(psamples_3_discr$nc2)
nc3_m <- mean(psamples_3_discr$nc3)
nc4_m <- mean(psamples_3_discr$nc4)
ng0_m <- mean(psamples_3_discr$ng0)
ng1_m <- mean(psamples_3_discr$ng1)
ng2_m <- mean(psamples_3_discr$ng2)
ng3_m <- mean(psamples_3_discr$ng3)
ng4_m <- mean(psamples_3_discr$ng4)

m <- c(bc0_m,bc1_m,bc2_m,bc3_m,bc4_m,
       bg0_m,bg1_m,bg2_m,bg3_m,bg4_m,
       bn0_m,bn1_m,bn2_m,bn3_m,bn4_m,
       cb0_m,cb1_m,cb2_m,cb3_m,cb4_m,
       cg0_m,cg1_m,cg2_m,cg3_m,cg4_m,
       cn0_m,cn1_m,cn2_m,cn3_m,cn4_m,
       gb0_m,gb1_m,gb2_m,gb3_m,gb4_m,
       gc0_m,gc1_m,gc2_m,gc3_m,gc4_m,
       gn0_m,gn1_m,gn2_m,gn3_m,gn4_m,
       nb0_m,nb1_m,nb2_m,nb3_m,nb4_m,
       nc0_m,nc1_m,nc2_m,nc3_m,nc4_m,
       ng0_m,ng1_m,ng2_m,ng3_m,ng4_m)

# probability that coeff. is bigger than 0 (i.e. choice is above chance)
bc0_above_chance <- length(which(psamples_3_discr$bc0 > 0.2)) / length(psamples_3_discr$bc0)
bg0_above_chance <- length(which(psamples_3_discr$bg0 > 0.2)) / length(psamples_3_discr$bg0)
bn0_above_chance <- length(which(psamples_3_discr$bn0 > 0.2)) / length(psamples_3_discr$bn0)
cb0_above_chance <- length(which(psamples_3_discr$cb0 > 0.2)) / length(psamples_3_discr$cb0)
cg0_above_chance <- length(which(psamples_3_discr$cg0 > 0.2)) / length(psamples_3_discr$cg0)
cn0_above_chance <- length(which(psamples_3_discr$cn0 > 0.2)) / length(psamples_3_discr$cn0)
gb0_above_chance <- length(which(psamples_3_discr$gb0 > 0.2)) / length(psamples_3_discr$gb0)
gc0_above_chance <- length(which(psamples_3_discr$gc0 > 0.2)) / length(psamples_3_discr$gc0)
gn0_above_chance <- length(which(psamples_3_discr$gn0 > 0.2)) / length(psamples_3_discr$gn0)
nb0_above_chance <- length(which(psamples_3_discr$nb0 > 0.2)) / length(psamples_3_discr$nb0)
nc0_above_chance <- length(which(psamples_3_discr$nc0 > 0.2)) / length(psamples_3_discr$nc0)
ng0_above_chance <- length(which(psamples_3_discr$ng0 > 0.2)) / length(psamples_3_discr$ng0)

bc1_above_chance <- length(which(psamples_3_discr$bc1 > 0.2)) / length(psamples_3_discr$bc1)
bg1_above_chance <- length(which(psamples_3_discr$bg1 > 0.2)) / length(psamples_3_discr$bg1)
bn1_above_chance <- length(which(psamples_3_discr$bn1 > 0.2)) / length(psamples_3_discr$bn1)
cb1_above_chance <- length(which(psamples_3_discr$cb1 > 0.2)) / length(psamples_3_discr$cb1)
cg1_above_chance <- length(which(psamples_3_discr$cg1 > 0.2)) / length(psamples_3_discr$cg1)
cn1_above_chance <- length(which(psamples_3_discr$cn1 > 0.2)) / length(psamples_3_discr$cn1)
gb1_above_chance <- length(which(psamples_3_discr$gb1 > 0.2)) / length(psamples_3_discr$gb1)
gc1_above_chance <- length(which(psamples_3_discr$gc1 > 0.2)) / length(psamples_3_discr$gc1)
gn1_above_chance <- length(which(psamples_3_discr$gn1 > 0.2)) / length(psamples_3_discr$gn1)
nb1_above_chance <- length(which(psamples_3_discr$nb1 > 0.2)) / length(psamples_3_discr$nb1)
nc1_above_chance <- length(which(psamples_3_discr$nc1 > 0.2)) / length(psamples_3_discr$nc1)
ng1_above_chance <- length(which(psamples_3_discr$ng1 > 0.2)) / length(psamples_3_discr$ng1)
bc2_above_chance <- length(which(psamples_3_discr$bc2 > 0.2)) / length(psamples_3_discr$bc2)
bg2_above_chance <- length(which(psamples_3_discr$bg2 > 0.2)) / length(psamples_3_discr$bg2)
bn2_above_chance <- length(which(psamples_3_discr$bn2 > 0.2)) / length(psamples_3_discr$bn2)
cb2_above_chance <- length(which(psamples_3_discr$cb2 > 0.2)) / length(psamples_3_discr$cb2)
cg2_above_chance <- length(which(psamples_3_discr$cg2 > 0.2)) / length(psamples_3_discr$cg2)
cn2_above_chance <- length(which(psamples_3_discr$cn2 > 0.2)) / length(psamples_3_discr$cn2)
gb2_above_chance <- length(which(psamples_3_discr$gb2 > 0.2)) / length(psamples_3_discr$gb2)
gc2_above_chance <- length(which(psamples_3_discr$gc2 > 0.2)) / length(psamples_3_discr$gc2)
gn2_above_chance <- length(which(psamples_3_discr$gn2 > 0.2)) / length(psamples_3_discr$gn2)
nb2_above_chance <- length(which(psamples_3_discr$nb2 > 0.2)) / length(psamples_3_discr$nb2)
nc2_above_chance <- length(which(psamples_3_discr$nc2 > 0.2)) / length(psamples_3_discr$nc2)
ng2_above_chance <- length(which(psamples_3_discr$ng2 > 0.2)) / length(psamples_3_discr$ng2)
bc3_above_chance <- length(which(psamples_3_discr$bc3 > 0.2)) / length(psamples_3_discr$bc3)
bg3_above_chance <- length(which(psamples_3_discr$bg3 > 0.2)) / length(psamples_3_discr$bg3)
bn3_above_chance <- length(which(psamples_3_discr$bn3 > 0.2)) / length(psamples_3_discr$bn3)
cb3_above_chance <- length(which(psamples_3_discr$cb3 > 0.2)) / length(psamples_3_discr$cb3)
cg3_above_chance <- length(which(psamples_3_discr$cg3 > 0.2)) / length(psamples_3_discr$cg3)
cn3_above_chance <- length(which(psamples_3_discr$cn3 > 0.2)) / length(psamples_3_discr$cn3)
gb3_above_chance <- length(which(psamples_3_discr$gb3 > 0.2)) / length(psamples_3_discr$gb3)
gc3_above_chance <- length(which(psamples_3_discr$gc3 > 0.2)) / length(psamples_3_discr$gc3)
gn3_above_chance <- length(which(psamples_3_discr$gn3 > 0.2)) / length(psamples_3_discr$gn3)
nb3_above_chance <- length(which(psamples_3_discr$nb3 > 0.2)) / length(psamples_3_discr$nb3)
nc3_above_chance <- length(which(psamples_3_discr$nc3 > 0.2)) / length(psamples_3_discr$nc3)
ng3_above_chance <- length(which(psamples_3_discr$ng3 > 0.2)) / length(psamples_3_discr$ng3)
bc4_above_chance <- length(which(psamples_3_discr$bc4 > 0.2)) / length(psamples_3_discr$bc4)
bg4_above_chance <- length(which(psamples_3_discr$bg4 > 0.2)) / length(psamples_3_discr$bg4)
bn4_above_chance <- length(which(psamples_3_discr$bn4 > 0.2)) / length(psamples_3_discr$bn4)
cb4_above_chance <- length(which(psamples_3_discr$cb4 > 0.2)) / length(psamples_3_discr$cb4)
cg4_above_chance <- length(which(psamples_3_discr$cg4 > 0.2)) / length(psamples_3_discr$cg4)
cn4_above_chance <- length(which(psamples_3_discr$cn4 > 0.2)) / length(psamples_3_discr$cn4)
gb4_above_chance <- length(which(psamples_3_discr$gb4 > 0.2)) / length(psamples_3_discr$gb4)
gc4_above_chance <- length(which(psamples_3_discr$gc4 > 0.2)) / length(psamples_3_discr$gc4)
gn4_above_chance <- length(which(psamples_3_discr$gn4 > 0.2)) / length(psamples_3_discr$gn4)
nb4_above_chance <- length(which(psamples_3_discr$nb4 > 0.2)) / length(psamples_3_discr$nb4)
nc4_above_chance <- length(which(psamples_3_discr$nc4 > 0.2)) / length(psamples_3_discr$nc4)
ng4_above_chance <- length(which(psamples_3_discr$ng4 > 0.2)) / length(psamples_3_discr$ng4)

above_chance <- c(bc0_above_chance,bc1_above_chance,bc2_above_chance,bc3_above_chance,bc4_above_chance,
                  bg0_above_chance,bg1_above_chance,bg2_above_chance,bg3_above_chance,bg4_above_chance,
                  bn0_above_chance,bn1_above_chance,bn2_above_chance,bn3_above_chance,bn4_above_chance,
                  cb0_above_chance,cb1_above_chance,cb2_above_chance,cb3_above_chance,cb4_above_chance,
                  cg0_above_chance,cg1_above_chance,cg2_above_chance,cg3_above_chance,cg4_above_chance,
                  cn0_above_chance,cn1_above_chance,cn2_above_chance,cn3_above_chance,cn4_above_chance,
                  gb0_above_chance,gb1_above_chance,gb2_above_chance,gb3_above_chance,gb4_above_chance,
                  gc0_above_chance,gc1_above_chance,gc2_above_chance,gc3_above_chance,gc4_above_chance,
                  gn0_above_chance,gn1_above_chance,gn2_above_chance,gn3_above_chance,gn4_above_chance,
                  nb0_above_chance,nb1_above_chance,nb2_above_chance,nb3_above_chance,nb4_above_chance,
                  nc0_above_chance,nc1_above_chance,nc2_above_chance,nc3_above_chance,nc4_above_chance,
                  ng0_above_chance,ng1_above_chance,ng2_above_chance,ng3_above_chance,ng4_above_chance)

post_data_discr <- data.frame(correct_ans,incorrect_ans,lci,uci,m,fpc,fpc_n, above_chance)

#write.table(post_data_discr, "exp3_posterior_discr.csv")


###################################################
#### plotting stcked barplot for discr - Fig 5 ####
###################################################

df <- post_data_discr

# calculate halves of the neutral category
df.split <-
  df %>% 
  filter(fpc_n == 3) %>% 
  mutate(m = as.numeric(m/2)) 

# replace old neutral-category
df <- 
  df %>% 
  filter(!fpc_n == 3)

df <- full_join(df,df.split) %>% 
  arrange(fpc_n) %>% 
  arrange(desc(fpc_n)) 

#split dataframe
df1 <- df %>% 
  filter(fpc_n == 3 | fpc_n== 2 | fpc_n==1) %>% 
  mutate(m = m *-1)

df2 <- df %>% 
  filter(fpc_n == 5 | fpc_n== 4 | fpc_n==3)

# reorder factor "fpc_n"
df1$fpc_n  <- factor(df1$fpc_n, levels=rev(unique(df1$fpc_n)))
df2$fpc_n  <- factor(df2$fpc_n, levels=rev(unique(df2$fpc_n)))

df2$fpc <- factor(df2$fpc, levels = rev(levels(df2$fpc)))
levels(df2$incorrect_ans) <- c("Broad", "Contrastive", "Given", "Narrow")
levels(df2$correct_ans) <- c("Broad", "Contrastive", "Given", "Narrow")


#Plot  
x3_stacked_discr <-
  ggplot() +
  geom_hline(yintercept = 0, color =c("black"))+
  geom_bar(data = df1, aes(x = incorrect_ans, y = m * 100, fill = fpc), position = "stack", stat = "identity") +
  geom_bar(data = df2, aes(x = incorrect_ans, y = m * 100, fill = fpc), position = "stack", stat = "identity") +
  #geom_text(data = df1, aes(x = incorrect_ans, label = round(m *-100,0), y = round(m * 100,0)), size = 5, colour = "black", hjust = 0.5, position = position_stack(vjust = 0.5)) +
  #geom_text(data = df2, aes(x = incorrect_ans, label = round(m *100,0), y = round(m * 100,0)), size = 5, colour = "black", hjust = 0.5, position = position_stack(vjust = 0.5)) +
  facet_wrap( ~ correct_ans, ncol = 2) +
  theme_bw() + 
  #coord_flip() +
  guides(fill=guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette="PiYG", name="",labels = c("never","dispreferred","equal","preferred","always")) +
  ylab("Predicted percentage of response\n") +
  xlab("\n Focus competitor") +
  labs(title = "Result for Experiment 3: one context - two prosodic forms",
       subtitle = "posterior means centered around the 'equal' category\n") +
  scale_y_continuous(expand = c(0, 0), breaks = (c(-80,-50,0,50,100)), labels = c("80%","50%","0%","50%","100%"), limits = c(-80,100)) +
  theme_classic() + 
  theme(legend.position = "right",
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
ggsave(filename = "x3_stacked_discr.jpeg",
       plot = x3_stacked_discr,
       device = "jpeg",
       width = 220, 
       height = 147,
       units = "mm", 
       bg = "transparent",
       dpi = 1200)



############################################
#### extract posterior values for ident ####
############################################

psamples_3_ident = posterior_samples(xmdl_x3_ident) %>% 
  mutate(bc0 = plogis(`b_Intercept[1]` - b_incorrect_anscontrastive), 
         bc_01 = plogis(`b_Intercept[2]` - b_incorrect_anscontrastive),
         bc_02 = plogis(`b_Intercept[3]` - b_incorrect_anscontrastive),
         bc_03 = plogis(`b_Intercept[4]` - b_incorrect_anscontrastive),
         bc1  = bc_01 - bc0,
         bc2  =  bc_02 - bc_01,
         bc3  =  bc_03 - bc_02,
         bc4  =  1 - bc_03,
         bg0  = plogis(`b_Intercept[1]` - b_incorrect_ansgiven), 
         bg_01 = plogis(`b_Intercept[2]` - b_incorrect_ansgiven),
         bg_02 = plogis(`b_Intercept[3]` - b_incorrect_ansgiven),
         bg_03 = plogis(`b_Intercept[4]` - b_incorrect_ansgiven),
         bg1  = bg_01 - bg0,
         bg2  =  bg_02 - bg_01,
         bg3  =  bg_03 - bg_02,
         bg4  =  1 - bg_03,
         bn0  = plogis(`b_Intercept[1]` - b_incorrect_ansnarrow), 
         bn_01 = plogis(`b_Intercept[2]` - b_incorrect_ansnarrow),
         bn_02 = plogis(`b_Intercept[3]` - b_incorrect_ansnarrow),
         bn_03 = plogis(`b_Intercept[4]` - b_incorrect_ansnarrow),
         bn1  = bn_01 - bn0,
         bn2  =  bn_02 - bn_01,
         bn3  =  bn_03 - bn_02,
         bn4  =  1 - bn_03,
         cb0  = plogis(`b_Intercept[1]` - b_correct_anscontrastive), 
         cb_01 = plogis(`b_Intercept[2]` - b_correct_anscontrastive),
         cb_02 = plogis(`b_Intercept[3]` - b_correct_anscontrastive),
         cb_03 = plogis(`b_Intercept[4]` - b_correct_anscontrastive),
         cb1  = cb_01 - cb0,
         cb2  =  cb_02 - cb_01,
         cb3  =  cb_03 - cb_02,
         cb4  =  1 - cb_03,
         cg0  = plogis(`b_Intercept[1]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`), 
         cg_01 = plogis(`b_Intercept[2]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`),
         cg_02 = plogis(`b_Intercept[3]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`),
         cg_03 = plogis(`b_Intercept[4]` - b_correct_anscontrastive - b_incorrect_ansgiven - `b_correct_anscontrastive:incorrect_ansgiven`),
         cg1  = cg_01 - cg0,
         cg2  =  cg_02 - cg_01,
         cg3  =  cg_03 - cg_02,
         cg4  =  1 - cg_03,
         cn0  = plogis(`b_Intercept[1]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`), 
         cn_01 = plogis(`b_Intercept[2]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`),
         cn_02 = plogis(`b_Intercept[3]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`),
         cn_03 = plogis(`b_Intercept[4]` - b_correct_anscontrastive - b_incorrect_ansnarrow - `b_correct_anscontrastive:incorrect_ansnarrow`),
         cn1  = cn_01 - cn0,
         cn2  =  cn_02 - cn_01,
         cn3  =  cn_03 - cn_02,
         cn4  =  1 - cn_03,
         gb0  = plogis(`b_Intercept[1]` - b_correct_ansgiven), 
         gb_01 = plogis(`b_Intercept[2]` - b_correct_ansgiven),
         gb_02 = plogis(`b_Intercept[3]` - b_correct_ansgiven),
         gb_03 = plogis(`b_Intercept[4]` - b_correct_ansgiven),
         gb1  = gb_01 - gb0,
         gb2  =  gb_02 - gb_01,
         gb3  =  gb_03 - gb_02,
         gb4  =  1 - gb_03,
         gc0  = plogis(`b_Intercept[1]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`), 
         gc_01 = plogis(`b_Intercept[2]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`),
         gc_02 = plogis(`b_Intercept[3]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`),
         gc_03 = plogis(`b_Intercept[4]` - b_correct_ansgiven - b_incorrect_anscontrastive - `b_correct_ansgiven:incorrect_anscontrastive`),
         gc1  = gc_01 - gc0,
         gc2  =  gc_02 - gc_01,
         gc3  =  gc_03 - gc_02,
         gc4  =  1 - gc_03,
         gn0  = plogis(`b_Intercept[1]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`), 
         gn_01 = plogis(`b_Intercept[2]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`),
         gn_02 = plogis(`b_Intercept[3]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`),
         gn_03 = plogis(`b_Intercept[4]` - b_correct_ansgiven - b_incorrect_ansnarrow - `b_correct_ansgiven:incorrect_ansnarrow`),
         gn1  = gn_01 - gn0,
         gn2  =  gn_02 - gn_01,
         gn3  =  gn_03 - gn_02,
         gn4  =  1 - gn_03,
         nb0  = plogis(`b_Intercept[1]` - b_correct_ansnarrow), 
         nb_01 = plogis(`b_Intercept[2]` - b_correct_ansnarrow),
         nb_02 = plogis(`b_Intercept[3]` - b_correct_ansnarrow),
         nb_03 = plogis(`b_Intercept[4]` - b_correct_ansnarrow),
         nb1  = nb_01 - nb0,
         nb2  =  nb_02 - nb_01,
         nb3  =  nb_03 - nb_02,
         nb4  =  1 - nb_03,
         nc0  = plogis(`b_Intercept[1]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`), 
         nc_01 = plogis(`b_Intercept[2]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`),
         nc_02 = plogis(`b_Intercept[3]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`),
         nc_03 = plogis(`b_Intercept[4]` - b_correct_ansnarrow - b_incorrect_anscontrastive - `b_correct_ansnarrow:incorrect_anscontrastive`),
         nc1  = nc_01 - nc0,
         nc2  =  nc_02 - nc_01,
         nc3  =  nc_03 - nc_02,
         nc4  =  1 - nc_03,
         ng0  = plogis(`b_Intercept[1]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`), 
         ng_01 = plogis(`b_Intercept[2]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`),
         ng_02 = plogis(`b_Intercept[3]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`),
         ng_03 = plogis(`b_Intercept[4]` - b_correct_ansnarrow - b_incorrect_ansgiven - `b_correct_ansnarrow:incorrect_ansgiven`),
         ng1  = ng_01 - ng0,
         ng2  =  ng_02 - ng_01,
         ng3  =  ng_03 - ng_02,
         ng4  =  1 - ng_03)


# 95% HDI of difference

correct_ans <- c(rep("Broad",15),rep("Contrastive",15),rep("Given",15),rep("Narrow",15))

incorrect_ans <- c(rep("Contrastive",5), rep("Given",5), rep("Narrow",5), 
                   rep("Broad",5), rep("Given",5), rep("Narrow",5),
                   rep("Broad",5), rep("Contrastive",5), rep("Narrow",5),
                   rep("Broad",5), rep("Contrastive",5), rep("Given",5))

fpc <- rep(c("1: never",
             "2: dispreferred",
             "3: equal",
             "4: preferred",
             "5: always"), 12)

fpc_n <- rep(c(1:5), 12)

bc0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc0))[1]
bc1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc1))[1]
bc2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc2))[1]
bc3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc3))[1]
bc4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc4))[1]
bg0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg0))[1]
bg1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg1))[1]
bg2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg2))[1]
bg3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg3))[1]
bg4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg4))[1]
bn0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn0))[1]
bn1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn1))[1]
bn2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn2))[1]
bn3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn3))[1]
bn4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn4))[1]
cb0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb0))[1]
cb1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb1))[1]
cb2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb2))[1]
cb3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb3))[1]
cb4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb4))[1]
cg0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg0))[1]
cg1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg1))[1]
cg2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg2))[1]
cg3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg3))[1]
cg4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg4))[1]
cn0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn0))[1]
cn1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn1))[1]
cn2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn2))[1]
cn3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn3))[1]
cn4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn4))[1]
gb0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb0))[1]
gb1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb1))[1]
gb2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb2))[1]
gb3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb3))[1]
gb4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb4))[1]
gc0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc0))[1]
gc1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc1))[1]
gc2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc2))[1]
gc3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc3))[1]
gc4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc4))[1]
gn0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn0))[1]
gn1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn1))[1]
gn2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn2))[1]
gn3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn3))[1]
gn4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn4))[1]
nb0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb0))[1]
nb1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb1))[1]
nb2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb2))[1]
nb3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb3))[1]
nb4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb4))[1]
nc0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc0))[1]
nc1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc1))[1]
nc2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc2))[1]
nc3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc3))[1]
nc4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc4))[1]
ng0_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng0))[1]
ng1_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng1))[1]
ng2_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng2))[1]
ng3_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng3))[1]
ng4_lci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng4))[1]

lci <- c(bc0_lci,bc1_lci,bc2_lci,bc3_lci,bc4_lci,
         bg0_lci,bg1_lci,bg2_lci,bg3_lci,bg4_lci,
         bn0_lci,bn1_lci,bn2_lci,bn3_lci,bn4_lci,
         cb0_lci,cb1_lci,cb2_lci,cb3_lci,cb4_lci,
         cg0_lci,cg1_lci,cg2_lci,cg3_lci,cg4_lci,
         cn0_lci,cn1_lci,cn2_lci,cn3_lci,cn4_lci,
         gb0_lci,gb1_lci,gb2_lci,gb3_lci,gb4_lci,
         gc0_lci,gc1_lci,gc2_lci,gc3_lci,gc4_lci,
         gn0_lci,gn1_lci,gn2_lci,gn3_lci,gn4_lci,
         nb0_lci,nb1_lci,nb2_lci,nb3_lci,nb4_lci,
         nc0_lci,nc1_lci,nc2_lci,nc3_lci,nc4_lci,
         ng0_lci,ng1_lci,ng2_lci,ng3_lci,ng4_lci
)

bc0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc0))[2]
bc1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc1))[2]
bc2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc2))[2]
bc3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc3))[2]
bc4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bc4))[2]
bg0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg0))[2]
bg1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg1))[2]
bg2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg2))[2]
bg3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg3))[2]
bg4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bg4))[2]
bn0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn0))[2]
bn1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn1))[2]
bn2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn2))[2]
bn3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn3))[2]
bn4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$bn4))[2]
cb0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb0))[2]
cb1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb1))[2]
cb2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb2))[2]
cb3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb3))[2]
cb4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cb4))[2]
cg0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg0))[2]
cg1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg1))[2]
cg2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg2))[2]
cg3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg3))[2]
cg4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cg4))[2]
cn0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn0))[2]
cn1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn1))[2]
cn2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn2))[2]
cn3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn3))[2]
cn4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$cn4))[2]
gb0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb0))[2]
gb1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb1))[2]
gb2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb2))[2]
gb3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb3))[2]
gb4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gb4))[2]
gc0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc0))[2]
gc1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc1))[2]
gc2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc2))[2]
gc3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc3))[2]
gc4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gc4))[2]
gn0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn0))[2]
gn1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn1))[2]
gn2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn2))[2]
gn3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn3))[2]
gn4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$gn4))[2]
nb0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb0))[2]
nb1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb1))[2]
nb2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb2))[2]
nb3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb3))[2]
nb4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nb4))[2]
nc0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc0))[2]
nc1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc1))[2]
nc2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc2))[2]
nc3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc3))[2]
nc4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$nc4))[2]
ng0_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng0))[2]
ng1_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng1))[2]
ng2_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng2))[2]
ng3_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng3))[2]
ng4_uci <- coda::HPDinterval(as.mcmc(psamples_3_ident$ng4))[2]

uci <- c(bc0_uci,bc1_uci,bc2_uci,bc3_uci,bc4_uci,
         bg0_uci,bg1_uci,bg2_uci,bg3_uci,bg4_uci,
         bn0_uci,bn1_uci,bn2_uci,bn3_uci,bn4_uci,
         cb0_uci,cb1_uci,cb2_uci,cb3_uci,cb4_uci,
         cg0_uci,cg1_uci,cg2_uci,cg3_uci,cg4_uci,
         cn0_uci,cn1_uci,cn2_uci,cn3_uci,cn4_uci,
         gb0_uci,gb1_uci,gb2_uci,gb3_uci,gb4_uci,
         gc0_uci,gc1_uci,gc2_uci,gc3_uci,gc4_uci,
         gn0_uci,gn1_uci,gn2_uci,gn3_uci,gn4_uci,
         nb0_uci,nb1_uci,nb2_uci,nb3_uci,nb4_uci,
         nc0_uci,nc1_uci,nc2_uci,nc3_uci,nc4_uci,
         ng0_uci,ng1_uci,ng2_uci,ng3_uci,ng4_uci
)

bc0_m <- mean(psamples_3_ident$bc0)
bc1_m <- mean(psamples_3_ident$bc1)
bc2_m <- mean(psamples_3_ident$bc2)
bc3_m <- mean(psamples_3_ident$bc3)
bc4_m <- mean(psamples_3_ident$bc4)
bg0_m <- mean(psamples_3_ident$bg0)
bg1_m <- mean(psamples_3_ident$bg1)
bg2_m <- mean(psamples_3_ident$bg2)
bg3_m <- mean(psamples_3_ident$bg3)
bg4_m <- mean(psamples_3_ident$bg4)
bn0_m <- mean(psamples_3_ident$bn0)
bn1_m <- mean(psamples_3_ident$bn1)
bn2_m <- mean(psamples_3_ident$bn2)
bn3_m <- mean(psamples_3_ident$bn3)
bn4_m <- mean(psamples_3_ident$bn4)
cb0_m <- mean(psamples_3_ident$cb0)
cb1_m <- mean(psamples_3_ident$cb1)
cb2_m <- mean(psamples_3_ident$cb2)
cb3_m <- mean(psamples_3_ident$cb3)
cb4_m <- mean(psamples_3_ident$cb4)
cg0_m <- mean(psamples_3_ident$cg0)
cg1_m <- mean(psamples_3_ident$cg1)
cg2_m <- mean(psamples_3_ident$cg2)
cg3_m <- mean(psamples_3_ident$cg3)
cg4_m <- mean(psamples_3_ident$cg4)
cn0_m <- mean(psamples_3_ident$cn0)
cn1_m <- mean(psamples_3_ident$cn1)
cn2_m <- mean(psamples_3_ident$cn2)
cn3_m <- mean(psamples_3_ident$cn3)
cn4_m <- mean(psamples_3_ident$cn4)
gb0_m <- mean(psamples_3_ident$gb0)
gb1_m <- mean(psamples_3_ident$gb1)
gb2_m <- mean(psamples_3_ident$gb2)
gb3_m <- mean(psamples_3_ident$gb3)
gb4_m <- mean(psamples_3_ident$gb4)
gc0_m <- mean(psamples_3_ident$gc0)
gc1_m <- mean(psamples_3_ident$gc1)
gc2_m <- mean(psamples_3_ident$gc2)
gc3_m <- mean(psamples_3_ident$gc3)
gc4_m <- mean(psamples_3_ident$gc4)
gn0_m <- mean(psamples_3_ident$gn0)
gn1_m <- mean(psamples_3_ident$gn1)
gn2_m <- mean(psamples_3_ident$gn2)
gn3_m <- mean(psamples_3_ident$gn3)
gn4_m <- mean(psamples_3_ident$gn4)
nb0_m <- mean(psamples_3_ident$nb0)
nb1_m <- mean(psamples_3_ident$nb1)
nb2_m <- mean(psamples_3_ident$nb2)
nb3_m <- mean(psamples_3_ident$nb3)
nb4_m <- mean(psamples_3_ident$nb4)
nc0_m <- mean(psamples_3_ident$nc0)
nc1_m <- mean(psamples_3_ident$nc1)
nc2_m <- mean(psamples_3_ident$nc2)
nc3_m <- mean(psamples_3_ident$nc3)
nc4_m <- mean(psamples_3_ident$nc4)
ng0_m <- mean(psamples_3_ident$ng0)
ng1_m <- mean(psamples_3_ident$ng1)
ng2_m <- mean(psamples_3_ident$ng2)
ng3_m <- mean(psamples_3_ident$ng3)
ng4_m <- mean(psamples_3_ident$ng4)

m <- c(bc0_m,bc1_m,bc2_m,bc3_m,bc4_m,
       bg0_m,bg1_m,bg2_m,bg3_m,bg4_m,
       bn0_m,bn1_m,bn2_m,bn3_m,bn4_m,
       cb0_m,cb1_m,cb2_m,cb3_m,cb4_m,
       cg0_m,cg1_m,cg2_m,cg3_m,cg4_m,
       cn0_m,cn1_m,cn2_m,cn3_m,cn4_m,
       gb0_m,gb1_m,gb2_m,gb3_m,gb4_m,
       gc0_m,gc1_m,gc2_m,gc3_m,gc4_m,
       gn0_m,gn1_m,gn2_m,gn3_m,gn4_m,
       nb0_m,nb1_m,nb2_m,nb3_m,nb4_m,
       nc0_m,nc1_m,nc2_m,nc3_m,nc4_m,
       ng0_m,ng1_m,ng2_m,ng3_m,ng4_m)


# probability that coeff. is bigger than 0 (meaning is above chance)
bc0_above_chance <- length(which(psamples_3_ident$bc0 > 0.2)) / length(psamples_3_ident$bc0)
bg0_above_chance <- length(which(psamples_3_ident$bg0 > 0.2)) / length(psamples_3_ident$bg0)
bn0_above_chance <- length(which(psamples_3_ident$bn0 > 0.2)) / length(psamples_3_ident$bn0)
cb0_above_chance <- length(which(psamples_3_ident$cb0 > 0.2)) / length(psamples_3_ident$cb0)
cg0_above_chance <- length(which(psamples_3_ident$cg0 > 0.2)) / length(psamples_3_ident$cg0)
cn0_above_chance <- length(which(psamples_3_ident$cn0 > 0.2)) / length(psamples_3_ident$cn0)
gb0_above_chance <- length(which(psamples_3_ident$gb0 > 0.2)) / length(psamples_3_ident$gb0)
gc0_above_chance <- length(which(psamples_3_ident$gc0 > 0.2)) / length(psamples_3_ident$gc0)
gn0_above_chance <- length(which(psamples_3_ident$gn0 > 0.2)) / length(psamples_3_ident$gn0)
nb0_above_chance <- length(which(psamples_3_ident$nb0 > 0.2)) / length(psamples_3_ident$nb0)
nc0_above_chance <- length(which(psamples_3_ident$nc0 > 0.2)) / length(psamples_3_ident$nc0)
ng0_above_chance <- length(which(psamples_3_ident$ng0 > 0.2)) / length(psamples_3_ident$ng0)

bc1_above_chance <- length(which(psamples_3_ident$bc1 > 0.2)) / length(psamples_3_ident$bc1)
bg1_above_chance <- length(which(psamples_3_ident$bg1 > 0.2)) / length(psamples_3_ident$bg1)
bn1_above_chance <- length(which(psamples_3_ident$bn1 > 0.2)) / length(psamples_3_ident$bn1)
cb1_above_chance <- length(which(psamples_3_ident$cb1 > 0.2)) / length(psamples_3_ident$cb1)
cg1_above_chance <- length(which(psamples_3_ident$cg1 > 0.2)) / length(psamples_3_ident$cg1)
cn1_above_chance <- length(which(psamples_3_ident$cn1 > 0.2)) / length(psamples_3_ident$cn1)
gb1_above_chance <- length(which(psamples_3_ident$gb1 > 0.2)) / length(psamples_3_ident$gb1)
gc1_above_chance <- length(which(psamples_3_ident$gc1 > 0.2)) / length(psamples_3_ident$gc1)
gn1_above_chance <- length(which(psamples_3_ident$gn1 > 0.2)) / length(psamples_3_ident$gn1)
nb1_above_chance <- length(which(psamples_3_ident$nb1 > 0.2)) / length(psamples_3_ident$nb1)
nc1_above_chance <- length(which(psamples_3_ident$nc1 > 0.2)) / length(psamples_3_ident$nc1)
ng1_above_chance <- length(which(psamples_3_ident$ng1 > 0.2)) / length(psamples_3_ident$ng1)
bc2_above_chance <- length(which(psamples_3_ident$bc2 > 0.2)) / length(psamples_3_ident$bc2)
bg2_above_chance <- length(which(psamples_3_ident$bg2 > 0.2)) / length(psamples_3_ident$bg2)
bn2_above_chance <- length(which(psamples_3_ident$bn2 > 0.2)) / length(psamples_3_ident$bn2)
cb2_above_chance <- length(which(psamples_3_ident$cb2 > 0.2)) / length(psamples_3_ident$cb2)
cg2_above_chance <- length(which(psamples_3_ident$cg2 > 0.2)) / length(psamples_3_ident$cg2)
cn2_above_chance <- length(which(psamples_3_ident$cn2 > 0.2)) / length(psamples_3_ident$cn2)
gb2_above_chance <- length(which(psamples_3_ident$gb2 > 0.2)) / length(psamples_3_ident$gb2)
gc2_above_chance <- length(which(psamples_3_ident$gc2 > 0.2)) / length(psamples_3_ident$gc2)
gn2_above_chance <- length(which(psamples_3_ident$gn2 > 0.2)) / length(psamples_3_ident$gn2)
nb2_above_chance <- length(which(psamples_3_ident$nb2 > 0.2)) / length(psamples_3_ident$nb2)
nc2_above_chance <- length(which(psamples_3_ident$nc2 > 0.2)) / length(psamples_3_ident$nc2)
ng2_above_chance <- length(which(psamples_3_ident$ng2 > 0.2)) / length(psamples_3_ident$ng2)
bc3_above_chance <- length(which(psamples_3_ident$bc3 > 0.2)) / length(psamples_3_ident$bc3)
bg3_above_chance <- length(which(psamples_3_ident$bg3 > 0.2)) / length(psamples_3_ident$bg3)
bn3_above_chance <- length(which(psamples_3_ident$bn3 > 0.2)) / length(psamples_3_ident$bn3)
cb3_above_chance <- length(which(psamples_3_ident$cb3 > 0.2)) / length(psamples_3_ident$cb3)
cg3_above_chance <- length(which(psamples_3_ident$cg3 > 0.2)) / length(psamples_3_ident$cg3)
cn3_above_chance <- length(which(psamples_3_ident$cn3 > 0.2)) / length(psamples_3_ident$cn3)
gb3_above_chance <- length(which(psamples_3_ident$gb3 > 0.2)) / length(psamples_3_ident$gb3)
gc3_above_chance <- length(which(psamples_3_ident$gc3 > 0.2)) / length(psamples_3_ident$gc3)
gn3_above_chance <- length(which(psamples_3_ident$gn3 > 0.2)) / length(psamples_3_ident$gn3)
nb3_above_chance <- length(which(psamples_3_ident$nb3 > 0.2)) / length(psamples_3_ident$nb3)
nc3_above_chance <- length(which(psamples_3_ident$nc3 > 0.2)) / length(psamples_3_ident$nc3)
ng3_above_chance <- length(which(psamples_3_ident$ng3 > 0.2)) / length(psamples_3_ident$ng3)
bc4_above_chance <- length(which(psamples_3_ident$bc4 > 0.2)) / length(psamples_3_ident$bc4)
bg4_above_chance <- length(which(psamples_3_ident$bg4 > 0.2)) / length(psamples_3_ident$bg4)
bn4_above_chance <- length(which(psamples_3_ident$bn4 > 0.2)) / length(psamples_3_ident$bn4)
cb4_above_chance <- length(which(psamples_3_ident$cb4 > 0.2)) / length(psamples_3_ident$cb4)
cg4_above_chance <- length(which(psamples_3_ident$cg4 > 0.2)) / length(psamples_3_ident$cg4)
cn4_above_chance <- length(which(psamples_3_ident$cn4 > 0.2)) / length(psamples_3_ident$cn4)
gb4_above_chance <- length(which(psamples_3_ident$gb4 > 0.2)) / length(psamples_3_ident$gb4)
gc4_above_chance <- length(which(psamples_3_ident$gc4 > 0.2)) / length(psamples_3_ident$gc4)
gn4_above_chance <- length(which(psamples_3_ident$gn4 > 0.2)) / length(psamples_3_ident$gn4)
nb4_above_chance <- length(which(psamples_3_ident$nb4 > 0.2)) / length(psamples_3_ident$nb4)
nc4_above_chance <- length(which(psamples_3_ident$nc4 > 0.2)) / length(psamples_3_ident$nc4)
ng4_above_chance <- length(which(psamples_3_ident$ng4 > 0.2)) / length(psamples_3_ident$ng4)

above_chance <- c(bc0_above_chance,bc1_above_chance,bc2_above_chance,bc3_above_chance,bc4_above_chance,
                  bg0_above_chance,bg1_above_chance,bg2_above_chance,bg3_above_chance,bg4_above_chance,
                  bn0_above_chance,bn1_above_chance,bn2_above_chance,bn3_above_chance,bn4_above_chance,
                  cb0_above_chance,cb1_above_chance,cb2_above_chance,cb3_above_chance,cb4_above_chance,
                  cg0_above_chance,cg1_above_chance,cg2_above_chance,cg3_above_chance,cg4_above_chance,
                  cn0_above_chance,cn1_above_chance,cn2_above_chance,cn3_above_chance,cn4_above_chance,
                  gb0_above_chance,gb1_above_chance,gb2_above_chance,gb3_above_chance,gb4_above_chance,
                  gc0_above_chance,gc1_above_chance,gc2_above_chance,gc3_above_chance,gc4_above_chance,
                  gn0_above_chance,gn1_above_chance,gn2_above_chance,gn3_above_chance,gn4_above_chance,
                  nb0_above_chance,nb1_above_chance,nb2_above_chance,nb3_above_chance,nb4_above_chance,
                  nc0_above_chance,nc1_above_chance,nc2_above_chance,nc3_above_chance,nc4_above_chance,
                  ng0_above_chance,ng1_above_chance,ng2_above_chance,ng3_above_chance,ng4_above_chance)

post_data_ident <- data.frame(correct_ans,incorrect_ans,lci,uci,m,fpc,fpc_n, above_chance)

#####################################################
#### plotting posterior values for ident - Fig 5 ####
#####################################################

df <- post_data_ident

# calculate halves of the neutral category
df.split <-
  df %>% 
  filter(fpc_n == 3) %>% 
  mutate(m = as.numeric(m/2)) 

# replace old neutral-category
df <- 
  df %>% 
  filter(!fpc_n == 3)

df <- full_join(df,df.split) %>% 
  arrange(fpc_n) %>% 
  arrange(desc(fpc_n)) 

#split dataframe
df1 <- df %>% 
  filter(fpc_n == 3 | fpc_n== 2 | fpc_n==1) %>% 
  mutate(m = m *-1)

df2 <- df %>% 
  filter(fpc_n == 5 | fpc_n== 4 | fpc_n==3)

# reorder factor "fpc_n"
df1$fpc_n  <- factor(df1$fpc_n, levels=rev(unique(df1$fpc_n)))
df2$fpc_n  <- factor(df2$fpc_n, levels=rev(unique(df2$fpc_n)))

df2$fpc <- factor(df2$fpc, levels = rev(levels(df2$fpc)))
levels(df2$incorrect_ans) <- c("Broad", "Contrastive", "Given", "Narrow")
levels(df2$correct_ans) <- c("Broad", "Contrastive", "Given", "Narrow")


###################################################
#### plotting stcked barplot for ident - Fig 6 ####
###################################################

x3_stacked_ident <-
  ggplot() +
  geom_hline(yintercept = 0, color =c("black"))+
  geom_bar(data = df1, aes(x = incorrect_ans, y = m * 100, fill = fpc), position = "stack", stat = "identity") +
  geom_bar(data = df2, aes(x = incorrect_ans, y = m * 100, fill = fpc), position = "stack", stat = "identity") +
  #geom_text(data = df1, aes(x = incorrect_ans, label = round(m *-100,0), y = round(m * 100,0)), size = 5, colour = "black", hjust = 0.5, position = position_stack(vjust = 0.5)) +
  #geom_text(data = df2, aes(x = incorrect_ans, label = round(m *100,0), y = round(m * 100,0)), size = 5, colour = "black", hjust = 0.5, position = position_stack(vjust = 0.5)) +
  facet_wrap( ~ correct_ans, ncol = 2) +
  theme_bw() + 
  #coord_flip() +
  guides(fill=guide_legend(title="",reverse=TRUE)) +
  scale_fill_brewer(palette="PiYG", name="",labels = c("never","dispreferred","equal","preferred","always")) +
  ylab("Predicted percentage of response\n") +
  xlab("\n Focus competitor") +
  labs(title = "Result for Experiment 3: two contexts - one prosodic form",
       subtitle = "posterior means centered around the 'equal' category\n") +
  scale_y_continuous(expand = c(0, 0), breaks = (c(-80,-50,0,50,100)), labels = c("80%","50%","0%","50%","100%"), limits = c(-80,100)) +
  theme_classic() + 
  theme(legend.position = "right",
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
ggsave(filename = "x3_stacked_ident.jpeg",
       plot = x3_stacked_ident,
       device = "png",
       width = 220, 
       height = 147,
       units = "mm", 
       #bg = "transparent",
       dpi = 500)


## save tables for later use
#save(xmdl_x3_discr, xmdl_x3_ident, post_data_discr, post_data_ident, file = "CMR_Bayesian_model_x3.RData")

## extract tables for paper
post_data_discr_t <- 
  post_data_discr %>%
  mutate(estimate = paste0(round(m,2), " (", round(lci,2), ",", round(uci,2), ")"))
  
post_data_discr_t <- 
  post_data_discr_t %>%
  select(correct_ans,incorrect_ans,fpc,estimate) %>%
  spread(fpc, estimate)

post_data_ident_t <- 
  post_data_ident %>%
  mutate(estimate = paste0(round(m,2), " (", round(lci,2), ",", round(uci,2), ")"))

post_data_ident_t <- 
  post_data_ident_t %>%
  select(correct_ans,incorrect_ans,fpc,estimate) %>%
  spread(fpc, estimate)

write.table(post_data_discr_t, "exp3_posterior_discr.csv")
write.table(post_data_ident_t, "exp3_posterior_ident.csv")

