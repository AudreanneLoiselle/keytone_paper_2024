# Audreanne Loiselle
# Created : 2023-04-13
# Last update : 2024-08-05
# keystone_paper_2024

# This script contains an example of the method we used in the following 
# manuscript :
## Adapting the concept of keystone species to diversity monitoring for plant, 
## bird, and fish communities of wetlands.

# In this paper, we proposed a modified version of the  Community Importance 
# index (CIi) to identify keystone species for any ecosystem service or function.

# We based our work on that of Avolio et al. (2019), but proposed a different
# mathematical framerwork that :
## 1) Allows the identification of species that positively or negatively
##    contribute to an ecosystem service or function
## 2) Uses a bootstrap method to identify species with significant contributions
##    and with an adaptive significance threshold.
## 3) Is adapted to the use of observational data, not only to data from removal 
##    experiments

# We tested our approach on diversity of plant, bird and fish communities in 3 
# types of wetlands, but our approach can be used for :
## 1) Any taxa
## 2) Any ecosystem type
## 3) Any ecpsystem service or function

# The present code showcases an example of our approach with one of the dataset
# used in our manuscrpti analysis : plant species.

# Plant dataset preparation ----------------------------------------------------------

# First, We need to load the species relative frequency data set, found in this 
# github repository. These values are the relative frequencies of each species
# at each sampled sites, as calculated using equation 1 in our manuscript.

plant <- read.csv("plant_DCi.csv", header=TRUE,sep=";")


# We also need the site_data dataset, also found in our repository
# It contains information about species richness and wetland type

site <- read.csv("site_data.csv", header=TRUE,sep=";")
site.plant <- site %>% 
  select(Site, wet_type, plant_rich)


# We will compute our analysis separately for the 3 types of wetlands we studied,
# as keystoneness can be affected by ecosystem characteristics. We suggest you
# identify the right ecosystem classification for your own analysis. We based
# ours on the Canadian wetland classification system, which uses plant, soil,
# and hydrology data to classify wetland classes. We also tested this
# classification using our own plant data, as explained in our paper.

# To do so we need a group vector

plant.wet <- plant %>% 
  mutate(wet_type = site.plant$wet_type)
ncol(plant.wet)


# Then, we will subset the plant dataset and the site_data dataset by type

# For alder swamps
DCi.alder.sub <- plant.wet[plant.wet$wet_type=="alder",][,-c(1,337)]
site.plant.alder <- site.plant[plant.wet$wet_type=="alder",]

# For peatlands
DCi.peat.sub <- plant.wet[plant.wet$wet_type=="peat",][,-c(1,337)]
site.plant.peat <- site.plant[plant.wet$wet_type=="peat",]

# For ash swamps
DCi.ash.sub <- plant.wet[plant.wet$wet_type=="swamp",][,-c(1,337)]
site.plant.ash <- site.plant[plant.wet$wet_type=="swamp",]


# As all species are not present in each wetland type, we need to remove columns 
# that sum to zero in each of the data subset.

# For alder swamps
DCi.alder <- DCi.alder.sub[,colSums(DCi.alder.sub)>0]

#For peatlands
DCi.peat <- DCi.peat.sub[,colSums(DCi.peat.sub)>0]

#xFor ash swamps
DCi.ash <- DCi.ash.sub[,colSums(DCi.ash.sub)>0]


# Function to compute CIi -------------------------------------------------

# Now that out datasets are ready, we need to create a function that will allow 
# us to compute the CIi values of each species for a wetland type.
# We will also integrate a plotting function to visualize our results.

# First, let's load the necessary packages
library(ggplot2)
library(ggrepel)

# Name the function
CIi.fun = function(
  
  #Name the objects to be used in the function
  ## esf = Ecosystem service or function of your choice
  ## data = the subset of relative frequencies to be used
  ## title.CIi = the title you want to give to your plot
  ## alpha = the alpha value of the confidence interval for the bootstrap
  ## yneg/ypos = the y axis scale for the plot, as it varies across datasets
  ## xmax = the maximum value of x axis, as you as it varies across datasets
  ## jump = the interval of breaks on the x axis,as it varies across datasets
  esf, data, title.CIi, alpha, yneg, ypos, xmax, jump)
  
  {
  
  # the equation we will use is as follows :
  ## CIi = [[tA-tB]/[ta+tb]]*[1/[ra-rb]]
  ## ta = esf value in site a
  ## tb = esf value in site b
  ## ra = relative frequency value of species i in site a
  ## ra = relative frequency value of species i in site b
  # So CIi gives us the impact of a species on the provisioning of an esf
  # To evaluate the global effect of the abundance of a species on an esf, 
  # we need to compute the CIi values for that species for EACH pair of 
  # sampled sites.
  # We need to repeat this for EACH species
  
  # To simplify coding, we used the following nomenclature :
  # diffesf_minus = ta-tb
  # diffesf_plus = ta+tb
  # t_i = the left side of the equation
  # p_i = the right side of the equation
  
  # The first step is thefore to calculate t_i for all sites
  set.seed(1)
  sr <- esf
  diffesf_minus <- outer(sr, sr, "-")
  diffesf_plus <- outer(sr, sr, "+")
  t_i <- (diffesf_minus/diffesf_plus)
  
  # Then, we calculate CI_i for each species at each site
  dimrc <- dim(data)
  #set the output matrix
  out_ci <- matrix(nrow = dimrc[2], ncol = 1)
  for (i in 1:dimrc[2]) {
    # calculate p_i for species "i"
    p_i <- outer(data[,i], data[,i], "-")
    # calculate all pairwise CI_i for species "i"
    CI_i <- t_i/p_i
    # remove all divisions by zero : when ra = rb
    CI_i <- CI_i[is.finite(CI_i)]
    # take the median CI_i for species "i" and store it in the output matrix
    out_ci[i] = median(CI_i) 
  }
  
  # Now we bootstrap the CI_i values for all species
  set.seed(123)
  # set the output matrix for bootsrapped values
  outs_ci <- matrix(nrow = dimrc[2], ncol = 999)
  for (j in 1: 999) { 
    sr_s <- sample(sr, replace=T); diffesf_minus <- outer(sr_s, sr_s, "-"); 
    diffesf_plus <- outer(sr_s, sr_s, "+"); ts_i <- (diffesf_minus/diffesf_plus)
    for (i in 1:dimrc[2]) {
      p_i <- outer(data[,i], data[,i], "-")
      CI_i <- ts_i/p_i
      CI_i <- CI_i[is.finite(CI_i)]
      outs_ci[i,j] = median(CI_i)
    }}
  
  # Calculate the confidence interval values for each species
  # here, you can change the confidence interval to make it more strict (99%)
  # or more flexible (95%)
  bootvalhigh <- apply(outs_ci, 1, quantile, probs = 1-alpha/2)
  bootvallow <- apply(outs_ci, 1, quantile, probs = 0+alpha/2)
  
  # Create a data frame with all the necessary data
  DC_i <- data.frame(apply(data, 2, max))
  
  all = data.frame(cbind(DC_i,out_ci,bootvalhigh,bootvallow, colnames(data)))
  colnames(all) <- c("DCi", "out_ci", "bootvalhigh ", "bootvallow","species" )
  CIi.col <- ifelse(all$out_ci>=1 & all$bootvalhigh<all$out_ci, "blue", 
                    ifelse(all$out_ci<=-1 & all$bootvallow>all$out_ci,"red","black"))
  CIi.pos <- ifelse(all$out_ci>=1 & all$bootvalhigh<all$out_ci, 1, 0)
  sum.pos <- sum(CIi.pos)
  Cii.neg <-ifelse(all$out_ci<=-1 & all$bootvallow>all$out_ci, 1, 0)
  sum.neg <- sum(Cii.neg)
  
  # Separate the positive and the negative values
  DCi.sub.pos <- all %>% 
    subset(out_ci>=1 & bootvalhigh<out_ci)
  
  DCi.sub.neg <- all %>% 
    subset(out_ci<=-1 & bootvallow>out_ci)
  
  # Generate a plot of the data
  # If the species significantly contributes to the esf tester, it's name will
  # appear in text. Otherwise, it will only be a point, representing its median
  # CIi value.
  gg = ggplot(data = all, aes(x = DCi, y = out_ci, label = species)) +
    geom_point(size = 1, col = CIi.col)+
    geom_text(data = DCi.sub.pos, aes(x = DCi, y = out_ci, label = species))+
    geom_text(data = DCi.sub.neg, aes(x = DCi, y = out_ci, label = species))+
    geom_abline(slope = 0, intercept = 0,linetype = "solid")+
    geom_abline(slope = 0, intercept = 1,linetype = "dashed", col = "blue")+
    geom_abline(slope = 0, intercept = -1,linetype = "dashed", col = "red")+
    labs(title = title.CIi, x = "Relative frequency", y = "Community Importance index")+
    ylim(yneg,ypos)+
    scale_x_continuous(limits=c(0,xmax), breaks=seq(0,xmax,jump))
  #return(all)
}

# Ok now, we can test the function. Let's try using the alder swamps data with 
# a 99% confidence interval

plant.alder.99 <- CIi.fun(site.plant.alder$plant_rich, DCi.alder, "Alder swamps",
                       0.01, -30, 30, 0.11, 0.01)

# As we only consider the top 5 most abundant species to be potential keystones,
# here, e can see that there are no KS, as the only species which median CIi 
# value falls outside of the significance range is juncs_effus, which has a low 
# abundance.

# Now let's try with the 95% confidence interval
plant.alder.95 <- CIi.fun(site.plant.alder$plant_rich, DCi.alder, "Alder swamps",
                       0.05, -30, 30, 0.11, 0.01)

# here, we can see that onocl_sensi is both significant AND among the top 5
# most abundant species. Therefore, it can be considered a KS in alder swamps.

# We see that with this threshold, we have a lot of potential keystone species
## 14 with positive contributions (blue)
## 7 with a negative contribution (red)

# You can adapt this code to your own dataset and adjust the alpha threshold 
# and explore potentiel KS in any species dataset!
# Happy coding!
