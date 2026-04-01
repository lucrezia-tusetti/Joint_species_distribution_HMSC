###########################
# Hmsc modeling
###########################

# Apply occurrence models to study species presence probabilities.
# Response variables (Y) are converted to binary:
#   0 = species not observed
#   1 = species observed

# Model overview
# Model 1: Environmental predictors only (fixed effects)
# Model 2: Environmental predictors + functional traits + phylogenetic correlation
# Model 3: Latent factor model, it captures co-occurrence patterns not explained
#          by environmental predictors or traits

library(Hmsc)       
library(dplyr)      
library(tidyr)      
library(ggplot2)    
library(ggrepel)    
library(coda)       
library(reshape2)   
library(pheatmap)
library(sf)         
library(scales) 

# Model 1: Environmental predictors only ---------------------------------------

Y <- as.matrix(Y %>% st_set_geometry(NULL) %>%      
                 select("Mulinia lateralis", "Pagurus bernhardus", "Chamelea striatula", "Bathyporeia elegans"))
rownames(Y) <- Y1$ID
Y[Y>= 1] <- 1

XData <- X %>% select(n_platforms, n_discharge, dist_coast) %>% 
  as.data.frame()
rownames(XData) <- X$ID
XFormula = ~ n_platforms + n_discharge + dist_coast

m_PRESENZA_ENV = Hmsc(Y = Y, 
                      XData = XData, 
                      XFormula = XFormula,              
                      distr = "probit")
nChains = 2
thin = 100
samples = 1000
transient = 500*thin
verbose = 1*thin

m_PRESENZA_ENV <- sampleMcmc(
  m_PRESENZA_ENV, 
  thin = thin,
  samples = samples, 
  transient = transient,
  nChains = nChains, 
  verbose = verbose
)
computeWAIC(m_PRESENZA_ENV)

# saveRDS(m_PRESENZA_ENV, file = "Presenza_ambiente.rds")
# m_PRESENZA_ENV <- readRDS("Presenza_ambiente.rds")

# Convert Hmsc model output into a coda-compatible object
# This allows MCMC diagnostics and posterior summaries using coda functions

mpost = convertToCodaObject(
  m_PRESENZA_ENV,
  spNamesNumbers = c(T,F),
  covNamesNumbers = c(T,F)
)

effectiveSize(mpost$Beta)                          # Number of effectively independent samples 
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf # Values close to 1 indicate good convergence

plot(mpost$Beta) # Traceplots

summary(mpost$Beta) # Obtained estimates

predY = computePredictedValues(m_PRESENZA_ENV)
MF = evaluateModelFit(m_PRESENZA_ENV, predY) 
MF # model fit metrics

# Map showing how species occurrence probability varies with distance from the coast
# Visualizes predicted probabilities from the Hmsc model for each species

Grad_disc = constructGradient(m_PRESENZA_ENV, focalVariable = "dist_coast", 
                              non.focalVariables=list(n_discharge=list(1), n_platforms = list(1)), 
                              ngrid = 200)
pred_disc = predict(m_PRESENZA_ENV, Gradient = Grad_disc, expected = TRUE)
plotGradient(m_PRESENZA_ENV, Grad_disc, pred = pred_disc, measure = "Y", index = 1,  xaxt = "n", 
             showData = F, main = " ", xlab="Distanza dalla costa", cex.lab = 0.8)
plotGradient(m_PRESENZA_ENV, Grad_disc, pred = pred_disc, measure = "Y", index = 2,  xaxt = "n", 
             showData = F, main = " ", xlab="Distanza dalla costa", cex.lab = 0.8)
plotGradient(m_PRESENZA_ENV, Grad_disc, pred = pred_disc, measure = "Y", index = 3,  xaxt = "n", 
             showData = F, main = " ", xlab="Distanza dalla costa", cex.lab = 0.8)
plotGradient(m_PRESENZA_ENV, Grad_disc, pred = pred_disc, measure = "Y", index = 4,  xaxt = "n",
             showData = F, main = " ", xlab="Distanza dalla costa", cex.lab = 0.8)
x_ticks_real <- c(0e5, 1e5, 2e5, 3e5, 4e5, 5e5)  
x_labels <- 0:5
axis(1, at = x_ticks_real, labels = x_labels)


# Model 2: Env +  traits + phylo correlation -----------------------------------

m_PA_1 = Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
              TrData = TrData, TrFormula = TrFormula, phyloTree = phyloTree,   distr = "probit")

nChains = 2
thin = 100
samples = 1000
transient = 500 * thin
verbose = 1 * thin
m_PA_1 <- sampleMcmc(m_PA_1, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

# saveRDS(m_PA_1, file = "Presenza_ambienteTratti.rds")
# m_PA_1 <- readRDS("Presenza_ambienteTratti.rds")

mpost = convertToCodaObject(
  m_PA_1,
  spNamesNumbers = c(T,F),
  covNamesNumbers = c(T,F)
)
effectiveSize(mpost$Beta) 
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

effectiveSize(mpost$Rho) 
gelman.diag(mpost$Rho, multivariate = FALSE)$psrf

summary(mpost$Beta)
summary(mpost$Rho)

predY = computePredictedValues(m_PA_1)
MF = evaluateModelFit(m_PA_1, predY)
MF

computeWAIC(m_PA_1)
# Model 3: Latent factor model -------------------------------------------------

Y <- as.matrix(Y %>% st_set_geometry(NULL) %>%      
                 select("Mulinia lateralis", "Pagurus bernhardus", "Chamelea striatula", "Bathyporeia elegans"))
rownames(Y) <- Y1$ID
Y[Y>= 1] <- 1

XData <- X %>% select(n_platforms, n_discharge, dist_coast) %>% 
  as.data.frame()
rownames(XData) <- X$ID
XFormula = ~ n_platforms + n_discharge + dist_coast

# Prepare spatial random effects for Hmsc latent factor model
# xy: matrix of geographic coordinates (longitude and latitude)
# studyDesign: data frame defining sampling units (here, each site is a factor)
# rownames(xy) set to match site IDs for consistency
# rL: spatial random level using Nearest Neighbor Gaussian Process (NNGP)
#      nNeighbours = 10 defines the number of neighbors considered for spatial correlation

xy = as.matrix(cbind(x=S$lon, y=S$lat))       
studyDesign = data.frame(sito = as.factor(X$ID))
rownames(xy) = studyDesign[,1]
rL = HmscRandomLevel(sData = xy, sMethod = "NNGP", nNeighbours = 10)

m_effCASUALI_4 = Hmsc(Y = Y, XData = XData, XFormula = XFormula, distr = "probit",
                      studyDesign = studyDesign, ranLevels = list(sito = rL))
nChains = 2
thin = 1000
samples = 1000
transient = 500*thin
verbose = thin

m_effCASUALI_4 <- sampleMcmc(m_effCASUALI_4, thin = thin, samples = samples, 
                             transient = transient, nChains = nChains, verbose = verbose)

# saveRDS(m_effCASUALI_4, file = "PROBIT_effCASUALI_4.rds")
# m_effCASUALI_4 <- readRDS("PROBIT_effCASUALI_4.rds")

computeWAIC(m_effCASUALI_4)

mpost = convertToCodaObject(
  m_effCASUALI_4,
  spNamesNumbers = c(T,F),
  covNamesNumbers = c(T,F)
)

effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

effectiveSize(mpost$Omega[[1]]) 
gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf

effectiveSize(mpost$Lambda[[1]]) 
gelman.diag(mpost$Lambda[[1]], multivariate = FALSE)$psrf

effectiveSize(mpost$Alpha[[1]]) 
gelman.diag(mpost$Alpha[[1]], multivariate = FALSE)$psrf

postLambda = getPostEstimate(m_effCASUALI_4, parName = "Lambda")
factor_importance = apply(postLambda$mean, 1, var)
factor_importance   # Variance of latent factors
                    # It is progressively shrunk towards zero thanks to the shrinking prior

summary(mpost$Beta)
summary(mpost$Omega[[1]])
summary(mpost$Alpha[[1]])
summary(mpost$Lambda[[1]])

predY = computePredictedValues(m_effCASUALI_4)
MF = evaluateModelFit(m_effCASUALI_4, predY)
MF

# Residual association matrix (standardized)
# Positive values: species co-occur more than expected based on environmental predictors
# Negative values: species co-occur less than expected based on environmental predictors

Omega1 <- matrix(c(
   0.042, -0.005, -0.002, -0.019,
  -0.005,  0.068,  0.009, -0.003,
  -0.002,  0.009,  0.041, -0.001,
  -0.019, -0.003, -0.001,  0.059
), 
nrow = 4, byrow = TRUE)
Omega2 <- cov2cor(Omega1) # standardized

species_names <- c("Mulinia lateralis", "Pagurus bernhardus",
                   "Chamelea striatula", "Bathyporeia elegans")
species_names <- gsub(" ", "\n", species_names)
rownames(Omega2) <- species_names
colnames(Omega2) <- species_names

df_cor <- melt(Omega2)
colnames(df_cor) <- c("Species1","Species2","Correlation")
df_cor$Species2 <- factor(df_cor$Species2, levels=rev(species_names))

ggplot(df_cor, aes(x=Species1, y=Species2, fill=Correlation)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(Correlation,3)), size=4) +
  scale_fill_gradient2(low="darkblue", mid="white", high="red3",
                       midpoint=0, limits=c(-max_val, max_val)) +
  theme_minimal(base_family = "serif") + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = .5, face="italic", size=12),
    axis.text.y = element_text(face="italic", size=12),
    legend.position = "bottom",
    legend.title = element_text(face="plain", size = 12, margin = margin(r = 20, b=20)),
    legend.text  = element_text(size = 12),
    legend.key.width  = unit(1, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  labs(fill="Correlazione residua", x=NULL, y=NULL, title="")


### Latent factor plot of sites and species ###

# This plot visualizes the first two latent factors from the Hmsc model:
# - Each point represents a sampling site positioned according to its latent factor values (Eta1, Eta2)
# - Color of points indicates distance from the coast
# - Triangles represent species, positioned according to their factor loadings (Lambda)
# - Species names are added with non-overlapping labels
# Interpretation:
# - Sites close to a species point have a higher predicted probability of that species occurring
# - The plot shows how sites and species relate in the latent factor space and how environmental gradients (e.g., distance from the coast) influence these patterns

postEta = getPostEstimate(m_effCASUALI_4, parName = "Eta")
Eta1 = postEta$mean[,1]
Eta2 = postEta$mean[,2]

postLambda = getPostEstimate(m_effCASUALI_4, parName = "Lambda")
species_coords = as.data.frame(t(postLambda$mean))
species_coords$Species = m_effCASUALI_4$spNames 

df_points <- data.frame(
  x = Eta1, y = Eta2,
  dist_coast = XData$dist_coast,  # color
  n_platforms = XData$n_platforms, n_discharge = XData$n_discharge
)

# Define species positions in the latent factor space for plotting
# x and y: coordinates based on the species' factor loadings (Lambda) from the Hmsc model
# Values are scaled by 2.5 to improve visual separation in the plot

df_species <- data.frame(
  x = c(-0.0820299497 ,-0.100224130 ,-0.0002582960 , 0.1412065093),
  y = c(0.0108565761 , 0.012981015 ,-0.0165607251 ,-0.0086433129),
  name = c("Mulinia lateralis", "Pagurus bernhardus", "Chamelea striatula", "Bathyporeia elegans")
)
df_species$x <- df_species$x*2.5
df_species$y <- df_species$y*2.5

ggplot() +
  geom_point(data = df_points, aes(x = x, y = y, color = dist_coast), size = 3) +
  geom_point(data = df_species, aes(x = x, y = y), shape = 17, size = 4, color = "gold1") +
  geom_text_repel(data = df_species, aes(x = x, y = y, label = name),
                  size = 6,
                  nudge_y = 0.005,      
                  nudge_x = 0.005,      
                  segment.color = "grey95") +
  theme_minimal() +
  labs(x = "Fattore latente 1", y = "Fattore latente 2", color = "Distanza dalla costa") +
  scale_color_gradientn(
    colors = c("darkblue", "grey90", "red3"),
    values = scales::rescale(c(min(df_points$dist_coast),
                               median(df_points$dist_coast),
                               max(df_points$dist_coast))),
    labels = function(x) round(x/100000,0)  # scala in migliaia
  )  +
  guides(color = guide_colorbar(direction = "horizontal", title.position = "top", barwidth = 15, barheight = 1)) +
  theme(
    legend.position = "bottom",
    legend.title.align = 0.5,
    legend.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12)
  )

# Simple version of the latent factor biplot

biPlot(m_effCASUALI_4,
       etaPost = postEta,
       lambdaPost = postLambda,
       colVar = "dist_coast") 
