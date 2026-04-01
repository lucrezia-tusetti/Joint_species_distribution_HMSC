# Joint_species_distribution_HMSC
Application of Bayesian models (HMSC) for the joint distribution of North Sea macrobenthic species

## 🛠️ Keywords & Skills

**Keywords:** `Joint Species Distribution` `Bayesian Modeling` `HMSC` `Macrobenthos` `North Sea` `Species Occurrence` `Latent Factors` `Residual Associations`

**Skills Applied:**  
- `Bayesian theory and statistical modeling`  
- `Species distribution modeling (HMSC)`  
- `R programming` `dplyr` `tidyr`  
- `Spatial analysis` `sf` `raster` `ggplot2`  
- `Visualization` `ggplot2` `ggrepel`  
- `Model evaluation` `coda` `Hmsc`  
- `Environmental data handling` `EMODnet datasets`  
- `Functional traits` `Phylogenetic analysis`  

**Main Libraries Used in R:**  
- `Hmsc`  
- `sf`  
- `dplyr` `tidyr`  
- `ggplot2` `ggrepel`  
- `coda`  
- `ape`  
- `reshape2`  
- `sp`  
- `rnaturalearth`  
- `units`

## 🌊 Introduction

This project builds on the methodological framework presented in Ovaskainen, O. & Abrego, N. (2020). *Joint Species Distribution Modelling: With Applications in R*. Cambridge University Press., applying advanced hierarchical Bayesian models (HMSC) to real-world ecological data.

The analysis is based on data from the GEANS project, which focuses on assessing the quality of seafloor habitats—an essential indicator of marine ecosystem health. GEANS promotes the use of fast, accurate, and cost-effective DNA-based approaches, enabling more efficient and standardized environmental assessments across the North Sea.

This work was developed as part of my Master's thesis and focuses on the joint distribution of macrobenthic species. The aim is to investigate how environmental conditions, species traits, and spatial processes interact to shape ecological communities, using a fully Bayesian modeling framework.

## 📊 Data

After data cleaning and preprocessing, the dataset includes **649 distinct sampling sites** and a total of **570 species**.

Each site is considered as a sampling unit (i.e. a visited location), and species information is used to construct presence/absence data. For the purpose of this analysis, the focus is restricted to four target species:

- *Mulinia lateralis*  
- *Pagurus bernhardus*  
- *Chamelea striatula*  
- *Bathyporeia elegans*  

For each site, the response variable is defined as **binary occurrence**:
- **1** = species observed  
- **0** = species not observed  

This allows the application of occurrence-based models to investigate how environmental conditions and spatial processes influence species distribution patterns.

### 🌍 Environmental variables

The analysis includes the following environmental predictors:

- **Oil and gas extraction platforms**  
  Data on installations related to oil and gas exploration and exploitation were obtained from EMODnet:  
  [EMODnet – Oil and Gas Platforms](https://ows.emodnet-humanactivities.eu/geonetwork/srv/api/records/ddbe3597-4e3f-4e74-8d31-947c4efef2e9)

- **Urban wastewater discharge points**  
  Data on wastewater discharge locations were obtained from EMODnet:  
  [EMODnet – Wastewater Discharge Points](https://emodnet.ec.europa.eu/geonetwork/srv/api/records/0bd23b5e-b288-4273-b8b7-d073538ada52)

- **Distance from the coast**  
  For each sampling site, the distance from the coastline was computed and used as a continuous environmental predictor.

To reduce spatial noise and account for local aggregation, a coarse spatial grid was applied, and the number of platforms and discharge points was counted within each grid cell.

## 🧠 Models

The analysis is based on the **HMSC (Hierarchical Modelling of Species Communities)** framework, a Bayesian approach for joint species distribution modeling. HMSC allows the simultaneous analysis of multiple species, accounting for environmental responses, species traits, phylogenetic relationships, and residual species associations.

In this project, three models of increasing complexity are applied:

- **Model 1 – Environmental model**  
  Includes only environmental predictors (platforms, discharge points, distance from the coast) to explain species occurrence.

- **Model 2 – Traits and phylogeny model**  
  Extends the environmental model by incorporating functional traits and phylogenetic relationships among species.

- **Model 3 – Latent factor model**  
  Introduces spatial random effects and latent factors to capture residual variation and species associations not explained by environmental variables.

All models are fitted using a Bayesian framework with MCMC sampling, and model performance is evaluated through convergence diagnostics and predictive metrics.

## 📈 Main Results

- **Distance from the coast** emerged as the most important environmental predictor, acting as a key ecological gradient for the studied species.  
- Species show **negative responses to wastewater discharge points and extraction platforms**, suggesting that human activities impact their distribution and that these species are sensitive to pollution.  
- **Residual associations** captured by the latent factor model reveal patterns of co-occurrence not explained by measured environmental variables.  
- Overall, the **latent factor model (Model 3)** provided the best predictive performance, achieving the highest AUC values among the models tested.
