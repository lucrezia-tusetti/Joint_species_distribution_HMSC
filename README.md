# Joint_species_distribution_HMSC
Application of Bayesian models (HMSC) for the joint distribution of North Sea macrobenthic species

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
