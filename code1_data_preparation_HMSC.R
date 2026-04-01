###########################
# Data preparation
###########################

# This script performs the following tasks:
# 1. Loads GEANS species observation data, platforms, and discharge points.
# 2. Creates spatial objects and grids for environmental predictors.
# 3. Prepares matrices for Hmsc modeling: response (Y), environmental predictors (X), 
#    geographic coordinates (S), functional traits (Tr), and taxonomic tree (C).

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggspatial)
library(grid)
library(ape)
library(units)

# GEANS data -------------------------------------------------------------------

geans <- read_tsv("DS-GEANS1.tsv")

# # Observations data: contains species and coordinates

geans <- geans %>%
  distinct(sampleid, .keep_all = TRUE)
geans <- geans %>% filter(!is.na(coord), !is.na(species),
         collection_date_start >= as.Date("2000-01-01"))
geans <- geans %>%
  select(coord, phylum, class, order, family, species)

# Map of observations in the North Sea 

geans <- geans %>%
  mutate(
    coord_clean = str_remove_all(coord, "\\[|\\]"),
    lat = as.numeric(str_trim(str_split_fixed(coord_clean, ",", 2)[,1])),
    lon = as.numeric(str_trim(str_split_fixed(coord_clean, ",", 2)[,2]))
  ) %>%
  select(-coord_clean)

geans_sf <- st_as_sf(geans, coords = c("lon", "lat"), crs = 4326)

world <- ne_countries(scale = "medium", returnclass = "sf")

geans <- geans %>%
  filter(lat >= 50, lat <= 60,
         lon >= -3, lon <= 13)
geans <- geans %>%
  filter(
    !(
      (lon >= -2.44 & lon <= -2.36 & lat >= 55.59 & lat <= 55.63) |  
        (lon >= 6     & lon <= 11     & lat >= 50    & lat <= 53)      
    )
  )
geans_sf <- st_as_sf(geans, coords = c("lon", "lat"), crs = 4326)

g1 <- ggplot() +
  geom_sf(data = world, fill = "grey98", color = "black") +
  geom_sf(data = geans_sf, color = "springgreen3", size = 1, alpha = 0.4) +
  geom_point(aes(x=Inf, y=Inf, color="Area d'interesse")) +
  geom_point(aes(x=Inf, y=Inf, color="Organismi osservati")) +
  scale_color_manual(name = "", values = c("Organismi osservati" = "springgreen3")) +
  guides(color = guide_legend(ncol = 2, override.aes = list(shape = c(19), 
                                                            size = c(3)))) +
  coord_sf(xlim = c(-5.5, 14), ylim = c(48, 64)) +
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_orienteering(fill=c("black", "black"))) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 14))







# Environmental predictors -----------------------------------------------------

# Platforms and discharge: create a coarse grid and count number of platforms
# This approach applies spatial smoothing to the data

# Pred 1) Extraction platforms

platforms <- read.csv("platforms.csv")
platforms$the_geom <- str_replace_all(platforms$the_geom,
  "(\\s*)([-0-9\\.]+) ([-0-9\\.]+)",
  "\\1\\3 \\2")                         # Latitude and longitude were swapped
platforms_sf <- st_as_sf(platforms, wkt = "the_geom", crs = 4326) 
platforms_sf <- platforms_sf %>% select(the_geom)

# Pred 2) Wastewater discharge points

discharge <- read.csv("dischargepoints.csv")
discharge_sf <- st_as_sf(discharge, coords = c("longitude", "latitude"), crs = 4326)
discharge_sf <- discharge_sf %>% select(geometry)

# Bounding box 

bbox_sf <- st_as_sfc(
  st_bbox(c(xmin = -3, xmax = 13, ymin = 49.9, ymax = 60), crs = 4326))
bbox_utm <- st_transform(bbox_sf, crs = 32632)
bin_size_m <- 300000 # Cell size 

# Grid from UTM(meters) -> long/lat

grid_utm <- st_make_grid(
  bbox_utm,
  cellsize = bin_size_m,
  what = "polygons"
)
grid_utm <- st_sf(grid_id = seq_along(grid_utm), geometry = grid_utm) # UTM
grid_sf <- st_transform(grid_utm, 4326)                               # long/lat
grid_sf <- grid_sf[lengths(st_intersects(grid_sf, bbox_sf)) > 0, ]    # Remove extra
grid_sf <- grid_sf[ !st_within(grid_sf, st_union(world), sparse = FALSE)[,1], ]
# grid_utm <- st_transform(grid_sf, crs = 32632)

# For each point (platform or discharge), determine which grid cell it falls into
# This adds a 'grid_id' column to each dataset for later aggregation

platforms_sf <- st_join(platforms_sf, grid_sf %>% select(grid_id), join = st_within)
discharge_sf <- st_join(discharge_sf, grid_sf %>% select(grid_id), join = st_within)

# Map of all points: platforms and discharge points in the North Sea

g2 <- ggplot() +
  geom_sf(data = world, fill = "grey98", color = "black") +
  geom_sf(data = platforms_sf, col="brown3", size = 0.7) +
  geom_sf(data = discharge_sf, col="steelblue", size=0.7) +
  geom_point(aes(x=Inf, y=Inf, color="Installazioni per attività di sfruttamento ed esplorazione")) +
  geom_point(aes(x=Inf, y=Inf, color="Punti di scarico delle acque reflue urbane")) +
  scale_color_manual(name = "", values = c("Installazioni per attività di sfruttamento ed esplorazione" = "brown3",
                                           "Punti di scarico delle acque reflue urbane" = "steelblue")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(shape = c(19, 19), 
                                                            size = c(3,3)))) +
  coord_sf(xlim = c(-6, 12), ylim = c(48, 64)) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_orienteering(fill = c("black", "black"))) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 14))

# Count platforms and discharge points in each grid cell

grid_sf <- grid_sf %>%
  left_join(
    platforms_sf %>%
      st_drop_geometry() %>%          
      group_by(grid_id) %>%
      summarise(n_platforms = n(), .groups = "drop"),
    by = "grid_id"
  ) %>%
  left_join(
    discharge_sf %>%
      st_drop_geometry() %>%          
      group_by(grid_id) %>%
      summarise(n_discharge = n(), .groups = "drop"),
    by = "grid_id"
  ) %>%
  mutate(
    n_platforms = replace_na(n_platforms, 0),
    n_discharge = replace_na(n_discharge, 0)
  )

# Point and grid maps

ggplot() +
  geom_sf(data = world, fill = "grey98", color = "black") +
  geom_sf(data = grid_sf, aes(fill = n_discharge), color = NA, alpha = 0.8) +
  scale_fill_gradient(name = "n_discharge", low = "white", high = "darkblue", na.value = NA) +
  geom_sf(data = discharge_sf, aes(color = "Punti di scarico delle acque reflue urbane"), size = 0.7) +
  scale_color_manual(name = NULL, values = c("Punti di scarico delle acque reflue urbane" = "steelblue")) +
  guides(
    fill = guide_colorbar(title = "n_discharge", title.size = 10, barwidth = 15, barheight = 0.8),
    color = guide_legend(ncol = 1, override.aes = list(size = 4), title.theme = element_text(size = 30), order=1)
  ) +
  coord_sf(xlim = c(-6, 12), ylim = c(48, 64)) +
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"), style = north_arrow_orienteering(fill = c("black", "black"))) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_minimal() +
  theme(axis.title = element_blank(), legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 10))


ggplot() +
  geom_sf(data = world, fill = "grey98", color = "black") +
  geom_sf(data = grid_sf, aes(fill = n_platforms), color = NA, alpha = 0.8) +
  scale_fill_gradient(name = "n_platforms", low = "white", high = "brown", na.value = NA) +
  geom_sf(data = platforms_sf, aes(color = "Installazioni per attività di sfruttamento ed esplorazione di petrolio e gas"), size = 0.7) +
  scale_color_manual(name = NULL, values = c("Installazioni per attività di sfruttamento ed esplorazione di petrolio e gas" = "brown3")) +
  guides(
    fill = guide_colorbar(title = "n_platforms", title.size = 10, barwidth = 15, barheight = 0.8),
    color = guide_legend(ncol = 1, override.aes = list(size = 4), title.theme = element_text(size = 30), order=1)
  ) +
  coord_sf(xlim = c(-6, 12), ylim = c(48, 64)) +
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"), style = north_arrow_orienteering(fill = c("black", "black"))) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme_minimal() +
  theme(axis.title = element_blank(), legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 10))



# Matrix creation for Hmsc functions -------------------------------------------

### Y Response variables ###

Y <- geans_sf %>%
  st_drop_geometry() %>%            
  mutate(geometry = geans_sf$geometry) %>%  
  count(geometry, species) %>%       
  pivot_wider(
    names_from = species,
    values_from = n,
    values_fill = 0
  ) %>%
  st_as_sf() %>%
  mutate(ID = row_number())

Y <- Y %>%
  select(ID, geometry,
    `Mulinia lateralis`, `Pagurus bernhardus`,
    `Chamelea striatula`,`Bathyporeia elegans`)

### X Environmental variables ###

X <- Y %>%
  select(ID, geometry)

X <- st_join(
  X,           
  grid_sf[, c("n_platforms", "n_discharge")],  
  join = st_intersects  
)

# Add a third environmental variable: distance from the coast

costa <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(3857)  
X <- st_transform(X, 3857)
X <- X %>%
  mutate(dist_coast = st_distance(geometry, st_union(costa)) %>% as.numeric())

X <- st_drop_geometry(X)

### S Geographic coordinates ###

S <- Y %>% select(ID, geometry) %>%
  st_centroid() %>%                     
  mutate(
    coords = st_coordinates(.)          
  ) %>%
  st_drop_geometry() %>%                
  mutate(
    lon = coords[,1],                   
    lat = coords[,2]                    
  ) %>%
  select(ID, lon, lat)

### Tr Functional traits ###

Tratti <- read.csv("Tratti.csv", sep=";")

# life_habit: infaunal, epifaunal/epibenthic
# skeleton:   exoskeleton, endoskeleton, shell
# min_depth
# max_depth

### C Taxonomic tree ###

# Representation of phylogenetic correlation: simplification using taxonomy

tassonomia <- geans_sf %>% 
  filter(species %in% c("Mulinia lateralis", "Pagurus bernhardus",
                        "Chamelea striatula", "Bathyporeia elegans")) %>%
  st_set_geometry(NULL) %>%  
  select(species, family, order, class, phylum) %>%
  distinct()

livelli <- c("phylum", "class", "order", "family", "species")
tassonomia[livelli] <- lapply(tassonomia[livelli], as.factor)

tree <- as.phylo(~ phylum/class/order/family/species, data=tassonomia, collapse = F)
tree$edge.length = rep(1, length(tree$edge))

plot(tree, cex=1.25)





