rm(list=ls())
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(janitor)
library(viridis)
library(scales)
library(RColorBrewer)
library(stringr)
library(spdep)
library(tidyverse)

df <- fread("/distn_duffy.csv") %>%
      clean_names() %>% mutate(region = str_trim(region))

world <- ne_countries(scale = "medium", returnclass = "sf")
duffy_df <- world %>% left_join(df, by = c("name_long" = "region"))

########################
nb <- poly2nb(duffy_df, queen = TRUE)

duffy_df <- duffy_df %>% mutate(
      neighbor_null_mean = sapply(1:n(), function(i) {
      neighbors <- nb[[i]]
      if(length(neighbors) == 1 && neighbors == 0) return(NA)
            mean(null[neighbors], na.rm = TRUE)
    })
  )

duffy_df <- duffy_df %>% 
  mutate(null2 = case_when(
    subregion == "Antarctica" ~ NA,
    !is.na(null) ~ null,
    is.na(null) & !is.na(neighbor_null_mean) ~ neighbor_null_mean,
    is.na(null) ~ 0
  ))

# Add bold line for EUR
pal_fun <- colorRampPalette(brewer.pal(9, "Reds"))

p_map <- ggplot() + geom_sf(data = duffy_df, aes(fill = null2), color = "gray50", linewidth = 0.1) +
                  scale_fill_gradientn(colours = pal_fun(100), 
                                      name = "Mean frequency", na.value = "gray90") +
                  #line for EUR
                  geom_sf(data = subset(duffy_df, continent %in% c("North America", "Europe") &
                                        name_long != "Russian Federation" & iso_a3 != "NIU" & iso_a3 != "GRL"),
                  fill = NA, color = "black", linewidth = 0.5) +
                  theme_minimal() +
                  theme(legend.position = "none", panel.grid = element_blank(),
                        axis.title = element_blank(),  axis.text  = element_blank(),
                        axis.ticks = element_blank())

# Global map highlighting different groups
global <- ne_countries(scale = "medium", returnclass = "sf")
global <- global %>%
  mutate(
    ancestry = case_when(
      continent == "Europe"  ~ "European",
      continent == "South America"  ~ "European2",
      continent == "North America"  ~ "European",
      continent == "Africa"  ~ "African",
      continent == "Asia" & subregion == "Southern Asia" ~ "South Asian",
      continent == "Asia" & subregion == "Western Asia" ~ "South Asian",
      continent == "Asia" & subregion == "Eastern Asia"  ~ "East Asian",
      continent == "Asia" & subregion == "South-Eastern Asia"  ~ "East Asian",
      continent == "Asia" & subregion == "Central Asia"  ~ "East Asian",
      continent == "Oceania"  ~ "European",
      TRUE ~ NA_character_
    )
  )

global$ancestry <- factor(global$ancestry, levels=c("East Asian", "European2", 
                                                    "European", "African", "South Asian"))

library(RColorBrewer)
cols <- c("#66C2A5","#EFC000FF", "#003C67FF", "#CD534CFF","#A6D854")

p_map2 <- ggplot() + geom_sf(data = subset(global, is.na(ancestry)), 
                   aes(fill = ancestry), linewidth = 0.1, color = "gray50") + # line for Antarctica
          geom_sf(data = subset(global, !is.na(ancestry)), # line for others
                  aes(fill = ancestry), color = "black", linewidth = 0.5, alpha = 0.5) + 
          scale_fill_manual(values = cols, na.value = "gray90")+
          #scale_fill_brewer(palette = cols, na.value = "gray90") +
          theme_minimal() +
          theme(legend.position = "none",
                panel.grid = element_blank(),
                axis.title = element_blank(),
                axis.text  = element_blank(),
                axis.ticks = element_blank()) 
 
##############################################
# Global map colored with EUR
global2 <- global %>% mutate(
           ancestry2 = case_when(
                       ancestry == "European" ~ "EUR",
                       TRUE ~ NA_character_))

ggplot(global2) +
  geom_sf(aes(fill = ancestry2), color = "gray50", linewidth = 0.1, alpha=0.5) +  #gray50
  scale_fill_brewer(
    palette = "Set2",
    na.value = "gray90"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) 

##############################################
# Make border thicker
ggplot() + geom_sf(
          data = global,
          aes(fill = ancestry),
          color = "gray50",
          linewidth = 0.1,
          alpha = 0.5
        ) +
        geom_sf(data = subset(global,
            continent %in% c("North America", "Europe") &
              name_long != "Russian Federation" & 
              iso_a3 != "NIU" &
              iso_a3 != "GRL"),
          fill = NA,
          color = "black",
          linewidth = 0.4
        ) +
        
        scale_fill_brewer(
          palette = "Set2",
          na.value = "gray90"
        ) +
        theme_minimal() +
        theme(
          legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank()
        )


##################################################################

df <- fread("/Users/ums3/Documents/PRS_CR/data/geo.csv")
select <- dplyr::select

df_f <- df %>% filter(if_any(starts_with("phe"), ~ !is.na(.))) %>% 
  select(lat, lon, starts_with("phe"))

df_f <- df_f %>%
  select(lat, lon, starts_with("phe")) %>%
  mutate(across(starts_with("phe"), ~ ntile(., 10)))

df_long <- df_f %>%
  pivot_longer(
    cols = starts_with("phe"),
    names_to = "phenotype",
    values_to = "value"
  )

# Define custom palettes per phenotype
palettes <- list(
  pheab = colorRampPalette(brewer.pal(9, "Purples")),
  phea  = colorRampPalette(brewer.pal(9, "Blues")),
  pheb  = colorRampPalette(brewer.pal(9, "Greens")),
  phe0  = colorRampPalette(brewer.pal(9, "Reds"))
)

plots <- list()
ind <- unique(df_long$phenotype)

# World polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

for(k in ind){
  temp <- df_long %>% filter(phenotype == k)
  
  # Convert to sf points
  temp_sf <- st_as_sf(temp, coords = c("lon","lat"), crs = 4326)
  
  # Spatial join: assign each point to a country
  temp_with_country <- st_join(temp_sf, world["iso_a3"])
  
  # Aggregate: mean per country
  allele_freq_df <- temp_with_country %>%
    st_drop_geometry() %>%
    group_by(iso_a3) %>%
    summarise(freq = mean(value, na.rm = TRUE), .groups = "drop")
  
  # Join back to polygons
  world_data <- world %>%
    left_join(allele_freq_df, by = "iso_a3")
  
  # Pick palette for this phenotype
  pal_fun <- palettes[[k]]
  
  # Plot
  plot <- ggplot(world_data) +
    geom_sf(aes(fill = freq)) +
    scale_fill_gradientn(
      colours = pal_fun(100), 
      name = "Mean frequency", 
      na.value = "gray90"
    ) +
    theme_minimal() +
    labs(title = k) +
    theme(legend.position = "none")
  
  plots[[k]] <- plot
}

# Arrange plots into one figure
plot <- grid.arrange(plots$pheab, plots$phea, plots$pheb, 
             plots$phe0, ncol = 2)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(gridExtra)

df <- fread("~/frequencies_by_country.csv")

df <- df %>%
  mutate(
    Country = case_when(
      Country == "United States (Utah)" ~ "United States of America",
      Country == "Congo (DRC)" ~ "Dem. Rep. Congo",
      grepl("China", Country) ~ "China",
      TRUE ~ Country
    )
  )

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::select(name, iso_a3, geometry)

coords <- st_centroid(world) %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(name, lon, lat)

df_geo <- df %>%
  left_join(coords, by = c("Country" = "name"))

df_geo <- df_geo %>%
  rename(
    pheab = `Fy(a+b−)`,   # Fy(a+b−)
    pheb  = `Fy(a−b+)`,   # Fy(a−b+)
    phea  = `Fy(a+b+)`,   # Fy(a+b+)
    phe0  = `Fy(a−b−)`    # Fy(a−b−)
  )

df_f <- df_geo %>%
  filter(if_any(starts_with("phe"), ~ !is.na(.))) %>%
  select(lat, lon, starts_with("phe"))

df_long <- df_f %>%
  pivot_longer(
    cols = starts_with("phe"),
    names_to = "phenotype",
    values_to = "value"
  )

palettes <- list(
  pheab = colorRampPalette(brewer.pal(9, "Purples")),
  phea  = colorRampPalette(brewer.pal(9, "Blues")),
  pheb  = colorRampPalette(brewer.pal(9, "Greens")),
  phe0  = colorRampPalette(brewer.pal(9, "Reds"))
)

world <- ne_countries(scale = "medium", returnclass = "sf")

plots <- list()
phenos <- unique(df_long$phenotype)

for(k in phenos){
  temp <- df_long %>% filter(phenotype == k)
  
  # Convert to sf points
  temp_sf <- st_as_sf(temp, coords = c("lon", "lat"), crs = 4326)
  
  # Spatial join: assign to countries
  temp_with_country <- st_join(temp_sf, world["iso_a3"])
  
  # Aggregate: mean per country
  allele_freq_df <- temp_with_country %>%
    st_drop_geometry() %>%
    group_by(iso_a3) %>%
    summarise(freq = mean(value, na.rm = TRUE), .groups = "drop")
  
  # Join back to polygons
  world_data <- world %>%
    left_join(allele_freq_df, by = "iso_a3")
  
  # Pick palette
  pal_fun <- palettes[[k]]
  
  # Make map
  p <- ggplot(world_data) +
    geom_sf(aes(fill = freq)) +
    scale_fill_gradientn(
      colours = pal_fun(100),
      name = "Frequency (%)",
      limits = c(0, 100),
      na.value = "gray90"
    ) +
    theme_minimal() +
    labs(
      title = case_when(
        k == "pheab" ~ "Fy(a+b−)",
        k == "phea"  ~ "Fy(a+b+)",
        k == "pheb"  ~ "Fy(a−b+)",
        k == "phe0"  ~ "Fy(a−b−)"
      ),
      fill = "Frequency (%)"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  plots[[k]] <- p
}

grid.arrange(
  plots$pheab, plots$phea,
  plots$pheb, plots$phe0,
  ncol = 2
)

