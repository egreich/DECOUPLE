


library(tidyverse)
library(raster)
library(terra)
library(cowplot)

# Create necessary folders if they do not already exist
if(!file.exists("plots")) { dir.create("plots")}

path_out = "./plots/" # set save path


# Load site data
df_sites <- read.csv("./data_input/sitedata.csv")

#df_sites <- df_sites[-c(32:47),] # get rid of extra rows


# Extract WorldClim normals 1970-2020 for all sites
# WorldClim 2.1 data at 30s (~1 km) spatial resolution
# Means for 1970-2000 
# https://www.worldclim.org/data/worldclim21.html#google_vignette
# Citation: Fick, S.E. and R.J. Hijmans, 2017. 
# WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. 
# International Journal of Climatology 37 (12): 4302-4315.

r <- getData("worldclim",var="bio",res=10)

r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")

lats <- df_sites$lat
lons <- df_sites$lon

coords <- data.frame(x=lons,y=lats)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r,points)

df_world <- cbind.data.frame(coordinates(points),values)

# Don't forget that WorldClim data has a scale factor of 10, so Temp = -37 is -3.7 ºC.
df_world$Temp <- df_world$Temp/10

df_sites <- df_sites %>%
  mutate(MAP = ifelse(is.na(MAP_mm), df_world$Prec, MAP_mm), 
         MAT = ifelse(is.na(MAT_C), df_world$Temp, MAT_C))


plot(r[[1]])
plot(points,add=T)


# mapping
#sm=rast("./data_input/aridity_36km.tif") # aridity
sm=rast("./data_input/SMAP_SM.tif") # SMAP soil moisture data
print(sm)
# Make custom color palette
library(unikn)
mypal2 = rev(cetcolor::cet_pal(20, name = "r2", alpha = 0.5)) 
unikn::seecol(mypal2)


# plot using ggplot and tidyterra

library(tidyterra)
p <- ggplot() +
  geom_spatraster(data = sm) +
  geom_point(data = df_sites, aes(x=lon, y=lat), color = "black", size=1) +
  scale_fill_gradientn(colors=mypal2,                               # Use user-defined colormap
                       name = "Soil Moisture",                                 # Colorbar name
                       na.value = "transparent",                    # transparent NA cells
                       labels=(c("0", "0.2", "0.4", "0.6", "0.8")), # Colorbar labels
                       breaks=seq(0,0.8,by=0.2),                    # Set breaks of colorbar
                       limits=c(0,0.8))+                            # Z-axis limits
  theme_void()                                 # Try different themes: theme_bw(), theme_gray(), theme_minimal
  
  p
  ggsave2("testmap.png", plot = p, path = path_out)

#~~~Projection 2: Robinson projection

#NE_places - SpatialPointsDataFrame with city and town points

WorldSHP=terra::vect(spData::world)

#--- convert to a SpatRaster ---#
points2 <- vect(points)

RobinsonPlot <- ggplot() +
  geom_spatraster(data = sm)+                   # Plot SpatRaster layer               
  geom_spatvector(data = WorldSHP, 
                  fill = "transparent") +       # Add world political map
  geom_spatvector(data = points2, # add points
                  fill   = "transparent",
                  colour = "black",
                  #stroke = 1,
                  alpha  = 0.5) +
  ggtitle("Site Locations") +              # Add title
  scale_fill_gradientn(colors=mypal2,           # Use user-defined colormap
                       name = "Soil Moisture",  # Name of the colorbar
                       na.value = "transparent",# Set color for NA values
                       lim=c(0,0.8))+           # Z axis limit
  theme_minimal()+                              # Select theme. Try 'theme_void'
  theme(plot.title = element_text(hjust =0.5),  # Place title in the middle of the plot
        text = element_text(size = 12))+        # Adjust plot text size for visibility
  coord_sf(crs = "ESRI:54030",                  # Reproject to World Robinson
           xlim = c(-152,152)*100000,    
           ylim = c(-55,90)*100000)

print(RobinsonPlot)
ggsave2("testmap2.png", plot = RobinsonPlot, path = path_out)


#### plot whittiker biome plot
#devtools::install_github("valentinitnelav/plotbiomes")
library(plotbiomes)


df_sites$MAP <- df_sites$MAP/10 # convert to cm

whittaker_base_plot()


plot_1 <- ggplot() +
  # add biome polygons
  geom_polygon(data = Whittaker_biomes,
               aes(x    = temp_c,
                   y    = precp_cm,
                   fill = biome),
               # adjust polygon borders
               colour = "gray98",
               size   = 1) +
  theme_bw()
plot_1
plot_2 <- plot_1 +
  # fill the polygons with predefined colors
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors)
plot_2

p <- whittaker_base_plot() +
  # add the temperature - precipitation data points
  geom_point(data = df_sites, 
             aes(x = MAT, 
                 y = MAP), 
             size   = 3,
             shape  = 21,
             colour = "gray95", 
             fill   = "black",
             stroke = 1,
             alpha  = 0.5) +
  theme_bw()
p

ggsave2("whittaker.png", plot = p, path = path_out)
