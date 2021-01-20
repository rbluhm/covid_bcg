# required
lop <- c("sf","haven", "tidyverse", "lubridate", "RColorBrewer", "units")
newp <- lop[!(lop %in% installed.packages()[,"Package"])]
if(length(newp)) install.packages(newp)
lapply(lop, require, character.only = TRUE)

# get od flows
input <- read_dta("./data/commuter_panel_filled_dec2019.dta")
input <- input %>% select(origin= id_kreis_a, destination = id_kreis_b, total=outgoing)
input

# read shape
kreise.df <- st_read('./data/RKI_Corona_Landkreise.gpkg', stringsAsFactors=F) %>% st_transform(4326)
kreise.df <- kreise.df %>% select(id=RS, county=GEN, type = BEZ, nuts_id = NUTS, state = BL, pop_latest = EWZ, area = KFL)
names(kreise.df)

## now aggregate over berlin
kreise.df <- kreise.df %>% mutate(id = replace(id, state=="Berlin", 11000))
kreise.df <- kreise.df %>% group_by(id,state) %>% summarize(county = first(county), pop = sum(pop_latest))

# define ddr
kreise.df <- kreise.df %>% mutate(ddr = case_when(
  state == "Mecklenburg-Vorpommern" ~ 1,
  state == "Brandenburg" ~ 1,
  state == "Berlin" ~ 1,
  state == "Sachsen-Anhalt" ~ 1,
  state == "Sachsen" ~ 1,
  state == "Thüringen" ~ 1,
  TRUE ~ 0
))

xy <- kreise.df %>%  st_centroid() %>% st_coordinates()
kreise.df$lon <- xy[,1]
kreise.df$lat <- xy[,2]
centroids.df <- kreise.df %>%  st_set_geometry(NULL) %>% select(id, ddr, county,  lat, lon, pop)

#Lots of joining to get the xy coordinates joined to the origin and then the destination points.
or.xy <- merge(input, centroids.df, by.x="origin", by.y="id")
names(or.xy)<- c("origin", "destination", "trips", "ddrX", "o_name", "oX", "oY", "popX")

dest.xy <-  merge(or.xy, centroids.df, by.x="destination", by.y="id")
names(dest.xy)<- c("origin", "destination", "trips", "ddrX", "o_name", "oX", "oY", 
                   "popX","ddrY", "d_name", "dX", "dY", "popY")


x1 <- dest.xy %>% select(oY, oX, trips, origin) %>%  st_as_sf(coords= c("oY","oX"))
x2 <- dest.xy %>% select(dY, dX, trips, origin) %>%  st_as_sf(coords= c("dY","dX"))
x1x2 <- st_nearest_points(x1,x2, pairwise = T)
st_geometry(dest.xy) <- x1x2
dest.xy <- dest.xy %>% st_set_crs(4326) %>% mutate(length_km = st_length(.), length_km = set_units(length_km, "km")) 

ddr_border <- st_read("./data/ddr_border.gpkg", stringsAsFactors=F)

png("./output/Supplement_Figure_S1_a_west.png", width = 5, height = 6, unit="in", res=300)
plot(kreise.df[1], lwd=0.5, col="grey95", reset=F, main="")
plot(ddr_border, col="red", lwd=2, add=T)
plot(dest.xy %>% filter(trips>1000, ddrY==0) %>% select(trips) %>%  st_geometry(), col="blue", add=T)
dev.off()

dest.xy <- dest.xy %>% filter(length_km > set_units(50, "km"))

png("./output/Supplement_Figure_S1_b_east_50k.png", width = 5, height = 6, unit="in", res=300)
plot(kreise.df[1], lwd=0.5, col="grey95", reset=F, main="")
plot(ddr_border, col="red", lwd=2, add=T)
plot(dest.xy %>% filter(trips>1000, ddrY==1) %>% select(trips) %>%  st_geometry(), col="blue", add=T)
dev.off()


png("./output/Supplement_Figure_S1_a_west_50k.png", width = 5, height = 6, unit="in", res=300)
plot(kreise.df[1], lwd=0.5, col="grey95", reset=F, main="")
plot(ddr_border, col="red", lwd=2, add=T)
plot(dest.xy %>% filter(trips>1000, ddrY==0) %>% select(trips) %>%  st_geometry(), col="blue", add=T)
dev.off()

png("./output/Supplement_Figure_S1_b_east.png", width = 5, height = 6, unit="in", res=300)
plot(kreise.df[1], lwd=0.5, col="grey95", reset=F, main="")
plot(ddr_border, col="red", lwd=2, add=T)
plot(dest.xy %>% filter(trips>1000, ddrY==1) %>% select(trips) %>%  st_geometry(), col="blue", add=T)
dev.off()
