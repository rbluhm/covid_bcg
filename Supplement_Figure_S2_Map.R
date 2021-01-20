# required
lop <- c("sf", "tidyverse", "haven", "lubridate", "RColorBrewer")
newp <- lop[!(lop %in% installed.packages()[,"Package"])]
if(length(newp)) install.packages(newp)
lapply(lop, require, character.only = TRUE)

# read shape
kreise.df <- st_read('./data/RKI_Corona_Landkreise.gpkg', stringsAsFactors=F)
table(kreise.df$BL)

# only use geoinformation and take data from daily scrapes
kreise.df <- kreise.df %>% select(id=RS, county=GEN, type = BEZ, nuts_id = NUTS, state = BL, pop_latest = EWZ, area = KFL)
names(kreise.df)

## now aggregate over berlin
kreise.df <- kreise.df %>% mutate(id = replace(id, state=="Berlin", 11000))
kreise.df <- kreise.df %>% group_by(id) %>% summarize()

## add simulated data
input.df <- read_dta("./intermediate/cases_simulation_donut_50.dta") 
kreise.df <- kreise.df %>%  mutate(id = as.numeric(id)) %>% left_join(input.df)


# add pop
input.df <- read_csv("./data/pop_2017_2018.csv")  %>%  mutate(id = as.numeric(id))
kreise.df <- kreise.df %>% left_join(input.df)

# add border
eastwest.borders.df <- st_read("./data/ddr_border.gpkg")

kreise.df <- kreise.df %>% mutate(log_cases_per_mil = log(1+ 1e6*cases_simulation/pop_2018))

png("./output/Supplement_Figure_S2_Map.png", width=6, height=6, units="in", res=300)
plot(kreise.df["log_cases_per_mil"], breaks = "quantile", nbreaks=9, lwd= 0.5,
     pal=brewer.pal(9,"Blues"), main = "", #main = "Log(1 + Cases/Million)", 
     key.pos = 4,
     key.width = lcm(1.3), key.length = .8, reset = F)
plot(st_geometry(eastwest.borders.df), col="darkred", lwd=2, add=T)
dev.off()

pdf("./output/Supplement_Figure_S2_Map.pdf", width=6, height=6)
plot(kreise.df["log_cases_per_mil"], breaks = "quantile", nbreaks=9, lwd= 0.5,
     pal=brewer.pal(9,"Blues"), main = "", #main = "Log(1 + Cases/Million)", 
     key.pos = 4,
     key.width = lcm(1.3), key.length = .8, reset = F)
plot(st_geometry(eastwest.borders.df), col="darkred", lwd=2, add=T)
dev.off()
