# required
lop <- c("sf", "tidyverse", "lubridate", "RColorBrewer")
newp <- lop[!(lop %in% installed.packages()[,"Package"])]
if(length(newp)) install.packages(newp)
lapply(lop, require, character.only = TRUE)

# read shape
kreise.df <- st_read('./data/RKI_Corona_Landkreise.gpkg', stringsAsFactors=F)
table(kreise.df$BL)

# only use geoinformation and take data from daily scrapes
kreise.df <- kreise.df %>% select(id=RS, county=GEN, type = BEZ, nuts_id = NUTS, state = BL, pop_latest = EWZ, area = KFL)
names(kreise.df)

## add daily rki data, use apr 26
input.df <- read_csv("./data/kreise_counts_panel_daily.csv") %>% select(-state, -pop_latest,-area, -type)
input.df <- input.df %>%  group_by(id) %>% filter(date==as_date("2020-04-26")) %>%  select(-starts_with("new"))
kreise.df <- kreise.df %>% left_join(input.df)

## now aggregate over berlin
kreise.df <- kreise.df %>% mutate(id = replace(id, state=="Berlin", 11000))
kreise.df <- kreise.df %>% group_by(id) %>% summarize(cum_cases = sum(cum_cases), cum_deaths = sum(cum_deaths))
kreise.df %>% summarize(cum_cases = sum(cum_cases), cum_deaths = sum(cum_deaths))

# add pop
input.df <- read_csv("./data/pop_2017_2018.csv")
kreise.df <- kreise.df %>% left_join(input.df)

# add border
eastwest.borders.df <- st_read("./data/ddr_border.gpkg")

kreise.df <- kreise.df %>% mutate(log_cases_per_mil = log(1+ 1e6*cum_cases/pop_2018))

png("./output/Figure_1.png", width=6, height=6, units="in", res=300)
plot(kreise.df["log_cases_per_mil"], breaks = "quantile", nbreaks=9, lwd= 0.5,
     pal=brewer.pal(9,"Blues"), main = "", #main = "Log(1 + Cases/Million)", 
     key.pos = 4,
     key.width = lcm(1.3), key.length = .8, reset = F)
plot(st_geometry(eastwest.borders.df), col="darkred", lwd=2, add=T)
dev.off()

pdf("./output/Figure_1.pdf", width=6, height=6)
plot(kreise.df["log_cases_per_mil"], breaks = "quantile", nbreaks=9, lwd= 0.5,
     pal=brewer.pal(9,"Blues"), main = "", #main = "Log(1 + Cases/Million)", 
     key.pos = 4,
     key.width = lcm(1.3), key.length = .8, reset = F)
plot(st_geometry(eastwest.borders.df), col="darkred", lwd=2, add=T)
dev.off()
