# Set the working directory
setwd("/Users/nozo/Library/CloudStorage/OneDrive-UBC/Quantitative_Ecology_Lab/Nozomu_Hirama_Undergraduate_Honours/Mammalian_Sperm_Abnormality/")
set.seed(88210729)

# Load packages
packages <- 
  list("ggplot2", "ggmcmc", "ggthemes", "ggridges", "tidyverse", "viridis",
       "sp", "sf", "rnaturalearth", "rnaturalearthdata", "patchwork", "brms", 
       "gridExtra", "bayesplot", "phangorn", "litsearchr", "ggtree", "shinystan", 
       "png", "ape", "phytools", "plotrix", "pals", "RColorBrewer")
sapply(packages, library, character=TRUE)

#Lit Review Randomized List of Species
Mammalian_Species <- read.csv("Data/MDD/MDD_v1.11_6649species.csv", na.strings=c("","NA"))
#Remove bats and rodents
table(Mammalian_Species$order)
Mammalian_Species <- Mammalian_Species[!Mammalian_Species$order=="CHIROPTERA",]
Mammalian_Species <- Mammalian_Species[!Mammalian_Species$order=="RODENTIA",]
#Generate randomized list of species
Ordered_Mammalian_Species <- Mammalian_Species[sample(1:nrow(Mammalian_Species)), ]
Ordered_Mammalian_Species$scientific_name <- gsub("_", " ", Ordered_Mammalian_Species$sciName)
Ordered_Mammalian_Species$order_id <- 1:nrow(Ordered_Mammalian_Species)
List_Species <- Ordered_Mammalian_Species[c("order_id", "mainCommonName", "scientific_name", "otherCommonNames")]
write.csv(List_Species, file="./Data/MDD/Randomized_Mammalian_Species.csv", row.names=FALSE)

#sample papers from search results
#create function
sample.paper <- function(results){
  #subtract 0.01 from min to ensure min is included in sample
  search_year_min <- min(results$year)-0.01
  #max is inclusive in function
  search_year_max <- max(results$year)
  year_bins <- seq(search_year_min, search_year_max, by = (search_year_max - search_year_min)/5)
  results$year_group <- cut(results$year, breaks = year_bins, labels = FALSE)
  articles <- results %>% group_by(year_group) %>% slice_sample(n=4)
  while(nrow(articles)<10){
    articles <- rbind(articles, results[!results$title %in% articles$title, ] %>% sample_n(1))
  }
  return(articles)
}

#import results
search_results <- import_results(file="./Data/Litsearch/Domestic_Horse.bib"); search_results$year <- as.numeric(search_results$year)
articles <- sample.paper(search_results)


# Load data
data <- read.csv("Data/Sperm_Abnormality_Database.csv", na.strings=c("","NA"))
supp_data <- read.csv("Data/Sperm_Abnormality_Supplementary.csv", na.strings=c("","NA"))


#Modify and clean data
data$normal_sperm <- data$X._Normal_Morphology/100
data$abnormal_sperm <- 1-data$normal_sperm
data$intact_acrosome <- data$X._Intact_Acrosome/100
data$head_abnormal <- data$X._Head_Abnormalities/100
data$midpiece_abnormal <- data$X._Midpiece_Abnormalities/100
data$tail_abnormal <- data$X._Tail_Abnormalities/100
data$motility <- data$Mean_Motility/100
data <- data[!is.na(data$Species),]
data$Species <- as.factor(data$Species)
#Round coordinates to 1 decimal place
data <- data %>%
  mutate(Lat = ifelse(is.na(Lat) | Lat == "multiple", Lat, round(as.numeric(Lat), 1)))
data <- data %>%
  mutate(Long = ifelse(is.na(Long) | Long == "multiple", Long, round(as.numeric(Long), 1)))
supp_data$Lat <- round(as.numeric(supp_data$Lat),1)
supp_data$Long <- round(as.numeric(supp_data$Long),1)
data$Publish_Date.date <- as.Date(data$Publish_Date)
data$Publish_Date <- as.numeric(as.Date(data$Publish_Date))


#Aggregate mutliple-entry data
data$original_ID <- data$ID
full_data <- data %>%
  left_join(supp_data, by = c("original_ID" = "ID"))
# for (supp_variable in noquote(ls(full_data)[ls(full_data) %in% ls(supp_data)|ls(full_data) %in% paste(ls(supp_data), ".y", sep="")])){
#   full_data %>%
#     mutate(supp_variable=coalesce(supp_variable, B))
# }

#Delete any species that only appeared once
# data <- data[(data$Species%in%((data%>%count(Species))[-which((data%>%count(Species))$n<=1),])$Species),]

length(unique(full_data$Binomial))

#Extract country data from rnaturalearth package
world <- ne_countries(scale = "medium", returnclass = "sf")

# #With multiple omitted
# countries <- data %>% count(Country, name="country_freq")
# countries <- countries[-which(countries$Country=="multiple"),]
# countries <- countries[!is.na(countries$Country),]
# world <- merge(world, countries, all = TRUE, by.x="iso_a3", by.y="Country")
# 
# #Plot country frequencies onto global map
# ggplot() +
#   geom_sf(data = world, aes(fill = country_freq), color = "black", size = 0.2) +
#   scale_fill_gradient(low = alpha("blue",0.2), high = alpha("blue",0.9), na.value = "grey") +
#   scale_alpha_continuous(range = c(0, 1), guide = "legend") +
#   labs(title = "Frequency Distribution of Countries") +
#   theme_minimal()
# 
# #Plot with variation of colors, original data
# for (color in c("blue", "#9E2EAF", "#A50F15")) {
#   maps_test <- list()
#   for (a in list(c(0.1, 0.7), c(0.1, 0.9), c(0.3, 0.7), c(0.3,0.9))) {
#     map_test <- ggplot() +
#       geom_sf(data = world, aes(fill = country_freq), color = "black", size = 0.2) +
#       scale_fill_gradient(low = alpha(color ,a[1]), high = alpha(color ,a[2]), na.value = "grey") +
#       scale_alpha_continuous(range = c(0, 1), guide = "legend") +
#       labs(title = paste("Frequency Distribution of Countries,", "low=", as.character(a[1]), "high=", as.character(a[2])), subtitle=paste0("(", length(unique(countries$Country)), " countries)")) +
#       xlab("Longitude") + ylab("Latitude") +
#       theme_minimal() +
#       theme(
#         legend.key.height = unit(0.3, "cm"),  # Adjust the height of legend keys
#         legend.key.width = unit(0.3, "cm"),   # Adjust the width of legend keys
#         plot.title = element_text(size = 10),  # Adjust title text size
#         plot.subtitle = element_text(size = 5),
#         axis.title=element_text(size=5),
#         legend.text = element_text(size = 5),  # Adjust legend text size
#         legend.title = element_text(size = 8)  # Adjust legend title size
#       )
#     maps_test[[length(maps_test) + 1]] <- map_test
#   }
#   wrap_plots(maps_test, ncol=1)
#   ggsave(last_plot(),
#          width = 10, height = 8,
#          dpi = 1200,
#          path="./Images",
#          bg = "transparent",
#          file=paste(color, ".png", sep=""))
# }
# 
# maps_test <- list()
# for (color in c("blue", "#9E2EAF", "#A50F15")) {
#   map_test <- ggplot() +
#     geom_sf(data = world, aes(fill = country_freq), color = "black", size = 0.2) +
#     scale_fill_gradient(low = alpha(color ,a[1]), high = alpha(color ,a[2]), na.value = "grey") +
#     scale_alpha_continuous(range = c(0, 1), guide = "legend") +
#     labs(title = paste("Frequency Distribution of Countries,", "low=", as.character(a[1]), "high=", as.character(a[2])), subtitle=paste0("(", length(unique(countries$Country)), " countries)")) +
#     xlab("Longitude") + ylab("Latitude") +
#     theme_minimal() +
#     theme(
#       legend.key.height = unit(0.3, "cm"),  # Adjust the height of legend keys
#       legend.key.width = unit(0.3, "cm"),   # Adjust the width of legend keys
#       plot.title = element_text(size = 10),  # Adjust title text size
#       plot.subtitle = element_text(size = 5),
#       axis.title=element_text(size=5),
#       legend.text = element_text(size = 5),
#       
#       legend.title = element_text(size = 8)
#     )
#     maps_test[[length(maps_test) + 1]] <- map_test
#   }
# wrap_plots(maps_test, ncol=3)
# ggsave(last_plot(),
#        width = 12, height = 8,
#        dpi = 1200,
#        path="./Images",
#        bg = "transparent",
#        file=paste("map_colors.png"))


coords <- data[c("ID", "Lat", "Long", "original_ID")][!(is.na(data$Lat)|(is.na(data$Long))), ] %>%
  left_join(supp_data[c("ID", "Lat", "Long")][!(is.na(supp_data$Lat)|(is.na(supp_data$Long))), ], by = c("original_ID" = "ID"))
coords <- coords %>%
  mutate(Lat = ifelse(Lat.x == "multiple", Lat.y, Lat.x)) %>%
  mutate(Long = ifelse(Long.x == "multiple", Long.y, Long.x)) %>%
  mutate(coord_vector = paste(Lat, Long, sep = ", "))
coord_freq <- data.frame(
  coord_vector = unique(coords$coord_vector),
  freq = tabulate(match(coords$coord_vector, unique(coords$coord_vector)))
)

unique_coords <- coords %>% 
  distinct(Lat, Long, .keep_all = TRUE)

# Calculate frequency within unique_coords
unique_coords <- coords %>%
  count(Lat, Long, name = "freq")

unique_coords$Lat <- as.numeric(unique_coords$Lat)
unique_coords$Long <- as.numeric(unique_coords$Long)
unique_coords$freq <- as.numeric(unique_coords$freq)

ggplot() +
  geom_sf(data = world, color = "black", size = 0.2, fill = "grey") +
  geom_point(data = unique_coords, aes(x = Long, y = Lat, size = freq), color = "red") +
  labs(title = "Frequency Distribution of Countries") +
  xlab("Longitude") + ylab("Latitude") +
  theme_minimal() +
  coord_sf(lims_method = "geometry_bbox") +
  scale_size_continuous(name = "Frequency")
ggsave(last_plot(),
       width = 12, height = 8,
       dpi = 1200,
       path="./Images",
       bg = "transparent",
       file=paste("map_circles_red.png"))


maps_test <- list()
for (color in c("blue", "#9E2EAF", "#A50F15")) {
  map_test <- ggplot() +
    geom_sf(data = world, color = "black", size = 0.2, fill = "grey") +
    geom_point(data = unique_coords, aes(x = Long, y = Lat, size = freq), color = color) +
    labs(title = "Frequency Distribution of Countries") +
    theme_minimal() +
    coord_sf(lims_method = "geometry_bbox")
  maps_test[[length(maps_test) + 1]] <- map_test
}
wrap_plots(maps_test, ncol=3)
ggsave(last_plot(),
       width = 12, height = 8,
       dpi = 1200,
       path="./Images",
       bg = "transparent",
       file="map_circles.png")



#Visualizing Data
ggplot(data,aes(x=Publish_Date.date,y=normal_sperm)) + 
  geom_point()

#Incorporate phylogeny data
source("./Scripts/Phylogeny.R")

data$Binomial <- gsub(" ", "_", data$Binomial)

#Fit the model for normal sperm
norm_sperm_species <- brm(
  normal_sperm ~ Publish_Date + (1+Publish_Date|Binomial), 
  data = data, 
  init=0, chains = 3, iter = 3000, warmup = 1000, 
  family = brmsfamily("Beta"), 
  set_prior("normal(0, 0.0001)", class = "b"))

species_ranef <- as.data.frame(ranef(norm_sperm_species))
phylo.ultra <- force.ultrametric(phylogeny, method="extend", message=FALSE)

# Add entries of "NA" for those not included in the dataset, but is in the phylogeny
missing_ranef <- data.frame(matrix(NA, ncol = 8, nrow = length(phylo.ultra$tip.label[!(phylo.ultra$tip.label %in% rownames(species_ranef))])),
                            row.names = phylo.ultra$tip.label[!(phylo.ultra$tip.label %in% rownames(species_ranef))])
names(missing_ranef) <- colnames(species_ranef)
species_full_ranef <- rbind(species_ranef, missing_ranef)

# Create vector of slopes
species_slopes <- setNames(species_full_ranef[,"Binomial.Estimate.Publish_Date"], rownames(species_full_ranef))
species_slopes <- species_slopes*365.25

slopes_available <- species_slopes[!is.na(species_slopes)]

# Plot slopes on tree (anc.ML method allows for NA entries)
lims <- c(min(species_slopes[!is.na(species_slopes)]), max(species_slopes[!is.na(species_slopes)]))
lims <- c(-max(abs(lims)), max(abs(lims)))
slope_tree <- contMap(phylo.ultra, slopes_available, lims=lims, method="anc.ML", plot=FALSE, 
                      sig=8)

slope_tree <- setMap(slope_tree, brewer.pal(n=11, name="RdYlBu"))

tips.na <- sapply(names(which(is.na(species_slopes))), function(x,y) which(y==x), y=slope_tree$tree$tip.label)
ancestors.na <- Ancestors(phylo.ultra, tips.na)
ancestors.available <- Ancestors(phylo.ultra, setdiff(1:Ntip(phylo.ultra), tips.na))
slopes.na <- setdiff(c(tips.na, unique(unlist(ancestors.na))),
                     c(setdiff(1:Ntip(phylo.ultra),tips.na),
                       unique(unlist(ancestors.available))))

slope_tree$tree <- paintBranches(slope_tree$tree, slopes.na, "NA")

cols <- setNames(c("grey", slope_tree$cols), c("NA",names(slope_tree$cols)))


species.info <- Mammalian_Species[Mammalian_Species$sciName %in% names(species_slopes),c("sciName", "order")]

species.info <- species.info %>%
  group_by(order) %>%
  mutate(num_in_order = length(sciName))
species.info <- species.info[-which(species.info$order=="CARNIVORA"),]
# Remove CARNIVORA as it spans throughout the tree which causes issues when adding lines in a fan plot.
# They will be added in manually based on the tip later.

species.info <- species.info %>%
  group_by(order) %>%
  mutate(common_ancestor = ifelse(num_in_order>1, findMRCA(phylo.ultra, sciName), which(phylo.ultra$tip.label==sciName)))

orders.info <- unique(species.info[c("order", "num_in_order", "common_ancestor")])
orders.name <- orders.info$order
orders.node <- orders.info$common_ancestor
orders.col <- unname(alphabet())
orders.text.size <- 0.3 + orders.info$num_in_order/50

line.offset <- c(1.2, 1.2, 1.28, rep(1.2, 14))
label.offset <- c(1.24, 1.24, 1.3, 1.24, 1.24, 1.24, 1.24, 1.28, 
                  1.24, 1.24, 1.28, rep(1.24, 6))

# Plot phylogeny with slopes
png("./Images/slopes.png", pointsize=10, width=6000, height=8000, res=600)
plot(slope_tree$tree, cols, type="fan", outline=TRUE, lwd=6,
     ftype="i", fsize=0.4, offset=4)
orders.labels <- mapply(arc.cladelabels, text=orders.name, node=orders.node, col=orders.col,
                        ln.offset=line.offset, lab.offset=label.offset, MoreArgs=list(mark.node=FALSE, lwd=6), 
                        cex=orders.text.size)
# Change radius and degrees depending on plot size and other factors
draw.arc(radius=376.5, deg1=64.2, deg2=-3, lwd=6, col=orders.col[length(orders.info$order)-1])
arctext("CARNIVORA", radius=390,
        middle=mean(c((64.2*pi/180), (-3*pi/180)), cex=2))
add.color.bar(600, slope_tree$cols, title="Estimated Annual Change in Normal Sperm Proportion \n", lims=lims,
              digits=3, subtitle="", x=-300, y=-475, prompt=FALSE, fsize=1.5, lwd=10)
legend(x=-312, y=-475, legend="missing", pch=22,
       pt.bg="grey", bty="n", pt.cex=3, cex=1.5)
dev.off()


# plot(obj$tree,cols,split.vertical=TRUE,outline=TRUE,lwd=6,
#      ftype="i",fsize=0.7)
# add.color.bar(40,obj$cols,title="trait",lims=range(x,na.rm=TRUE),
#               digits=2,subtitle="",x=0,y=0.97*Ntip(mammal.tree),prompt=FALSE)
# legend(x=-3,y=0.95*Ntip(mammal.tree),legend="missing",pch=22,
#        pt.bg="white",bty="n",pt.cex=2)


ggtree(phylogeny, branch.length = "none") + 
  theme_tree2() +
  geom_tiplab(size=2)

ggsave(last_plot(),
       width = 10, height = 12,
       dpi = 1200,
       path="./Images",
       file="phylo.png")

# norm_sperm_species_ce_95 <- conditional_effects(norm_sperm_species, effects = "Publish_Date:Binomial", prob=0.95)
# plot(norm_sperm_species_ce_95, points=TRUE) [[1]] + 
#   scale_x_continuous(labels = function(x) as.Date(x, origin = "1970-01-01")) + 
#   labs(x = "Year", y = "Normal Sperm with 95% CI") +
#   theme(legend.position="bottom", legend.key.size = unit(2, "mm")) +
#   guides(fill=guide_legend(nrow=2, byrow=TRUE))
# ggsave(last_plot(),
#        width = 12, height = 8,
#        dpi = 1200,
#        path="./Images",
#        bg = "transparent",
#        file="Species_Effect_95_Legend.png")

#All Measures
measures <- c("normal_sperm", "intact_acrosome", "head_abnormal", "midpiece_abnormal", "tail_abnormal", "motility")
labels <- c("Normal Sperm", "Intact Acrosomes", "Head Abnormalities", "Mid-piece Abnormalities", "Tail Abnormalities", "Motile Sperm")

models <- list()
for (measure in measures) {
  formula <- as.formula(paste(measure, "~ Publish_Date + (1|gr(Binomial, cov = phylo.cov))"))
  model <- brm(
    formula,
    init=0, chains = 3, iter = 3000, warmup = 1000,
    set_prior("normal(0, 0.0001)", class = "b"), 
    data = data,
    family = brmsfamily("Beta"),
    data2 = list(phylo.cov = phylo.cov), 
    silent = TRUE
  )
  models[[measure]] <- model
}

conditional_effects_list <- list()
for (measure in measures){
  conditional_effects_list[[measure]] <- conditional_effects(models[[measure]], effects = "Publish_Date", prob = 0.95)
}

plots <- list()
tags=c("A", "B", "C", "D", "E", "F")
if (length(measures)==length(labels)){
  for (i in 1:length(measures)){
    plots[[measures[i]]] <- plot(conditional_effects_list[[measures[i]]], plot=FALSE, points=TRUE) [[1]] + 
      scale_x_continuous(labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%Y"), 
                         breaks = c(seq(as.Date("1975-01-01"), as.Date("2025-06-15"), by = "5 years"))) +
      labs(x = "Year", y = labels[i], tag = tags[i]) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.tag = element_text(size = 12, family = "sans", face = "bold"),
            axis.title.y = element_text(size = 25, family = "sans", face = "bold", vjust = 1.5),
            axis.title.x = element_text(size = 25, family = "sans", face = "bold"),
            axis.text.y = element_text(size = 10, family = "sans"),
            axis.text.x = element_text(size = 10, family = "sans"), 
            panel.background = element_blank(),
            panel.border = element_rect(linewidth = 0.2, fill = NA),
            plot.background = element_blank(),
            plot.margin = unit(c(0.2,0.4,0,0.2), "cm"))
  }
} else {
  print("Measures count and Labels count do not match up")
}

Traits.time <- 
  arrangeGrob(plots[["normal_sperm"]], 
              plots[["intact_acrosome"]], 
              plots[["head_abnormal"]], 
              plots[["midpiece_abnormal"]], 
              plots[["tail_abnormal"]], 
              plots[["motility"]], 
              ncol=3)

ggsave(Traits.time,
       width = 15, height = 10,
       dpi = 1200,
       path="./Images",
       file="Sperm_Traits.png")

#Change font to ggplot recognized font
plot(models[["normal_sperm"]], ask = FALSE)
bayesplot_theme_set(theme_default() + theme(text=element_text(family="Arial")))
pp_check(models[["normal_sperm"]], ndraws = 100)

# launch_shinystan(models[["normal_sperm"]])


