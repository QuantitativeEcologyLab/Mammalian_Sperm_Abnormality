# Set the working directory
setwd("/Users/nozo/Library/CloudStorage/OneDrive-UBC/Quantitative_Ecology_Lab/Nozomu_Hirama_Undergraduate_Honours/Mammalian_Sperm_Abnormality/")
set.seed(88210729)

# Load packages
packages <- 
  list("ggplot2", "ggmcmc", "ggthemes", "ggridges", "tidyverse", "viridis",
       "sp", "sf", "rnaturalearth", "rnaturalearthdata", "patchwork", "brms", 
       "gridExtra", "grid", "bayesplot", "phangorn", "litsearchr", "ggtree", "shinystan", 
       "png", "ape", "phytools", "plotrix", "pals", "RColorBrewer")
sapply(packages, library, character=TRUE)

#Lit Review Randomized List of Species
Mammalian_Species <- read.csv("Data/MDD/MDD_v1.11_6649species.csv", na.strings=c("", NA))
#Remove bats and rodents
table(Mammalian_Species$order)
Mammalian_Species <- Mammalian_Species[!Mammalian_Species$order=="CHIROPTERA",]
Mammalian_Species <- Mammalian_Species[!Mammalian_Species$order=="RODENTIA",]

#Generate randomized list of species
Ordered_Mammalian_Species <- Mammalian_Species[sample(1:nrow(Mammalian_Species)), ]
Ordered_Mammalian_Species$scientific_name <- gsub("_", " ", Ordered_Mammalian_Species$sciName)
Ordered_Mammalian_Species$order_id <- 1:nrow(Ordered_Mammalian_Species)

#Remove infraorder cetacea (whales, dolphins)
#If conducting fresh literature search, remove these before generating randomized list with the bats and rodents.
#If adding on to pre-existing database, continue with initial randomized order, since the initial list used
#for the study was generated with cetaceans included.
Mammalian_Species <- subset(Mammalian_Species, !(infraorder=="CETACEA"&!is.na(infraorder)))
Ordered_Mammalian_Species <- subset(Ordered_Mammalian_Species, !(infraorder=="CETACEA"&!is.na(infraorder)))
List_Species <- Ordered_Mammalian_Species[c("order_id", "mainCommonName", "scientific_name", "otherCommonNames")]
write.csv(List_Species, file="./Data/MDD/Randomized_Mammalian_Species.csv", row.names=FALSE)

#Optional function to sample 4 papers from each of the 5 evenly divided time frames in search result
#Useful when there are too many papers. Uncomment to use.
# #sample papers from search results
# #create function
# sample.paper <- function(results){
#   #subtract 0.01 from min to ensure min is included in sample
#   search_year_min <- min(results$year)-0.01
#   #max is inclusive in function
#   search_year_max <- max(results$year)
#   year_bins <- seq(search_year_min, search_year_max, by = (search_year_max - search_year_min)/5)
#   results$year_group <- cut(results$year, breaks = year_bins, labels = FALSE)
#   articles <- results %>% group_by(year_group) %>% slice_sample(n=4)
#   while(nrow(articles)<20){
#     articles <- rbind(articles, results[!results$title %in% articles$title, ] %>% sample_n(1))
#   }
#   return(articles)
# }
# 
# #import results
# search_results <- import_results(file="./Data/Litsearch/Domestic_Horse.bib"); search_results$year <- as.numeric(search_results$year)
# articles <- sample.paper(search_results)


# Load data
data <- read.csv("Data/Sperm_Abnormality_Database.csv", na.strings=c("",NA))
supp_data <- read.csv("Data/Sperm_Abnormality_Supplementary.csv", na.strings=c("",NA))


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


length(unique(full_data$Binomial))

#Extract country data from rnaturalearth package
world <- ne_countries(scale = "medium", returnclass = "sf")

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

map.plot <- ggplot() +
  geom_sf(data = world, color = "#6f6f6f", size = 0.2, fill = "#d3d3d3") +
  geom_point(data = unique_coords, aes(x = Long, y = Lat, size = freq), color = "#0651f0") +
  # labs(title = "Frequency Distribution of Countries") +
  # xlab("Longitude") + ylab("Latitude") +
  theme_minimal() +
  coord_sf(lims_method = "geometry_bbox") +
  scale_size_continuous(name = "Frequency") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        # legend.position = "none",
        panel.background = element_blank(),
        plot.background = element_blank())

map.plot
ggsave(last_plot(),
       width = 12, height = 8,
       dpi = 600,
       path="./Images",
       bg = "transparent",
       file=paste("map_circles_blue.png"))

#Visualizing Data
ggplot(data,aes(x=Publish_Date.date,y=normal_sperm)) + 
  geom_point()

#Create mosaic plot of collection method and housing conditions
data.present <- subset(data, Presence_of_Data=="Yes")
collection.prop <- data.frame(prop.table(table(data.present$Collection_Method, exclude=NULL)))
colnames(collection.prop) <- c("Collection_Method", "Proportion")
housing.prop <- data.frame(prop.table(table(data.present$Housing, exclude=NULL)))
colnames(housing.prop) <- c("Housing_Condition", "Proportion")

collection.prop$Collection_Method <- factor(collection.prop$Collection_Method,
                                            levels = collection.prop$Collection_Method[order(-collection.prop$Proportion, decreasing = TRUE)])
housing.prop$Housing_Condition <- factor(housing.prop$Housing_Condition,
                                         levels = housing.prop$Housing_Condition[order(-housing.prop$Proportion, decreasing=TRUE)])

collection_plot <- ggplot(collection.prop, aes(x = "", y = Proportion, fill = Collection_Method)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Set3", na.value="grey") +
  xlab("Collection Method") +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) 

# Create stacked bar plot for Housing
housing_plot <- ggplot(housing.prop, aes(x = "", y = Proportion, fill = Housing_Condition)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Set2", na.value="grey") +
  xlab("Housing Condition") +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) 

combined_plot <- grid.arrange(arrangeGrob(map.plot,
                                          arrangeGrob(collection_plot, housing_plot, ncol = 2), 
                                          heights = c(2, 1)), nrow = 1)

ggsave(combined_plot,
       width = 12, height = 10,
       dpi = 600,
       path="./Images",
       bg = "transparent",
       file=paste("data_distribution.png"))



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


species.info <- Mammalian_Species[Mammalian_Species$sciName %in% names(species_slopes), c("sciName", "order", "iucnStatus")]

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
pics <- c("Crocidura_russula","Equus_ferus","Felis_catus","Gorilla_gorilla","Homo_sapiens", 
          "Loxodonta_africana", "Myrmecophaga_tridactyla", "Ovis_orientalis", "Panthera_leo", "Vulpes_lagopus")

# Plot phylogeny with slopes
png("./Images/slopes.png", pointsize=10, width=6000, height=8000, res=600)
plot(slope_tree$tree, cols, type="fan", outline=TRUE, lwd=6,
     ftype="i", fsize=0.4, offset=4)
orders.labels <- mapply(arc.cladelabels, text=orders.name, node=orders.node, col=orders.col,
                        ln.offset=line.offset, lab.offset=label.offset, MoreArgs=list(mark.node=FALSE, lwd=6), 
                        cex=orders.text.size)
# Change radius and degrees depending on plot size and other factors
draw.arc(radius=376.5, deg1=66.5, deg2=-3, lwd=6, col=orders.col[length(orders.info$order)-1])
arctext("CARNIVORA", radius=390,
        middle=mean(c((66.5*pi/180), (-3*pi/180)), cex=2))
add.color.bar(600, slope_tree$cols, title="Estimated Annual Change in Normal Sperm Proportion \n", lims=lims,
              digits=3, subtitle="", x=-300, y=-475, prompt=FALSE, fsize=1.5, lwd=10)
legend(x=-312, y=-475, legend="missing", pch=22,
       pt.bg="grey", bty="n", pt.cex=3, cex=1.5)

dev.off()

# IUCN and Slopes
species.iucn <- Mammalian_Species[Mammalian_Species$sciName %in% names(slopes_available), c("sciName", "iucnStatus")]
slopes.iucn <- data.frame(sciName=names(slopes_available), 
                          slope=slopes_available)
slopes.iucn <- merge(slopes.iucn, species.iucn, by.x = "sciName", by.y = "sciName")
slopes.iucn$iucnStatus <- factor(slopes.iucn$iucnStatus, 
                                    levels=c("DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX", "NE"))

ggplot(slopes.iucn, aes(x=iucnStatus, y=slope, fill=iucnStatus)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, na.rm = TRUE) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  labs(x = "IUCN Status", y = "Annual Change in Normal Sperm Proportion")


ggsave(last_plot(),
       width = 8, height = 8,
       dpi = 600,
       path="./Images",
       file="IUCN_Slopes.png")

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
            plot.tag = element_text(size = 18, face = "bold"),
            axis.title.y = element_text(size = 25, face = "bold", vjust = 1.5),
            axis.title.x = element_text(size = 18, face = "bold"),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10), 
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
       dpi = 600,
       path="./Images",
       file="Sperm_Traits.png")


plot(models[["normal_sperm"]])
plot(models[["intact_acrosome"]])
plot(models[["head_abnormal"]])
plot(models[["midpiece_abnormal"]])
plot(models[["tail_abnormal"]])
plot(models[["motility"]])



