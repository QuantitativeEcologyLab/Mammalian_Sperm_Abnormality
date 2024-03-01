# Set the working directory
setwd("/Users/nozo/Library/CloudStorage/OneDrive-UBC/Quantitative_Ecology_Lab/Nozomu_Hirama_Undergraduate_Honours/Mammalian_Sperm_Abnormality/")
set.seed(88210729)

# Load any packages
library(ggplot2)
library(ggmcmc)
library(ggthemes)
library(ggridges)
library(tidyverse)
library(viridis)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(lubridate)
library(brms)
library(glmmTMB)
library(dplyr)
library(litsearchr)

#Lit Review Randomized List of Species
Mammalian_Species<-read.csv("Data/MDD/MDD_v1.11_6649species.csv", na.strings=c("","NA"))
#Remove bats and rodents
table(Mammalian_Species$order)
Mammalian_Species<-Mammalian_Species[!Mammalian_Species$order=="CHIROPTERA",]
Mammalian_Species<-Mammalian_Species[!Mammalian_Species$order=="RODENTIA",]
#Generate randomized list of species
Ordered_Mammalian_Species<-Mammalian_Species[sample(1:nrow(Mammalian_Species)), ]
Ordered_Mammalian_Species$scientific_name<-gsub("_", " ", Ordered_Mammalian_Species$sciName)
Ordered_Mammalian_Species$order_id<-1:nrow(Ordered_Mammalian_Species)
List_Species<-Ordered_Mammalian_Species[c("order_id", "mainCommonName", "scientific_name", "otherCommonNames")]
write.csv(List_Species, file="./Data/MDD/Randomized_Mammalian_Species.csv", row.names=FALSE)

#sample papers from search results
#create function
sample.paper<-function(results){
  #subtract 0.01 from min to ensure min is included in sample
  search_year_min<-min(results$year)-0.01
  #max is inclusive in function
  search_year_max<-max(results$year)
  year_bins<-seq(search_year_min, search_year_max, by = (search_year_max - search_year_min)/5)
  results$year_group <- cut(results$year, breaks = year_bins, labels = FALSE)
  articles <- results %>% group_by(year_group) %>% slice_sample(n=2)
  while(nrow(articles)<10){
    articles<-rbind(articles, results[!results$title %in% articles$title, ] %>% sample_n(1))
  }
  return(articles)
}

#import results
search_results<-import_results(file="./Data/Litsearch/Common_Brushtailed_Possum.bib"); search_results$year<-as.numeric(search_results$year)
articles<-sample.paper(search_results)
# View(articles)
# View(articles[,"doi"])

#remove the articles specified at the given indices, and re-run paper sample to obtain new ones
renew.papers<-function(search_results, articles, paper.indices){
  search_results<-search_results[!search_results$unique_id %in% articles$unique_id[paper.indices],]
  articles<-articles[-paper.indices,]
  while(nrow(articles)<10){
    articles<-rbind(articles, search_results[!search_results$title %in% articles$title, ] %>% sample_n(1))
  }
  return(articles)
}

articles<-renew.papers(search_results, articles, 1:10)
# View(articles)
# View(articles[,"doi"])



# Load data
data<-read.csv("Data/Sperm_Abnormality_Database.csv", na.strings=c("","NA"))
supp_data<-read.csv("Data/Sperm_Abnormality_Supplementary.csv", na.strings=c("","NA"))


#Modify and clean data
data$normal_sperm<-data$X._Normal_Morphology/100
data$abnormal_sperm<-1-data$normal_sperm
data$intact_acrosome<-data$X._Intact_Acrosome/100
data$head_abnormal<-data$X._Head_Abnormalities/100
data$midpiece_abnormal<-data$X._Midpiece_Abnormalities/100
data$tail_abnormal<-data$X._Tail_Abnormalities/100
data$motility<-data$Mean_Motility/100
data<-data[!is.na(data$Species),]
data$Species<-as.factor(data$Species)
data <- data %>%
  mutate(Lat = ifelse(is.na(Lat) | Lat == "multiple", Lat, round(as.numeric(Lat), 1)))
data <- data %>%
  mutate(Long = ifelse(is.na(Long) | Long == "multiple", Long, round(as.numeric(Long), 1)))
supp_data$Lat<-round(as.numeric(supp_data$Lat),1)
supp_data$Long<-round(as.numeric(supp_data$Long),1)
data$Publish_Date.date<-as.Date(data$Publish_Date)
data$Publish_Date<-as.numeric(as.Date(data$Publish_Date))


#Aggregate mutliple-entry data
data$original_ID<-data$ID
full_data <- data %>%
  left_join(supp_data, by = c("original_ID" = "ID"))
# for (supp_variable in noquote(ls(full_data)[ls(full_data) %in% ls(supp_data)|ls(full_data) %in% paste(ls(supp_data), ".y", sep="")])){
#   full_data %>%
#     mutate(supp_variable=coalesce(supp_variable, B))
# }

#Delete any species that only appeared once
# data<-data[(data$Species%in%((data%>%count(Species))[-which((data%>%count(Species))$n<=1),])$Species),]

length(unique(full_data$Binomial))

#Extract country data from rnaturalearth package
world <- ne_countries(scale = "medium", returnclass = "sf")

#With multiple omitted
countries<-data %>% count(Country, name="country_freq")
countries<-countries[-which(countries$Country=="multiple"),]
countries<-countries[!is.na(countries$Country),]
world <- merge(world, countries, all = TRUE, by.x="iso_a3", by.y="Country")

#Plot country frequencies onto global map
ggplot() +
  geom_sf(data = world, aes(fill = country_freq), color = "black", size = 0.2) +
  scale_fill_gradient(low = alpha("blue",0.2), high = alpha("blue",0.9), na.value = "grey") +
  scale_alpha_continuous(range = c(0, 1), guide = "legend") +
  labs(title = "Frequency Distribution of Countries") +
  theme_minimal()

#Plot with variation of colors, original data

for (color in c("blue", "#9E2EAF", "#A50F15")) {
  maps_test<-list()
  for (a in list(c(0.1, 0.7), c(0.1, 0.9), c(0.3, 0.7), c(0.3,0.9))) {
    map_test<- ggplot() +
      geom_sf(data = world, aes(fill = country_freq), color = "black", size = 0.2) +
      scale_fill_gradient(low = alpha(color ,a[1]), high = alpha(color ,a[2]), na.value = "grey") +
      scale_alpha_continuous(range = c(0, 1), guide = "legend") +
      labs(title = paste("Frequency Distribution of Countries,", "low=", as.character(a[1]), "high=", as.character(a[2])), subtitle=paste0("(", length(unique(countries$Country)), " countries)")) +
      xlab("Longitude") + ylab("Latitude") +
      theme_minimal() +
      theme(
        legend.key.height = unit(0.3, "cm"),  # Adjust the height of legend keys
        legend.key.width = unit(0.3, "cm"),   # Adjust the width of legend keys
        plot.title = element_text(size = 10),  # Adjust title text size
        plot.subtitle = element_text(size = 5),
        axis.title=element_text(size=5),
        legend.text = element_text(size = 5),  # Adjust legend text size
        legend.title = element_text(size = 8)  # Adjust legend title size
      )
    maps_test[[length(maps_test) + 1]]<-map_test
  }
  wrap_plots(maps_test, ncol=1)
  ggsave(last_plot(),
         width = 10, height = 8,
         dpi = 1200,
         path="./Images",
         bg = "transparent",
         file=paste(color, ".png", sep=""))
}

maps_test<-list()
for (color in c("blue", "#9E2EAF", "#A50F15")) {
  map_test<- ggplot() +
    geom_sf(data = world, aes(fill = country_freq), color = "black", size = 0.2) +
    scale_fill_gradient(low = alpha(color ,a[1]), high = alpha(color ,a[2]), na.value = "grey") +
    scale_alpha_continuous(range = c(0, 1), guide = "legend") +
    labs(title = paste("Frequency Distribution of Countries,", "low=", as.character(a[1]), "high=", as.character(a[2])), subtitle=paste0("(", length(unique(countries$Country)), " countries)")) +
    xlab("Longitude") + ylab("Latitude") +
    theme_minimal() +
    theme(
      legend.key.height = unit(0.3, "cm"),  # Adjust the height of legend keys
      legend.key.width = unit(0.3, "cm"),   # Adjust the width of legend keys
      plot.title = element_text(size = 10),  # Adjust title text size
      plot.subtitle = element_text(size = 5),
      axis.title=element_text(size=5),
      legend.text = element_text(size = 5),  # Adjust legend text size
      legend.title = element_text(size = 8)  # Adjust legend title size
    )
    maps_test[[length(maps_test) + 1]]<-map_test
  }
wrap_plots(maps_test, ncol=3)
ggsave(last_plot(),
       width = 12, height = 8,
       dpi = 1200,
       path="./Images",
       bg = "transparent",
       file=paste("map_colors.png"))


coords <- data[c("ID", "Lat", "Long", "original_ID")][!(is.na(data$Lat)|(is.na(data$Long))), ] %>%
  left_join(supp_data[c("ID", "Lat", "Long")][!(is.na(supp_data$Lat)|(is.na(supp_data$Long))), ], by = c("original_ID" = "ID"))
coords <- coords %>%
  mutate(Lat = ifelse(Lat.x == "multiple", Lat.y, Lat.x)) %>%
  mutate(Long = ifelse(Long.x == "multiple", Long.y, Long.x)) %>%
  mutate(coord_vector = paste(Lat, Long, sep = ", "))
coord_freq <- data.frame(
  coord_vector = unique(coords$coord_vector),  # Extract unique values from coord_vector
  freq = tabulate(match(coords$coord_vector, unique(coords$coord_vector)))  # Count occurrences of each unique value
)

unique_coords <- coords %>% 
  distinct(Lat, Long, .keep_all = TRUE)

# Calculate frequency within unique_coords
unique_coords <- coords %>%
  count(Lat, Long, name = "freq")

unique_coords$Lat<-as.numeric(unique_coords$Lat)
unique_coords$Long<-as.numeric(unique_coords$Long)
unique_coords$freq<-as.numeric(unique_coords$freq)

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


maps_test<-list()
for (color in c("blue", "#9E2EAF", "#A50F15")) {
  map_test<- ggplot() +
    geom_sf(data = world, color = "black", size = 0.2, fill = "grey") +
    geom_point(data = unique_coords, aes(x = Long, y = Lat, size = freq), color = color) +
    labs(title = "Frequency Distribution of Countries") +
    theme_minimal() +
    coord_sf(lims_method = "geometry_bbox")
  maps_test[[length(maps_test) + 1]]<-map_test
}
wrap_plots(maps_test, ncol=3)
ggsave(last_plot(),
       width = 12, height = 8,
       dpi = 1200,
       path="./Images",
       bg = "transparent",
       file="map_circles.png")

#Fit the model for normal sperm
norm_sperm_brm<-brm(normal_sperm ~ Publish_Date*Species + (1+Publish_Date|Species), data = data, init=0, chains = 3, iter = 3000, warmup = 1000, family = brmsfamily("Beta", link="logit"))
norm_sperm_ce<-conditional_effects(norm_sperm_brm, effects = "Publish_Date")
norm_sperm_plot<-plot(norm_sperm_ce, plot=FALSE)[[1]]

ggplot()+
  theme_void() +
  norm_sperm_plot

ggsave(last_plot(),
       width = 12, height = 8,
       dpi = 1200,
       path="./Images",
       file="Normal_Sperm_Mixed_95.png")

  

  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE) +  # Add species-specific trend lines
  geom_smooth(data = overall_df, aes(y = normal_sperm),  # Add overall trend line
              color = "black", linetype = "dashed", size = 1.5) +  # Change size to make it thicker
  labs(x = "Year", y = "% Normal Sperm",  # Axes labels
       title = "Temporal Trends in Mammalian Sperm Morphology") +  # Title
  theme_minimal() +  # Optional: Change plot theme
  theme(legend.position = "none",  # Remove legend
        plot.title = element_text(hjust = 0.5)) +  # Center title
  scale_color_viridis_d(option = "C")  # Use colorblind-friendly palette

ggsave(last_plot(),
       width = 12, height = 8,
       dpi = 1200,
       path="./Images",
       file="Normal_Sperm.png")



newdata=data.frame(Publish_Date = as.Date(data$Publish_Date))





norm_sperm_mixed <- lmer(normal_sperm ~ Publish_Date + (1|Species), data=data, na.action=na.omit)

overall_lm <- lm(normal_sperm ~ Publish_Date, data = data, na.action=na.omit)

# Create a new dataframe for predicted values from overall_lm
overall_df <- data.frame(
  Publish_Date = data$Publish_Date, 
  normal_sperm = predict(overall_lm, newdata = data),  # Use newdata argument to ensure the same length
  Species = rep("Overall", nrow(data))  # Adding a constant "Overall" for species
)

# Plot
ggplot(data, aes(x = Publish_Date, y = normal_sperm, color = Species)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE) +  # Add species-specific trend lines
  geom_smooth(data = overall_df, aes(y = normal_sperm),  # Add overall trend line
              color = "black", linetype = "dashed", size = 1.5) +  # Change size to make it thicker
  labs(x = "Year", y = "% Normal Sperm",  # Axes labels
       title = "Temporal Trends in Mammalian Sperm Morphology") +  # Title
  theme_minimal() +  # Optional: Change plot theme
  theme(legend.position = "none",  # Remove legend
        plot.title = element_text(hjust = 0.5)) +  # Center title
  scale_color_viridis_d(option = "C")  # Use colorblind-friendly palette

ggsave(last_plot(),
       width = 5, height = 5,
       dpi = 1200,
       path="./Images",
       file="Normal_Sperm.png")






# Plot
ggplot(data, aes(x = Publish_Date, y = normal_sperm, color = Species)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE) +  # Add species-specific trend lines
  geom_smooth(data = overall_df, aes(y = normal_sperm),  # Add overall trend line
              color = "black", linetype = "dashed", size = 1.5) +  # Change size to make it thicker
  labs(x = "Publish Date", y = "% Normal Sperm",  # Axes labels
       title = "Relationship between Publish Date and % Normal Sperm") +  # Title
  theme_minimal() +  # Optional: Change plot theme
  theme(legend.position = "none") +  # Remove legend
  scale_color_viridis_d(option = "C")  # Use colorblind-friendly palette
















ggplot(data, aes(x = Publish_Date, y = normal_sperm, color = Species)) +
  geom_point() +
  geom_line(aes(y = predicted)) +
  labs(x = "Publish Date", y = "% Normal Sperm") +
  ggtitle("Relationship between Publish Date and Normal Sperm by Species") +
  theme_minimal()








norm_sperm_glme <- glmer(normal_sperm ~ Publish_date + (date|Species), 
                         family = binomial, 
                         data = data[!is.na(data$normal_sperm), ])

summary(norm_sperm_glme)

test <- data.frame(slope = as.vector(getME(norm_sperm_glme, "b")),
                   species = unique(data[!is.na(data$normal_sperm), ]$Species))

coefficients_norm_sperm_glmer <- fixef(norm_sperm_glme)

norm_sperm_logit_fit <- function(x) {
  B_0 <- coefficients_norm_sperm_glmer[1]
  B_1 <- coefficients_norm_sperm_glmer[2]
  mu = exp(B_0 + B_1*x)/(1+exp(B_0 + B_1*x))
  mu
}

norm_sperm_data <- data.frame(Year = seq(1977,2022, 0.1),
                   normal_sperm = norm_sperm_logit_fit(seq(1977,2022, 0.1))*100)

norm_sperm_fig <- 
  ggplot() +
  geom_point(data = data[!is.na(data$X._Normal_Morphology), ], aes(y = X._Normal_Morphology, x = decimal_date(as.POSIXlt(Publish_Date)), col = Species),
             alpha = 0.5, stroke = 0, shape=16) +
  geom_smooth(data = norm_sperm_data, aes(x = Year, y = normal_sperm), se = F, span = 1.5, col = "black") +
  scale_y_continuous(limits = c(0,100), expand = c(0,1)) +
  #scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015)) +
  scale_color_viridis(discrete = T) +
  ylab("Percent normal sperm (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=12, family = "sans", face = "bold"),
        axis.title.x = element_text(size=12, family = "sans", face = "bold"),
        axis.text.y = element_text(size=10, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))

ggsave(norm_sperm_fig,
       width = 3.23, height = 3, units = "in",
       dpi = 600,
       bg = "white",
       path="./Images",
       file="Fertility_Declines.png")

#####
measures<-c("normal_sperm", "abnormal_sperm", "intact_acrosome", "head_abnormal", "midpiece_abnormal", "tail_abnormal", "motility")
for (trait in measures){
  paste(trait, "_glme", sep="") <- glmer(trait ~ decimal_date(as.POSIXlt(Publish_Date)) + (1|Species), 
                           family = binomial, 
                           data = data[!is.na(data$trait), ])
  
  summary(paste(trait, "_glme", sep=""))
  
  test <- data.frame(slope = as.vector(getME(paste(trait, "_glme", sep=""), "b")),
                     species = unique(data[!is.na(data$trait), ]$Species))
  
  paste("coefficients_", trait, "_glme", sep="") <- fixef(paste(trait, "_glme", sep=""))
  
  paste(trait, "_logit_fit", sep="") <- function(x) {
    B_0 <- paste("coefficients_", trait, "_glme", sep="")[1]
    B_1 <- paste("coefficients_", trait, "_glme", sep="")[2]
    mu = exp(B_0 + B_1*x)/(1+exp(B_0 + B_1*x))
    mu
  }
}


#####

norm_sperm_glme <- glmer(normal_sperm ~ decimal_date(as.POSIXlt(Publish_Date)) + (date|Species), 
                         family = binomial, 
                         data = data[!is.na(data$normal_sperm), ])

summary(norm_sperm_glme)

test <- data.frame(slope = as.vector(getME(norm_sperm_glme, "b")),
                   species = unique(data[!is.na(data$normal_sperm), ]$Species))

coefficients_norm_sperm_glmer <- fixef(norm_sperm_glme)


norm_sperm_brm <- brm(normal_sperm | trials(decimal_date(as.POSIXlt(Publish_Date))) ~
                             applicant.gender + (applicant.gender|dept),
                           family = binomial,
                           data = UCBadmit,
                           iter = 5000, chains = 3) 

initial_values <- list(
  b_Intercept = rnorm(nlevels(data$Species), 0, 1),  # Random intercepts for each Species
  b_Publish_Date = 0,  # Initialize slope parameter
  shape_param_alpha = 2,  # Example value for beta distribution parameter
  shape_param_beta = 2    # Example value for beta distribution parameter
)


norm_sperm_brm <- brm(
  data = data[!is.na(data$normal_sperm), ],
  family = brmsfamily("beta", link = "logit"),  
  normal_sperm ~ decimal_date(as.POSIXlt(Publish_Date)) + (1 + decimal_date(as.POSIXlt(Publish_Date)) | Species),
  iter =  5000,    
  warmup = 2000,   
  chains = 4,        
  cores = 4
)

norm_sperm_glm<-glmmTMB(normal_sperm ~ fixed_var1 + fixed_var2 + (random_var1 | random_var2),
        family = betabinomial(link = "logit"), data = data[!is.na(data$normal_sperm), ])

norm_sperm_glm <- glmmTMB(
  formula = normal_sperm ~ decimal_date(as.POSIXlt(Publish_Date)) + (1 + decimal_date(as.POSIXlt(Publish_Date)) | Species),
  data = data[!is.na(data$normal_sperm), ],
  family = betabinomial(link = "logit")
)

fixed_effects <- fixef(norm_sperm_glm)

# Extract random effects coefficients for each Species
random_effects <- ranef(norm_sperm_glm)$Species

# Combine fixed and random effects for each species
species_coefs <- coef(norm_sperm_glm)$Species

# Creating a sequence of dates from the dataset
seq_dates <- seq(min(as.numeric(as.POSIXlt(data_clean$Publish_Date))), 
                 max(as.numeric(as.POSIXlt(data_clean$Publish_Date))), 
                 by = 1)

# Creating a grid of dates for each species
pred_df <- expand.grid(Species = unique(data_clean$Species), decimal_date = seq_dates)

# Predicting values using the model for each species and date
pred_df$normal_sperm <- predict(norm_sperm_glm, newdata = pred_df, type = "response")

# Plotting normal_sperm over time for each species
ggplot(pred_df, aes(x = decimal_date, y = normal_sperm, color = Species)) +
  geom_line() +
  labs(title = "Normal Sperm Over Time by Species",
       x = "Decimal Date",
       y = "Normal Sperm") +
  theme_minimal()




library(ggplot2)
library(gridExtra)

# List of variables
variables <- c("normal_sperm", "abnormal_sperm", "intact_acrosome", "head_abnormal", "midpiece_abnormal", "tail_abnormal", "motility")

# List to store individual plots
plots <- list()

# Iterate over each variable
for (variable in variables) {
  # Fit linear model
  lm_model <- lm(paste(variable, "~ Publish_Date"), data = data, na.action = na.omit)
  
  # Create a new dataframe for predicted values
  predicted_df <- data.frame(
    Publish_Date = data$Publish_Date, 
    value = predict(lm_model, newdata = data),  # Use newdata argument to ensure the same length
    Variable = variable
  )
  
  # Create individual plot
  plot <- ggplot(data, aes(x = Publish_Date, y = !!sym(variable), color = Species)) +
    geom_point() +  # Add points
    geom_smooth(method = "lm", se = FALSE) +  # Add species-specific trend lines
    geom_smooth(data = predicted_df, aes(y = value),  # Add overall trend line
                color = "black", linetype = "dashed", size = 1.5) +  # Change size to make it thicker
    labs(x = "Publish Date", y = paste("%", variable, sep = " "),  # Axes labels
         title = paste("Relationship between Publish Date and", variable)) +  # Title
    theme_minimal() +  # Optional: Change plot theme
    theme(legend.position = "none") +  # Remove legend
    scale_color_viridis_d(option = "C")  # Use colorblind-friendly palette
  
  # Store the plot in the list
  plots[[variable]] <- plot
}

# Calculate the number of rows and columns based on the number of plots
num_plots <- length(plots)
num_cols <- ifelse(num_plots %% 2 == 0, 2, 3)
num_rows <- ceiling(num_plots / num_cols)

# Set a base height for each row
base_height <- rep(2, num_rows)  # Adjust this value as needed

# Calculate the total height
total_height <- sum(base_height)

# Arrange plots in a grid with adjusted aspect ratio
grid_arrange <- grid.arrange(grobs = plots, ncol = num_cols, heights = base_height)

# Save the grid of plots as an image
ggsave("Parameter_Trends.png", grid_arrange, width = 10, height = total_height, dpi = 300, path = "./Images", bg = "transparent")


