# Set the working directory
setwd("~/Dropbox/UBC/Grants/Microplastics_Synergy/Supporting_Files")

# Load any packages
library(nlme)
library(ggplot2)
library(lme4)
library(viridis)

# Load in the data and carry out some basic data carpentry
data <- read.csv("sperm_database.csv")
data$normal_sperm <- data$X..normal.sperm/100
data$motility <- data$Mean.motility/100
data <- data[!is.na(data$normal_sperm),]

length(unique(data$Specie))

#Fit the model
FIT <- glmer(normal_sperm ~ Year + (1|Specie),
             family = binomial,
             data = data)

summary(FIT)

test <- data.frame(slope = as.vector(getME(FIT, "b")),
                   species = unique(data$Specie))


Logit_Fit <- function(x) {
  B_0 <- 40.099779
  B_1 <- -0.019753
  mu = exp(B_0 + B_1*x)/(1+exp(B_0 + B_1*x))
  mu
}

DATA <- data.frame(Year = seq(1977,2022, 0.1),
                   normal_sprem = Logit_Fit(seq(1977,2022, 0.1))*100)

FIG <- 
  ggplot() +
  geom_point(data = data, aes(y = X..normal.sperm, x = Year, col = Specie),
             alpha = 0.5, stroke = 0, shape=16) +
  geom_smooth(data = DATA, aes(x = Year, y = normal_sprem), se = F, span = 1.5, col = "black") +
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

ggsave(FIG,
       width = 3.23, height = 3, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Fertility_Declines.png")


#################

data$abnormal_sprem <- 1 - data$normal_sperm

#Fit the model
FIT <- glmer(abnormal_sprem ~ Year + (1|Specie),
             family = binomial,
             data = data)

summary(FIT)
coef(FIT)

row.names(coef(FIT)$Specie)[order(coef(FIT)$Specie[,1], decreasing = T)]

Logit_Fit <- function(x) {
  B_0 <- -40.099779
  B_1 <- 0.019753
  mu = exp(B_0 + B_1*x)/(1+exp(B_0 + B_1*x))
  mu
}

DATA <- data.frame(Year = seq(1977,2022, 0.1),
                   abnormal_sprem = Logit_Fit(seq(1977,2022, 0.1))*100)


FIG <- 
  ggplot() +
  geom_point(data = data, aes(y = abnormal_sprem*100, x = Year, col = Specie),
             alpha = 0.5, stroke = 0, shape=16) +
  geom_smooth(data = DATA, aes(x = Year, y = abnormal_sprem), se = F, span = 1.5, col = "black") +
  scale_y_continuous(limits = c(0,100), expand = c(0,1)) +
  #scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015)) +
  scale_color_viridis(discrete = T) +
  ylab("Percent abnormal sperm (%)") +
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

ggsave(FIG,
       width = 3.23, height = 3, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Abnormality_Increase.png")
