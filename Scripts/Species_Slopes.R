setwd("/Users/nozo/Library/CloudStorage/OneDrive-UBC/Quantitative_Ecology_Lab/Nozomu_Hirama_Undergraduate_Honours/Mammalian_Sperm_Abnormality/")
set.seed(88210729)

packages <- 
  list("ggplot2", "tidyverse", "brms", "ape", "phytools")
sapply(packages, library, character=TRUE)

# Load data
data <- read.csv("Data/Sperm_Abnormality_Database.csv", na.strings=c("","NA"))

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
data$Publish_Date.date <- as.Date(data$Publish_Date)
data$Publish_Date <- as.numeric(as.Date(data$Publish_Date))

#Incorporate phylogeny data
phylo.trees <- read.nexus("Data/Phylogeny/output.nex")
phylogeny <- ls.consensus(phylo.trees)

# Rename species to match entries for analyses
phylogeny$tip.label[which(phylogeny$tip.label=="Callicebus_caligatus")] <- "Plecturocebus_caligatus"
phylogeny$tip.label[which(phylogeny$tip.label=="Macropus_eugenii")] <- "Notamacropus_eugenii"
phylogeny$tip.label[which(phylogeny$tip.label=="Vicugna_vicugna")] <- "Lama_vicugna"
phylogeny$tip.label[which(phylogeny$tip.label=="Herpestes_fuscus")] <- "Urva_fusca"
phylogeny$tip.label[which(phylogeny$tip.label=="Mazama_pandora")] <- "Odocoileus_pandora"
phylogeny$tip.label[which(phylogeny$tip.label=="Mustela_africana")] <- "Neogale_africana"
phylogeny$tip.label[which(phylogeny$tip.label=="Lutra_maculicollis")] <- "Hydrictis_maculicollis"
phylogeny$tip.label[which(phylogeny$tip.label=="Galagoides_demidovii")] <- "Galagoides_demidoff"
phylogeny$tip.label[which(phylogeny$tip.label=="Cephalophus_nigrifrons")] <- "Cephalophorus_nigrifrons"
phylogeny$tip.label[which(phylogeny$tip.label=="Mustela_felipei")] <- "Neogale_felipei"
phylogeny$tip.label[which(phylogeny$tip.label=="Cryptotis_medellinia")] <- "Cryptotis_medellinius"
phylogeny$tip.label[which(phylogeny$tip.label=="Leopardus_colocolo")] <- "Leopardus_pajeros"
phylogeny$tip.label[which(phylogeny$tip.label=="Leontocebus_weddelli")] <- "Saguinus_weddelli"
phylogeny$tip.label[which(phylogeny$tip.label=="Callicebus_urubambensis")] <- "Plecturocebus_urubambensis"
phylogeny$tip.label[which(phylogeny$tip.label=="Parascaptor_leucura")] <- "Parascaptor_leucurus"
phylogeny$tip.label[which(phylogeny$tip.label=="Callicebus_regulus")] <- "Cheracebus_regulus"
phylogeny$tip.label[which(phylogeny$tip.label=="Herpestes_brachyurus")] <- "Urva_brachyura"
phylogeny$tip.label[which(phylogeny$tip.label=="Ictonyx_libyca")] <- "Poecilictis_libyca"
phylogeny$tip.label[which(phylogeny$tip.label=="Chimarrogale_platycephalus")] <- "Chimarrogale_platycephala"

# Remove species with lacking phylogenetic data
remove.species<-c("Petrogale wilkinsi", "Chodsigoa dabieshanensis", "Cheirogaleus andysabini", "Crocidura caudipilosa", "Cryptotis evaristoi", 
                  "Sorex obscurus", "Episoriculus umbrinus", "Cyclopes thomasi", "Cephalophorus brookei", "Lutra cinerea", "Ailurus styani", 
                  "Episoriculus soluensis", "Cercopithecus wolfi", "Ochotona vizier", "Philander nigratus", "Solenodon arredondoi", 
                  "Sorex madrensis", "Xanthonycticebus intermedius", "Rhynchocyon stuhlmanni", "Nycticebus hilleri", "Euroscaptor malayanus", 
                  "Cebuella niveiventris", "Sylvilagus daulensis", "Tupaia ferruginea", "Marmosa germana", "Madoqua cavendishi", 
                  "Saguinus leucogenys", "Sylvilagus tapetillus")

data <- data[!data$Binomial %in% remove.species,]

# Add missing species
# Add the European Mouflon with divergence time 378,000 years
# https://doi.org/10.1093/jhered/esi127
node <- which(phylogeny$tip.label=="Ovis_aries")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Ovis_musimon", 
                      where = node,
                      edge.length = 0.378,
                      position = 0.378)

# Add the Black Langur with divergence time 330,000 years
# https://doi.org/10.24272/j.issn.2095-8137.2020.254
node <- which(phylogeny$tip.label=="Trachypithecus_hatinhensis")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Trachypithecus_ebenus", 
                      where = node,
                      edge.length = 0.33,
                      position = 0.33)

# Add the Vietnamese Water Shrew with divergence time 2,310,000 years
# https://doi.org/10.1371/journal.pone.0077156
node <- which(phylogeny$tip.label=="Chimarrogale_himalayica")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Chimarrogale_varennei",
                      where = node,
                      edge.length = 2.31,
                      position = 2.31)

# Add the Groves's Titi with divergence time 1,300,000 years
# https://doi.org/10.1016/j.ympev.2018.11.012
node <- which(phylogeny$tip.label=="Callicebus_moloch")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Plecturocebus_grovesi",
                      where = node,
                      edge.length = 1.3,
                      position = 1.3)

# Add the Mitered Langur with divergence time 2,600,000 years
# https://doi.org/10.1038/s41598-020-66007-8
node <- which(phylogeny$tip.label=="Presbytis_femoralis")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Presbytis_mitrata",
                      where = node,
                      edge.length = 3.3,
                      position = 3.3)

# Add the East Sumatran Banded Langur with divergence time 2,600,000 years
# https://doi.org/10.1038/s41598-020-66007-8
node <- which(phylogeny$tip.label=="Presbytis_femoralis")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Presbytis_percura",
                      where = node,
                      edge.length = 2.6,
                      position = 2.6)

# Add the Cayman Nesophontes with divergence time 10,000,000 years
# https://doi.org/10.1093/molbev/msaa137
node <- which(phylogeny$tip.label=="Nesophontes_major")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Nesophontes_hemicingulus",
                      where = node,
                      edge.length = 10,
                      position = 10)

# Add the Buton Macaque with divergence time 1,000,000 years
# https://doi.org/10.7554/eLife.78169
node <- which(phylogeny$tip.label=="Macaca_ochreata")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Macaca_brunnescens",
                      where = node,
                      edge.length = 1,
                      position = 1)

# Add the Bryde's Whale with divergence time 6,300,000 years
# https://doi.org/10.1016/j.ympev.2006.03.032
node <- which(phylogeny$tip.label=="Balaenoptera_edeni")
node <- getMRCA(phylogeny, tip = c("Balaenoptera_edeni", "Balaenoptera_physalus"))
phylogeny <- bind.tip(phylogeny,
                      tip.label="Balaenoptera_brydei",
                      where = node,
                      edge.length = 6.3,
                      position = 6.3)

# Add the Indochinese Gray Langur with divergence time 360,000 years
# https://doi.org/10.1007/s10764-017-0008-4
node <- which(phylogeny$tip.label=="Trachypithecus_barbei")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Trachypithecus_crepusculus",
                      where = node,
                      edge.length = 0.36,
                      position = 0.36)

# Add the Makassar Tarsier with divergence time 350,000 years
# https://doi.org/10.1098/rsbl.2021.0642
node <- which(phylogeny$tip.label=="Tarsius_dentatus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Tarsius_fuscus",
                      where = node,
                      edge.length = 0.35,
                      position = 0.35)

# Add the East African Potto with divergence time 8,250,000 years
# and the Central African Potto with divergence time 5,500,000 years from the East African Potto
# https://doi.org/10.1111/zoj.12286
node <- which(phylogeny$tip.label=="Perodicticus_potto")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Perodicticus_ibeanus",
                      where = node,
                      edge.length = 8.25,
                      position = 8.25)
node <- which(phylogeny$tip.label=="Perodicticus_ibeanus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Perodicticus_edwardsi",
                      where = node,
                      edge.length = 5.5,
                      position = 5.5)

# Add the Urial with divergence time 1,260,000 years
# https://doi.org/10.1016/j.ympev.2009.10.037
node <- which(phylogeny$tip.label=="Ovis_orientalis")
node <- getMRCA(phylogeny, tip = c("Ovis_orientalis", "Ovis_aries", "Ovis_musimon"))
phylogeny <- bind.tip(phylogeny,
                      tip.label="Ovis_vignei",
                      where = node,
                      edge.length = 1.26,
                      position = 1.26)

# Add the Red Bald Uacari with divergence time 460,000 years
# https://doi.org/10.1016/j.ympev.2022.107509
node <- which(phylogeny$tip.label=="Cacajao_calvus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Cacajao_rubicundus",
                      where = node,
                      edge.length = 0.46,
                      position = 0.46)

# Add the Dabie Mountains Shrew Mole with divergence time 2,410,000 years
# https://doi.org/10.24272/j.issn.2095-8137.2020.266
node <- which(phylogeny$tip.label=="Uropsilus_gracilis")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Uropsilus_dabieshanensis",
                      where = node,
                      edge.length = 2.41,
                      position = 2.41)

# Add the Domestic Ass with divergence time 400,000 years
# https://doi.org/10.1073/pnas.1412627111
node <- which(phylogeny$tip.label=="Equus_africanus")
phylogeny <- bind.tip(phylogeny,
                      tip.label="Equus_asinus",
                      where = node,
                      edge.length = 0.4,
                      position = 0.4)

data$Binomial <- gsub(" ", "_", data$Binomial)


#Fit the model for normal sperm
norm_sperm_species <- brm(
  normal_sperm ~ Publish_Date + (1+Publish_Date|Binomial), 
  data = data, 
  init=0, chains = 3, iter = 3000, warmup = 1000, 
  family = brmsfamily("Beta"), 
  set_prior("normal(0, 0.0001)", class = "b"))

species_ranef <- as.data.frame(ranef(norm_sperm_species))

# Add entries of "NA" for those not included in the dataset, but is in the phylogeny
missing_ranef <- data.frame(matrix(NA, ncol = 8, nrow = length(phylogeny$tip.label[!(phylogeny$tip.label %in% rownames(species_ranef))])),
                            row.names = phylogeny$tip.label[!(phylogeny$tip.label %in% rownames(species_ranef))])
names(missing_ranef) <- colnames(species_ranef)
species_full_ranef <- rbind(species_ranef, missing_ranef)

# Create vector of slopes
species_slopes <- setNames(species_full_ranef[,"Binomial.Estimate.Publish_Date"], rownames(species_full_ranef))

slopes_available <- species_slopes[!is.na(species_slopes)]

# Plot slopes on tree (anc.ML method allows for NA entries)
slope_tree <- contMap(phylogeny, slopes_available, method="anc.ML", plot=FALSE)

slope_tree <- setMap(slope_tree, viridis(n=10))

tips.na <- sapply(names(which(is.na(species_slopes))), function(x,y) which(y==x), y=slope_tree$tree$tip.label)
ancestors.na <- Ancestors(phylogeny, tips.na)
ancestors.available <- Ancestors(phylogeny, setdiff(1:Ntip(phylogeny), tips.na))
slopes.na <- setdiff(c(tips.na, unique(unlist(ancestors.na))),
            c(setdiff(1:Ntip(phylogeny),tips.na),
              unique(unlist(ancestors.available))))

slope_tree$tree <- paintBranches(slope_tree$tree, slopes.na, "NA")

cols <- setNames(c("grey", slope_tree$cols), c("NA",names(slope_tree$cols)))
plot(slope_tree$tree, cols, split.vertical=TRUE, outline=TRUE, lwd=6,
     ftype="i",fsize=0.7)
add.color.bar(40, slope_tree$cols, title="trait", lims=range(species_slopes ,na.rm=TRUE),
              digits=2, subtitle="", prompt=FALSE)
