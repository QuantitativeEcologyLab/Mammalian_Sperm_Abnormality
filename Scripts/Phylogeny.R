# This script assembles the mammalian phylogeny used to support the analyses
# presented in the main text. The base phylogeny was obtained from the VertLife
# repository. Additional information on missing species and subspecies
# Are then inserted manually with divergence times taken from the peer-reviewed
# literature when available.
setwd("/Users/nozo/Library/CloudStorage/OneDrive-UBC/Quantitative_Ecology_Lab/Nozomu_Hirama_Undergraduate_Honours/Mammalian_Sperm_Abnormality/")
set.seed(88210729)


#Load in the necessary packages
library(ape)
library(phytools)


#Phylogenetic Data
#Create list of all species
unique_species <- data %>%
  distinct(Binomial)
write.csv(unique_species, file="./Data/Unique_Species_List.csv", row.names=FALSE)

#Add other species needed for node references
unique_species <- rbind(unique_species, "Ovis aries")
unique_species <- rbind(unique_species, "Mustela felipei")
unique_species <- rbind(unique_species, "Trachypithecus hatinhensis")
unique_species <- rbind(unique_species, "Callicebus caligatus")
unique_species <- rbind(unique_species, "Macropus eugenii")
unique_species <- rbind(unique_species, "Vicugna vicugna")
unique_species <- rbind(unique_species, "Herpestes fuscus")
unique_species <- rbind(unique_species, "Mazama pandora")
unique_species <- rbind(unique_species, "Chimarrogale himalayica")
unique_species <- rbind(unique_species, "Mustela africana")
unique_species <- rbind(unique_species, "Lutra maculicollis")
unique_species <- rbind(unique_species, "Galagoides demidovii")
unique_species <- rbind(unique_species, "Cephalophus nigrifrons")
unique_species <- rbind(unique_species, "Callicebus moloch")
unique_species <- rbind(unique_species, "Cryptotis medellinia")
unique_species <- rbind(unique_species, "Presbytis femoralis")
unique_species <- rbind(unique_species, "Leopardus colocolo")
unique_species <- rbind(unique_species, "Leontocebus weddelli")
unique_species <- rbind(unique_species, "Callicebus urubambensis")
unique_species <- rbind(unique_species, "Nesophontes major")
unique_species <- rbind(unique_species, "Perodicticus potto")
unique_species <- rbind(unique_species, "Macaca ochreata")
unique_species <- rbind(unique_species, "Parascaptor leucura")
unique_species <- rbind(unique_species, "Callicebus regulus")
unique_species <- rbind(unique_species, "Herpestes brachyurus")
unique_species <- rbind(unique_species, "Balaenoptera edeni")
unique_species <- rbind(unique_species, "Trachypithecus barbei")
unique_species <- rbind(unique_species, "Tarsius dentatus")
unique_species <- rbind(unique_species, "Ovis orientalis")
unique_species <- rbind(unique_species, "Cacajao calvus")
unique_species <- rbind(unique_species, "Uropsilus gracilis")
unique_species <- rbind(unique_species, "Ictonyx libyca")
unique_species <- rbind(unique_species, "Equus africanus")
unique_species <- rbind(unique_species, "Chimarrogale platycephalus")


#Update species list
write.csv(unique_species, file="./Data/Unique_Species_List.csv", row.names=FALSE)

#----------------------------------------------------------------------
# Import the phylogeny
#----------------------------------------------------------------------

#Base phylogeny
phylo.trees <- read.nexus("Data/Phylogeny/output.nex")

#Estimate the consensus tree
phylogeny <- ls.consensus(phylo.trees)

unique_species$Binomial <- gsub(" ", "_", unique_species$Binomial)

missing_species <- subset(unique_species, !(unique_species$Binomial %in% phylogeny[["tip.label"]]))

branch_lengths <- setNames(phylogeny$edge.length[sapply(1:length(phylogeny$tip.label),function(x,y) which (y==x),y=phylogeny$edge[,2])],phylogeny$tip.label)

#----------------------------------------------------------------------
#Rename a few species to match the subspecies used in the analyses
#----------------------------------------------------------------------

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


#----------------------------------------------------------------------
# Remove species with insufficient phylogenetic data
#----------------------------------------------------------------------

# The following species have been removed due to the lack of phylogenetic data/divergence times

# The Wilkins's Rock Wallaby (Petrogale wilkinsi)
# The Dabieshan Brown-toothed Shrew	(Chodsigoa dabieshanensis)
# The Montagne d'Ambre Dwarf Lemur (Cheirogaleus andysabini)
# The Hairy-tailed White-toothed Shrew (Crocidura caudipilosa)
# The Evaristo's Small-eared Shrew (Cryptotis evaristoi)
# The Northern Montane Shrew (Sorex obscurus)
# The Hidden Brown-toothed Shrew (Episoriculus umbrinus)
# The Thomas's Silky Anteater	(Cyclopes thomasi)
# The Brooke's Duiker	(Cephalophorus brookei)
# The Asian Small-clawed Otter (Lutra cinerea)
# The Eastern Red Panda	(Ailurus styani)
# The Nepalese Brown-toothed Shrew (Episoriculus soluensis)
# The Wolf's Monkey	(Cercopithecus wolfi)
# The Iranian Pika (Ochotona vizier)
# The Black Four-eyed Opossum	(Philander nigratus)
# The Giant Solenodon	(Solenodon arredondoi)
# The Sierra Shrew (Sorex madrensis)
# The Northern Pygmy Slow Loris	(Xanthonycticebus intermedius)
# The White-tailed Sengi (Rhynchocyon stuhlmanni)
# The Hiller's Slow Loris (Nycticebus hilleri)
# The Malaysian Mole (Euroscaptor malayanus)
# The Southern Pygmy Marmoset (Cebuella niveiventris)
# The Ecuadorean Tapeti	(Sylvilagus daulensis)
# The Sumatran Treeshrew (Tupaia ferruginea)
# The North-western Woolly Mouse Opossum (Marmosa germana)
# The Naivasha Dik-dik (Madoqua cavendishi)
# The Andean Saddle-back Tamarin (Saguinus leucogenys)
# The Coastal Tapeti (Sylvilagus tapetillus)

remove.species<-c("Petrogale wilkinsi", "Chodsigoa dabieshanensis", "Cheirogaleus andysabini", "Crocidura caudipilosa", "Cryptotis evaristoi", 
                  "Sorex obscurus", "Episoriculus umbrinus", "Cyclopes thomasi", "Cephalophorus brookei", "Lutra cinerea", "Ailurus styani", 
                  "Episoriculus soluensis", "Cercopithecus wolfi", "Ochotona vizier", "Philander nigratus", "Solenodon arredondoi", 
                  "Sorex madrensis", "Xanthonycticebus intermedius", "Rhynchocyon stuhlmanni", "Nycticebus hilleri", "Euroscaptor malayanus", 
                  "Cebuella niveiventris", "Sylvilagus daulensis", "Tupaia ferruginea", "Marmosa germana", "Madoqua cavendishi", 
                  "Saguinus leucogenys", "Sylvilagus tapetillus")

data <- data[!data$Binomial %in% remove.species,]

#----------------------------------------------------------------------
# Add in the missing species/sub-species
#----------------------------------------------------------------------

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


#----------------------------------------------------------------------
# Check final tree
#----------------------------------------------------------------------

unique_species <- data %>%
  distinct(Binomial)

unique_species$Binomial <- gsub(" ", "_", unique_species$Binomial)

missing_species <- subset(unique_species, !(unique_species$Binomial %in% phylogeny[["tip.label"]]))

branch_lengths <- setNames(phylogeny$edge.length[sapply(1:length(phylogeny$tip.label),function(x,y) which (y==x),y=phylogeny$edge[,2])],phylogeny$tip.label)

#----------------------------------------------------------------------
# Construct covariance matrix from phyolgeny
#----------------------------------------------------------------------

phylo.cov <- vcv.phylo(phylogeny)
