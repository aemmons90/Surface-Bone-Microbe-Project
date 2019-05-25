

###############################SIMPER
#use rarefied data
otus<-otu_table(Bones_scale1)
#otu<-t(otus)
A<-otu_table(A_NMDS)
B<-otu_table(B_NMDS)
C<-otu_table(C_NMDS)

#use simper scripts; load into R (from https://github.com/asteinberger9/seq_scripts)
source("G:/My Drive/BoneSurfaceProject/16S_Data/Phyloseq2018/simper_pretty.r")
source("G:/My Drive/BoneSurfaceProject/16S_Data/Phyloseq2018/R_krusk.r")

#############################Look at SIMPER by body region in each individual
#same code was used for each individual with object and designators changed each time
simper.pretty(C, sample_data(C_NMDS), c('Body.region'), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, 'C_BR')


simper.results = data.frame(read.csv("C_BR_clean_simper.csv"))


kruskal.pretty(C, sample_data(C_NMDS), simper.results, c('Body.region'), 'C_BR')


krusk<-read.table("region_krusk_simper1.csv", header=TRUE, sep=",")
krusk1<-krusk%>% filter(fdr_krusk_p.val<=0.05)

#get taxonomic info to match to the table
sigtaxa<- c(as.character(krusk1$OTU))


simpdat<-Bones_scale1 %>% subset_taxa(rownames(tax_table(Bones_scale1)) %in% sigtaxa)

a<-as.data.frame(tax_table(simpdat))

a$OTU<-rownames(a)

krusk2<-left_join(krusk1, a, by="OTU" )


#want to see how all members of these genera, not just these OTUs, compare
a$Genus<-as.character(a$Genus)
sigtaxa1<-c(a$Genus)
simptax<-Bones_scale1 %>% subset_taxa(Genus %in% sigtaxa1) %>%
  tax_glom(taxrank = "Genus")


n <- 18
palette <- distinctColorPalette(n)


plot_bar(simptax, "BONE.TYPE", fill="Genus")+ facet_grid(Individual~Body.region, space="free", scales="free")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1))+
  scale_fill_manual(values = palette)


###################taxa differentiating individuals
simper.pretty(otus, sample_data(Bones_scale1), c('Individual'), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, 'IndSimper')

simper.results2 = data.frame(read.csv("IndSimper_clean_simper.csv"))



kruskal.pretty(otus, sample_data(Bones_scale1), simper.results2, c('Individual'), 'Indv', tax)

#keep only significant taxa
kruskind<-read.table("Indv_krusk_simper.csv", header=TRUE, sep=",")
kruskind1<-kruskind%>% filter(fdr_krusk_p.val<=0.05)


#get taxonomic info to match to the table
sigtaxa1<- c(as.character(kruskind1$OTU))

#merge data sets so that i can see which taxa I'm looking at
krusk2<-merge(kruskind1, as.data.frame(tax_table(simpdat1)))


simpdat1<-Bones_scale1 %>% subset_taxa(rownames(tax_table(Bones_scale1)) %in% sigtaxa1)

plot_bar(simpdat1, "Family", fill="Genus", facet_grid=~Individual)+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1))



#plot relative abundance
Simtrans <- Bones_scale1  %>%  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() 


#Average relative abundance by bodysite. 
Indsimper <- Simtrans %>%
  dplyr::arrange(Individual, OTU) %>%
  dplyr::group_by(Individual, Phylum,Class, Order, Family, Genus, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::filter(OTU %in% sigtaxa1 & avg.rel.abund > 0)

Indsimper$Genus<-as.character(as.factor(Indsimper$Genus))
Indsimper$Family<-as.character(as.factor(Indsimper$Family))
Indsimper$Genus[is.na(Indsimper$Genus)] <- Indsimper$Family[is.na(Indsimper$Genus)]



library(randomcoloR)
n <- 19
palette <- distinctColorPalette(n)

# Plot #need to fix grid facet headers so that they fit
ind<-ggplot(Indsimper, aes(x = Family, y = avg.rel.abund, fill = Genus)) + 
  facet_grid(Individual~., scales="free", space="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_bw()+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Average Relative Abundance") +
  # ggtitle("Class Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=8)) 

#ggsave("figS6.png", width = 5, height = 4, dpi=300)

#look at deinococcus thermus
ftaxa<-subset_taxa(Bones_scale1, Phylum == "Deinococcus-Thermus") 

#look at derma
dtaxa<-subset_taxa(Bones_scale1, Family == "Dermacoccaceae") 
#plot

plot_bar(dtaxa, "Family", fill="Genus", facet_grid=Individual~Body.region)+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1))
##################################################Create dataframes to be used in python for random forests
##look at genus level for python
#example code

scale1<-tax_glom(Bones_scale1, "Genus") 

scalenoteeth<-subset_samples(scale1, Body.region != "tooth")

newOTU<-as.data.frame(t(otu_table(scalenoteeth)))


newOTU$HDNA <- sample_data(scalenoteeth)$Human.DNA 
newOTU$HDNA <- log10((newOTU$HDNA+1))
write.csv(newOTU, "scalenoteeth.csv", sep = ",", col.names=TRUE, row.names=TRUE)

###############################
#look at important features to log transformed human DNA 
#read in important features resulting from python script
imp<-read.csv("impOTU16noteeth.csv")
imp<-imp[,-1]

imp<-imp%>% filter(X1>0)

colnames(imp)<-c("OTU", "Importance")
imp$OTU<- as.character(imp$OTU)
IOTU<-c(imp$OTU)


#subset physeq object by important taxa
my_subset <- subset(otu_table(Bones_scale1), rownames(otu_table(Bones_scale1)) %in% IOTU)
new_physeq <- merge_phyloseq(my_subset, tax_table(Bones_scale1), sample_data(Bones_scale1))

#merge taxonomic info with importance values
row.names(imp)<-imp$OTU
df<-as.data.frame(tax_table(new_physeq))
df$OTU<- row.names(df)
imp2<-inner_join(imp, df, by="OTU")

#plot taxa sums by human DNA
rfimp <-Bones_scale1 %>% 
  tax_glom(taxrank = "Genus") %>%  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02)  %>%                       # Filter out low abundance taxa
  dplyr::arrange(Individual,Body.region, BONE.TYPE, Sample.ID) %>%
  dplyr::group_by(Individual,Body.region, BONE.TYPE, Genus) %>%
  dplyr::select(Human.DNA, Abundance) 

rfimp2<- rfimp%>%
  filter(Genus=="Clostridium_sensu_stricto_15" |Genus=="Paracoccus" | Genus=="Actinotalea" | Genus=="Dermacoccaceae_unclassified") 

#drop non-important levels
rfimp2<-droplevels(rfimp2)

#make an important column
rfimp2$importance<-rfimp2$Genus 
levels(rfimp2$importance)<-c("0.04", "0.23", "0.07", "0.04" )
levels(rfimp2$Genus)<-c("Actinotalea", "Clostridium", "Dermacoccaceae", "Paracoccus" )
#fix BONE.TYPE 
rfimp2$BONE.TYPE <-as.character(rfimp2$BONE.TYPE)
rfimp2$Body.region <-as.character(rfimp2$Body.region)


rfimp2<- rfimp2 %>% dplyr::arrange(Genus) %>% filter(Abundance >0)

library(ggrepel)
a<-ggplot(rfimp2,aes(x=Human.DNA, y=Abundance, group=Genus))+
  geom_point(aes(color=importance), size=1)+
  facet_grid(Genus~.)+
  #geom_line(aes(color=Body.region))+
  xlab("Human DNA (ng/gbp)")+
  theme_classic()+
  scale_color_viridis(discrete=TRUE)+
  geom_text_repel(aes(label = ifelse(Human.DNA > 400,Body.region, "")), size = 2, segment.size = 0.1)+
  theme(legend.title= element_blank())+
  #geom_hline(yintercept=0, color="red", linetype="dotted")+
  theme(legend.position = "bottom")+
  ggtitle("A")

#ggsave("rf_Surface.png", width = 3.0, height = 5.5, dpi=300)

#plot outcome of random forests using OTUs at the family level
rfmod<-read.csv("test16noteeth.csv")
#file with teeth = testscale.csv
b<-ggplot(rfmod, aes(x=x, y=y))+
  geom_point()+
  theme_bw()+
  geom_smooth(method='lm',formula=y~x, color="red")+
  xlab("Test values")+
  ylab("Predicted values")+
  geom_abline(slope=1, linetype="dotted", color="red")+
  ggtitle("B")

rf<-lm(y~x, data=rfmod)

######################################################18S
#looking at random forest and human DNA quantity/quality
#in this case Bones_scale1 is from the eukaryotic file 
scale18<-tax_glom(Bones_scale1, "Genus") 

scalenoteeth18<-subset_samples(scale18, Body.region != "tooth")

newOTU18<-as.data.frame((otu_table(scalenoteeth18)))
newOTU18<-as.data.frame(t(otu_table(scalenoteeth18)))


newOTU18$HDNA <- sample_data(scalenoteeth18)$Human.DNA..Qf. 

#newOTU$HDNA <- log10((newOTU$HDNA+1))
write.csv(newOTU, "scale18.csv", sep = ",", col.names=TRUE, row.names=TRUE)


#################################Combine eukaryote and bacteria datasets
#merge 16S and 18S datasets 
#get the data frame to work with
newOTU18<-as.data.frame((otu_table(scale18)))

#put an E before the O in "OTU" in the rownames
row.names(newOTU18)<-gsub("O","EO", row.names(newOTU18))

#transpose data frame
newOTU18<-as.data.frame(t(newOTU18))

#add human DNA
newOTU18$HDNA <- sample_data(scale18)$Human.DNA..Qf. 

#merge 16S and 18S data frames using dplyr
newOTU18$sample<-as.character(row.names(newOTU18))
newOTU$sample <-as.character(row.names(newOTU))
combinedrf<-inner_join(newOTU18, newOTU, by="sample")

#remove weird columns
combinedrf<-subset(combinedrf, select=-c(HDNA.x, sample))

#rename HDNA.y column
names(combinedrf)[names(combinedrf) == 'HDNA.y'] <- 'HDNA'

#write to csv file
write.csv(combinedrf, "combinedrf.csv", col.names=TRUE, row.names=TRUE)



#get predictions vs test values
rfmod<-read.csv("testcombinedteeth.csv")

ggplot(rfmod, aes(x=x, y=y))+
  geom_point()+
  theme_bw()+
  geom_smooth(method='lm',formula=y~x, color="red")+
  xlab("Test values")+
  ylab("Predicted values")+
  geom_abline(slope=1, linetype="dotted", color="red")
# scale_y_continuous(breaks = c(0,50,100,150,200,250,300)) +
#scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950)) 

rf<-lm(y~x, data=rfmod)

#extract important OTUs
imp<-read.csv("impOTUcombinedteeth.csv")
imp<-imp[,-1]

imp<-imp%>% filter(X1>0)

colnames(imp)<-c("OTU", "Importance")
imp$OTU<- as.character(imp$OTU)
IOTU<-c(imp$OTU)


#subset physeq object by important taxa
my_subset <- subset(otu_table(scale1), rownames(otu_table(scale1)) %in% IOTU)
new_physeq16 <- merge_phyloseq(my_subset, tax_table(scale1), sample_data(scale1))



#remove E for eukaryotic set
impE<-imp[c(1,10,11),]
impE$OTU<-gsub("EO","O", impE$OTU)
IOTU<-c(impE$OTU)

my_subset <- subset(otu_table(Bones_scale1), rownames(otu_table(Bones_scale1)) %in% IOTU)
new_physeq18 <- merge_phyloseq(my_subset, tax_table(Bones_scale1),
                               sample_data(Bones_scale1))

rfimp <-Bones_scale1 %>% 
  tax_glom(taxrank = "Genus") %>%  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02)  %>%                       # Filter out low abundance taxa
  dplyr::arrange(Individual,Body.region, BONE.TYPE, Sample.ID) %>%
  dplyr::group_by(Individual,Body.region, BONE.TYPE, Genus) %>%
  dplyr::select(Human.DNA..Qf., Abundance) %>%
  filter(Genus == "Saccharomycetales_unclassified")

rfimp$Body.region<-as.character(rfimp$Body.region)
rfimp$Individual<-as.character(rfimp$Individual)

ggplot(rfimp,aes(x=Human.DNA..Qf., y=Abundance, group=Genus))+
  geom_point(aes(color=Body.region), size=3)+
  #facet_grid(Individual~.)+
  xlab("Human DNA (ng/gbp)")+
  ylab("Relative Abundance")+
  theme_classic()+
  scale_color_viridis(discrete=TRUE)+
  geom_text_repel(aes(label = ifelse(Body.region=="skull"&Individual=="B",Individual, "")), size = 2, segment.size = 0.1)+
  theme(legend.title= element_blank())+
  #geom_hline(yintercept=0, color="red", linetype="dotted")+
  theme(legend.position = "bottom")


imp<-read.csv("impOTUscale18.csv")
imp<-imp[,-1]

imp<-imp%>% filter(X1>0)

colnames(imp)<-c("OTU", "Importance")
imp$OTU<- as.character(imp$OTU)
IOTU<-c(imp$OTU)


#subset physeq object by important taxa
my_subset <- subset(otu_table(Bones_scale1), rownames(otu_table(Bones_scale1)) %in% IOTU)
new_physeq <- merge_phyloseq(my_subset, tax_table(Bones_scale1), sample_data(Bones_scale1))

#merge taxonomic info with importance values
row.names(imp)<-imp$OTU
df<-as.data.frame(tax_table(new_physeq))
df$OTU<- row.names(df)
imp2<-inner_join(imp, df, by="OTU")


#plot outcome of random forests using OTUs at the genus level
rfmod<-read.csv("scaletest18.csv")

b<-ggplot(rfmod, aes(x=x, y=y))+
  geom_point()+
  theme_bw()+
  geom_smooth(method='lm',formula=y~x, color="red")+
  xlab("Test values")+
  ylab("Predicted values")+
  geom_abline(slope=1, linetype="dotted", color="red")+
  ggtitle("B")

rf<-lm(y~x, data=rfmod)


