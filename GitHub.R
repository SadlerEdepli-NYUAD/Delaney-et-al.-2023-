#All codes below were used to generate figures
#Codes are appear in order of the figure

#Figure S1 ----
#Single liver correlation and as3mt expression
library(readr)
Table_S1_endogenous_single_livers <- read_csv("Table S1 endogenous single livers.csv", 
                                              col_types = cols(...1 = col_skip()))
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 16))

ggplot(Table_S1_endogenous_single_livers , aes(x=YM, y=PM)) + 
  geom_point(alpha=0.3)

highlight <- subset(Table_S1_endogenous_single_livers [4736,])


ggplot(Table_S1_endogenous_single_livers , aes(x=YM, y=PM)) + 
  geom_point(alpha=0.3) +
  geom_point(data=highlight, 
             aes(x=YM,y=PM), 
             color='red',
             size=3) + 
  xlab("Mean Read Count (Yelena)") + ylab("Mean Read Count (Patrice)")

#genes less than as3mt
lower = subset.data.frame(Table_S1_endogenous_single_livers , Table_S1_endogenous_single_livers $YM < 8.33) #total 32105 / lower 30802 / 0.96
higher = subset.data.frame(Table_S1_endogenous_single_livers , Table_S1_endogenous_single_livers $YM > 8.33)  #total 32105 / higher 1303  / 0.04






#Figure 2 Venn & Cross plot & GO -----
#Upload datasets
#Whole Larvae as3mt mut
WL <- read.csv("Table S2 as3mtmut liver.csv")
#Pooled Larval Livers as3mt mut
LL <- read.csv("Table S3 as3mtmut whole larvae.csv")

#Subset Sig
WLSig <- subset(WL, WL$padj < 0.05)
LLSig <- subset(LL, LL$padj < 0.05)

#Venn 
library(eulerr)
make_2way_eulerr = function(list1,list2,name,colours,scale){
  n1 = length(list1)
  n2 = length(list2)
  
  n12 = length(intersect(list1,list2))
  
  
  comb = c('1'=n1,'2'=n2,
           '1&2'=n12)
  venn = euler(comb, input = 'union', shape = 'ellipse')
  save_venn_plots(venn,name,colours,scale)
}
save_venn_plots = function(eulerr,name,colours,scale){
  labelname = paste(name,'_labels.png',sep='')
  nolabelname = paste(name,'.png',sep='')
  #to see labels to add manually in illustrator
  png(labelname)
  print({plot(eulerr,
              labels=FALSE,
              quantities = TRUE,
              fills = colours)
  })
  dev.off()
  
  png(nolabelname,width = 200*scale,height = 200*scale)
  print({plot(eulerr,edges = list(lwd=1*scale),
              labels=FALSE,
              quantities = FALSE,
              fills = colours)
  })
  dev.off()
}

Liver = '#C779D0' 
WholeLarvae =  '#FEAC5E' 

#Found colors using https://hexcolorspicker.com/google-color-picker/

As3mtMUT_venn = make_2way_eulerr(WLSig$gene.ids, LLSig$gene.ids, 'Figures/Venn',
                                 c(WholeLarvae, Liver),7)


#Crossplot
#Subset and rename columns (WL) X
WT_1mM_forMerge <- WL[, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(WT_1mM_forMerge)[names(WT_1mM_forMerge) == "log2FoldChange"] <- "log2FoldChange.x"
names(WT_1mM_forMerge)[names(WT_1mM_forMerge) == "padj"] <- "padj.x"
WT_1mM_Sig <- subset.data.frame(WL, WL$padj<0.05)

#Subset and rename columns (LL) Y
A_1mM_forMerge  <- LL [, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(A_1mM_forMerge)[names(A_1mM_forMerge) == "log2FoldChange"] <- "log2FoldChange.y"
names(A_1mM_forMerge)[names(A_1mM_forMerge) == "padj"] <- "padj.y"
As3mt_1mM_Sig <- subset.data.frame(LL, LL$padj<0.05)

#Mergedatasets
Mergedfile = merge(WT_1mM_forMerge,A_1mM_forMerge, by="gene.ids", all=T)

#Clean up merged (remove duplicates and NA) 
Mergefile_clean <- na.omit(Mergedfile) #removes NA
#install.packages("tidyverse")
library(tidyverse)
Mergefile_final <- Mergefile_clean[!duplicated(Mergefile_clean$gene.ids), ] #removes duplicates


Mergefile_final$Legend <- NA
#assign attributes for color based on significance
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "WL < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "LL < 0.05"
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "Both < 0.05"
Mergefile_final[Mergefile_final$log2FoldChange.x>0&Mergefile_final$log2FoldChange.y<0, ]$Legend <- "NEG Cor"
Mergefile_final[Mergefile_final$log2FoldChange.x<0&Mergefile_final$log2FoldChange.y>0, ]$Legend <- "NEG Cor"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "Neither"


ActsUP <- Mergefile_final[Mergefile_final$log2FoldChange.x < 0 & Mergefile_final$log2FoldChange.y < 0, ]    
ActsDOWN <- Mergefile_final[Mergefile_final$log2FoldChange.x > 0 & Mergefile_final$log2FoldChange.y > 0, ]
Acts_All_WL_LL_PosCor <- rbind(ActsUP,ActsDOWN)

#Remove neither
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Neither"), ]
#Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Both < 0.05"), ]#Add a extra column to the dataset

library(ggplot2)
#graph
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=1, shape=19)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Larval Liver", limits=c(-13, 13)) + #play with margin
  scale_x_continuous(name="Whole Larvae", limits=c(-13, 13)) + #play with margin
  scale_color_manual(values=c('#4BC0C8','#C779D0', 'grey', '#FEAC5E')) +
  geom_smooth(method='lm', color="magenta", size=0.5, formula= y~x) +
  geom_point(data = Acts_All_WL_LL_PosCor %>% filter(gene.ids == "ENSDARG00000027572"), color = "cyan") #as3mt
all_plot

r <- cor(Acts_All_WL_LL_PosCor$log2FoldChange.x, Acts_All_WL_LL_PosCor$log2FoldChange.y, method = "pearson")
r #0.640

#Go on all Positively Correlated only
#subset all POS_COR UP OR DOWN
Mergefile_final_1 <-subset(Mergefile_final, Mergefile_final$log2FoldChange.x>0 & Mergefile_final$log2FoldChange.y>0)
#write_csv(Mergefile_final_1, file="PosCor_LL_WL_UP_1457")
Mergefile_final_2 <-subset(Mergefile_final, Mergefile_final$log2FoldChange.x<0 & Mergefile_final$log2FoldChange.y<0)
#write_csv(Mergefile_final_2, file="PosCor_LL_WL_DOWN_1823")
Mergefile_final_all <- rbind(Mergefile_final_1,Mergefile_final_2)

#Remove NEGCOR going in different directions

ActsUP <- Mergefile_final[Mergefile_final$log2FoldChange.x < 0 & Mergefile_final$log2FoldChange.y < 0, ]    
ActsDOWN <- Mergefile_final[Mergefile_final$log2FoldChange.x > 0 & Mergefile_final$log2FoldChange.y > 0, ]
Acts_All_WL_LL_PosCor <- rbind(ActsUP,ActsDOWN)
LL_Acts <- subset(LL, LL$gene.ids %in% Acts_All_WL_LL_PosCor$gene.ids)
#write.csv(LL_Acts_3280, file="LL_Acts_608.csv")


Mergefile_final <- Mergefile_final[, -c(6)]
Mergefile_final$Legend <- NA

#assign attributes for color based on significance
Mergefile_final[Mergefile_final$log2FoldChange.x<0&Mergefile_final$log2FoldChange.y<0, ]$Legend <- "Pos Cor"
Mergefile_final[Mergefile_final$log2FoldChange.x>0&Mergefile_final$log2FoldChange.y>0, ]$Legend <- "Pos Cor"
Mergefile_final[Mergefile_final$log2FoldChange.x<0&Mergefile_final$log2FoldChange.y>0, ]$Legend <- "Neg Cor"                
Mergefile_final[Mergefile_final$log2FoldChange.x>0&Mergefile_final$log2FoldChange.y<0, ]$Legend <- "Neg Cor"  

library(ggrepel)

#Crossplot
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=1, shape=19)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Whole Larvae DEGs", limits=c(-13, 13)) + #play with margin
  scale_x_continuous(name="Larval Liver DEGs", limits=c(-13, 13)) + #play with margin
  scale_color_manual(values=c('grey','#4BC0C8')) +
  geom_smooth(method='lm', color="magenta", size=0.5, formula= y~x) +
  geom_point(data = Acts_All_WL_LL_PosCor %>% filter(gene.ids == "ENSDARG00000027572"), color = "magenta")
all_plot

#Now Revigo on 608
library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Dr.eg.db

# Curate GO terms
UP_ACT <- subset(Acts_All_WL_LL_PosCor, Acts_All_WL_LL_PosCor$log2FoldChange.x > 0 & Acts_All_WL_LL_PosCor$log2FoldChange.y > 0) 
DOWN_ACT <- subset(Acts_All_WL_LL_PosCor, Acts_All_WL_LL_PosCor$log2FoldChange.x < 0 & Acts_All_WL_LL_PosCor$log2FoldChange.y < 0) 

#Plotting poscor, DEGs in one, other or both
fig = list(UP=UP_ACT$gene.ids, DOWN=DOWN_ACT$gene.ids)

#Functions
for_plotting_comparecluster = function(list,finalgo,no_top_IDs,text_wrapsize){
  finaldb = as.data.frame(finalgo)
  topids = vector()
  names = names(list)
  
  #extract top ids for each group
  for (i in 1:length(names)){
    name = names[i]
    top = finaldb[finaldb$Cluster==name,]
    top = top[order(top$p.adjust),]
    top = top$ID[1:no_top_IDs]
    topids = c(topids,top)
  }
  
  dropids = unique(finaldb$ID)
  dropids = dropids[!dropids %in% topids]
  #subset GO to only include top IDs
  plotdb = dropGO(finalgo,term=dropids)
  
  library(dplyr)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  
  plot_table = as.data.frame(plotdb)
  plot_table$GeneRatio = as.numeric(sub('/.*','',plot_table$GeneRatio))/as.numeric(sub('.*/','',plot_table$GeneRatio))
  plot_table$BgRatio = as.numeric(sub('/.*','',plot_table$BgRatio))/as.numeric(sub('.*/','',plot_table$BgRatio))
  plot_table$Description = str_wrap(plot_table$Description, width = text_wrapsize)
  plot_table$p.adjust = as.numeric(plot_table$p.adjust)
  plot_table$Description = factor(plot_table$Description,levels=rev(unique(plot_table$Description)))
  
  return(plot_table)
}

# create table for uploading on REVIGO
create_revigo_table = function(go, filename){
  db = as.data.frame(go)
  db = db[order(db$p.adjust),]
  keepgos = data.frame(ID=vector(),p.adjust=vector())
  #create list of unique GO IDs with lowest p-value of all groups
  for (i in 1:nrow(db)){
    if (!db$ID[i] %in% keepgos$ID){
      keepgos = rbind(keepgos,db[i,])
    }
  }
  write.table(keepgos[,c(2,6)],quote = F,sep = ' ',row.names = F,col.names=F,file = filename)
}

drop_from_revigo = function(current_go,revigo_filename){
  revi = read.csv("Revigo_BP_Table.csv")
  db = as.data.frame(current_go)
  ids = db$ID[!db$ID %in% revi$TermID[revi$Representative==-1]]
  dropped_go = dropGO(current_go,term = ids)
  return(dropped_go)
}

fig_go = compareCluster(fig, fun = "enrichGO", keyType = "ENSEMBL", OrgDb = "org.Dr.eg.db" , ont = "BP",readable = TRUE) #change MF to CC or BP depending on what you want
create_revigo_table(fig_go, "fig_goidswithp.txt")

# Copy and paste this file to REVIGO (http://revigo.irb.hr) (add Danio rerio species) and retrieve summary of GO terms as csv file and save it in the working directory as REVIGO_fig.csv
revi = read.delim("Revigo_BP_Table_FORPAPER.csv")
revi=as.data.frame(revi)
db = as.data.frame(fig_go)

ids = db$ID[!db$ID %in% revi$TermID[revi$Dispensability<=0.5]]
dropped_go = dropGO(fig_go,term = ids)

scale = 1.5
wrapsize = 42

#Plot terms - 8 terms
fig_plottable = for_plotting_comparecluster(fig,dropped_go,10,wrapsize) # change number depending on how many top terms you would like to be shown
write.table(fig_plottable, file = "GO_BP_10terms.txt", quote = F, sep = "\t ", row.names = F)

plot_fig = ggplot(fig_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="magenta", high="cyan", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c("Down"))+
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size = 22,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18),
        text = element_text(size=15),
        axis.title.x =element_text(size=22))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 18*scale),breaks=c(0.02,0.04,0.08,0.14))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.002,fontface='bold'),colour='white',show.legend = F)

plot_fig

postscript("3280_BP.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 11, height = 13)
plot_fig
dev.off()

#Code for Figure 2F & G after Figure S4

#Repeated and save for BP and CC for figure 2 to get both plots


#Figure S3 Cross plot & GO ----
WL <- read.csv("Table S3 as3mtmut whole larvae.csv")
LL <- read.csv("Table S2 as3mtmut liver.csv")

#Subset and rename columns (WL) X
WT_1mM_forMerge <- WL[, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(WT_1mM_forMerge)[names(WT_1mM_forMerge) == "log2FoldChange"] <- "log2FoldChange.x"
names(WT_1mM_forMerge)[names(WT_1mM_forMerge) == "padj"] <- "padj.x"
WT_1mM_Sig <- subset.data.frame(WL, WL$padj<0.05)

#Subset and rename columns (LL) Y
A_1mM_forMerge  <- LL [, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(A_1mM_forMerge)[names(A_1mM_forMerge) == "log2FoldChange"] <- "log2FoldChange.y"
names(A_1mM_forMerge)[names(A_1mM_forMerge) == "padj"] <- "padj.y"
As3mt_1mM_Sig <- subset.data.frame(LL, LL$padj<0.05)

#Mergedatasets
Mergedfile = merge(WT_1mM_forMerge,A_1mM_forMerge, by="gene.ids", all=T)

#Clean up merged (remove duplicates and NA) 
Mergefile_clean <- na.omit(Mergedfile) #removes NA
#install.packages("tidyverse")
library(tidyverse)
Mergefile_final <- Mergefile_clean[!duplicated(Mergefile_clean$gene.ids), ] #removes duplicates


Mergefile_final$Legend <- NA
#assign attributes for color based on significance
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "WL < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "LL < 0.05"
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "Both < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "Neither"

#Remove neither
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Neither"), ]
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "WL < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "LL < 0.05"
#Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Both < 0.05"), ]#Add a extra column to the dataset

library(ggplot2)
#graph
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=1, shape=19)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Larval Liver", limits=c(-13, 13)) + #play with margin
  scale_x_continuous(name="Whole Larvae", limits=c(-13, 13)) + #play with margin
  scale_color_manual(values=c('#4BC0C8','#C779D0', '#FEAC5E'))
all_plot

#write.csv(ACTS_ALL_701_gene.ids, file="ACTS_ALL_701_gene.ids.csv")

#Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Both < 0.05"), ]#Add a extra column to the dataset


#now plot ONLY ACTS
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "WL < 0.05"), ]
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "LL < 0.05"), ]

library(ggplot2)
library(tidyverse)
library(ggrepel)
#graph
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=1, shape=19)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Whole Larvae DEGs", limits=c(-13, 13)) + #play with margin
  scale_x_continuous(name="Larval Liver DEGs", limits=c(-13, 13)) + #play with margin
  scale_color_manual(values=c('#4BC0C8')) +
  geom_smooth(method='lm', color="magenta", size=0.5, formula= y~x)
all_plot


#Remove ACTs going in different directions

ActsUP <- Mergefile_final[Mergefile_final$log2FoldChange.x < 0 & Mergefile_final$log2FoldChange.y < 0, ]    
ActsDOWN <- Mergefile_final[Mergefile_final$log2FoldChange.x > 0 & Mergefile_final$log2FoldChange.y > 0, ]
Acts_All_WL_LL_PosCor <- rbind(ActsUP,ActsDOWN)
LL_Acts_608 <- subset(LL, LL$gene.ids %in% Acts_All_WL_LL_PosCor$gene.ids)
#write.csv(LL_Acts_608, file="LL_Acts_608.csv")


#plot for paper:

Mergefile_final <- Mergefile_final[, -c(6)]
Mergefile_final$Legend <- NA

#assign attributes for color based on significance
Mergefile_final[Mergefile_final$log2FoldChange.x<0&Mergefile_final$log2FoldChange.y<0, ]$Legend <- "Pos Cor"
Mergefile_final[Mergefile_final$log2FoldChange.x>0&Mergefile_final$log2FoldChange.y>0, ]$Legend <- "Pos Cor"
Mergefile_final[Mergefile_final$log2FoldChange.x<0&Mergefile_final$log2FoldChange.y>0, ]$Legend <- "Neg Cor"                
Mergefile_final[Mergefile_final$log2FoldChange.x>0&Mergefile_final$log2FoldChange.y<0, ]$Legend <- "Neg Cor"  

library(ggrepel)

#Crossplot
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=1, shape=19)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Whole Larvae DEGs", limits=c(-13, 13)) + #play with margin
  scale_x_continuous(name="Larval Liver DEGs", limits=c(-13, 13)) + #play with margin
  scale_color_manual(values=c('grey','#4BC0C8')) +
  geom_smooth(method='lm', color="magenta", size=0.5, formula= y~x) +
  geom_point(data = Acts_All_WL_LL_PosCor %>% filter(gene.ids == "ENSDARG00000027572"), color = "magenta")
all_plot

r <- cor(Acts_All_WL_LL_PosCor$log2FoldChange.x, Acts_All_WL_LL_PosCor$log2FoldChange.y, method = "pearson")
r #0.832

#GO of UP and down of ACTS
Acts <- subset(Mergefile_final, Mergefile_final$Legend == "Pos Cor")

DOWN_ACTs <- subset(Acts, Acts$log2FoldChange.x < 0)
UP_ACTs <- subset(Acts, Acts$log2FoldChange.x > 0)

library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Dr.eg.db

# Curate GO terms
#Functions
for_plotting_comparecluster = function(list,finalgo,no_top_IDs,text_wrapsize){
  finaldb = as.data.frame(finalgo)
  topids = vector()
  names = names(list)
  
  #extract top ids for each group
  for (i in 1:length(names)){
    name = names[i]
    top = finaldb[finaldb$Cluster==name,]
    top = top[order(top$p.adjust),]
    top = top$ID[1:no_top_IDs]
    topids = c(topids,top)
  }
  
  dropids = unique(finaldb$ID)
  dropids = dropids[!dropids %in% topids]
  #subset GO to only include top IDs
  plotdb = dropGO(finalgo,term=dropids)
  
  library(dplyr)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  
  plot_table = as.data.frame(plotdb)
  plot_table$GeneRatio = as.numeric(sub('/.*','',plot_table$GeneRatio))/as.numeric(sub('.*/','',plot_table$GeneRatio))
  plot_table$BgRatio = as.numeric(sub('/.*','',plot_table$BgRatio))/as.numeric(sub('.*/','',plot_table$BgRatio))
  plot_table$Description = str_wrap(plot_table$Description, width = text_wrapsize)
  plot_table$p.adjust = as.numeric(plot_table$p.adjust)
  plot_table$Description = factor(plot_table$Description,levels=rev(unique(plot_table$Description)))
  
  return(plot_table)
}

# create table for uploading on REVIGO
create_revigo_table = function(go, filename){
  db = as.data.frame(go)
  db = db[order(db$p.adjust),]
  keepgos = data.frame(ID=vector(),p.adjust=vector())
  #create list of unique GO IDs with lowest p-value of all groups
  for (i in 1:nrow(db)){
    if (!db$ID[i] %in% keepgos$ID){
      keepgos = rbind(keepgos,db[i,])
    }
  }
  write.table(keepgos[,c(2,6)],quote = F,sep = ' ',row.names = F,col.names=F,file = filename)
}

drop_from_revigo = function(current_go,revigo_filename){
  revi = read.csv("Revigo_BP_Table.csv")
  db = as.data.frame(current_go)
  ids = db$ID[!db$ID %in% revi$TermID[revi$Representative==-1]]
  dropped_go = dropGO(current_go,term = ids)
  return(dropped_go)
}

# Perform GO on these clusters
fig = list(UP=Mergefile_final_1$gene.ids, DOWN=Mergefile_final_2$gene.ids)

fig_go = compareCluster(fig, fun = "enrichGO", keyType = "ENSEMBL", OrgDb = "org.Dr.eg.db" , ont = "BP",readable = TRUE)
create_revigo_table(fig_go, "fig_goidswithp.txt")

# Copy and paste this file to REVIGO (http://revigo.irb.hr) (add Danio rerio species) and retrieve summary of GO terms as csv file and save it in the working directory as REVIGO_fig.csv
revi = read.delim("Revigo_BP_Table.csv")
revi=as.data.frame(revi)
db = as.data.frame(fig_go)

ids = db$ID[!db$ID %in% revi$TermID[revi$Dispensability<=0.5]]
dropped_go = dropGO(fig_go,term = ids)

scale = 1.5
wrapsize = 42

#Plot terms - 8 terms
fig_plottable = for_plotting_comparecluster(fig,dropped_go,10,wrapsize) # change number depending on how many top terms you would like to be shown
write.table(fig_plottable, file = "GO_BP_10terms.txt", quote = F, sep = "\t ", row.names = F)

plot_fig = ggplot(fig_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="magenta", high="cyan", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c('UP', "Down"))+
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size = 22,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18),
        text = element_text(size=15),
        axis.title.x =element_text(size=22))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 18*scale),breaks=c(0.02,0.04,0.08,0.14))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.002,fontface='bold'),colour='white',show.legend = F)

plot_fig

postscript("AS3MTWLvsLL_CC_PosCOR_608_BP.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 11, height = 13)
plot_fig
dev.off()


#Figure S4 ----
Larvae <- read.csv("Table S2 as3mtmut liver.csv")
LarSig <- subset(Larvae, Larvae$padj < 0.05)

Female <- read.csv("Table S5 Adult Female.csv")
FemaleSig <- subset(Female, Female$padj < 0.05)

Male <- read.csv("Table S4 Adult Male.csv")
MaleSig <- subset(Male, Male$padj < 0.05)

#Crossplot of Male vs Female
#Subset and rename columns (FEMALE) X
AdultFemale <- read.csv("Table S5 Adult Female.csv")
#Create dataframe with ONLY columns I need: log2FC, padj and gene names
AdultFemale <- AdultFemale [, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(AdultFemale)[names(AdultFemale) == "log2FoldChange"] <- "log2FoldChange.x"
names(AdultFemale)[names(AdultFemale) == "padj"] <- "padj.x"

#Subset and rename columns (MALE) Y
AdultMale <- read.csv("Table S4 Adult Male.csv")
AdultMale <- AdultMale [, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(AdultMale)[names(AdultMale) == "log2FoldChange"] <- "log2FoldChange.y"
names(AdultMale)[names(AdultMale) == "padj"] <- "padj.y"

#Mergedatasets
Mergedfile = merge(AdultFemale,AdultMale, by="gene.ids", all=T)

#Clean up merged (remove duplicates and NA if necessary - this merge was fine and did not need to remove dups or NAs but useful when combining more complex datasets) 
Mergefile_clean <- na.omit(Mergedfile) #removes NA
#install.packages("tidyverse")
library(tidyverse)
Mergefile_final <- Mergefile_clean[!duplicated(Mergefile_clean$gene.ids), ] #removes duplicates

#Now let's color things --> lets make a blank legend column to add to DF

Mergefile_final$Legend <- NA

#assign attributes for color based on padj
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "Female < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "Male < 0.05"
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "Both < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "Neither"

#Remove neither
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Neither"), ]
#Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Both < 0.05"), ] #I wanted to see both for this analysis but you can remove

library(ggplot2)
#graph
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=2, shape=19)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Male Adult Liver as3mt-/-", limits=c(-15, 15)) + #play with margin
  scale_x_continuous(name="Female Adult Liver as3mt-/-", limits=c(-18, 18)) + #play with margin
  scale_color_manual(values=c('#660066','#33ffcc', '#FF6633'))
all_plot

# Venn
make_3way_eulerr = function(list1,list2,list3,name,colours,scale){
  n1 = length(list1)
  n2 = length(list2)
  n3 = length(list3)
  
  n12 = length(intersect(list1,list2))
  n23 = length(intersect(list2,list3))
  n13 = length(intersect(list1,list3))
  
  n123 = length(intersect(intersect(list1,list2),list3))
  
  comb = c('1'=n1,'2'=n2,'3'=n3,
           '1&2'=n12,'2&3'=n23,'1&3'=n13,
           '1&2&3'=n123)
  venn = euler(comb, input = 'union', shape = 'ellipse')
  save_venn_plots(venn,name,colours,scale)
}

save_venn_plots = function(eulerr,name,colours,scale){
  labelname = paste(name,'_labels.png',sep='')
  nolabelname = paste(name,'.png',sep='')
  #to see labels to add manually in illustrator
  png(labelname)
  print({plot(eulerr,
              labels=FALSE,
              quantities = TRUE,
              fills = colours)
  })
  dev.off()
  
  png(nolabelname,width = 400*scale,height = 400*scale)
  print({plot(eulerr,edges = list(lwd=1*scale),
              labels=FALSE,
              quantities = FALSE,
              fills = colours)
  })
  dev.off()
}

larLiver = '#C779D0' 
as3mtfemale = '#33ffcc' 
as3mtmale =  '#FF6633' 

#Found colors using https://hexcolorspicker.com/google-color-picker/

As3mtMUT_venn = make_3way_eulerr(LarSig$gene.ids, FemaleSig$gene.ids, MaleSig$gene.ids, 'Figures/Venn',
                                 c(larLiver, as3mtfemale, as3mtmale),7)

# Heatmap
# Let's make heatmaps of subsets

#Upload files
LarSig <- subset(Larvae, Larvae$padj < 0.05)
FemaleSig <- subset(Female, Female$padj < 0.05)
MaleSig <- subset(Male, Male$padj < 0.05)

#I just need column 4 (log2fc) and 9 (gene.ids)
LarH <- Larvae[, -c(1,3,5,6,7,8,9,10,11,12,13,14,15,16,17,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(LarH)[names(LarH) == "log2FoldChange"] <- "Lar"

#I just need column 4 (log2fc) and 9 (gene.names)
FemaleH <- Female[, -c(1,3,5,6,7,8,9,10,11,12,13,14,15,16)]
names(FemaleH)[names(FemaleH) == "log2FoldChange"] <- "Female"

#I just need column 4 (log2fc) and 9 (gene.names)
MaleH <- Male[, -c(1,3,5,6,7,8,9,10,11,12,13,14,15,16)]
names(MaleH)[names(MaleH) == "log2FoldChange"] <- "Male"

# Subset - 1st all 50 shared
Fifty <- merge.data.frame(LarH, MaleH, by='gene.ids', all=T)
Fifty <- merge.data.frame(Fifty, FemaleH, by='gene.ids', all=T)
Fifty <- subset(Fifty, Fifty$gene.ids %in% LarSig$gene.ids)
Fifty <- subset(Fifty, Fifty$gene.ids %in% FemaleSig$gene.ids)
Fifty <- subset(Fifty, Fifty$gene.ids %in% MaleSig$gene.ids)


#fifty in lar

LL_50 <- subset(LL, LL$gene.ids %in% Fifty$gene.ids)

#end with 50 genes, yay! 

#make matrix
dat_matrix <- as.matrix(Fifty[1: nrow(Fifty), 2: ncol(Fifty)])
rownames(dat_matrix) <- Fifty$gene.ids
colnames(dat_matrix) <- colnames(Fifty)[2: ncol(Fifty)]
head(dat_matrix) #looks good: gene.names are row names and log2FC are the row values


#better legend
library(RColorBrewer)
library(pheatmap)
col <- brewer.pal(n=10,"PRGn") #I used 8 because my scale has 10 breakpoints but you can increase or decrease this number based on your needs. Look at RColorBrewer options online to change the coloring if you would like! 

#make better legend 
breaksList = seq(-10, 10, by = 2) #even out axis so that colors are centered at 0 - easier to understand the heatmap! IMPORTANT: Edit the upper and lower limit based on log2fc range in your dataset
pheatmap(dat_matrix, cluster_cols=F, fontsize_row=2, cellheight=2, cellwidth = 100, cutree_rows =2, display_numbers = FALSE, breaks = breaksList, color=col, border_color = "black")

# Subset - 2nd 85
EightyFive <- merge.data.frame(LarH, MaleH, by='gene.ids', all=T)
EightyFive <- merge.data.frame(EightyFive, FemaleH, by='gene.ids', all=T)
EightyFive <- subset(EightyFive, EightyFive$gene.ids %in% LarSig$gene.ids)
EightyFive <- subset(EightyFive, EightyFive$gene.ids %in% FemaleSig$gene.ids)
EightyFive <- EightyFive[!(EightyFive$gene.ids %in% Fifty$gene.ids),]
#end with 85 genes, yay! 

#make matrix
dat_matrix <- as.matrix(EightyFive[1: nrow(EightyFive), 2: ncol(EightyFive)])
rownames(dat_matrix) <- EightyFive$gene.ids
colnames(dat_matrix) <- colnames(EightyFive)[2: ncol(EightyFive)]
head(dat_matrix) #looks good: gene.names are row names and log2FC are the row values

#better legend
library(RColorBrewer)
col <- brewer.pal(n=10,"PRGn") #I used 8 because my scale has 10 breakpoints but you can increase or decrease this number based on your needs. Look at RColorBrewer options online to change the coloring if you would like! 

#make better legend 
breaksList = seq(-10, 10, by = 2) #even out axis so that colors are centered at 0 - easier to understand the heatmap! IMPORTANT: Edit the upper and lower limit based on log2fc range in your dataset
pheatmap(dat_matrix, cluster_cols=F, fontsize_row=2, cellheight=2, cellwidth = 100, cutree_rows =2, display_numbers = FALSE, breaks = breaksList, color=col, border_color = "black")

# Subset - 3rd 162
OneHunSixtTwo <- merge.data.frame(LarH, MaleH, by='gene.ids', all=T)
OneHunSixtTwo <- merge.data.frame(OneHunSixtTwo, FemaleH, by='gene.ids', all=T)
OneHunSixtTwo <- subset(OneHunSixtTwo, OneHunSixtTwo$gene.ids %in% LarSig$gene.ids)
OneHunSixtTwo <- subset(OneHunSixtTwo, OneHunSixtTwo$gene.ids %in% MaleSig$gene.ids)
OneHunSixtTwo <- OneHunSixtTwo[!(OneHunSixtTwo$gene.ids %in% Fifty$gene.ids),]
#end with 162 genes, yay! 

#make matrix
dat_matrix <- as.matrix(OneHunSixtTwo[1: nrow(OneHunSixtTwo), 2: ncol(OneHunSixtTwo)])
rownames(dat_matrix) <- OneHunSixtTwo$gene.ids
colnames(dat_matrix) <- colnames(OneHunSixtTwo)[2: ncol(OneHunSixtTwo)]
head(dat_matrix) #looks good: gene.names are row names and log2FC are the row values

#better legend
library(RColorBrewer)
col <- brewer.pal(n=10,"PRGn") #I used 8 because my scale has 10 breakpoints but you can increase or decrease this number based on your needs. Look at RColorBrewer options online to change the coloring if you would like! 

#make better legend 
breaksList = seq(-10, 10, by = 2) #even out axis so that colors are centered at 0 - easier to understand the heatmap! IMPORTANT: Edit the upper and lower limit based on log2fc range in your dataset
pheatmap(dat_matrix, cluster_cols=F, fontsize_row=2, cellheight=2, cellwidth = 100, cutree_rows =2, display_numbers = FALSE, breaks = breaksList, color=col, border_color = "black")

#Cross plot of larvae and male positively correlated genes

Male <- read.csv("Table S4 Adult Male.csv")
LL <- read.csv("Table S2 as3mtmut liver.csv")

#Subset and rename columns (WL) X
WT_1mM_forMerge <- Male[, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(WT_1mM_forMerge)[names(WT_1mM_forMerge) == "log2FoldChange"] <- "log2FoldChange.x"
names(WT_1mM_forMerge)[names(WT_1mM_forMerge) == "padj"] <- "padj.x"
WT_1mM_Sig <- subset.data.frame(Male, Male$padj<0.05)

#Subset and rename columns (LL) Y
A_1mM_forMerge  <- LL [, -c(1,3,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(A_1mM_forMerge)[names(A_1mM_forMerge) == "log2FoldChange"] <- "log2FoldChange.y"
names(A_1mM_forMerge)[names(A_1mM_forMerge) == "padj"] <- "padj.y"
As3mt_1mM_Sig <- subset.data.frame(LL, LL$padj<0.05)

#Mergedatasets
Mergedfile = merge(WT_1mM_forMerge,A_1mM_forMerge, by="gene.ids", all=T)

#Clean up merged (remove duplicates and NA) 
Mergefile_clean <- na.omit(Mergedfile) #removes NA
#install.packages("tidyverse")
library(tidyverse)
Mergefile_final <- Mergefile_clean[!duplicated(Mergefile_clean$gene.ids), ] #removes duplicates

Mergefile_final$Legend <- NA
#assign attributes for color based on significance
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "WL < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "LL < 0.05"
Mergefile_final[Mergefile_final$padj.x<0.05&Mergefile_final$padj.y<0.05, ]$Legend <- "Both < 0.05"
Mergefile_final[Mergefile_final$padj.x>0.05&Mergefile_final$padj.y>0.05, ]$Legend <- "Neither"

#Remove terms
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Neither"), ]
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "WL < 0.05"), ]
Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "LL < 0.05"), ]
#Mergefile_final <- Mergefile_final[!(Mergefile_final$Legend == "Both < 0.05"), ]#Add a extra column to the dataset

library(ggplot2)
#graph
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=1, shape=19)  +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Larval Liver", limits=c(-13, 13)) + #play with margin
  scale_x_continuous(name="Whole Larvae", limits=c(-13, 13)) + #play with margin
  scale_color_manual(values=c('#4BC0C8','#C779D0', '#FEAC5E'))
all_plot

#plot with pos and cor colored 
Mergefile_final <- Mergefile_final[, -c(6)]
Mergefile_final$Legend <- NA

#assign attributes for color based on significance
Mergefile_final[Mergefile_final$log2FoldChange.x<0&Mergefile_final$log2FoldChange.y<0, ]$Legend <- "Pos Cor"
Mergefile_final[Mergefile_final$log2FoldChange.x>0&Mergefile_final$log2FoldChange.y>0, ]$Legend <- "Pos Cor"
Mergefile_final[Mergefile_final$log2FoldChange.x<0&Mergefile_final$log2FoldChange.y>0, ]$Legend <- "Neg Cor"                
Mergefile_final[Mergefile_final$log2FoldChange.x>0&Mergefile_final$log2FoldChange.y<0, ]$Legend <- "Neg Cor"  

library(ggrepel)

#Crossplot
all_plot <-ggplot(Mergefile_final, aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) +
  geom_point(size=1, shape=19) +
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  scale_y_continuous(name="Whole Larvae DEGs", limits=c(-13, 13)) + #play with margin
  scale_x_continuous(name="Larval Liver DEGs", limits=c(-13, 13)) + #play with margin
  scale_color_manual(values=c('grey','#4BC0C8')) +
  geom_smooth(method='lm', color="magenta", size=0.5, formula= y~x) +
  geom_point(data = Mergefile_final %>% filter(gene.ids == "ENSDARG00000103614"), color = "magenta") +
  geom_point(data = Mergefile_final %>% filter(gene.ids == "ENSDARG00000045132"), color = "magenta")
all_plot
 
Acts_All_Male_LL_PosCor <- subset(Mergefile_final, Mergefile_final$Legend == "Pos Cor")

r <- cor(Acts_All_Male_LL_PosCor$log2FoldChange.x, Acts_All_Male_LL_PosCor$log2FoldChange.y, method = "pearson")
r #0.894


#Figure 2F ----
Male <- read.csv("Table S4 Adult Male.csv")
LL <- read.csv("Table S2 as3mtmut liver.csv")
Female <- read.csv("Table S5 Adult Female.csv")
library(readr)
fig_plottable_UPSET_TERMS <- read_csv("fig_plottable_UPSET_TERMS.csv")

#Mitp
Ribo_Lar <- subset(LL,LL$gene.names %in% fig_plottable_UPSET_TERMS$Mito)
Ribo_Adult_Male <- subset(Male, Male$gene.names %in% fig_plottable_UPSET_TERMS$Mito)
Ribo_Adult_Female <- subset(Female, Female$gene.names %in% fig_plottable_UPSET_TERMS$Mito)

#I just need column 4 (log2fc) and 9 (gene.ids)
LarH <- Ribo_Lar[, -c(1,2,3,5,6,7,8,10,11,12,13,14,15,16,17,17,18,19,20,21,22,23,24,25,26,27,28,29)]
names(LarH)[names(LarH) == "log2FoldChange"] <- "Lar"

#I just need column 4 (log2fc) and 9 (gene.names)
FemaleH <- Ribo_Adult_Female[, -c(1,2,3,5,6,7,8,10,11,12,13,14,15,16)]
names(FemaleH)[names(FemaleH) == "log2FoldChange"] <- "Female"

#I just need column 4 (log2fc) and 9 (gene.names)
MaleH <- Ribo_Adult_Male[, -c(1,2,3,5,6,7,8,10,11,12,13,14,15,16)]
names(MaleH)[names(MaleH) == "log2FoldChange"] <- "Male"

# Subset - 1st all 50 shared
Ribo <- merge.data.frame(FemaleH, MaleH, by='gene.names', all=T)
Ribo <- merge.data.frame(Ribo, LarH, by='gene.names', all=T)
Ribo <- Ribo[order(Ribo$Lar),]

#make matrix
dat_matrix <- as.matrix(Ribo[1: nrow(Ribo), 2: ncol(Ribo)])
rownames(dat_matrix) <- Ribo$gene.names
colnames(dat_matrix) <- colnames(Ribo)[2: ncol(Ribo)]
head(dat_matrix) #looks good: gene.names are row names and log2FC are the row values

#better legend
library(RColorBrewer)
library(pheatmap)
col <- brewer.pal(n=10,"PRGn") #I used 8 because my scale has 10 breakpoints but you can increase or decrease this number based on your needs. Look at RColorBrewer options online to change the coloring if you would like! 

#make better legend 
breaksList = seq(-5, 5, by = 1) #even out axis so that colors are centered at 0 - easier to understand the heatmap! IMPORTANT: Edit the upper and lower limit based on log2fc range in your dataset
pheatmap(dat_matrix,cluster_cols=FALSE, cluster_rows=FALSE, fontsize_row=11, cellheight=15, cellwidth = 15, cutree_rows =1, display_numbers = FALSE, breaks = breaksList, color=col, border_color = "black")

#Figure 3

#Figure 5 














#Figure 3