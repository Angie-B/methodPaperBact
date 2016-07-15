library(ggplot2)
library(reshape2)
require(graphics); require(grDevices)
library(Hmisc)
library(gtools)
library(cowplot)
require(gplots)
require(RColorBrewer)
library(xlsx)

## Get the averages ------------------------
setwd("M:/Summer_2016_metabolites/Methods_Paper_Samples/qc")

mydata.overloaded <- read.csv("M:/Summer_2016_metabolites/Methods_Paper_Samples/qc/QC.output.new.Bacteria_HILIC_160613_1to2_update_akb.csv" ,header=TRUE, comment.char="", as.is=TRUE)
mydata <- mydata.overloaded
# mydata$Area[grep("overloaded?",mydata$Notes)]<-2e10

s <- split(mydata$Area,list(mydata$Compound.Name,mydata$group))
means <- sapply(s,mean,na.rm=T)

## make into a useful dataframe of averages --------------
mydf<-as.data.frame(means)
cmpd <- c()
group <- c()
for (i in 1:nrow(mydf)){
     this.name <- row.names(mydf)[i]
     this.cmpd <- strsplit(this.name,"[.]")[[1]][1]
     this.group <- strsplit(this.name,"[.]")[[1]][2]
     cmpd <- c(cmpd, this.cmpd)
     group <- c(group,this.group)
}
mydf$group <- group
mydf$cmpd <- cmpd

wide <- reshape(mydf, v.names = "means", idvar = "group",
                timevar = "cmpd", direction = "wide")
dat <- wide[,2:141]
row.names(dat) <- wide$group

dat.2 <- dat
dat.2[is.na(dat.2)] <- 0

dat.2 <- as.matrix(dat.2)


## normalize to D-Methionine ------------------------

s.norm <- split(mydata$Area,mydata$Compound.Name)
test.norm <- c()
for (i in 1:length(s.norm)){
     test.norm <- cbind(test.norm,s.norm[[i]])
     print(paste(i,length(s.norm[[i]])))
}
colnames(test.norm) <- names(s.norm)
rownames(test.norm) <- mydata$Replicate.Name[c(1:(nrow(test.norm)-1),nrow(mydata))]
norm.df <- as.data.frame(test.norm)
new.df <- apply(norm.df,2,function(x) x/norm.df[,"D-Methionine"])

## remove overloaded compounds
OL.cmpds <- unique(mydata$Compound.Name[grep("overloaded?",mydata$Notes)])
data.no.overload.noOD <- new.df[,-match(OL.cmpds,colnames(new.df))]

## remove internal standards
IS.cmpds <- unique(mydata$Compound.Name[grep("13C",mydata$Compound.Name)])
IS.cmpds <- c(IS.cmpds, unique(mydata$Compound.Name[grep("D-",mydata$Compound.Name)]) )
IS.cmpds <- c(IS.cmpds, unique(mydata$Compound.Name[grep("Heavy",mydata$Compound.Name)]) )
IS.cmpds <- c(IS.cmpds, unique(mydata$Compound.Name[grep("d3",mydata$Compound.Name)]) )
IS.cmpds <- c(IS.cmpds, unique(mydata$Compound.Name[grep("d4",mydata$Compound.Name)]) )

data.no.overload.noOD <- data.no.overload.noOD[,-match(IS.cmpds,colnames(data.no.overload.noOD))]

## normalize to optical density -------------------------------
dens <- read.xlsx("M:/Summer_2016_metabolites/Methods_Paper_Samples/Cultures_Extraction_Log.xlsx", 1)  # read first sheet
dens.order <- order(dens$Vial.ID)[2:16]
new.df<- data.no.overload.noOD
new.df$OD <- as.numeric(as.character(dens$OD.1.10..at.Harvest.time[dens.order]))

OD.Norm.df <- apply(new.df,2,function(x) x/new.df[,"OD"])

## make some heatmap plots ------------------


## another version
require(graphics); require(grDevices)

x  <- as.matrix(OD.Norm.df)
x[is.na(x)] <- 0

ditch.list <- c()
for (i in 1:ncol(x)){
     if (colMeans(x)[i]==0){
          ditch.list <- c(ditch.list,i)
     }
}
x <- x[,-ditch.list]
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
## need to get rid of NAs here!!!!
hv <- heatmap(x, col = cm.colors(256), scale = "column",
               margins = c(5,5),
              xlab = "Compound", ylab =  "Sample type",
              main = "heatmap of bacteria 1:2 dilutions")
utils::str(hv) # the two re-ordering index vectors

## another version
## make my color bar
range(x)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 399)
colors.2 <- c(seq(-1000,-300.1,length=100),seq(-300,0,length=100),seq(0.1,300,length=100),seq(301,1000,length=100))
## make heatmap
hm <- heatmap.2(x, 
                breaks = colors.2, 
                col = my_palette, trace="none", 
                symm=F,symkey=F,symbreaks=T, scale="none",
                margins = c(5,10))
utils::str(hm) # the two re-ordering index vectors