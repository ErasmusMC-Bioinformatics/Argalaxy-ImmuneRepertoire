# ---------------------- load/install packages ----------------------

if (!("gridExtra" %in% rownames(installed.packages()))) {
  install.packages("gridExtra", repos="http://cran.xl-mirror.nl/") 
}
library(gridExtra)
if (!("ggplot2" %in% rownames(installed.packages()))) {
  install.packages("ggplot2", repos="http://cran.xl-mirror.nl/") 
}
library(ggplot2)
if (!("plyr" %in% rownames(installed.packages()))) {
  install.packages("plyr", repos="http://cran.xl-mirror.nl/") 
}
library(plyr)

if (!("data.table" %in% rownames(installed.packages()))) {
  install.packages("data.table", repos="http://cran.xl-mirror.nl/") 
}
library(data.table)

if (!("reshape2" %in% rownames(installed.packages()))) {
  install.packages("reshape2", repos="http://cran.xl-mirror.nl/")
}
library(reshape2)

if (!("lymphclon" %in% rownames(installed.packages()))) {
  install.packages("lymphclon", repos="http://cran.xl-mirror.nl/")
}
library(lymphclon)

# ---------------------- parameters ----------------------

args <- commandArgs(trailingOnly = TRUE)

infile = args[1] #path to input file
outfile = args[2] #path to output file
outdir = args[3] #path to output folder (html/images/data)
clonaltype = args[4] #clonaltype definition, or 'none' for no unique filtering
ct = unlist(strsplit(clonaltype, ","))
species = args[5] #human or mouse
locus = args[6] # IGH, IGK, IGL, TRB, TRA, TRG or TRD
filterproductive = ifelse(args[7] == "yes", T, F) #should unproductive sequences be filtered out? (yes/no)
clonality_method = args[8]


# ---------------------- Data preperation ----------------------

print("Report Clonality - Data preperation")

inputdata = read.table(infile, sep="\t", header=TRUE, fill=T, comment.char="", stringsAsFactors=F)

inputdata$Sample = as.character(inputdata$Sample)


print(paste("nrows: ", nrow(inputdata)))

setwd(outdir)

# remove weird rows
inputdata = inputdata[inputdata$Sample != "",]

print(paste("nrows: ", nrow(inputdata)))

#remove the allele from the V,D and J genes
inputdata$Top.V.Gene = gsub("[*]([0-9]+)", "", inputdata$Top.V.Gene)
inputdata$Top.D.Gene = gsub("[*]([0-9]+)", "", inputdata$Top.D.Gene)
inputdata$Top.J.Gene = gsub("[*]([0-9]+)", "", inputdata$Top.J.Gene)

print(paste("nrows: ", nrow(inputdata)))

#filter uniques
inputdata.removed = inputdata[NULL,]

print(paste("nrows: ", nrow(inputdata)))

inputdata$clonaltype = 1:nrow(inputdata)

#keep track of the count of sequences in samples or samples/replicates for the front page overview
input.sample.count = data.frame(data.table(inputdata)[, list(All=.N), by=c("Sample")])
input.rep.count = data.frame(data.table(inputdata)[, list(All=.N), by=c("Sample", "Replicate")])

PRODF = inputdata
UNPROD = inputdata
if(filterproductive){
  if("Functionality" %in% colnames(inputdata)) { # "Functionality" is an IMGT column
    #PRODF = inputdata[inputdata$Functionality == "productive" | inputdata$Functionality == "productive (see comment)", ]
    PRODF = inputdata[inputdata$Functionality %in% c("productive (see comment)","productive"),]
    
    PRODF.count = data.frame(data.table(PRODF)[, list(count=.N), by=c("Sample")])
    
    UNPROD = inputdata[inputdata$Functionality %in% c("unproductive (see comment)","unproductive"), ]
  } else {
    PRODF = inputdata[inputdata$VDJ.Frame != "In-frame with stop codon" & inputdata$VDJ.Frame != "Out-of-frame" & inputdata$CDR3.Found.How != "NOT_FOUND" , ]
    UNPROD = inputdata[!(inputdata$VDJ.Frame != "In-frame with stop codon" & inputdata$VDJ.Frame != "Out-of-frame" & inputdata$CDR3.Found.How != "NOT_FOUND" ), ]
  }
}

for(i in 1:nrow(UNPROD)){
    if(!is.numeric(UNPROD[i,"CDR3.Length"])){
        UNPROD[i,"CDR3.Length"] = 0
    }
}

prod.sample.count = data.frame(data.table(PRODF)[, list(Productive=.N), by=c("Sample")])
prod.rep.count = data.frame(data.table(PRODF)[, list(Productive=.N), by=c("Sample", "Replicate")])

unprod.sample.count = data.frame(data.table(UNPROD)[, list(Unproductive=.N), by=c("Sample")])
unprod.rep.count = data.frame(data.table(UNPROD)[, list(Unproductive=.N), by=c("Sample", "Replicate")])

clonalityFrame = PRODF

#remove duplicates based on the clonaltype
if(clonaltype != "none"){
  clonaltype = paste(clonaltype, ",Sample", sep="") #add sample column to clonaltype, unique within samples
  PRODF$clonaltype = do.call(paste, c(PRODF[unlist(strsplit(clonaltype, ","))], sep = ":"))
  PRODF = PRODF[!duplicated(PRODF$clonaltype), ]
    
  UNPROD$clonaltype = do.call(paste, c(UNPROD[unlist(strsplit(clonaltype, ","))], sep = ":"))
  UNPROD = UNPROD[!duplicated(UNPROD$clonaltype), ]
  
  #again for clonalityFrame but with sample+replicate
  clonalityFrame$clonaltype = do.call(paste, c(clonalityFrame[unlist(strsplit(clonaltype, ","))], sep = ":"))
  clonalityFrame$clonality_clonaltype = do.call(paste, c(clonalityFrame[unlist(strsplit(paste(clonaltype, ",Replicate", sep=""), ","))], sep = ":"))
  clonalityFrame = clonalityFrame[!duplicated(clonalityFrame$clonality_clonaltype), ]
}

if(nrow(PRODF) == 0){
	stop("No sequences left after filtering")
}

prod.unique.sample.count = data.frame(data.table(PRODF)[, list(Productive_unique=.N), by=c("Sample")])
prod.unique.rep.count = data.frame(data.table(PRODF)[, list(Productive_unique=.N), by=c("Sample", "Replicate")])

unprod.unique.sample.count = data.frame(data.table(UNPROD)[, list(Unproductive_unique=.N), by=c("Sample")])
unprod.unique.rep.count = data.frame(data.table(UNPROD)[, list(Unproductive_unique=.N), by=c("Sample", "Replicate")])


PRODF$freq = 1

if(any(grepl(pattern="_", x=PRODF$ID))){ #the frequency can be stored in the ID with the pattern ".*_freq_.*"
  PRODF$freq = gsub("^[0-9]+_", "", PRODF$ID)
  PRODF$freq = gsub("_.*", "", PRODF$freq)
  PRODF$freq = as.numeric(PRODF$freq)
  if(any(is.na(PRODF$freq))){ #if there was an "_" in the ID, but not the frequency, go back to frequency of 1 for every sequence
    PRODF$freq = 1
  }
}

#make a names list with sample -> color
naive.colors = c('blue4', 'darkred', 'olivedrab3', 'red', 'gray74', 'darkviolet', 'lightblue1', 'gold', 'chartreuse2', 'pink', 'Paleturquoise3', 'Chocolate1', 'Yellow', 'Deeppink3', 'Mediumorchid1', 'Darkgreen', 'Blue', 'Gray36', 'Hotpink', 'Yellow4')
unique.samples = unique(PRODF$Sample)

if(length(unique.samples) <= length(naive.colors)){
	sample.colors = naive.colors[1:length(unique.samples)]
} else {
	sample.colors = rainbow(length(unique.samples))
}

names(sample.colors) = unique.samples

print("Sample.colors")
print(sample.colors)


#write the complete dataset that is left over, will be the input if 'none' for clonaltype and 'no' for filterproductive
write.table(PRODF, "allUnique.txt", sep="\t",quote=F,row.names=F,col.names=T)
#write.table(PRODF, "allUnique.csv", sep=",",quote=F,row.names=F,col.names=T)
write.table(UNPROD, "allUnproductive.txt", sep="\t",quote=F,row.names=F,col.names=T)

print("SAMPLE TABLE:")
print(table(PRODF$Sample))

#write the samples to a file
sampleFile <- file("samples.txt")
un = unique(inputdata$Sample)
un = paste(un, sep="\n")
writeLines(un, sampleFile)
close(sampleFile)

# ---------------------- Counting the productive/unproductive and unique sequences ----------------------

print("Report Clonality - counting productive/unproductive/unique")

#create the table on the overview page with the productive/unique counts per sample/replicate
#first for sample

sample.count = merge(input.sample.count, prod.sample.count, by="Sample", all.x=T)
sample.count$perc_prod = round(sample.count$Productive / sample.count$All * 100)
sample.count = merge(sample.count, prod.unique.sample.count, by="Sample", all.x=T)
sample.count$perc_prod_un = round(sample.count$Productive_unique / sample.count$All * 100)

sample.count = merge(sample.count , unprod.sample.count, by="Sample", all.x=T)
sample.count$perc_unprod = round(sample.count$Unproductive / sample.count$All * 100)
sample.count = merge(sample.count, unprod.unique.sample.count, by="Sample", all.x=T)
sample.count$perc_unprod_un = round(sample.count$Unproductive_unique / sample.count$All * 100)

#then sample/replicate
rep.count = merge(input.rep.count, prod.rep.count, by=c("Sample", "Replicate"), all.x=T)

print(rep.count)

fltr = is.na(rep.count$Productive)
if(any(fltr)){
	rep.count[fltr,"Productive"] = 0
}

print(rep.count)

rep.count$perc_prod = round(rep.count$Productive / rep.count$All * 100)
rep.count = merge(rep.count, prod.unique.rep.count, by=c("Sample", "Replicate"), all.x=T)
rep.count$perc_prod_un = round(rep.count$Productive_unique / rep.count$All * 100)

rep.count = merge(rep.count, unprod.rep.count, by=c("Sample", "Replicate"), all.x=T)
rep.count$perc_unprod = round(rep.count$Unproductive / rep.count$All * 100)
rep.count = merge(rep.count, unprod.unique.rep.count, by=c("Sample", "Replicate"), all.x=T)
rep.count$perc_unprod_un = round(rep.count$Unproductive_unique / rep.count$All * 100)

rep.count$Sample = paste(rep.count$Sample, rep.count$Replicate, sep="_")
rep.count = rep.count[,names(rep.count) != "Replicate"]

count = rbind(sample.count, rep.count)



write.table(x=count, file="productive_counting.txt", sep=",",quote=F,row.names=F,col.names=F)

# ---------------------- V+J+CDR3 sequence count ----------------------

VJCDR3.count = data.frame(table(clonalityFrame$Top.V.Gene, clonalityFrame$Top.J.Gene, clonalityFrame$CDR3.Seq.DNA))
names(VJCDR3.count) = c("Top.V.Gene", "Top.J.Gene", "CDR3.Seq.DNA", "Count")

VJCDR3.count = VJCDR3.count[VJCDR3.count$Count > 0,]
VJCDR3.count = VJCDR3.count[order(-VJCDR3.count$Count),]

write.table(x=VJCDR3.count, file="VJCDR3_count.txt", sep="\t",quote=F,row.names=F,col.names=T)

# ---------------------- Frequency calculation for V, D and J ----------------------

print("Report Clonality - frequency calculation V, D and J")

PRODFV = data.frame(data.table(PRODF)[, list(Length=sum(freq)), by=c("Sample", "Top.V.Gene")])
Total = ddply(PRODFV, .(Sample), function(x) data.frame(Total = sum(x$Length)))
PRODFV = merge(PRODFV, Total, by.x='Sample', by.y='Sample', all.x=TRUE)
PRODFV = ddply(PRODFV, c("Sample", "Top.V.Gene"), summarise, relFreq= (Length*100 / Total))

PRODFD = data.frame(data.table(PRODF)[, list(Length=sum(freq)), by=c("Sample", "Top.D.Gene")])
Total = ddply(PRODFD, .(Sample), function(x) data.frame(Total = sum(x$Length)))
PRODFD = merge(PRODFD, Total, by.x='Sample', by.y='Sample', all.x=TRUE)
PRODFD = ddply(PRODFD, c("Sample", "Top.D.Gene"), summarise, relFreq= (Length*100 / Total))

PRODFJ = data.frame(data.table(PRODF)[, list(Length=sum(freq)), by=c("Sample", "Top.J.Gene")])
Total = ddply(PRODFJ, .(Sample), function(x) data.frame(Total = sum(x$Length)))
PRODFJ = merge(PRODFJ, Total, by.x='Sample', by.y='Sample', all.x=TRUE)
PRODFJ = ddply(PRODFJ, c("Sample", "Top.J.Gene"), summarise, relFreq= (Length*100 / Total))

# ---------------------- Setting up the gene names for the different species/loci ----------------------

print("Report Clonality - getting genes for species/loci")

Vchain = ""
Dchain = ""
Jchain = ""

if(species == "custom"){
	print("Custom genes: ")
	splt = unlist(strsplit(locus, ";"))
	print(paste("V:", splt[1]))
	print(paste("D:", splt[2]))
	print(paste("J:", splt[3]))
	
	Vchain = unlist(strsplit(splt[1], ","))
	Vchain = data.frame(v.name = Vchain, chr.orderV = 1:length(Vchain))
	
	Dchain = unlist(strsplit(splt[2], ","))
	if(length(Dchain) > 0){
		Dchain = data.frame(v.name = Dchain, chr.orderD = 1:length(Dchain))
	} else {
		Dchain = data.frame(v.name = character(0), chr.orderD = numeric(0))
	}
	
	Jchain = unlist(strsplit(splt[3], ","))
	Jchain = data.frame(v.name = Jchain, chr.orderJ = 1:length(Jchain))

} else {
	genes = read.table("genes.txt", sep="\t", header=TRUE, fill=T, comment.char="")

	Vchain = genes[grepl(species, genes$Species) & genes$locus == locus & genes$region == "V",c("IMGT.GENE.DB", "chr.order")]
	colnames(Vchain) = c("v.name", "chr.orderV")
	Dchain = genes[grepl(species, genes$Species) & genes$locus == locus & genes$region == "D",c("IMGT.GENE.DB", "chr.order")]
	colnames(Dchain) = c("v.name", "chr.orderD")
	Jchain = genes[grepl(species, genes$Species) & genes$locus == locus & genes$region == "J",c("IMGT.GENE.DB", "chr.order")]
	colnames(Jchain) = c("v.name", "chr.orderJ")
}
useD = TRUE
if(nrow(Dchain) == 0){
  useD = FALSE
  cat("No D Genes in this species/locus")
}
print(paste(nrow(Vchain), "genes in V"))
print(paste(nrow(Dchain), "genes in D"))
print(paste(nrow(Jchain), "genes in J"))

# ---------------------- merge with the frequency count ----------------------

PRODFV = merge(PRODFV, Vchain, by.x='Top.V.Gene', by.y='v.name', all.x=TRUE)

PRODFD = merge(PRODFD, Dchain, by.x='Top.D.Gene', by.y='v.name', all.x=TRUE)

PRODFJ = merge(PRODFJ, Jchain, by.x='Top.J.Gene', by.y='v.name', all.x=TRUE)

# ---------------------- Create the V, D and J frequency plots and write the data.frame for every plot to a file ----------------------

print("Report Clonality - V, D and J frequency plots")

pV = ggplot(PRODFV)
pV = pV + geom_bar( aes( x=factor(reorder(Top.V.Gene, chr.orderV)), y=relFreq, fill=Sample), stat='identity', position="dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pV = pV + xlab("Summary of V gene") + ylab("Frequency") + ggtitle("Relative frequency of V gene usage") + scale_fill_manual(values=sample.colors)
pV = pV + theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())
write.table(x=PRODFV, file="VFrequency.txt", sep="\t",quote=F,row.names=F,col.names=T)

png("VPlot.png",width = 1280, height = 720)
pV
dev.off()

ggsave("VPlot.pdf", pV, width=13, height=7)

if(useD){
  pD = ggplot(PRODFD)
  pD = pD + geom_bar( aes( x=factor(reorder(Top.D.Gene, chr.orderD)), y=relFreq, fill=Sample), stat='identity', position="dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  pD = pD + xlab("Summary of D gene") + ylab("Frequency") + ggtitle("Relative frequency of D gene usage") + scale_fill_manual(values=sample.colors)
  pD = pD + theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())
  write.table(x=PRODFD, file="DFrequency.txt", sep="\t",quote=F,row.names=F,col.names=T)
  
  png("DPlot.png",width = 800, height = 600)
  print(pD)
  dev.off()
  
  ggsave("DPlot.pdf", pD, width=10, height=7)
}

pJ = ggplot(PRODFJ)
pJ = pJ + geom_bar( aes( x=factor(reorder(Top.J.Gene, chr.orderJ)), y=relFreq, fill=Sample), stat='identity', position="dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pJ = pJ + xlab("Summary of J gene") + ylab("Frequency") + ggtitle("Relative frequency of J gene usage") + scale_fill_manual(values=sample.colors)
pJ = pJ + theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())
write.table(x=PRODFJ, file="JFrequency.txt", sep="\t",quote=F,row.names=F,col.names=T)

png("JPlot.png",width = 800, height = 600)
pJ
dev.off()

ggsave("JPlot.pdf", pJ)

# ---------------------- Now the frequency plots of the V, D and J families ----------------------

print("Report Clonality - V, D and J family plots")

VGenes = PRODF[,c("Sample", "Top.V.Gene")]
VGenes$Top.V.Gene = gsub("-.*", "", VGenes$Top.V.Gene)
VGenes = data.frame(data.table(VGenes)[, list(Count=.N), by=c("Sample", "Top.V.Gene")])
TotalPerSample = data.frame(data.table(VGenes)[, list(total=sum(.SD$Count)), by=Sample])
VGenes = merge(VGenes, TotalPerSample, by="Sample")
VGenes$Frequency = VGenes$Count * 100 / VGenes$total
VPlot = ggplot(VGenes)
VPlot = VPlot + geom_bar(aes( x = Top.V.Gene, y = Frequency, fill = Sample), stat='identity', position='dodge' ) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Distribution of V gene families") + 
  ylab("Percentage of sequences") +
  scale_fill_manual(values=sample.colors) +
  theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())
png("VFPlot.png")
VPlot
dev.off()
ggsave("VFPlot.pdf", VPlot)

write.table(x=VGenes, file="VFFrequency.txt", sep="\t",quote=F,row.names=F,col.names=T)

if(useD){
  DGenes = PRODF[,c("Sample", "Top.D.Gene")]
  DGenes$Top.D.Gene = gsub("-.*", "", DGenes$Top.D.Gene)
  DGenes = data.frame(data.table(DGenes)[, list(Count=.N), by=c("Sample", "Top.D.Gene")])
  TotalPerSample = data.frame(data.table(DGenes)[, list(total=sum(.SD$Count)), by=Sample])
  DGenes = merge(DGenes, TotalPerSample, by="Sample")
  DGenes$Frequency = DGenes$Count * 100 / DGenes$total
  DPlot = ggplot(DGenes)
  DPlot = DPlot + geom_bar(aes( x = Top.D.Gene, y = Frequency, fill = Sample), stat='identity', position='dodge' ) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggtitle("Distribution of D gene families") + 
    ylab("Percentage of sequences") + 
    scale_fill_manual(values=sample.colors) +
    theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())
  png("DFPlot.png")
  print(DPlot)
  dev.off()
  
  ggsave("DFPlot.pdf", DPlot)
  write.table(x=DGenes, file="DFFrequency.txt", sep="\t",quote=F,row.names=F,col.names=T)
}

# ---------------------- Plotting the cdr3 length ----------------------

print("Report Clonality - CDR3 length plot")

CDR3Length = data.frame(data.table(PRODF)[, list(Count=.N), by=c("Sample", "CDR3.Length")])
TotalPerSample = data.frame(data.table(CDR3Length)[, list(total=sum(.SD$Count)), by=Sample])
CDR3Length = merge(CDR3Length, TotalPerSample, by="Sample")
CDR3Length$Frequency = CDR3Length$Count * 100 / CDR3Length$total
CDR3LengthPlot = ggplot(CDR3Length)
CDR3LengthPlot = CDR3LengthPlot + geom_bar(aes( x = factor(reorder(CDR3.Length, as.numeric(CDR3.Length))), y = Frequency, fill = Sample), stat='identity', position='dodge' ) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Length distribution of CDR3") + 
  xlab("CDR3 Length") + 
  ylab("Percentage of sequences") +
  scale_fill_manual(values=sample.colors) +
  theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())
png("CDR3LengthPlot.png",width = 1280, height = 720)
CDR3LengthPlot
dev.off()

ggsave("CDR3LengthPlot.pdf", CDR3LengthPlot, width=12, height=7)

write.table(x=CDR3Length, file="CDR3LengthPlot.txt", sep="\t",quote=F,row.names=F,col.names=T)

# ---------------------- Plot the heatmaps ----------------------

#get the reverse order for the V and D genes
revVchain = Vchain
revDchain = Dchain
revVchain$chr.orderV = rev(revVchain$chr.orderV)
revDchain$chr.orderD = rev(revDchain$chr.orderD)

if(useD){
  print("Report Clonality - Heatmaps VD")
  plotVD <- function(dat){
    if(length(dat[,1]) == 0){
      return()
    }
    
    img = ggplot() + 
      geom_tile(data=dat, aes(x=factor(reorder(Top.D.Gene, chr.orderD)), y=factor(reorder(Top.V.Gene, chr.orderV)), fill=relLength)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      scale_fill_gradient(low="gold", high="blue", na.value="white") + 
      ggtitle(paste(unique(dat$Sample), " (N=" , sum(dat$Length, na.rm=T) ,")", sep="")) + 
      xlab("D genes") + 
      ylab("V Genes") +
      theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), panel.grid.major = element_line(colour = "gainsboro"))
    
    png(paste("HeatmapVD_", unique(dat[3])[1,1] , ".png", sep=""), width=200+(15*length(Dchain$v.name)), height=100+(15*length(Vchain$v.name)))
    print(img)
    dev.off()
    
    ggsave(paste("HeatmapVD_", unique(dat[3])[1,1] , ".pdf", sep=""), img, height=13, width=8)
    
    write.table(x=acast(dat, Top.V.Gene~Top.D.Gene, value.var="Length"), file=paste("HeatmapVD_", unique(dat[3])[1,1], ".txt", sep=""), sep="\t",quote=F,row.names=T,col.names=NA)
  }
  
  VandDCount = data.frame(data.table(PRODF)[, list(Length=.N), by=c("Top.V.Gene", "Top.D.Gene", "Sample")])
  
  VandDCount$l = log(VandDCount$Length)
  maxVD = data.frame(data.table(VandDCount)[, list(max=max(l)), by=c("Sample")])
  VandDCount = merge(VandDCount, maxVD, by.x="Sample", by.y="Sample", all.x=T)
  VandDCount$relLength = VandDCount$l / VandDCount$max
  check = is.nan(VandDCount$relLength)
  if(any(check)){
	  VandDCount[check,"relLength"] = 0
  }

  completeVD = merge(VandDCount, revVchain, by.x="Top.V.Gene", by.y="v.name", all=TRUE)
  completeVD = merge(completeVD, Dchain, by.x="Top.D.Gene", by.y="v.name", all=TRUE)
  
  fltr = is.nan(completeVD$relLength)
  if(all(fltr)){
	  completeVD[fltr,"relLength"] = 0
  }
  
  VDList = split(completeVD, f=completeVD[,"Sample"])
  lapply(VDList, FUN=plotVD)
}

print("Report Clonality - Heatmaps VJ")

plotVJ <- function(dat){
  if(length(dat[,1]) == 0){
    return()
  }
  cat(paste(unique(dat[3])[1,1]))
  img = ggplot() + 
    geom_tile(data=dat, aes(x=factor(reorder(Top.J.Gene, chr.orderJ)), y=factor(reorder(Top.V.Gene, chr.orderV)), fill=relLength)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_gradient(low="gold", high="blue", na.value="white") + 
    ggtitle(paste(unique(dat$Sample), " (N=" , sum(dat$Length, na.rm=T) ,")", sep="")) + 
    xlab("J genes") + 
    ylab("V Genes") +
    theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), panel.grid.major = element_line(colour = "gainsboro"))
  
  png(paste("HeatmapVJ_", unique(dat[3])[1,1] , ".png", sep=""), width=200+(15*length(Jchain$v.name)), height=100+(15*length(Vchain$v.name)))
  print(img)
  dev.off()
  
  ggsave(paste("HeatmapVJ_", unique(dat[3])[1,1] , ".pdf", sep=""), img, height=11, width=4)
  
  write.table(x=acast(dat, Top.V.Gene~Top.J.Gene, value.var="Length"), file=paste("HeatmapVJ_", unique(dat[3])[1,1], ".txt", sep=""), sep="\t",quote=F,row.names=T,col.names=NA)
}



VandJCount = data.frame(data.table(PRODF)[, list(Length=.N), by=c("Top.V.Gene", "Top.J.Gene", "Sample")])

VandJCount$l = log(VandJCount$Length)
maxVJ = data.frame(data.table(VandJCount)[, list(max=max(l)), by=c("Sample")])
VandJCount = merge(VandJCount, maxVJ, by.x="Sample", by.y="Sample", all.x=T)
VandJCount$relLength = VandJCount$l / VandJCount$max

check = is.nan(VandJCount$relLength)
if(any(check)){
	VandJCount[check,"relLength"] = 0
}

completeVJ = merge(VandJCount, revVchain, by.x="Top.V.Gene", by.y="v.name", all=TRUE)
completeVJ = merge(completeVJ, Jchain, by.x="Top.J.Gene", by.y="v.name", all=TRUE)

fltr = is.nan(completeVJ$relLength)
if(any(fltr)){
	completeVJ[fltr,"relLength"] = 1
}

VJList = split(completeVJ, f=completeVJ[,"Sample"])
lapply(VJList, FUN=plotVJ)



if(useD){
  print("Report Clonality - Heatmaps DJ")	
  plotDJ <- function(dat){
    if(length(dat[,1]) == 0){
      return()
    }
    img = ggplot() + 
      geom_tile(data=dat, aes(x=factor(reorder(Top.J.Gene, chr.orderJ)), y=factor(reorder(Top.D.Gene, chr.orderD)), fill=relLength)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      scale_fill_gradient(low="gold", high="blue", na.value="white") + 
      ggtitle(paste(unique(dat$Sample), " (N=" , sum(dat$Length, na.rm=T) ,")", sep="")) + 
      xlab("J genes") + 
      ylab("D Genes") +
      theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), panel.grid.major = element_line(colour = "gainsboro"))
    
    png(paste("HeatmapDJ_", unique(dat[3])[1,1] , ".png", sep=""), width=200+(15*length(Jchain$v.name)), height=100+(15*length(Dchain$v.name)))
    print(img)
    dev.off()
    
    ggsave(paste("HeatmapDJ_", unique(dat[3])[1,1] , ".pdf", sep=""), img, width=4, height=7)
    
    write.table(x=acast(dat, Top.D.Gene~Top.J.Gene, value.var="Length"), file=paste("HeatmapDJ_", unique(dat[3])[1,1], ".txt", sep=""), sep="\t",quote=F,row.names=T,col.names=NA)
  }
  
  
  DandJCount = data.frame(data.table(PRODF)[, list(Length=.N), by=c("Top.D.Gene", "Top.J.Gene", "Sample")])
  
  DandJCount$l = log(DandJCount$Length)
  maxDJ = data.frame(data.table(DandJCount)[, list(max=max(l)), by=c("Sample")])
  DandJCount = merge(DandJCount, maxDJ, by.x="Sample", by.y="Sample", all.x=T)
  DandJCount$relLength = DandJCount$l / DandJCount$max
  
  check = is.nan(DandJCount$relLength)
  if(any(check)){
    DandJCount[check,"relLength"] = 0
  }
  
  cartegianProductDJ = expand.grid(Top.D.Gene = Dchain$v.name, Top.J.Gene = Jchain$v.name)
  
  completeDJ = merge(DandJCount, revDchain, by.x="Top.D.Gene", by.y="v.name", all=TRUE)
  completeDJ = merge(completeDJ, Jchain, by.x="Top.J.Gene", by.y="v.name", all=TRUE)
  
  fltr = is.nan(completeDJ$relLength)
  if(any(fltr)){
	  completeDJ[fltr, "relLength"] = 1
  }
  
  DJList = split(completeDJ, f=completeDJ[,"Sample"])
  lapply(DJList, FUN=plotDJ)
}


# ---------------------- output tables for the circos plots ----------------------

print("Report Clonality - Circos data")

for(smpl in unique(PRODF$Sample)){
	PRODF.sample = PRODF[PRODF$Sample == smpl,]
	
	fltr = PRODF.sample$Top.V.Gene == ""
	if(any(fltr, na.rm=T)){
	  PRODF.sample[fltr, "Top.V.Gene"] = "NA"
	}
	
	fltr = PRODF.sample$Top.D.Gene == ""
	if(any(fltr, na.rm=T)){
	  PRODF.sample[fltr, "Top.D.Gene"] = "NA"
	}

	fltr = PRODF.sample$Top.J.Gene == ""
	if(any(fltr, na.rm=T)){
	  PRODF.sample[fltr, "Top.J.Gene"] = "NA"
	}
	
	v.d = table(PRODF.sample$Top.V.Gene, PRODF.sample$Top.D.Gene)
	v.j = table(PRODF.sample$Top.V.Gene, PRODF.sample$Top.J.Gene)
	d.j = table(PRODF.sample$Top.D.Gene, PRODF.sample$Top.J.Gene)

	write.table(v.d, file=paste(smpl, "_VD_circos.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)
	write.table(v.j, file=paste(smpl, "_VJ_circos.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)
	write.table(d.j, file=paste(smpl, "_DJ_circos.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)
}

# ---------------------- calculating the clonality score ----------------------

if("Replicate" %in% colnames(inputdata)) #can only calculate clonality score when replicate information is available
{
  print("Report Clonality - Clonality")
  write.table(clonalityFrame, "clonalityComplete.txt", sep="\t",quote=F,row.names=F,col.names=T)
  if(clonality_method == "boyd"){
    samples = split(clonalityFrame, clonalityFrame$Sample, drop=T)
   
    for (sample in samples){
      res = data.frame(paste=character(0))
      sample_id = unique(sample$Sample)[[1]]
      for(replicate in unique(sample$Replicate)){
        tmp = sample[sample$Replicate == replicate,]
        clone_table = data.frame(table(tmp$clonaltype))
        clone_col_name = paste("V", replicate, sep="")
        colnames(clone_table) = c("paste", clone_col_name)
        res = merge(res, clone_table, by="paste", all=T)
      }
      
      res[is.na(res)] = 0
      
      write.table(res, file=paste("raw_clonality_", sample_id, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=F)
      write.table(as.matrix(res[,2:ncol(res)]), file=paste("raw_clonality2_", sample_id, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=F)
      
      res = read.table(paste("raw_clonality_", sample_id, ".txt", sep=""), header=F, sep="\t", quote="", stringsAsFactors=F, fill=T, comment.char="")
      
      infer.result = infer.clonality(as.matrix(res[,2:ncol(res)]))
      
      #print(infer.result)
      
      write.table(data.table(infer.result[[12]]), file=paste("lymphclon_clonality_", sample_id, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=F)
      
      res$type = rowSums(res[,2:ncol(res)])
      
      coincidence.table = data.frame(table(res$type))
      colnames(coincidence.table) = c("Coincidence Type",  "Raw Coincidence Freq")
      write.table(coincidence.table, file=paste("lymphclon_coincidences_", sample_id, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=T)
    }
  }
  clonalFreq = data.frame(data.table(clonalityFrame)[, list(Type=.N), by=c("Sample", "clonaltype")])

  #write files for every coincidence group of >1
  samples = unique(clonalFreq$Sample)
  for(sample in samples){
      clonalFreqSample = clonalFreq[clonalFreq$Sample == sample,]
      if(max(clonalFreqSample$Type) > 1){
          for(i in 2:max(clonalFreqSample$Type)){
              clonalFreqSampleType = clonalFreqSample[clonalFreqSample$Type == i,]
              clonalityFrame.sub = clonalityFrame[clonalityFrame$clonaltype %in% clonalFreqSampleType$clonaltype,]
              clonalityFrame.sub = clonalityFrame.sub[order(clonalityFrame.sub$clonaltype),]
              write.table(clonalityFrame.sub, file=paste("coincidences_", sample, "_", i, ".txt", sep=""), sep="\t",quote=F,row.names=F,col.names=T)
          }
      }
  }

  clonalFreqCount = data.frame(data.table(clonalFreq)[, list(Count=.N), by=c("Sample", "Type")])
  clonalFreqCount$realCount = clonalFreqCount$Type * clonalFreqCount$Count
  clonalSum = data.frame(data.table(clonalFreqCount)[, list(Reads=sum(realCount)), by=c("Sample")])
  clonalFreqCount = merge(clonalFreqCount, clonalSum, by.x="Sample", by.y="Sample")

  ct = c('Type\tWeight\n2\t1\n3\t3\n4\t6\n5\t10\n6\t15')
  tcct = textConnection(ct)
  CT  = read.table(tcct, sep="\t", header=TRUE)
  close(tcct)
  clonalFreqCount = merge(clonalFreqCount, CT, by.x="Type", by.y="Type", all.x=T)
  clonalFreqCount$WeightedCount = clonalFreqCount$Count * clonalFreqCount$Weight

  ReplicateReads = data.frame(data.table(clonalityFrame)[, list(Type=.N), by=c("Sample", "Replicate", "clonaltype")])
  ReplicateReads = data.frame(data.table(ReplicateReads)[, list(Reads=.N), by=c("Sample", "Replicate")])
  clonalFreqCount$Reads = as.numeric(clonalFreqCount$Reads)
  ReplicateReads$Reads = as.numeric(ReplicateReads$Reads)
  ReplicateReads$squared = as.numeric(ReplicateReads$Reads * ReplicateReads$Reads)

  ReplicatePrint <- function(dat){
      write.table(dat[-1], paste("ReplicateReads_", unique(dat[1])[1,1] , ".txt", sep=""), sep="\t",quote=F,na="-",row.names=F,col.names=F)
  }

  ReplicateSplit = split(ReplicateReads, f=ReplicateReads[,"Sample"])
  lapply(ReplicateSplit, FUN=ReplicatePrint)

  ReplicateReads = data.frame(data.table(ReplicateReads)[, list(ReadsSum=sum(as.numeric(Reads)), ReadsSquaredSum=sum(as.numeric(squared))), by=c("Sample")])
  clonalFreqCount = merge(clonalFreqCount, ReplicateReads, by.x="Sample", by.y="Sample", all.x=T)

  ReplicateSumPrint <- function(dat){
      write.table(dat[-1], paste("ReplicateSumReads_", unique(dat[1])[1,1] , ".txt", sep=""), sep="\t",quote=F,na="-",row.names=F,col.names=F)
  }

  ReplicateSumSplit = split(ReplicateReads, f=ReplicateReads[,"Sample"])
  lapply(ReplicateSumSplit, FUN=ReplicateSumPrint)

  clonalFreqCountSum = data.frame(data.table(clonalFreqCount)[, list(Numerator=sum(WeightedCount, na.rm=T)), by=c("Sample")])
  clonalFreqCount = merge(clonalFreqCount, clonalFreqCountSum, by.x="Sample", by.y="Sample", all.x=T)
  clonalFreqCount$ReadsSum = as.numeric(clonalFreqCount$ReadsSum) #prevent integer overflow
  clonalFreqCount$Denominator = (((clonalFreqCount$ReadsSum * clonalFreqCount$ReadsSum) - clonalFreqCount$ReadsSquaredSum) / 2)
  clonalFreqCount$Result = (clonalFreqCount$Numerator + 1) / (clonalFreqCount$Denominator + 1)

  ClonalityScorePrint <- function(dat){
      write.table(dat$Result, paste("ClonalityScore_", unique(dat[1])[1,1] , ".txt", sep=""), sep="\t",quote=F,na="-",row.names=F,col.names=F)
  }

  clonalityScore = clonalFreqCount[c("Sample", "Result")]
  clonalityScore = unique(clonalityScore)

  clonalityScoreSplit = split(clonalityScore, f=clonalityScore[,"Sample"])
  lapply(clonalityScoreSplit, FUN=ClonalityScorePrint)

  clonalityOverview = clonalFreqCount[c("Sample", "Type", "Count", "Weight", "WeightedCount")]



  ClonalityOverviewPrint <- function(dat){
      dat = dat[order(dat[,2]),]
      write.table(dat[-1], paste("ClonalityOverView_", unique(dat[1])[1,1] , ".txt", sep=""), sep="\t",quote=F,na="-",row.names=F,col.names=F)
  }

  clonalityOverviewSplit = split(clonalityOverview, f=clonalityOverview$Sample)
  lapply(clonalityOverviewSplit, FUN=ClonalityOverviewPrint)
  
}

bak = PRODF
bakun = UNPROD

imgtcolumns = c("X3V.REGION.trimmed.nt.nb","P3V.nt.nb", "N1.REGION.nt.nb", "P5D.nt.nb", "X5D.REGION.trimmed.nt.nb", "X3D.REGION.trimmed.nt.nb", "P3D.nt.nb", "N2.REGION.nt.nb", "P5J.nt.nb", "X5J.REGION.trimmed.nt.nb", "X3V.REGION.trimmed.nt.nb", "X5D.REGION.trimmed.nt.nb", "X3D.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb", "N1.REGION.nt.nb", "N2.REGION.nt.nb", "P3V.nt.nb", "P5D.nt.nb", "P3D.nt.nb", "P5J.nt.nb")
if(all(imgtcolumns %in% colnames(inputdata)))
{
  print("found IMGT columns, running junction analysis")
    
  #ensure certain columns are in the data (files generated with older versions of IMGT Loader)
  col.checks = c("N.REGION.nt.nb", "N1.REGION.nt.nb", "N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb")
  for(col.check in col.checks){
	  if(!(col.check %in% names(PRODF))){
		  print(paste(col.check, "not found adding new column"))
		  if(nrow(PRODF) > 0){ #because R is anoying...
			PRODF[,col.check] = 0
		  } else {
			PRODF = cbind(PRODF, data.frame(N3.REGION.nt.nb=numeric(0), N4.REGION.nt.nb=numeric(0)))
		  }
		  if(nrow(UNPROD) > 0){
			UNPROD[,col.check] = 0
		  } else {
			UNPROD = cbind(UNPROD, data.frame(N3.REGION.nt.nb=numeric(0), N4.REGION.nt.nb=numeric(0)))
		  }
	  }
  }
  
  PRODF.with.D = PRODF[nchar(PRODF$Top.D.Gene, keepNA=F) > 2,]
  PRODF.no.D = PRODF[nchar(PRODF$Top.D.Gene, keepNA=F) < 4,]
  write.table(PRODF.no.D, "productive_no_D.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=T)
  
  UNPROD.with.D = UNPROD[nchar(UNPROD$Top.D.Gene, keepNA=F) > 2,]
  UNPROD.no.D = UNPROD[nchar(UNPROD$Top.D.Gene, keepNA=F) < 4,]
  write.table(UNPROD.no.D, "unproductive_no_D.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=T)
  
  num_median = function(x, na.rm=T) { as.numeric(median(x, na.rm=na.rm)) }

  newData = data.frame(data.table(PRODF.with.D)[,list(unique=.N, 
                                               VH.DEL=mean(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                               P1=mean(.SD$P3V.nt.nb, na.rm=T),
                                               N1=mean(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb"), with=F], na.rm=T)),
                                               P2=mean(.SD$P5D.nt.nb, na.rm=T),
                                               DEL.DH=mean(.SD$X5D.REGION.trimmed.nt.nb, na.rm=T),
                                               DH.DEL=mean(.SD$X3D.REGION.trimmed.nt.nb, na.rm=T),
                                               P3=mean(.SD$P3D.nt.nb, na.rm=T),
                                               N2=mean(rowSums(.SD[,c("N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
                                               P4=mean(.SD$P5J.nt.nb, na.rm=T),
                                               DEL.JH=mean(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
                                               Total.Del=mean(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5D.REGION.trimmed.nt.nb", "X3D.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
                                               Total.N=mean(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb", "N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
                                               Total.P=mean(rowSums(.SD[,c("P3V.nt.nb", "P5D.nt.nb", "P3D.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
                                               Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
                                         by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisProd_mean_wD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
  
  newData = data.frame(data.table(PRODF.with.D)[,list(unique=.N, 
                                               VH.DEL=num_median(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                               P1=num_median(.SD$P3V.nt.nb, na.rm=T),
                                               N1=num_median(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb"), with=F], na.rm=T)),
                                               P2=num_median(.SD$P5D.nt.nb, na.rm=T),
                                               DEL.DH=num_median(.SD$X5D.REGION.trimmed.nt.nb, na.rm=T),
                                               DH.DEL=num_median(.SD$X3D.REGION.trimmed.nt.nb, na.rm=T),
                                               P3=num_median(.SD$P3D.nt.nb, na.rm=T),
                                               N2=num_median(rowSums(.SD[,c("N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
                                               P4=num_median(.SD$P5J.nt.nb, na.rm=T),
                                               DEL.JH=num_median(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
											   Total.Del=num_median(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5D.REGION.trimmed.nt.nb", "X3D.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
											   Total.N=num_median(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb", "N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
											   Total.P=num_median(rowSums(.SD[,c("P3V.nt.nb", "P5D.nt.nb", "P3D.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
											   Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
                                         by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisProd_median_wD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
  
  newData = data.frame(data.table(UNPROD.with.D)[,list(unique=.N, 
                                                VH.DEL=mean(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                                P1=mean(.SD$P3V.nt.nb, na.rm=T),
                                                N1=mean(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb"), with=F], na.rm=T)),
                                                P2=mean(.SD$P5D.nt.nb, na.rm=T),
                                                DEL.DH=mean(.SD$X5D.REGION.trimmed.nt.nb, na.rm=T),
                                                DH.DEL=mean(.SD$X3D.REGION.trimmed.nt.nb, na.rm=T),
                                                P3=mean(.SD$P3D.nt.nb, na.rm=T),
                                                N2=mean(rowSums(.SD[,c("N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
                                                P4=mean(.SD$P5J.nt.nb, na.rm=T),
                                                DEL.JH=mean(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
                                                Total.Del=mean(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5D.REGION.trimmed.nt.nb", "X3D.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
                                                Total.N=mean(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb", "N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
                                                Total.P=mean(rowSums(.SD[,c("P3V.nt.nb", "P5D.nt.nb", "P3D.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
                                                Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
                                          by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisUnProd_mean_wD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
  
    newData = data.frame(data.table(UNPROD.with.D)[,list(unique=.N, 
                                                VH.DEL=num_median(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                                P1=num_median(.SD$P3V.nt.nb, na.rm=T),
                                                N1=num_median(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb"), with=F], na.rm=T)),
                                                P2=num_median(.SD$P5D.nt.nb, na.rm=T),
                                                DEL.DH=num_median(.SD$X5D.REGION.trimmed.nt.nb, na.rm=T),
                                                DH.DEL=num_median(.SD$X3D.REGION.trimmed.nt.nb, na.rm=T),
                                                P3=num_median(.SD$P3D.nt.nb, na.rm=T),
                                                N2=num_median(rowSums(.SD[,c("N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
                                                P4=num_median(.SD$P5J.nt.nb, na.rm=T),
                                                DEL.JH=num_median(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
                                                Total.Del=num_median(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5D.REGION.trimmed.nt.nb", "X3D.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
                                                Total.N=num_median(rowSums(.SD[,c("N.REGION.nt.nb", "N1.REGION.nt.nb", "N2.REGION.nt.nb", "N3.REGION.nt.nb", "N4.REGION.nt.nb"), with=F], na.rm=T)),
                                                Total.P=num_median(rowSums(.SD[,c("P3V.nt.nb", "P5D.nt.nb", "P3D.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
                                                Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
															by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisUnProd_median_wD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
  
  #---------------- again for no-D
  
  newData = data.frame(data.table(PRODF.no.D)[,list(unique=.N, 
                                               VH.DEL=mean(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                               P1=mean(.SD$P3V.nt.nb, na.rm=T),
                                               N1=mean(.SD$N.REGION.nt.nb, na.rm=T),
                                               P2=mean(.SD$P5J.nt.nb, na.rm=T),
                                               DEL.JH=mean(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
                                               Total.Del=mean(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
                                               Total.N=mean(.SD$N.REGION.nt.nb, na.rm=T),
                                               Total.P=mean(rowSums(.SD[,c("P3V.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
                                               Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
                                         by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisProd_mean_nD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
  
  newData = data.frame(data.table(PRODF.no.D)[,list(unique=.N, 
                                               VH.DEL=num_median(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                               P1=num_median(.SD$P3V.nt.nb, na.rm=T),
                                               N1=num_median(.SD$N.REGION.nt.nb, na.rm=T),
                                               P2=num_median(.SD$P5J.nt.nb, na.rm=T),
                                               DEL.JH=num_median(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
											   Total.Del=num_median(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
											   Total.N=num_median(.SD$N.REGION.nt.nb, na.rm=T),
											   Total.P=num_median(rowSums(.SD[,c("P3V.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
											   Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
                                         by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisProd_median_nD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
  
  newData = data.frame(data.table(UNPROD.no.D)[,list(unique=.N, 
                                                VH.DEL=mean(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                                P1=mean(.SD$P3V.nt.nb, na.rm=T),
                                                N1=mean(.SD$N.REGION.nt.nb, na.rm=T),
                                                P2=mean(.SD$P5J.nt.nb, na.rm=T),
                                                DEL.JH=mean(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
                                                Total.Del=mean(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
                                                Total.N=mean(.SD$N.REGION.nt.nb, na.rm=T),
                                                Total.P=mean(rowSums(.SD[,c("P3V.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
                                                Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
                                          by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisUnProd_mean_nD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
  
  
    newData = data.frame(data.table(UNPROD.no.D)[,list(unique=.N, 
                                                VH.DEL=num_median(.SD$X3V.REGION.trimmed.nt.nb, na.rm=T),
                                                P1=num_median(.SD$P3V.nt.nb, na.rm=T),
                                                N1=num_median(.SD$N.REGION.nt.nb, na.rm=T),
                                                P2=num_median(.SD$P5J.nt.nb, na.rm=T),
                                                DEL.JH=num_median(.SD$X5J.REGION.trimmed.nt.nb, na.rm=T),
                                                Total.Del=num_median(rowSums(.SD[,c("X3V.REGION.trimmed.nt.nb", "X5J.REGION.trimmed.nt.nb"), with=F], na.rm=T)),
                                                Total.N=num_median(.SD$N.REGION.nt.nb, na.rm=T),
                                                Total.P=num_median(rowSums(.SD[,c("P3V.nt.nb", "P5J.nt.nb"), with=F], na.rm=T)),
                                                Median.CDR3.l=as.double(median(as.numeric(.SD$CDR3.Length), na.rm=T))),
															by=c("Sample")])
  newData[,sapply(newData, is.numeric)] = round(newData[,sapply(newData, is.numeric)],1)
  write.table(newData, "junctionAnalysisUnProd_median_nD.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)
}

PRODF = bak
UNPROD = bakun


# ---------------------- D reading frame ----------------------

D.REGION.reading.frame = PRODF[,c("Sample", "D.REGION.reading.frame")]

chck = is.na(D.REGION.reading.frame$D.REGION.reading.frame)
if(any(chck)){
	D.REGION.reading.frame[chck,"D.REGION.reading.frame"] = "No D"
}

D.REGION.reading.frame.1 = data.frame(data.table(D.REGION.reading.frame)[, list(Freq=.N), by=c("Sample", "D.REGION.reading.frame")])

D.REGION.reading.frame.2 = data.frame(data.table(D.REGION.reading.frame)[, list(sample.sum=sum(as.numeric(.SD$D.REGION.reading.frame), na.rm=T)), by=c("Sample")])

D.REGION.reading.frame = merge(D.REGION.reading.frame.1, D.REGION.reading.frame.2, by="Sample")

D.REGION.reading.frame$percentage = round(D.REGION.reading.frame$Freq / D.REGION.reading.frame$sample.sum * 100, 1)

write.table(D.REGION.reading.frame, "DReadingFrame.txt" , sep="\t",quote=F,row.names=F,col.names=T)

D.REGION.reading.frame = ggplot(D.REGION.reading.frame)
D.REGION.reading.frame = D.REGION.reading.frame + geom_bar(aes( x = D.REGION.reading.frame, y = percentage, fill=Sample), stat='identity', position='dodge' ) + ggtitle("D reading frame") + xlab("Frame") + ylab("Frequency")
D.REGION.reading.frame = D.REGION.reading.frame + scale_fill_manual(values=sample.colors)
D.REGION.reading.frame = D.REGION.reading.frame + theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())

png("DReadingFrame.png")
D.REGION.reading.frame
dev.off()

ggsave("DReadingFrame.pdf", D.REGION.reading.frame)

# ---------------------- AA composition in CDR3 ----------------------

AACDR3 = PRODF[,c("Sample", "CDR3.Seq")]

TotalPerSample = data.frame(data.table(AACDR3)[, list(total=sum(nchar(as.character(.SD$CDR3.Seq)))), by=Sample])

AAfreq = list()

for(i in 1:nrow(TotalPerSample)){
	sample = TotalPerSample$Sample[i]
  AAfreq[[i]] = data.frame(table(unlist(strsplit(as.character(AACDR3[AACDR3$Sample == sample,c("CDR3.Seq")]), ""))))
  AAfreq[[i]]$Sample = sample
}

AAfreq = ldply(AAfreq, data.frame)
AAfreq = merge(AAfreq, TotalPerSample, by="Sample", all.x = T)
AAfreq$freq_perc = as.numeric(AAfreq$Freq / AAfreq$total * 100)


AAorder = read.table(sep="\t", header=TRUE, text="order.aa\tAA\n1\tR\n2\tK\n3\tN\n4\tD\n5\tQ\n6\tE\n7\tH\n8\tP\n9\tY\n10\tW\n11\tS\n12\tT\n13\tG\n14\tA\n15\tM\n16\tC\n17\tF\n18\tL\n19\tV\n20\tI")
AAfreq = merge(AAfreq, AAorder, by.x='Var1', by.y='AA', all.x=TRUE)

AAfreq = AAfreq[!is.na(AAfreq$order.aa),]

AAfreqplot = ggplot(AAfreq)
AAfreqplot = AAfreqplot + geom_bar(aes( x=factor(reorder(Var1, order.aa)), y = freq_perc, fill = Sample), stat='identity', position='dodge' )
AAfreqplot = AAfreqplot + annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0, ymax = Inf, fill = "red", alpha = 0.2)
AAfreqplot = AAfreqplot + annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 0, ymax = Inf, fill = "blue", alpha = 0.2)
AAfreqplot = AAfreqplot + annotate("rect", xmin = 5.5, xmax = 6.5, ymin = 0, ymax = Inf, fill = "blue", alpha = 0.2)
AAfreqplot = AAfreqplot + annotate("rect", xmin = 6.5, xmax = 7.5, ymin = 0, ymax = Inf, fill = "red", alpha = 0.2)
AAfreqplot = AAfreqplot + ggtitle("Amino Acid Composition in the CDR3") + xlab("Amino Acid, from Hydrophilic (left) to Hydrophobic (right)") + ylab("Percentage") + scale_fill_manual(values=sample.colors)
AAfreqplot = AAfreqplot + theme(panel.background = element_rect(fill = "white", colour="black"),text = element_text(size=15, colour="black"), panel.grid.major.y = element_line(colour = "black"), panel.grid.major.x = element_blank())

png("AAComposition.png",width = 1280, height = 720)
AAfreqplot
dev.off()

ggsave("AAComposition.pdf", AAfreqplot, width=12, height=7)

write.table(AAfreq, "AAComposition.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=T)

# ---------------------- AA median CDR3 length ----------------------

median.aa.l = data.frame(data.table(PRODF)[, list(median=as.double(median(as.numeric(.SD$CDR3.Length, na.rm=T), na.rm=T))), by=c("Sample")])
write.table(median.aa.l, "AAMedianBySample.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=F)

if(clonaltype != "none"){
	#generate the "Sequences that are present in more than one replicate" dataset
	clonaltype.in.replicates = inputdata
	clonaltype.in.replicates = clonaltype.in.replicates[clonaltype.in.replicates$Functionality %in% c("productive (see comment)","productive"),]
	clonaltype.in.replicates = clonaltype.in.replicates[!(is.na(clonaltype.in.replicates$ID) | is.na(clonaltype.in.replicates$Top.V.Gene) | is.na(clonaltype.in.replicates$Top.J.Gene)),]
	clonaltype = unlist(strsplit(clonaltype, ","))

	clonaltype.in.replicates$clonaltype = do.call(paste, c(clonaltype.in.replicates[clonaltype], sep = ":"))

	clonaltype.in.replicates = clonaltype.in.replicates[!duplicated(clonaltype.in.replicates$clonaltype),]

	clonaltype = clonaltype[-which(clonaltype == "Sample")]

	clonaltype.in.replicates$clonaltype = do.call(paste, c(clonaltype.in.replicates[clonaltype], sep = ":"))
	clonaltype.in.replicates = clonaltype.in.replicates[,c("clonaltype","Replicate", "ID", "Sequence", "Sample")]


	write.table(clonaltype.in.replicates, "clonaltypes_replicates_before_table.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=T)

	clonaltype.counts = data.frame(table(clonaltype.in.replicates$clonaltype))

	write.table(clonaltype.counts, "clonaltypes_counts.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=T)

	names(clonaltype.counts) = c("clonaltype", "coincidence")

	clonaltype.counts = clonaltype.counts[clonaltype.counts$coincidence > 1,]

	clonaltype.in.replicates = clonaltype.in.replicates[clonaltype.in.replicates$clonaltype %in% clonaltype.counts$clonaltype,]
	clonaltype.in.replicates = merge(clonaltype.in.replicates, clonaltype.counts, by="clonaltype")
	clonaltype.in.replicates = clonaltype.in.replicates[order(-clonaltype.in.replicates$coincidence, clonaltype.in.replicates$clonaltype, clonaltype.in.replicates$Replicate),c("coincidence","clonaltype", "Sample", "Replicate", "ID", "Sequence")]


	write.table(clonaltype.in.replicates, "clonaltypes_replicates.txt" , sep="\t",quote=F,na="-",row.names=F,col.names=T)
} else {
	cat("No clonaltype", file="clonaltypes_replicates_before_table.txt")
	cat("No clonaltype", file="clonaltypes_counts.txt")
	cat("No clonaltype", file="clonaltypes_replicates.txt")
}























