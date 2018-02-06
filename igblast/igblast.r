args <- commandArgs(trailingOnly = TRUE)

infile=args[1]
outfile=args[2]

blasted = read.table(infile, header=T, sep="\t", fill=T, stringsAsFactors=F, comment.char="")

blasted$ID = 1:nrow(blasted)
blasted$VDJ.Frame = "Out-of-frame"

search = blasted$inFrame == "true" & blasted$noStop == "false"
if(sum(search) > 0){
  blasted[search ,]$VDJ.Frame = "In-frame with stop codon"
}

search = blasted$inFrame == "true" & blasted$noStop == "true"
if(sum(search) > 0){
  blasted[search ,]$VDJ.Frame = "In-frame"
}

blasted$Top.V.Gene = blasted$vSegment
blasted$Top.D.Gene = blasted$dSegment
blasted$Top.J.Gene = blasted$jSegment
blasted$CDR1.Seq = blasted$cdr1aa
blasted$CDR1.Length = nchar(blasted$CDR1.Seq)
blasted$CDR2.Seq = blasted$cdr2aa
blasted$CDR2.Length = nchar(blasted$CDR2.Seq)
blasted$CDR3.Seq = blasted$cdr3aa
blasted$CDR3.Length = nchar(blasted$CDR3.Seq)
blasted$CDR3.Seq.DNA = blasted$cdr3nt
blasted$CDR3.Length.DNA = nchar(blasted$CDR3.Seq.DNA)
blasted$Strand = "+/-"
blasted$CDR3.Found.How = "found"

search = blasted$cdr3nt == ""
if(sum(search) > 0){
  blasted[search,]$CDR3.Found.How = "NOT_FOUND"
}

blasted$AA.JUNCTION = blasted$CDR3.Seq

n = c("X.reads_count", "ID", "VDJ.Frame", "Top.V.Gene", "Top.D.Gene", "Top.J.Gene", "CDR1.Seq", "CDR1.Length", "CDR2.Seq", "CDR2.Length", "CDR3.Seq", "CDR3.Length", "CDR3.Seq.DNA", "CDR3.Length.DNA", "Strand", "CDR3.Found.How", "Functionality", "AA.JUNCTION")

n[!(n %in% names(blasted))]

blasted = blasted[,c("X.reads_count", "ID", "VDJ.Frame", "Top.V.Gene", "Top.D.Gene", "Top.J.Gene", "CDR1.Seq", "CDR1.Length", "CDR2.Seq", "CDR2.Length", "CDR3.Seq", "CDR3.Length", "CDR3.Seq.DNA", "CDR3.Length.DNA", "Strand", "CDR3.Found.How", "AA.JUNCTION")]

names(blasted) = c("frequency.count", "ID", "VDJ Frame", "Top V Gene", "Top D Gene", "Top J Gene", "CDR1 Seq", "CDR1 Length", "CDR2 Seq", "CDR2 Length", "CDR3 Seq", "CDR3 Length", "CDR3 Seq DNA", "CDR3 Length DNA", "Strand", "CDR3 Found How", "AA JUNCTION")

#duplicate rows based on frequency.count
blasted = blasted[rep(seq_len(nrow(blasted)), blasted$frequency.count),]
blasted$ID = 1:nrow(blasted)

blasted = blasted[,c("ID", "VDJ Frame", "Top V Gene", "Top D Gene", "Top J Gene", "CDR1 Seq", "CDR1 Length", "CDR2 Seq", "CDR2 Length", "CDR3 Seq", "CDR3 Length", "CDR3 Seq DNA", "CDR3 Length DNA", "Strand", "CDR3 Found How", "AA JUNCTION")]

write.table(blasted, outfile, quote=F, sep="\t", row.names=F, col.names=T)
