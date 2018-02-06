args <- commandArgs(trailingOnly = TRUE)

print(args)

inputs = args[1:(length(args) - 1)]
output = args[length(args)]

current.id = ""
counter = 1

result = NULL

for(current in inputs){
	if(grepl("/", current)){ #its a path to a file
		print(paste("Adding file", counter, "to", current.id))
		dat = read.table(current, sep="\t", header=T, quote="", fill=T)
		
		if(nrow(dat) == 0){
			print(paste(counter, "of", current.id, "has no sequences, skipping"))
			next
		}
		
		#IMGT check
		
		dat$Sample = current.id
		dat$Replicate = counter
		
		if(is.null(result)){
			result = dat[NULL,]
		}
		
		result = rbind(result, dat)
		
		counter = counter + 1
		
	} else { #its an ID of a patient
		print(paste("New patient", current))
		current.id = current
		counter = 1
	}
}

write.table(result, output, sep="\t", quote=F, row.names=F, col.names=T)
