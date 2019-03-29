
N50 <- function(lengths) {
	lengths.sorted <- as.numeric(sort(lengths, decreasing = T))
	target <- .5 * sum(lengths.sorted)
    running_total <- 0
	for (length in lengths.sorted) {
		running_total <- running_total + length
		if (running_total >= target) {
			return(length)
		}
	}
}

pacbio <- read.table("pb.length.txt")
print("PacBio Metrics:")
print(summary(pacbio))
cat(paste0(" N50\t: ",N50(pacbio[[1]]), "\n\n"))

ont <- read.table("ont.length.txt")
print("ONT metrics:")
print(summary(ont))
cat(paste0(" N50\t: ",N50(ont[[1]]), "\n"))
