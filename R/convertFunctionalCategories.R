`convertFunctionalCategories` <-   
function (functionalCategories, fromSpecies = "h", toSpecies = "r", homologTable = "./homologene.data")
{	
	#this function was adapted from convertSmc() of the PGSEA package
	
	if(fromSpecies != toSpecies) {
		if(!file.exists(homologTable)) {
			warning(paste("File", homologTable, "not found. Please download from ftp://ftp.ncbi.nih.gov (see ?convertFunctionalCategories for details)."))
			#download.file("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data",destfile="homologene.data",mode="wb")
			NULL
		} else {
			fromSpecies <- switch(fromSpecies, h = 9606, r = 10116, m = 10090)
			toSpecies <- switch(toSpecies, h = 9606, r = 10116, m = 10090)
			homologTable <- read.delim(homologTable, col.names = c("HID", "TaxID", "GeneID", "Symbol", "Protein", "Protein-accession"))[, 1:3]
			hs <- which(homologTable$TaxID %in% fromSpecies)
			rn <- which(homologTable$TaxID %in% toSpecies)
			from <- homologTable[hs, ]
			to <- homologTable[rn, ]
			for (z in 1:length(functionalCategories)) {
				matched <- functionalCategories[[z]]
				ids <- to[match(from[match(matched, from[, "GeneID"]), "HID"], to[, "HID"]), "GeneID"]
				functionalCategories[[z]] <- unique(ids[!is.na(ids)])
				if (z%%50 == 0)	cat("finished ", z, "of", length(functionalCategories), "\n")
			}
		}
	}
	functionalCategories
}

