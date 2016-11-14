`convertGeneTable` <-   
function (geneTable, fromSpecies = "h", toSpecies = "r", homologTable = "./homologene.data")
{	
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
			toIDs <- function(x) { 
				ids <- to[to[,"HID"]==from[x[1]==from[,"GeneID"],"HID"],"GeneID"]
				if(length(ids) < 1) NA else ids
			}
			geneTable2 <- matrix(unlist(apply(geneTable, 1, function(x) {
				toids <- toIDs(x)
				t(cbind(toIDs(x), matrix(rep(x[-1], each=length(toids)), nrow=length(toids))))})), byrow=TRUE, ncol=ncol(geneTable))
			colnames(geneTable2) <- colnames(geneTable)
			geneTable <- geneTable2[!is.na(geneTable2[,1]), ]
		}
	}
	geneTable
}

