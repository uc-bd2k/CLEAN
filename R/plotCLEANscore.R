`plotCLEANscore` <- 
function(fClustAnnotations, fCategoryName = "GO", maxCLEANscore = NULL, add = FALSE, ...) 
{
	index <- grep(fCategoryName, names(fClustAnnotations))
	if(is.null(maxCLEANscore)) { 
		mF <- round(max(fClustAnnotations[[index]]$CLEANscore)) 
	} else { 
		mF <- maxCLEANscore
	}
	if(add) {
		lines(1:mF, sapply(1:mF, function(x) sum(fClustAnnotations[[index]]$CLEANscore > x, na.rm = TRUE)), ...)
	} else {
		plot(1:mF, sapply(1:mF, function(x) sum(fClustAnnotations[[index]]$CLEANscore > x, na.rm = TRUE)), log = "y", xlab = "CLEAN score",
			ylab = "# genes with CLEAN score <= threshold", type = "l", ...)
	}
}
