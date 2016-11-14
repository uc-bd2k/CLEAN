`r2gtr` <-
function (hr, file = "cluster.gtr", dec = ".") {
	height <- hr$height
	if (range(height)[1] >= 0 & range(height)[2] <= 1) {
	##if (substr(distance, 1, 1) == "p") {
		height <- 1 - height
	} else {
		height <- height + 1
		height <- height[1]/height
	}
	##height <- signif(height, digits = digits)
	n <- length(height)
	node <- 1:n
	node <- paste("NODE", node, "X", sep = "")
	merge1 <- hr$merge[, 1]
	merge11 <- paste("NODE", merge1, "X", sep = "")
	merge12 <- paste("GENE", -merge1, "X", sep = "")
	merge1[hr$merge[, 1] > 0] <- merge11[hr$merge[, 1] > 0]
	merge1[hr$merge[, 1] < 0] <- merge12[hr$merge[, 1] < 0]
	merge2 <- hr$merge[, 2]
	merge11 <- paste("NODE", merge2, "X", sep = "")
	merge12 <- paste("GENE", -merge2, "X", sep = "")
	merge2[hr$merge[, 2] > 0] <- merge11[hr$merge[, 2] > 0]
	merge2[hr$merge[, 2] < 0] <- merge12[hr$merge[, 2] < 0]
	data <- data.frame(cbind(node, merge1, merge2))
	data <- cbind(data, height)
	write.table(data, file = file, row.name = FALSE, col.names = FALSE,
		quote = FALSE, sep = "\t", dec = dec)
}

