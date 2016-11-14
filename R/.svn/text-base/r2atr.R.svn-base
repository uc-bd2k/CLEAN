`r2atr` <-
function (hc, file = "cluster.atr", dec = ".") {
	height <- hc$height
	if (range(height)[1] >= 0 & range(height)[2] <= 1) {
		height <- 1 - height
	} else {
		height <- height + 1
		height <- height[1]/height
	}
	##height <- signif(height, digits = digits)
	n <- length(height)
	node <- 1:n
	node <- paste("NODE", node, "X", sep = "")
	merge1 <- hc$merge[, 1]
	merge11 <- paste("NODE", merge1, "X", sep = "")
	merge12 <- paste("ARRY", -merge1, "X", sep = "")
	merge1[hc$merge[, 1] > 0] <- merge11[hc$merge[, 1] > 0]
	merge1[hc$merge[, 1] < 0] <- merge12[hc$merge[, 1] < 0]
	merge2 <- hc$merge[, 2]
	merge11 <- paste("NODE", merge2, "X", sep = "")
	merge12 <- paste("ARRY", -merge2, "X", sep = "")
	merge2[hc$merge[, 2] > 0] <- merge11[hc$merge[, 2] > 0]
	merge2[hc$merge[, 2] < 0] <- merge12[hc$merge[, 2] < 0]
	data <- data.frame(cbind(node, merge1, merge2))
	data <- cbind(data, height)
	if (dec == ",") {
		data <- apply(data, 2, function(u) {chartr(".", ",", u)})
		data[data == "NA"] <- NA
	}
	write.table(data, file = file, row.name = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

