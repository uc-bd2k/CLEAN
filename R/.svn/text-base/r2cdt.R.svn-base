`r2cdt` <-
function (hr, hc, data, labels = TRUE, description = TRUE, file = "cluster.cdt", dec = ".") {
	n <- dim(data)[1]
	data <- as.data.frame(cbind(1, data))
	if (!description) {
		if (!labels) {
			data <- cbind(row.names(data), data)
		} else {
			data <- cbind(as.factor(data[, 2]), data)
		}
	} else {
		data <- cbind(as.factor(data[, 3]), data[, -3])
	}
	if (!labels) {
		data <- cbind(row.names(data), data)
	} else {
		data <- cbind(as.factor(data[, 3]), data[, -3])
	}
	GID <- paste("GENE", 1:n, "X", sep = "")
	data <- cbind(as.factor(GID), data)
	data <- data[hr$order, ]
	m <- dim(data)[2]
	data[, 5:m] <- data[, 5:m][hc$order]
	colnames(data)[5:m] <- colnames(data)[5:m][hc$order]
	##data[, 5:m] <- signif(data[, 5:m], digits = 4)
	levels(data[, 1]) <- c(levels(data[, 1]), "1", "GID", "AID", "EWEIGHT")
	levels(data[, 2]) <- c(levels(data[, 2]), "1", "UNIQID")
	levels(data[, 3]) <- c(levels(data[, 3]), "1", "NAME")
	data <- rbind(1, data)
	nom <- c("GID", "UNIQID", "NAME", "GWEIGHT", names(data)[-c(1,2, 3, 4)])
	names(data) <- nom
	data[1, 1] <- "EWEIGHT"
	data <- rbind(c("AID", NA, NA, NA, paste("ARRY", hc$order, "X", sep = "")), data)
	data[2, 2:4] <- NA
	if (dec == ",") {
		data <- apply(data, 2, function(u) { chartr(".", ",", u) })
		data[data == "NA"] <- NA
	}
	write.table(data, file = file, sep = "\t", row.names = FALSE, col.names = TRUE, na = "", quote = FALSE)
}

