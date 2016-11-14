`getBinarySubTree` <-
function(tree,right=FALSE) {
	getRow <- function(merge,component,rows) {
		toDoList <- component
		doneList <- NULL
		while(length(toDoList)>0) {
			doneList <- c(toDoList[1], doneList)
			component <- merge[toDoList[1], 2:1]
			toDoList <- c(component[component>0], toDoList[-1])
		}
		doneList
	}
	n<-2-right
	component<-tree$merge[length(tree$merge)/n]
	if(component < 0) {
		element <-which(tree$merge==component)
		if(element <= (length(tree$merge)/2)){
			row<-element
		}else {
			row<-element-length(tree$merge)/2
		}
		res <- list(merge = NULL, 
					heights = tree$height[row],
					order = 1,
					labels = tree$labels[-component],
					method = tree$method,
					call = tree$call,
					dist.method = tree$dist.method)
	} else {
		which.merge.rows<-sort(getRow(tree$merge,component,NULL))
		merge.list <- tree$merge[which.merge.rows,]
		merge.height <- tree$height[which.merge.rows]
		pos.loc<-which(merge.list > 0)
		pos.elements <- sort(merge.list[pos.loc], index.return=TRUE)
		rank<-sort(pos.elements[[2]], index.return=TRUE)
		merge.list[pos.loc]<-rank[[2]]

		neg.loc<-which(merge.list < 0)
		neg.elements <- merge.list[neg.loc]
		minus.neg.elements <- -neg.elements
		neg.elements<-sort(minus.neg.elements, index.return=TRUE)
		rank<-sort(neg.elements[[2]], index.return=TRUE)
		merge.list[neg.loc]<- -rank[[2]]

		which.nodes <-sort(minus.neg.elements)
		names(which.nodes) <- NULL
		subtree.labels <- tree$labels[which.nodes]
		subtree.order <- tree$order[tree$labels[tree$order] %in% subtree.labels]
		
		old.order <- subtree.order
		subtree.order <- match(tree$labels[old.order], subtree.labels)
		res <- list(merge = merge.list,
					height = merge.height,
					order = subtree.order,
					labels = subtree.labels,
					method = tree$method,
					call = list(match.call(), tree$call),
					dist.method = tree$dist.method)
	}
	class(res) <- "hclust"
	res
}

