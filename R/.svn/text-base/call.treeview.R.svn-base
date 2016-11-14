`call.treeview` <-
function(cdtFile = "cluster.cdt", treeViewPath = .fTreeViewPath, wait = FALSE, ilincs=F, removeConfig=F) {
	removeConfigString=NULL
	if(removeConfig) removeConfigString=" -overwrite_global_config"
	if(Sys.info()["sysname"] == "Windows") {
		if(substr(path.expand(cdtFile), 2, 2) == ":") {   #if expands to absolute path
			cdtFile <- path.expand(cdtFile)               #explicitly make path absolute
		} else {                                          #else is relative file path
			cdtFile <- paste(getwd(), cdtFile, sep = "/") #explicitly make absolute
		}
	} else {
		if(substr(path.expand(cdtFile), 1, 1) == "/") { 
			cdtFile <- path.expand(cdtFile)
		} else {
			cdtFile <- paste(getwd(), cdtFile, sep = "/") 
		}
	}
	if(file.exists(cdtFile)){
		if(ilincs)
			system(paste("java -jar ",  treeViewPath, "LFTreeView.jar -new -r ", cdtFile,removeConfigString," -ilincs", sep = ""), wait = wait)
		else
			system(paste("java -jar ",  treeViewPath, "TreeView.jar -new -r ", cdtFile, sep = ""), wait = wait)
	}
	else if(file.exists(paste(cdtFile, "cdt", sep = ".")))	{
		if(ilincs)
			system(paste("java -jar ",  treeViewPath, "LFTreeView.jar -new -r ", paste(cdtFile, "cdt", sep = "."),removeConfigString," -ilincs", sep = ""), wait = wait)
		else
			system(paste("java -jar ",  treeViewPath, "TreeView.jar -new -r ", paste(cdtFile, "cdt", sep = "."), sep = ""), wait = wait)
	}
	else
		warning("The cdt file was not found. Aborted call to fTreeView")
}

