#`.First.lib` <-
`.onLoad` <-
function(lib, pkg)
{
    #.fTreeViewPath <- paste(.Library, "/CLEAN/doc/", sep = "")
    #    .fTreeViewPath <- paste(system.file(package = "CLEAN"), "/doc/", sep = "")
    #    assign(".fTreeViewPath", paste(.Library, "/CLEAN/doc/", sep = ""), envir = sys.frame())
    for(cpath in .libPaths()){
        cfpath<-paste(cpath, "/CLEAN/doc/", sep="")
        if(file.exists(paste(cfpath,"LFTreeView.jar", sep="")) | file.exists(paste(cfpath,"TreeView.jar", sep="")))
            assign(".fTreeViewPath",cfpath,envir=sys.frame())
}

