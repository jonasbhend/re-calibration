`colkd.na` <-
function (param=NULL) { 

    if (is.null(param)) {
        return(rgb(red=192,green=192,blue=192,maxColorValue=255)) }

    if (param == "sun.anom") {
        return(rgb(red=236,green=208,blue=252,maxColorValue=255)) 
    } else {
        return(rgb(red=192,green=192,blue=192,maxColorValue=255)) 
    }

}

