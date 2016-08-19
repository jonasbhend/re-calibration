colkd.ind.cold.v1 <-
function (type="frost") {

    if (type=="frost" | type=="ice") {

      r <- c(  14,  43, 171, 232, 215, 215, 201, 138,  51,  28 )
      g <- c(  97, 172, 207, 246, 230, 230, 235, 207, 166,  84 )
      b <- c(  52, 102,  99, 158, 196, 227, 245, 255, 232, 171 )

    } else {

      stop("Argument <type> in colkd.ind.cold(type=?) must be either frost or ice.")

    }

    rgb(red=r,green=g,blue=b,maxColorValue=255)

}
