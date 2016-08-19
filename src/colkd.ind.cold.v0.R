colkd.ind.cold.v0 <-
function (type="frost") {

    if (type=="frost") {

      r <- c(  14,  43, 171, 232, 215, 181, 143, 127,   5,   2 )
      g <- c(  97, 172, 207, 246, 227, 202, 179, 151, 112,  56 )
      b <- c(  52, 102,  99, 158, 238, 255, 255, 255, 176,  88 )

    } else if (type=="ice") {

      r <- c( 230, 176, 128,  49,  17,  46,  71,  92, 135, 202)
      g <- c( 255, 255, 209, 162, 135, 105, 145, 161, 194, 238)
      b <- c( 140, 112,  54,  36,  63, 120, 181, 209, 217, 237)

    } else {

      stop("Argument <type> in colkd.ind.cold(type=?) must be either frost or ice.")

    }

    rgb(red=r,green=g,blue=b,maxColorValue=255)

}
