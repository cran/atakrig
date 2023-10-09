## atakrig
## Author: Maogui Hu.


## discretize raster to fine resolution points. ----
# Input:
#   x: SpatRaster
#   cellsize: discretized fine grid size.
#   type: "value", "nodata", "all", whether only valid pixels, or only NODATA pixles, or all pixels extracted.
#   psf: PSF type, "equal", "gau", or user defined PSF matrix (normalized).
#   sigma: standard deviation for Gaussian PSF.
# Output: list(areaValues, discretePoints)
#   areaValues: original pixel coordinates and value, data.frame(areaId,centx,centy,value).
#   discretePoints: discretized pixels points and weight, data.frame(areaId,ptx,pty,weight).
discretizeRaster <- function(x, cellsize, type="value", psf="equal", sigma=2) {
  if(!is(x, "SpatRaster")) {
    stop("x isn't SpatRaster object!")
  }

  p <- as.points(x, na.rm=FALSE)
  p <- cbind(crds(p), values(p))
  names(p) <- c("x","y","value")
  p$value[is.nan(p$value)] <- NA

  type <- tolower(type)
  if(type == "value") {
    p <- subset(p, !is.na(p$value))
  } else if(type == "nodata") {
    p <- subset(p, is.na(p$value))
  } else if(type == "all"){

  } else {
    error("invalid type.")
  }

  p <- data.frame(areaId=1:nrow(p), p)

  blkResoX <- xres(x)
  blkResoY <- yres(x)
  xResoNum <- round(blkResoX / cellsize)
  yResoNum <- round(blkResoY / cellsize)

  if(is.numeric(psf)) {
    if(!all(dim(psf) == c(yResoNum, xResoNum))) {
      stop("Incorrect psf dimension!")
    }
    w <- psf
  } else if(psf == "equal") {
    w <- 1/(xResoNum*yResoNum)
  } else if(psf == "gau") {
    xy <- expand.grid(x=(1:xResoNum)-ceiling(xResoNum/2), y=(1:yResoNum)-ceiling(yResoNum/2))
    w <- 1/(2 * pi * sigma^2) * exp(-(xy$x^2 + xy$y^2)/(2 * sigma^2))
    # w <- matrix(w/sum(w), ncol = xResoNum, byrow = TRUE)
    w <- w/sum(w)
  } else {
    stop("Incorrect psf parameter!")
  }

  discretePoints <- lapply(1:nrow(p), function(i) {
    xy <- expand.grid(list(ptx = p$x[i] - 0.5 * blkResoX + (1:xResoNum-0.5) * cellsize,
                           pty = p$y[i] - 0.5 * blkResoY + (1:yResoNum-0.5) * cellsize))
    data.frame(areaId = p$areaId[i], xy, weight=w)
  })
  discretePoints <- do.call(rbind, discretePoints)

  colnames(p)[2:4] <- c("centx","centy","value")
  rslt <- list(areaValues = p, discretePoints = discretePoints)
  class(rslt) <- c("list", "discreteArea")

  return(rslt)
}


discretizePolygon <- function(x, cellsize, id=NULL, value=NULL, showProgressBar=FALSE) {

  if(is.null(id)) {
    id <- "_id_"
    x$`_id_` <- 1:nrow(x)
  }
  if(is.null(value)) {
    value <- "_value_"
    x$`_value_` <- NA
  }

  regularsample <- function(x, cellsize, iter=0) {
    bb <- st_bbox(x)
    xs <- tryCatch(seq(bb[1] + 0.5*cellsize, bb[3], by=cellsize), error = function(e) mean(bb[c(1,3)]))
    ys <- tryCatch(seq(bb[2] + 0.5*cellsize, bb[4], by=cellsize), error = function(e) mean(bb[c(2,4)]))
    xy <- expand.grid(xs, ys)
    names(xy) <- c("x","y")

    crs <- st_crs(x)
    pts <- st_filter(st_as_sf(xy, coords = c("x","y"), crs=crs), x)
    if(nrow(pts) == 0 && iter <= 0) iter <- 10

    if(iter > 0) {
      frac <- 0.5/(iter + 1)
      n <- nrow(xy)
      for (i in 1:iter) {
        offsetxy <- xy + frac * i * data.frame(x=rep(cellsize,n), y=rep(cellsize,n))
        pts0 <- st_filter(st_as_sf(offsetxy, coords = c("x","y"), crs=crs), x)
        if(nrow(pts0) > nrow(pts)) pts <- pts0

        offsetxy <- xy + frac * i * data.frame(x=rep(cellsize,n), y=-rep(cellsize,n))
        pts0 <- st_filter(st_as_sf(offsetxy, coords = c("x","y"), crs=crs), x)
        if(nrow(pts0) > nrow(pts)) pts <- pts0

        offsetxy <- xy + frac * i * data.frame(x=-rep(cellsize,n), y=rep(cellsize,n))
        pts0 <- st_filter(st_as_sf(offsetxy, coords = c("x","y"), crs=crs), x)
        if(nrow(pts0) > nrow(pts)) pts <- pts0

        offsetxy <- xy + frac * i * data.frame(x=-rep(cellsize,n), y=-rep(cellsize,n))
        pts0 <- st_filter(st_as_sf(offsetxy, coords = c("x","y"), crs=crs), x)
        if(nrow(pts0) > nrow(pts)) pts <- pts0
      }
    }

    return(pts)
  }

  # if(class(x[[id]]) == "factor") {
    x[[id]] <- as.character(x[[id]])
  # }

  if(showProgressBar) pb <- txtProgressBar(max = nrow(x), width = 50, style = 3)
  discretePoints <- lapply(1:nrow(x), function(i) {
    if(showProgressBar) setTxtProgressBar(pb, i)
    pts <- regularsample(x[i,], cellsize = cellsize)
    xys <- cbind(data.frame(areaId=x[[i,id]]), st_coordinates(pts)[,c("X", "Y"), drop=FALSE])
    xys$weight <- 1/nrow(xys)
    return(xys)
  })
  if(showProgressBar) close(pb)
  discretePoints <- do.call(rbind, discretePoints)
  names(discretePoints)[2:3] <- c("ptx", "pty")

  areaValues <- data.frame(areaId=x[[id]], value=x[[value]])

  centxy <- st_centroid(x)
  areaValues <- cbind(areaValues, st_coordinates(centxy))
  names(areaValues)[3:4] <- c("centx", "centy")
  areaValues <- areaValues[,c("areaId","centx","centy","value")]

  if(is(areaValues$areaId, "factor")) areaValues$areaId <- as.character(areaValues$areaId)
  if(is(discretePoints$areaId, "factor")) discretePoints$areaId <- as.character(discretePoints$areaId)

  rslt <- list(areaValues = areaValues, discretePoints = discretePoints)
  class(rslt) <- c("list", "discreteArea")

  return(rslt)
}


isValidDiscreteAreaObj <- function(x){
  if(!all(sort(names(x))) == c("areaValues", "discretePoints")) return(FALSE)
  if(!all(names(x$areaValues) == c("areaId","centx","centy","value"))) return(FALSE)
  if(!all(names(x$discretePoints) == c("areaId","ptx","pty","weight"))) return(FALSE)
  return(TRUE)
}


updateDiscreteAreaValue <- function(x, newval) {
  for (i in 1:nrow(newval)) {
    indx <- x$areaValues$areaId == newval$areaId
    x$areaValues$value[indx] <- newval$value[i]
  }
  return(x)
}


## select discretized data by area id. ----
subsetDiscreteArea <- function(x, selAreaId, revSel=FALSE) {
  if(!revSel) {
    rslt <- list(areaValues = x$areaValues[x$areaValues[,1] %in% selAreaId,,drop=FALSE],
         discretePoints = x$discretePoints[x$discretePoints[,1] %in% selAreaId,,drop=FALSE])
  } else {
    rslt <- list(areaValues = x$areaValues[!x$areaValues[,1] %in% selAreaId,,drop=FALSE],
         discretePoints = x$discretePoints[!x$discretePoints[,1] %in% selAreaId,,drop=FALSE])
  }

  class(rslt) <- c("list", "discreteArea")
  return(rslt)
}


## rbindDiscreteArea ----
# Input:
#   x, y: discretized area, list(areaValues, discretePoints):
#       areaValues: sample values, data.frame(areaId,centx,centy,value).
#       discretePoints: discretized area-samples, data.frame(areaId,ptx,pty,weight), weight is normalized.
rbindDiscreteArea <- function(x, y) {
  if(any(match(unique(y$areaValues$areaId), unique(x$areaValues$areaId)))) {
    x$areaValues$areaId <- paste0("x_", x$areaValues$areaId)
    x$discretePoints$areaId <- paste0("x_", x$discretePoints$areaId)

    y$areaValues$areaId <- paste0("y_", y$areaValues$areaId)
    y$discretePoints$areaId <- paste0("y_", y$discretePoints$areaId)
  }

  xy <- list(areaValues = rbind(x$areaValues, y$areaValues),
             discretePoints = rbind(x$discretePoints, y$discretePoints))
  class(xy) <- c("list", "discreteArea")
  return(xy)
}


