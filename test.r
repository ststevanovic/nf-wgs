# convert list of 2 into hash
conver <- function(x) {
  h <- list()
  for (i in 1:length(x)) {
    h[[i]] <- x[[i]]
  }
  return(h)
}
convert <- function(x) {
  if (is.null(x)) {
    return(NULL)
  } else {
    return(list(x[,1] = x[,2]))
  }
}