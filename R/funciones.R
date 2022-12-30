dyn.load("weibull_nr.dll")

wnr <- function(y, d, init, niter = 30, tol = 0.0001, hist = TRUE) {
    n <- length(y)
    codigo <- 0
    xfin <- c(0.00001, 0.0001)
    vcov <- c(0.01, 0.01, 0.01, 0.01)
    salida <- .C("weibull_nr",
        x0 = as.double(init),
        ops = as.double(c(niter, tol)),
        n = as.integer(n),
        ti = as.double(y),
        di = as.double(d),
        cod = as.integer(codigo),
        xfin = as.double(xfin),
        vcov = as.double(vcov),
        hist = as.integer(hist)
    )
    if (salida$cod == 1) {
       return(salida)
    } else {
       cat("El algoritmo no ha convengido\n")
    }
}
dyn.unload("weibull_nr.dll")