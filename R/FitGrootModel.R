#' @title FitGrootModel
#' @description This function fits gas data to Groot Model: A/(1+(C/time)^B).
#' Please read:Groot, J. C., Cone, J. W., Williams, B. A., Debersaques, F. M., & Lantinga, E. A. (1996). Multiphasic analysis of gas production kinetics for in vitro fermentation of ruminant feeds. Animal Feed Science Technology, 64(1), 77-89.
#' I also added a LagT defined by the intercept of y = 0 and the tangent line passing RMax.
#'
#' @param x A numeric vector. Time
#' @param y A numeric vector. gas volume
#'
#' @return Returns a vactor with all parameters. A: asymptotiv value. C: half time to reach A. B: switching characteristics. Rmax: maximum gas production rate.
#' TRmax: time to reach RMax. LagT: lag phase time. FinalCV: final cumulative volume
#' @export
#' @import nls2
#' @examples Coming.
FitGrootModel <- function(x,y){
  stp <- data.frame(A = max(y),
                    B = seq(from = 1, to = 10, by = 0.1),
                    C = median(x))
  Eq <- function(A, C, B, time){
    A/(1+(C/time)^B)
  }
  # input the equation into fun
  fun = expression(A/(1+(C/time)^B))
  # get the derivative equation of the fun
  #D(fun,"time")
  # put the derivative equation into TRmaxfun, I copy the last command line into the function here
  TRmaxfun <- function(A, C, B, time){
    A * ((C/time)^(B - 1) * (B * (C/time^2)))/(1 + (C/time)^B)^2
  }

  tryCatch({
    tmp_mono_out <- nls2(cCV ~ Eq(A, C, B, time),
                         data = df,
                         start = stp)
    A <- summary(tmp_mono_out)$parameters["A","Estimate"]
    C <- summary(tmp_mono_out)$parameters["C","Estimate"]
    B <- summary(tmp_mono_out)$parameters["B","Estimate"]

    Max <- optimize(TRmaxfun, # the derivative function
                    interval = c(0,50000),
                    tol = 0.00001, # the accuracy
                    A = A,
                    C = C,
                    B = B,
                    maximum = T)
    TRMax <- Max$maximum
    RMax <- Max$objective
  }, error=function(e){})

  FinalCV <- max(y)

  GetIntercept <- function(x,y,k){
    Inter <- y - k*x
    return(Inter)
  }

  GetGrootY <- function(A,B,C,t){
    Y <- A/(1+(C/t)^B)
    return(Y)
  }

  GetLagT <- function(k,b){
    LagT <- -b/k
  }

  LagT <- GetLagT(k = RMax,
                  b = GetIntercept(k = RMax,
                                   x = TRMax,
                                   y = GetGrootY(A = A,
                                                 B = B,
                                                 C = C,
                                                 t = TRMax)))

  Parameter <- c(A,C,B,RMax,TRMax,LagT, FinalCV)
  names(Parameter) <- c("A","C","B","RMax","TRMax","LagT","FinalCV")
  return(Parameter)
}

