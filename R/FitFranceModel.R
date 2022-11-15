#' @title FitGrootModel
#' @description This function fits gas data to France monophasic Model: A*(1 - exp(-b*(time - T) - c*(sqrt(time)-sqrt(T))))
#' Dhanoa, M. S., López, S., Powell, C. D., Sanderson, R., Ellis, J. L., Murray, J.-A., Garber, A., Williams, B. A., & France, J. (2021). An illustrative analysis of atypical gas production profiles obtained from in vitro digestibility studies using fecal inoculum. Animals, 11(4), 1069.
#'
#' @param x A numeric vector. Time
#' @param y A numeric vector. gas volume
#'
#' @return Returns a vactor with all parameters. A: asymptotiv value. Rmax: maximum gas production rate.
#' TRmax: time to reach RMax. LagT: lag phase time. FinalCV: final cumulative volume
#' @export
#' @import nls2
#' @examples Coming.
FitFranceModel <- function(x,y){
  # setup the curve functions
  Eq <- function(A, b, T, c, time){
    A*(1 - exp(-b*(time - T) - c*(sqrt(time)-sqrt(T))))
  }

  # input the equation into fun
  fun = expression(A*(1 - exp(-b*(time - T) - c*(sqrt(time)-sqrt(T)))))
  # get the derivative equation of the fun
  #D(fun,"time")
  # put the derivative equation into TRmaxfun, I copy the last command line into the function here
  TRmaxfun <- function(A, b, T, c, time){
    A * (exp(-b * (time - T) - c * (sqrt(time) - sqrt(T))) * (b + c * (0.5 * time^-0.5)))
  }

  tryCatch({
    tmp_mono_out <- nls2(cCV ~ Eq1(A, b, T, c, time),
                         data = Wno_gf_list[[e]],
                         start = list(A = max(y), b = 0.5, c = -1, T = 3))
    A <-  summary(tmp_mono_out)$parameters["A","Estimate"]
    b <-  summary(tmp_mono_out)$parameters["b","Estimate"]
    c <-  summary(tmp_mono_out)$parameters["c","Estimate"]
    LagT <-  summary(tmp_mono_out)$parameters["T","Estimate"]
    Max <- optimize(TRmaxfun, # the derivative function
                    c(LagT, max(x)), # the time interval
                    tol = 0.00001, # the accuracy
                    A = A,
                    b = b,
                    c = c,
                    T = LagT,
                    maximum = T)
    TRMax <- Max$maximum
    RMax <- Max$objective
  }, error=function(e){})

  FinalCV <- max(y)

  Parameter <- c(A,b,c,RMax,TRMax,LagT, FinalCV)
  names(Parameter) <- c("A","b","c","RMax","TRMax","LagT","FinalCV")
  return(Parameter)
}
