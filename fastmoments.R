library(Rcpp)

Rcpp::sourceCpp(code = 
  readChar("onlinemoments.cpp", file.info("onlinemoments.cpp")$size)
)


FastMoments <- function(values_) {
  t = new(moving_moments)
  return (t)
}

t = new(online_moments)
