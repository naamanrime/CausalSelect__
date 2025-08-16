expit <-
function(x){
  pr = ( exp(x) / (1+exp(x)) )
  return(pr)
}
