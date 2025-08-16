DR =  function(ps,or1,or0,A,Y){
  (A*Y-(A-ps)*or1)/ps -  (   (1-A)*Y +(A-ps)*or0)/(1-ps)
}