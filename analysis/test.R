for (i in 1:30000){

  print(test[i])
  print(unlist(lapply(lapply(test[i],function(x)strsplit(x,"_")),"[[",3))  )

}
