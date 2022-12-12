  library(fatBVARS)

  load("MonthlyData_Chain1.RData")
  ML_chain1 <- marginalLL(Chain1, ndraws = ndraws, numCores = numCores)
  save(ML_chain1, file=paste("ML_C1_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain2.RData")
  ML_chain2 <- marginalLL(Chain2, ndraws = ndraws, numCores = numCores)
  save(ML_chain2, file=paste("ML_C2_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain3.RData")
  ML_chain3 <- marginalLL(Chain3, ndraws = ndraws, numCores = numCores)
  save(ML_chain3, file=paste("ML_C3_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain4.RData")
  ML_chain4 <- marginalLL(Chain4, ndraws = ndraws, numCores = numCores)
  save(ML_chain4, file=paste("ML_C4_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain5.RData")
  ML_chain5 <- marginalLL(Chain5, ndraws = ndraws, numCores = numCores)
  save(ML_chain5, file=paste("ML_C5_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain11.RData")
  ML_chain11 <- marginalLL(Chain11, ndraws = ndraws, numCores = numCores)
  save(ML_chain11, file=paste("ML_C11_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain12.RData")
  ML_chain12 <- marginalLL(Chain12, ndraws = ndraws, numCores = numCores)
  save(ML_chain12, file=paste("ML_C12_", Sys.time() , ".RData", sep = ""))

  str(ML_chain12)


  load("MonthlyData_Chain6.RData")
  ML_chain6 <- marginalLL(Chain6, ndraws = ndraws, numCores = numCores)
  save(ML_chain6, file=paste("ML_C6_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain7.RData")
  ML_chain7 <- marginalLL(Chain7, ndraws = ndraws, numCores = numCores)
  save(ML_chain7, file=paste("ML_C7_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain8.RData")
  ML_chain8 <- marginalLL(Chain8, ndraws = ndraws, numCores = numCores)
  save(ML_chain8, file=paste("ML_C8_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain9.RData")
  ML_chain9 <- marginalLL(Chain9, ndraws = ndraws, numCores = numCores)
  save(ML_chain9, file=paste("ML_C9_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain10.RData")
  ML_chain10 <- marginalLL(Chain10, ndraws = ndraws, numCores = numCores)
  save(ML_chain10, file=paste("ML_C10_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain13.RData")
  ML_chain13 <- marginalLL(Chain13, ndraws = ndraws, numCores = numCores)
  save(ML_chain13, file=paste("ML_C13_", Sys.time() , ".RData", sep = ""))

  load("MonthlyData_Chain14.RData")
  ML_chain14 <- marginalLL(Chain14, ndraws = ndraws, numCores = numCores)
  save(ML_chain14, file=paste("ML_C14_", Sys.time() , ".RData", sep = ""))

  # save(ML_chain1, ML_chain2, ML_chain3, ML_chain4, ML_chain5, ML_chain6, ML_chain7,
  #      ML_chain8, ML_chain9, ML_chain10, ML_chain11, ML_chain12, ML_chain13, ML_chain14,
  #      file = paste("ML_", Sys.time() , ".RData", sep = ""))
