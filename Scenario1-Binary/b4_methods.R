library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)
#library(doBy)

true_prob <- function(datA, datB, K, prev, error){
  
  comp_mat <- compare4(datA, datB, K=K)
  X.save = comp_mat
  e = error
  p = 1/nrow(datA)
  N = nrow(comp_mat)
  
  indM = (comp_mat[,K+1]==1)
  indU = (comp_mat[,K+1]==0)
  
  comp_mat = comp_mat[,1:K]
  # To compute the true, you replace it with true parameters
  p0M = (1-e)*(1-prev) #(0,0)
  p1M = e*(1-prev) #(0,1) 
  p2M = e*prev #(1,0)
  p3M = 1-p2M-p1M-p0M # (1,1)
  
  p0U = (1-prev)*((1-e)*(1-prev)+e*prev) #(0,0)
  p1U = (1-prev)*((1-e)*prev + e*(1-prev))  #(0,1)
  p2U = prev*((1-e)*(1-prev) + e*prev)  #(1,0)
  p3U = 1- p2U- p1U - p0U # (1,1)
  
  c0 <- as.numeric(comp_mat==0)
  c1 <- as.numeric(comp_mat==1)
  c2 <- as.numeric(comp_mat==2)
  c3 <- as.numeric(comp_mat==3)
  
  p0M_mat = matrix(rep(p0M,N), nrow =N, byrow = TRUE)
  p1M_mat = matrix(rep(p1M,N), nrow =N, byrow = TRUE)
  p2M_mat = matrix(rep(p2M,N), nrow =N, byrow = TRUE)
  p3M_mat = matrix(rep(p3M,N), nrow =N, byrow = TRUE)
  
  p0U_mat = matrix(rep(p0U,N), nrow =N, byrow = TRUE)
  p1U_mat = matrix(rep(p1U,N), nrow =N, byrow = TRUE)
  p2U_mat = matrix(rep(p2U,N), nrow =N, byrow = TRUE)
  p3U_mat = matrix(rep(p3U,N), nrow =N, byrow = TRUE)
  
  
  probM = rowProds(p0M_mat^(c0)*p1M_mat^(c1)*p2M_mat^(c2)*p3M_mat^(c3))
  probU = rowProds(p0U_mat^(c0)*p1U_mat^(c1)*p2U_mat^(c2)*p3U_mat^(c3))
  
  
  g = p*probM/(p*probM+(1-p)*probU)
  
  m0 = p0M
  m1 = p1M + p2M 
  m2 = 1- m0 - m1
  m_true = rbind(m0,m1,m2)
  
  u0 = p0U
  u1 = p1U + p2U 
  u2 = 1- u0 - u1
  u_true = rbind(u0,u1,u2)
  
  return(list(g_true = g, indM = indM, indU = indU, comp_mat4 = X.save, m2 = m2, m1 = m1, m0 = m0, u2 = u2, u1 = u1, u0 = u0))
}

# Fellegi-Sunter method with binary comparison
FS <- function(datA, datB, K,tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB =  nrow(datB)
  
  comp_mat <- compare_binary(datA, datB, K)
  
  fit <- EM_binary(comp_mat, datA, datB, K,tol = tol, maxits = maxits)
  g = fit$g
  converge = fit$converge
  
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 1
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  

  
  return(c(TPR_FS5, PPV_FS5, converge))
}



# Fellegi-Sunter with 3 categorical comparison
FS3 <- function(datA, datB, K, tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB = nrow(datB)
  comp_mat <- compare3(datA, datB, K=K)
  
  
  ## Using EM with the above estimated starting point
  fit = EM3(comp_mat, datA, datB, K, tol=tol, maxits = maxits)
  g = fit$g
  converge = fit$converge
  
  
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 1
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  
  
  return(c(TPR_FS5, PPV_FS5, converge))
}


bayesian <- function(datA, datB, K){
  nB = nrow(datB)
  nA = nrow(datA)
  bayes = recordLink(datB[,1:K], datA[,1:K], eps_plus =0.01, eps_minus = 0.01,use_diff = FALSE)
  g = as.vector(t(bayes))
  
  #############################
  threshold = 0.5
  upper = which(bayes >= threshold)
  indexB = upper%%nB
  indexB[which(indexB==0)] = nB
  indexA = ceiling(upper/nB)
  idB = datB[indexB, K+1]
  idA = datA[indexA, K+1]
  
  TPR_bayes5 = sum(idA==idB)/nB
  if (length(idB) == 0){
    PPV_bayes5 = 1
  }else{
    PPV_bayes5 = sum(idA==idB)/length(idB)
  }
  
  
 
 
  return(c(TPR_bayes5, PPV_bayes5, 1))
}

### Fellegi-Sunter with 4 categorical comparison
FS4 <- function(datA, datB, K, comp_mat4, tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB = nrow(datB)
  comp_mat <- compare4(datA, datB, K=K)
  
  
  ## Using EM
  fit = EM4(comp_mat, datA, datB, K, tol=tol, maxits = maxits)
  
  g = fit$g
  converge = fit$converge

  
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 0
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  
  
  return(c(TPR_FS5, PPV_FS5, converge))
}

