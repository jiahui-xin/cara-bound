true_values<-function(model_func,rand="RAR",target="New", TB=100,n=1e6,seed=12345){
  set.seed(seed)
  model_output=model_func(n)
  Ymat=model_output$Ymat
  Y1=Ymat[,2]
  Y0=Ymat[,1]
  true_tau=mean(Y1)-mean(Y0)
  if(rand=="RAR"){
    alloc_prob=target_alloc(mean(Y1), sd(Y1), 
                            mean(Y0), sd(Y0), target, TB)
    true_bound=var(Y1)/alloc_prob+var(Y0)/(1-alloc_prob)
  }
  if(rand=="CARA"){
    strata=model_output$strata
    strata_set=unique(strata)
    strata_len=length(strata_set)
    tau=numeric(strata_len)
    alloc_prob=numeric(strata_len)
    var1=numeric(strata_len)
    var0=numeric(strata_len)
    strata_prob=numeric(strata_len)
    for(s in strata_set){
      Y1s=Y1[strata==s]
      Y0s=Y0[strata==s]
      tau[s]=mean(Y1s)-mean(Y0s)
      alloc_prob[s]=target_alloc(mean(Y1s), sd(Y1s), 
                              mean(Y0s), sd(Y0s), target, TB)
      var1[s]=var(Y1s)
      var0[s]=var(Y0s)
      strata_prob[s]=sum(strata==s)/n
    }
    true_bound=sum(strata_prob*(var1/alloc_prob+var0/(1-alloc_prob)+
                                  (tau-true_tau)^2))
  }
  return(c(true_tau,true_bound))
}


burn_in_plot<-function(model_func, n=500, n0=10,target="Neyman",
                             rand="RAR",TB=30,
                             repli=1000){
  if(rand=="RAR"){
  true_pi=oracle2(model_func(1e5),target,TB)$target_alloc
  df <- data.frame(t(replicate(repli,burn_in_test(model_func(n),n0,target,rand,TB)
  ))
  )
  
  names(df)=c("col1","col2")
  
  # 将数据转换为长格式（stacked format）
  df_long <- data.frame(
    value = c(df$col1, df$col2),
    probability = rep(c("all-sample estimate", "burn-in estimate"), each = repli)
  )
  
  result<-ggplot(df_long, aes(x = probability, y = value, fill = probability)) +
    geom_boxplot(outlier.shape = NA,alpha = 0.5) +  # 设置透明度
    geom_point(aes(x = "all-sample estimate", y = true_pi), shape = 17, color = "red",
               size = 4,show.legend = FALSE) +  # 在第一列上标记五角星
    geom_point(aes(x = "burn-in estimate", y = true_pi), shape = 17, color = "red",
               size = 4,show.legend = FALSE) +  # 在第二列上标记五角星
    labs(#title = "Boxplot Comparison of Two Columns with Star Marker at 0.5",
      x = "",
      y = "") +
    coord_cartesian(ylim=c(0,0.75))+
    theme_minimal() +
    scale_fill_manual(values = c("lightgreen", "skyblue"))  # 手动设置颜色
  }
  if(rand=="CARA"){
    true_pi=oracle(model2(1e5),target,TB)$target_alloc_strata
    n_stratum=length(true_pi)
    data <- data.frame(
      stratum = rep(paste0("stratum",c(1:n_stratum,1:n_stratum)), each = repli),
      probability = rep(c("all-sample estimate", "burn-in estimate"), 
                        each = n_stratum*repli),
      value = c(t(replicate(repli,burn_in_test(model_func(n),n0,target,rand,TB))))
    )
    
    # Plot the boxplot
    result=ggplot(data, aes(y = value, x = stratum, fill = probability)) +
      geom_boxplot(outlier.shape = NA,alpha = 0.4) +  # Transparent boxplot
      #geom_jitter(position = position_jitter(height = 0.2), alpha = 0.5) +  # Jitter points
      geom_point(aes(y = true_pi[1], x = "stratum1"), shape = 17, color = "red", size = 4,
                 show.legend = F) +  # Overlay original value as blue triangle
      geom_point(aes(y = true_pi[2], x = "stratum2"), shape = 17, color = "red", size = 4,
                 show.legend = F) +  # Overlay original value as blue triangle
      geom_point(aes(y = true_pi[3], x = "stratum3"), shape = 17, color = "red", size = 4,
                 show.legend = F) +  # Overlay original value as blue triangle
      labs(title = "",
           x = "", y = "") +
      coord_cartesian(ylim=c(0,1))+
      theme_minimal() +
      scale_fill_manual(values = c("lightgreen", "skyblue", "darkblue", "black")) +  # Color for randomizations
      theme(legend.position = "right")
  }
  return(result)
}

burn_in_experiment<-function(model_func, n=500, n0=10,target="Neyman",
                             rand="CARA",TB=30,
                             repli=1000,round_number=3){
  result=
    rbind(colMeans(t(replicate(1e3,burn_in_test(model_func(n),n0,target,rand,TB)
                               ))),
  apply(t(replicate(1e3,burn_in_test(model_func(n),n0,target,rand,TB) ) ),2,sd))
  if(rand=="CARA")colnames(result)=c(paste0("pi",1:3),paste0("pihat",1:3))
  if(rand=="RAR")colnames(result)=c("pi","pihat")
  rownames(result)=c("est","sd")
  return(round(result,round_number))
}

burn_in_test<-function(model_output,n0=10,target="Neyman",rand="CARA",TB=30){
  Ymat=model_output$Ymat
  strata=model_output$strata
  n=nrow(Ymat)
  pi_hat=numeric(max(strata))
  pi=numeric(max(strata))
  #An=numeric(n)
  if(rand=="RAR"){
      An=sample(rep(0:1,n0),2*n0,replace=F)
      Y1=Ymat[,2]
      Y0=Ymat[,1]
      Y1bar=mean(Y1)
      Y0bar=mean(Y0)
      sigma1=sd(Y1)
      sigma0=sd(Y0)
      pi=target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB)
      
      Y1=(Ymat[1:(2*n0),2])[An==1]
      Y0=(Ymat[1:(2*n0),1])[An==0]
      Y1bar=mean(Y1)
      Y0bar=mean(Y0)
      sigma1=sd(Y1)
      sigma0=sd(Y0)
      pi_hat=target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB)
      return(c(pi,pi_hat))
    }
    
  
  if(rand=="CARA"){
    
    for(s in 1:max(strata)){
      An=sample(rep(0:1,n0),2*n0,replace=F)
      ind=(1:n)[strata==s]
      Ymat_strata=Ymat[ind,]
      Y1=Ymat_strata[,2]
      Y0=Ymat_strata[,1]
      Y1bar=mean(Y1)
      Y0bar=mean(Y0)
      sigma1=sd(Y1)
      sigma0=sd(Y0)
      pi[s]=target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB)
      
      Y1=(Ymat_strata[1:(2*n0),2])[An==1]
      Y0=(Ymat_strata[1:(2*n0),1])[An==0]
      Y1bar=mean(Y1)
      Y0bar=mean(Y0)
      sigma1=sd(Y1)
      sigma0=sd(Y0)
      pi_hat[s]=target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB)
    }
    
  }
  if(rand=="CADBCD"){
    pi_hat=numeric(max(strata))
    pi=numeric(max(strata))
    An=sample(rep(0:1,n0),2*n0,replace=F)
    for(s in 1:max(strata)){
      ind=(1:n)[strata==s]
      Ymat_strata=Ymat[ind,]
      Y1=Ymat_strata[,2]
      Y0=Ymat_strata[,1]
      Y1bar=mean(Y1)
      Y0bar=mean(Y0)
      sigma1=sd(Y1)
      sigma0=sd(Y0)
      pi[s]=target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB)
      
      
      Y1=(Ymat[1:(2*n0),2])[An==1&strata[1:(2*n0)]==s]
      Y0=(Ymat[1:(2*n0),1])[An==0&strata[1:(2*n0)]==s]
      Y1bar=mean(Y1)
      Y0bar=mean(Y0)
      sigma1=sd(Y1)
      sigma0=sd(Y0)
      pi_hat[s]=target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB)
    }
  }
  
  
  return (c(pi,pi_hat))
  }
    


BCD<-function(n,lambda=0.75){
  An=numeric(n)
  An[1]=sample(0:1,1)
  D=An[1]-1/2
  for(i in 2:n){
    if(D>0)prob=0.75
    else if(D<0)prob=0.25
    else prob=0.5
    An[i]=runif(1)>prob
    D=D+An[i]-1/2
  }
  return(An)
}

MIN<-function(model_output, lambda=0.75){
  Ymat=model_output$Ymat
  strata=model_output$strata
  n=nrow(Ymat)
  An=numeric(n)
  strata_set=unique(strata)
  for(s in strata_set){
    strata_ind=(1:n)[strata==s]
    size_strata=sum(strata==s)
    An[strata_ind]=BCD(size_strata, lambda)
  }
  return(An)
}

CADBCD_R<-function(model_output,  gamma = 2, n0 = 30, target="Neyman", TB = 30){
  Ymat=model_output$Ymat
  strata=model_output$strata
  n=nrow(Ymat)
  An=numeric(n)
  pi=numeric(max(strata))
  An[1:(2*n0)]=sample(rep(0:1,n0),2*n0,replace=F)
  for(i in (2*n0+1):n){
    rho=0
    for(s in 1:max(strata)){
      pre_seq=1:(i-1)
      ind=pre_seq[strata[pre_seq]==s]
      ind1=intersect(ind,pre_seq[An[pre_seq]==1])
      ind0=intersect(ind,pre_seq[An[pre_seq]==0])
      Y1=Ymat[ind1,2]
      Y0=Ymat[ind0,1]
      Y1bar=mean(Y1)
      Y0bar=mean(Y0)
      sigma1=sd(Y1)
      sigma0=sd(Y0)
      pi[s]=target_alloc(Y1bar, sigma1, Y0bar, sigma0, target, TB)
      rho=rho+mean(strata[pre_seq]==s)*pi[s]
    }
    N1=sum(An[pre_seq])
    j=strata[i]
    prob=pi[j]*(rho/N1*(i-1))^gamma/
      (pi[j]*(rho/N1*(i-1))^gamma+(1-pi[j])*((1-rho)/(1-N1/(i-1)))^gamma)
    #print(prob)
    if(prob<0.1)prob=0.1
    if(prob>0.9)prob=0.9
    An[i]=runif(1)<prob
  }
  return(An)
}

CARA_R<-function(model_output,  gamma = 2, n0 = 10, 
               target = "Neyman", TB = 1){
  Ymat=model_output$Ymat
  strata=model_output$strata
  n=nrow(Ymat)
  An=numeric(n)
  strata_set=unique(strata)
  for(s in strata_set){
    strata_ind=(1:n)[strata==s]
    Ymat_strata=Ymat[strata_ind,]
    An[strata_ind]=DBCD(Ymat_strata,gamma,n0,target,TB)
  }
  return(An)
}

DID<-function(Ymat, An){
  y0=Ymat[,1]
  y1=Ymat[,2]
  y0_mean=sum(y0*(1-An))/sum(1-An)
  y1_mean=sum(y1*An)/sum(An)
  return(y1_mean-y0_mean)
}

SDID<-function(model_output, An){
  Ymat=model_output$Ymat
  strata=model_output$strata
  n=nrow(Ymat)
  strata_set=unique(strata)
  tau=0
  for(s in strata_set){
    n_strata=sum(strata==s)
    strata_ind=(1:n)[strata==s]
    Ymat_strata=Ymat[strata_ind,]
    An_strata=An[strata_ind]
    tau_strata=DID(Ymat_strata, An_strata)
    tau=tau+n_strata/n*tau_strata
  }
  return(tau)
}

model1<-function(n){
  strata=rbinom(n,1,1/3)+1
  y1=rbinom(n,1,0.9)*(strata==1)+rbinom(n,1,0.5)*(strata==2)#+rbinom(n,1,0.9)*(strata==3)
  y0=rbinom(n,1,0.5)*(strata==1)+rbinom(n,1,0.5)*(strata==2)#+rbinom(n,1,0.5)*(strata==3)
  Ymat=cbind(y0,y1)
  return(list(Ymat=Ymat,strata=strata))
}

model2<-function(n){
  strata=sample(1:3,n,replace=T)
  y0=(strata+1)*rt(n,5,strata)
  y1=rt(n,5,strata)+20
  vec=numeric(sum(strata==2))
  vec=y0[strata==2]
  y0[strata==2]=y1[strata==2]-10
  y1[strata==2]=vec+20
  Ymat=cbind(y0,y1)
  return(list(Ymat=Ymat,strata=strata))
}

model3<-function(n){
  strata=sample(1:3,n,replace=T)
  y0=rnorm(n,20,strata)
  y1=rnorm(n,60-5*strata,2*(4-strata))
  Ymat=cbind(y0,y1)
  return(list(Ymat=Ymat,strata=strata))
}

oracle<-function(model_output, target, TB=1){
  Ymat=model_output$Ymat
  strata=model_output$strata
  n=nrow(Ymat)
  strata_set=unique(strata)
  num_of_strata=length(strata_set)
  prob_strata=numeric(num_of_strata)
  tau_strata=numeric(num_of_strata)
  Y0bar_strata=numeric(num_of_strata)
  Y1bar_strata=numeric(num_of_strata)
  sigma0_strata=numeric(num_of_strata)
  sigma1_strata=numeric(num_of_strata)
  target_alloc_strata=numeric(num_of_strata)
  for(s in strata_set){
    n_strata=sum(strata==s)
    strata_ind=(1:n)[strata==s]
    Ymat_strata=Ymat[strata_ind,]
    prob_strata[s]=n_strata/n
    Y0bar_strata[s]=mean(Ymat_strata[,1])
    Y1bar_strata[s]=mean(Ymat_strata[,2])
    tau_strata[s]=Y1bar_strata[s]-Y0bar_strata[s]
    sigma0_strata[s]=sd(Ymat_strata[,1])
    sigma1_strata[s]=sd(Ymat_strata[,2])
    target_alloc_strata[s]=target_alloc(
      Y1bar_strata[s], sigma1_strata[s], Y0bar_strata[s], sigma0_strata[s],
      target, TB)
  }
  return(list(prob_strata=prob_strata, tau_strata=tau_strata,
              Y0bar_strata=Y0bar_strata, Y1bar_strata=Y1bar_strata,
              sigma0_strata=sigma0_strata, sigma1_strata=sigma1_strata,
              target_alloc_strata=target_alloc_strata))
}

oracle2<-function(model_output, target, TB=1){
  Ymat=model_output$Ymat
  Y0bar=mean(Ymat[,1])
  Y1bar=mean(Ymat[,2])
  tau=Y1bar-Y0bar
  sigma0=sd(Ymat[,1])
  sigma1=sd(Ymat[,2])
  target_alloc_val=target_alloc(
      Y1bar, sigma1, Y0bar, sigma0, target, TB)
  return(list(tau=tau, Y0bar=Y0bar, Y1bar=Y1bar, sigma0=sigma0, sigma1=sigma1,
              target_alloc=target_alloc_val))
}

repli_func<-function(model_func, n=500, rand="CARA", target="Neyman", 
                      TB=1, n0=10, gamma=2){
  model_output=model_func(n)
  Ymat=model_output$Ymat
  #if(target=="New")oracle_val=oracle(model_output, "New",TB)
  #else oracle_val=oracle(model_output, "Neyman",TB)
  if(rand=="BCD")An=BCD(n)
  if(rand=="MIN")An=MIN(model_output)
  if(rand=="CR")An=CRand(Ymat)
  if(rand=="RAR")An=DBCD(Ymat, gamma = gamma, target = target, TB = TB, n0 = n0)
  if(rand=="CARA")An=CARA(model_output, target = target, TB = TB,n0 = n0,
                          gamma=gamma)
  if(rand=="CADBCD")An=CADBCD(model_output, 
                              gamma = gamma, target = target, 
                              TB = TB,n0 = n0)
  sdid=SDID(model_output, An)
  did=DID(Ymat, An)
  Y=An*Ymat[,2]+(1-An)*Ymat[,1]
  strata=model_output$strata
  strata_set=unique(strata)
  num_of_strata=length(strata_set)
  Ymean=numeric(num_of_strata)
  for(s in strata_set){
    n_strata=sum(strata==s)
    strata_ind=(1:n)[strata==s]
    Ymat_strata=Ymat[strata_ind,]
    An_strata=An[strata==s]
    Y_strata=An_strata*Ymat_strata[,2]+(1-An_strata)*Ymat_strata[,1]
    Ymean[s]=mean(Y_strata)
  }
  result=list(sdid=sdid,did=did,Ymean=Ymean)
  return(unlist(result))
}

repli_func2<-function(model_func, n=500, rand="RAR", 
                      target="Neyman", TB=1, n0=10, gamma = 2){
  model_output=model_func(n)
  Ymat=model_output$Ymat
  #if(target=="New")oracle_val=oracle2(model_output, "New", TB)
  #else oracle_val=oracle2(model_output, "Neyman", TB)
  if(rand=="BCD")An=BCD(n)
  if(rand=="MIN")An=MIN(model_output)
  if(rand=="CR")An=CRand(Ymat)
  if(rand=="RAR")An=DBCD(Ymat, target = target, TB = TB, n0 = n0, gamma=gamma)
  if(rand=="CARA")An=CARA(model_output, gamma = gamma,
                          target = target, TB = TB, n0 = n0)
  did=DID(Ymat, An)
  Y=An*Ymat[,2]+(1-An)*Ymat[,1]
  Ymean=mean(Y)
  result=list(did=did,Ymean=Ymean)
  return(unlist(result))
}



 
experiment<-function(model_func=model2, n=500, randomise="CARA",
                     target="New", TB=1, 
                     n0=10,  gamma =2, repli=1e4,seed=12345){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
    results <- mclapply(
      1:repli, 
      function(x) {
        repli_func(model_func, n, randomise, target, TB,n0,gamma)
      },
      mc.cores = 8,
      mc.set.seed = TRUE)
    final_result <- as.data.frame(do.call(rbind, results))
    return(final_result)
}

experiment2<-function(model_func=model2, n=500, randomise="CARA", target="New", TB=1, 
                     n0=10, gamma=2, repli=1e4, seed=123){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  results <- mclapply(
    1:repli, 
    function(x) {
      repli_func2(model_func, n, randomise, target, TB,n0,gamma)
    },
    mc.cores = 8,
    mc.set.seed = TRUE)
  final_result <- as.data.frame(do.call(rbind, results))
  return(final_result)
}

experiment_rng<-function(model_func=model2, n=500, randomise="CARA",
                     target="New", TB=1, 
                     n0=10,  gamma =2, repli=3e3){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  results <- mclapply(
    1:repli, 
    function(x) {
      repli_func(model_func, n, randomise, target, TB,n0,gamma)
    },
    mc.cores = 8,
    mc.set.seed = TRUE)
  final_result <- as.data.frame(do.call(rbind, results))
  return(final_result)
}

experiment2_rng<-function(model_func=model2, n=500, randomise="CARA", target="New", TB=1, 
                      n0=10, gamma=2, repli=1e3, seed=123){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  results <- mclapply(
    1:repli, 
    function(x) {
      repli_func2(model_func, n, randomise, target, TB,n0,gamma)
    },
    mc.cores = 8,
    mc.set.seed = TRUE)
  final_result <- as.data.frame(do.call(rbind, results))
  return(final_result)
}


extract2<-function(final_result, round_num=3){
  bias=with(final_result,mean(did-true_tau))
  #bound= (colMeans(final_result)[6]^2/colMeans(final_result)[7]+
  #               colMeans(final_result)[5]^2/(1-colMeans(final_result)[7]))/n 
  variance=var(final_result$did)
  Ymean=mean(final_result$Ymean)
  info = c(true_bound/n,Ymean,bias,variance)
  names(info) = NULL
  return (round(info,round_num) )
}

extract<-function(final_result,round_num=3){
  #bound= sum(colMeans(final_result)[3:5]*(
  #  (colMeans(final_result)[6:8]-true_tau)^2 +
  #    colMeans(final_result)[15:17]^2/
  #    (1-colMeans(final_result)[21:23])+
  #    colMeans(final_result)[18:20]^2/
  #    colMeans(final_result)[21:23])/n )
  #tau = sum(colMeans(final_result)[3:5]*colMeans(final_result)[6:8])
  bias_did=mean(final_result$did)-true_tau
  bias_sdid=mean(final_result$sdid)-true_tau
  variance_did=var(final_result$did)
  variance_sdid=var(final_result$sdid)
  Ymean=colMeans(final_result)[3:5]
  info = c(true_bound/n,Ymean,bias_did,variance_did,bias_sdid,variance_sdid)
  names(info) = NULL
  return (round(info,round_num) )
}


extract2.binary<-function(final_result, round_num=3){
  bias=with(final_result,mean(did-true_tau))
  #bound= (colMeans(final_result)[6]^2/colMeans(final_result)[7]+
  #         colMeans(final_result)[5]^2/(1-colMeans(final_result)[7]))/n 
  variance=var(final_result$did)
  Ymean=mean(final_result$Ymean)
  bound=final_result$true_bound
  info = c(true_bound/n,Ymean,bias,variance)
  names(info) = NULL
  return (round(info,round_num) )
}

extract.binary<-function(final_result,round_num=3){
  #bound= sum(colMeans(final_result)[3:4]*(
  #  (colMeans(final_result)[5:6]-true_tau)^2 +
  #    colMeans(final_result)[11:12]^2/
  #    (1-colMeans(final_result)[15:16])+
  #   colMeans(final_result)[13:14]^2/
  #   colMeans(final_result)[15:16])/n )
  #tau = sum(colMeans(final_result)[3:4]*colMeans(final_result)[5:6])
  bias_did=mean(final_result$did)-true_tau
  bias_sdid=mean(final_result$sdid)-true_tau
  variance_did=var(final_result$did)
  variance_sdid=var(final_result$sdid)
  Ymean=colMeans(final_result)[3:4]
  info = c(true_bound/n,Ymean,bias_did,variance_did,bias_sdid,variance_sdid)
  names(info) = NULL
  return (round(info,round_num) )
}
