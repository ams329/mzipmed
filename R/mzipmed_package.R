#MZIP Estimation Function




#' Marginalized Zero-Inflated Poisson Regression Model
#'
#' This function uses Long et. al's(2014) MZIP model to allow you to fit counts variables with excess zeroes
#'    while allowing for easy interpretations. This function assumes that
#'    the outcome and covariates are all the same sample size without missing
#'    data. Covariates must be numerical, so binary predictors such as
#'    gender or race need to be dummy coded with zeroes and ones. For more
#'    information about this model and interpretations see Long, D Leann et al. "A marginalized zero-inflated Poisson regression model with overall exposure effects." Statistics in medicine vol. 33,29 (2014): 5151-65. doi:10.1002/sim.6293.
#'    Note: BFGS likelihood optimization was used for this R package
#' @param y is the outcome variable
#' @param pred is a vector of covariates (use cbind for multiple)
#' @param print if =T or =TRUE will print beta coefficient estimates, standard errors, and p-values for the  into the console. If =F or FALSE nothing will be printed into the console. Default is FALSE
#' @return The function will return a list of 22 elements.
#'     In the list G(Gamma) refers to the excess zero/logistic part of the model. \cr
#'     and A(Alpha) refers to the Poisson/mean part of the model for example. \cr
#'     Gest are the gamma coefficients for the logistic part of the MZIP model. \cr
#'     Aest are the alpha coefficients for the Poisson part of the MZIP model. \cr
#'     _ModelSE are the standard errors for each coefficient in the model.\cr
#'     _RobustSE are the robust standard errors for each coefficient in the model. \cr
#'     _ModelUpper are the upper confidence limits for each coefficient. \cr
#'     _ModelLower are the lower confidence limits. \cr
#'     _RobustUpper are the upper confidence limits based on robust standard error. \cr
#'     _RobustLower are the lower confidence limits based on robust standard errors. \cr
#'     _ModelZ are the Z scores for each coefficient. \cr
#'     _RobustZ are the robust Z scores for each coefficient. \cr
#'     _ModelZpval are the p-values based on the Z scores for the model. \cr
#'     _RobustZpval are the p-values based on the robust z scores. \cr
#'     AlphaCov is the covariance matrix for the poisson coefficient estimates \cr
#'     Cov is the covariance matrix for the MZIP model
#' @examples
#'     mzip(y=CarriesCount,pred=sex,Print=F)
#'     mzip(y=CarriesCount,pred=cbind(sex,BMI,SBP),print=T)
#'     mzip(y=data$outcome,pred=cbind(data$exposure,data$confounder))
#' @export



mzip = function(y,pred,print=F){
  
  intercept=rep(1,length(y))
  Z=cbind(intercept,pred)
  X=Z
  
  
  like = function(theta) {
    
    
    
    
    zg = Z%*%c(theta[1:dim(Z)[2]])
    z0g = Z[y==0,]%*%c(theta[1:dim(Z)[2]])
    x0a = X[y==0,]%*%c(theta[(1+dim(Z)[2]):(dim(Z)[2]+dim(X)[2])])
    z1g = Z[y>0,]%*%c(theta[1:dim(Z)[2]])
    x1a = X[y>0,]%*%c(theta[(1+dim(Z)[2]):(dim(Z)[2]+dim(X)[2])])
    
    #log likelihood
    bill=-1*sum(log(1+exp(zg)))+
      sum(log(exp(z0g)+exp(-(1+exp(z0g))*exp(x0a))))+
      sum(-(1+exp(z1g))*exp(x1a))+
      sum(y[y>0]*log(1+exp(z1g)))+
      sum(x1a*y[y>0])-
      sum(lgamma(y[y>0]+1))
    
    return(-bill)
  }
  
  estimates=stats::optim(rep(0,dim(Z)[2]+dim(X)[2]), like,hessian=T,method="BFGS")
  
  gamma_hat = estimates$par[1:dim(Z)[2]]
  alpha_hat = estimates$par[(1+dim(Z)[2]):(dim(Z)[2]+dim(X)[2])]
  
  z_gamma_hat 	= Z%*%gamma_hat # n x 1 vector
  x_alpha_hat 	= X%*%alpha_hat # n x 1 vector*/
  
  outcome=y
  
  psi_hat = exp(z_gamma_hat)/(1+exp(z_gamma_hat)) #n x 1 vector
  nu_hat = exp(x_alpha_hat) #n x 1 vector
  psi_hat2 = 1/(1-psi_hat)
  
  diag_gg	= (psi_hat^2*(1-psi_hat)*(psi_hat*psi_hat2*nu_hat+1)*(exp(nu_hat*psi_hat2)-nu_hat*(psi_hat2)-1))/(psi_hat*exp(nu_hat*psi_hat2)+(1-psi_hat))
  
  diag_aa = (nu_hat*(psi_hat*(exp(nu_hat*psi_hat2)-nu_hat*psi_hat2-1)+1))/(psi_hat*exp(nu_hat*psi_hat2)+(1-psi_hat))
  
  diag_ga = nu_hat*psi_hat*(1-exp(-nu_hat*psi_hat2)+(-(1+nu_hat*psi_hat*psi_hat2+nu_hat*(psi_hat*psi_hat2)^2+exp(-nu_hat*psi_hat2))/((psi_hat*psi_hat2)*exp(nu_hat*psi_hat2)+1)))
  
  #test for error
  dd_gg=diag(as.vector(diag_gg))
  dd_aa=diag(as.vector(diag_aa))
  
  tx_dd_gg=t(Z)%*%diag(as.vector(diag_gg))
  dd_gg_x=diag(as.vector(diag_gg))%*%(Z)
  
  I_gg			= t(Z)%*%diag(as.vector(diag_gg))%*%(Z)
  I_aa			= t(X)%*%diag(as.vector(diag_aa))%*%(X)
  I_ga			= t(X)%*%diag(as.vector(diag_ga))%*%(Z)
  I_ag			= t(I_ga)
  
  Inform=estimates$hessian
  Inv_inform		= MASS::ginv(Inform)
  w=nrow(Inv_inform)/2+1
  k=nrow(Inv_inform)
  aCov=Inv_inform[w:k,w:k]
  
  M1 = matrix(0,nrow=dim(Z)[2]+dim(X)[2],ncol=dim(Z)[2]+dim(X)[2])
  score_g = matrix(0,nrow=dim(Z)[2],ncol=1)
  score_a = matrix(0,nrow=dim(X)[2],ncol=1)
  
  for(qq in 1:dim(Z)[1]){
    
    y = outcome[qq]
    ph= psi_hat[qq]
    ph2=psi_hat2[qq]
    nu=nu_hat[qq]
    
    score_g = ((y==0)*(ph*ph2*(exp(nu*ph2)-nu))/(ph*ph2*exp(nu*ph2)+1)+ph*(y - 1) - (y>0)*ph*ph2*nu)%*%t(Z[qq,])
    
    score_a = ((y-nu*ph2)*(y>0) - (y==0)*(nu*ph2)/(ph*ph2*exp(nu*ph2)+1))%*%t(X[qq,])
    
    score = cbind(score_g,score_a)
    
    M1 = M1 + t(score)%*%(score)
    
  }
  
  robust = Inv_inform%*%M1%*%Inv_inform
  
  m_se			= sqrt(diag(Inv_inform))
  
  r_se			= sqrt(diag(robust))
  
  mupper = c(gamma_hat,alpha_hat) + 1.96*m_se
  mlower = c(gamma_hat,alpha_hat) - 1.96*m_se
  rupper = c(gamma_hat,alpha_hat) + 1.96*r_se
  rlower = c(gamma_hat,alpha_hat) - 1.96*r_se
  
  GWald=(gamma_hat/m_se[1:dim(Z)[2]])^2
  GPval=ifelse(1-stats::pchisq(q=(gamma_hat/m_se[1:dim(Z)[2]])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(gamma_hat/m_se[1:dim(Z)[2]])^2,df=1),digits=5))
  GRobWald=(gamma_hat/r_se[1:dim(Z)[2]])^2
  GRobPval=ifelse(1-stats::pchisq(q=(gamma_hat/r_se[1:dim(Z)[2]])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(gamma_hat/r_se[1:dim(Z)[2]])^2,df=1),digits=5))
  AWald=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2
  APval=ifelse(1-stats::pchisq(q=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1),digits=5))
  ARobWald=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2
  ARobPval=ifelse(1-stats::pchisq(q=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1)<.0001,"<.0001",round(1-pchisq(q=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1),digits=5))
  
  
  if(print){cat("Gamma Estimates:",gamma_hat,'\n',"Alpha estimates:",alpha_hat,'\n',"M SE: ", m_se,'\n',"R SE:",r_se,'\n',"Gamma P-Value",GPval,'\n',"R Gamma P-Val",GRobPval,'\n',"Alpha P-Val",APval,'\n',"R Alpha P-Val",ARobPval,'\n')}
  
  output = list( Gest = gamma_hat,Aest = alpha_hat, GModelSE = m_se[1:dim(Z)[2]], GRobustSE = r_se[1:dim(Z)[2]],
                 AModelSE = m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])], ARobustSE = r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GModelUpper = mupper[1:dim(Z)[2]], AModelUpper = mupper[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GModelLower = mlower[1:dim(Z)[2]], AModelLower = mlower[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GRobustUpper = rupper[1:dim(Z)[2]], ARobustUpper = rupper[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GRobustLower = rlower[1:dim(Z)[2]], ARobustLower = rlower[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 GModelZ = gamma_hat/m_se[1:dim(Z)[2]], GRobustZ = gamma_hat/r_se[1:dim(Z)[2]],
                 GModelZpval = 2*(1-stats::pnorm(abs(gamma_hat/m_se[1:dim(Z)[2]]))),
                 GRobustZpval = 2*(1-stats::pnorm(abs(gamma_hat/r_se[1:dim(Z)[2]]))),
                 AModelZ = alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 ARobustZ = alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])],
                 AModelZpval = 2*(1-stats::pnorm(abs(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])]))),
                 ARobustZpval = 2*(1-stats::pnorm(abs(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])]))),
                 AlphaCov=aCov,Cov=Inv_inform)
  
  return(output)
}





#' Mediation Analysis for Zero-Inflated Count Mediators using MZIP (Continuous Outcome)
#'
#' This function incorporates the MZIP model into the counterfactual approach to mediation analysis
#' as proposed by Vanderweele when the mediator is a Zero-Inflated count variable. Standard Errors for
#' direct and indirect effects are computed using delta method or bootstrapping. Note: This function
#' assumes that the outcome is continuous and all exposure, mediator, outcome, and confounder variables
#' have the same sample size. Binary variables must be dummy coded prior.
#' @param outcome is the continuous outcome variable
#' @param mediator is the zero-inflated mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetition if bootstrapped errors are used
#' @param C is a vector for theoretical values of each confounder. By default each each value of C will be the mean value of each confounder.
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @return The function will return a list of 12 elements.
#'     LM is the results of regressing the mediator+exposure+confounder on the outcome using a linear model \cr
#'     MZIP is the results of regressing the exposure and confounders on the mediator using the MZIP model \cr
#'     NDE is the direct effect \cr
#'     NIE is the indirect effect. \cr
#'     NDEse is the standard error for the  direct effect \cr
#'     NDEci is the 95% confidence interval for the direct effect\cr
#'     NIEse is the standard error for  the indirect effect \cr
#'     NIEci is the 95% confidence interval for the indirect effect \cr
#'     TE is the total effect \cr
#'     TEse is the standard error for the total effect \cr
#'     TECI is the confidence interval for the total effect \cr
#'     PM is the proportion mediated
#' @examples
#'     lmoutzimed(outcome=ContinuousOutcome,mediator=ZICount,exposure=race,n=1000,error='Boot')
#'     lmoutzimed(outcome=data$outcome,mediator=data$mediator,exposure=data$exp,confounder=cbind(data$var1,data$var2),X=10,Xstar=0,C=c(1,3))
#'     lmoutzimed(outcome=outcome,mediator=mediator,exposure=race,confounder=sex,C=0,error='Delta')
#' @export



lmoutzimed=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,error='Delta'){
  lmout=data.frame(outcome)
  
  if (is.null(confounder)){
    lmpred=data.frame(exposure,mediator)
    medreg=mzip(y=mediator,pred=exposure,print=F)
  } else{
    lmpred=data.frame(exposure,mediator,confounder)
    medreg=mzip(y=mediator,pred=cbind(exposure,confounder),print=F)
  }
  #as.formula and lm part of stats package
  lmdata=data.frame(lmout,lmpred)
  f<-stats::as.formula(paste(colnames(lmout),paste(colnames(lmpred),collapse="+"),sep="~"))
  
  m=ncol(lmdata)-1
  r=ncol(lmdata)
  outreg=stats::lm(f,data=lmdata)
  
  if (!is.null(confounder)){
    if (is.null(C)){
      confounder=cbind(confounder)
      C=colMeans(confounder)
    }
  }
  
  #Direct Effect
  DE=outreg$coefficients[[2]]*(X-Xstar)
  
  
  if (is.null(confounder)){
    IE=outreg$coefficients[[3]]*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
  } else{
    IE=outreg$coefficients[[3]]*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
  }
  #Proportion Mediated
  PM=IE/(IE+DE)
  TE=IE+DE
  
  if (error=='Delta'){
    
    if (is.null(confounder)){
      GamDE=c(0,0,0,1,0)
      d1=outreg$coefficients[[3]]*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      d2=outreg$coefficients[[3]]*(X*exp(medreg$Aest[1]+medreg$Aest[2]*X)-Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      d4=0
      d5=0
      d6=(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      GamIE=c(d1,d2,d4,d5,d6)
    } else {
      
      GamDE=c(0,0,0*C,0,1,0,0*C)
      d1=outreg$coefficients[[3]]*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      d2=outreg$coefficients[[3]]*(X*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))-Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      d3=C*outreg$coefficients[[3]]*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      d4=0
      d5=0
      d6=(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      d7=C*0
      GamIE=c(d1,d2,d3,d4,d5,d6,d7)
    }
    
    GamTE=GamIE+GamDE
    #Extract Covariance matrices
    MZIPCov=medreg$AlphaCov
    lmCov=stats::vcov(outreg) #uses stats package
    
    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nMZI,nlm)
    Botmleft=matrix(0,nlm,nMZI)
    
    Top=cbind(MZIPCov,Topright)
    Btm=cbind(Botmleft,lmCov)
    
    CovM=rbind(Top,Btm)
    
    #Standard Errors
    DEse=sqrt(GamDE %*% CovM %*% GamDE)*abs(X-Xstar)
    IEse=sqrt(GamIE %*% CovM %*% GamIE)*abs(X-Xstar)
    TEse=sqrt(GamTE %*% CovM %*% GamTE)*abs(X-Xstar)
    
    #Confidence Intervals
    LDEci=DE-1.96*DEse
    UDEci=DE+1.96*DEse
    DECI=c(LDEci,UDEci)
    
    LIEci=IE-1.96*IEse
    UIEci=IE+1.96*IEse
    IECI=c(LIEci,UIEci)
    
    LTEci=TE-1.96*TEse
    UTEci=TE+1.96*TEse
    TECI=c(LTEci,UTEci)
  } 
  
  if (error=='Boot'){
    datab=list()
    datab2=list()
    outregb=list()
    medregb=list()
    DEb=list()
    IEb=list()
    PMb=list()
    confb=list()
    TEb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(lmdata),replace=T)
      datab2[[i]]=lmdata[datab[[i]],]
      outregb[[i]]=stats::lm(f,data=datab2[[i]])
      if (is.null(confounder)){
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=datab2[[i]][["exposure"]],print=F)
        IEb[[i]]=outregb[[i]]$coefficients[[3]]*(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X)-exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar))
      } else {
        confb[[i]]=data.matrix(datab2[[i]][c(4:r)])
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=cbind(datab2[[i]][["exposure"]],confb[[i]]),print=F)
        IEb[[i]]=outregb[[i]]$coefficients[[3]]*(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+sum(medregb[[i]][["Aest"]][c(3:m)]*C))-exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:m)]*C)))
      }
      DEb[[i]]=outregb[[i]]$coefficients[[2]]*(X-Xstar)
      # PMb[[i]]=IEb[[i]]/(IEb[[i]]+DEb[[i]])
      TEb[[i]]=IEb[[i]]+DEb[[i]]
    }
    
    #matrixStats and stats
    #Direct Effect
    DEse=matrixStats::colSds(do.call(rbind,DEb))
    DECI=stats::quantile(do.call(rbind,DEb),c(0.025,.975))
    
    #Indirect Effect #Did not include the 1-0 here
    
    IEse=matrixStats::colSds(do.call(rbind,IEb))
    IECI=stats::quantile(do.call(rbind,IEb),c(0.025,.975))
    #Proportion Mediated
    TEse=matrixStats::colSds(do.call(rbind,TEb))
    TECI=stats::quantile(do.call(rbind,TEb),c(0.025,0.975))
  }
  
  output=list(lm=outreg,mzip=medreg,NDE=DE,NIE=IE,NDEse=DEse,NIEse=IEse,NDEci=DECI,NIEci=IECI,
              TE=TE,TEse=TEse,TEci=TECI,PM=PM)
}




#' Mediation Analysis for Zero-Inflated Count Mediators using MZIP with Exposure-Mediator Interactions (Continuous Outcome)
#'
#' This function will do the same thing as the lmoutzimed function, but includes an exposure-mediator interaction.
#' 4-way decomposition of total effect (Vanderweele) are included in the output.
#' @param outcome is the continuous outcome variable
#' @param mediator is the zero-inflated mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetitions for bootstrapping. Default is 1000. Setting n when using delta method errors will have no effect on output.
#' @param C is a vector for theoretical values of each confounder. If left out the default will be set to the mean of each confounder giving marginal effects
#' @param M is a fixed value for the mediator, M. If M is not specified, M will be set to its mean value 
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @return The function will return a list of 30 elements.
#'     LM is the results of regressing the mediator+exposure+confounder on the outcome using a linear model. To assess interaction effect individually look in the lm statement at the 4th parameter estimate \cr
#'     MZIP is the results of regressing the exposure and confounders on the mediator using the MZIP model \cr
#'     CDE is the controlled direct effect \cr
#'     NDE is the natural direct effect \cr
#'     NIE is the indirect effect. \cr
#'     PM is the proportion mediated\cr
#'     CDEse is the standard error for the  controlled direct effect \cr
#'     CDEci is the 95% confidence interval for the controlled direct effect\cr
#'     NDEste is the standard error for the  natural direct effect \cr
#'     NDEci is the 95% confidence interval for the natural direct effect\cr
#'     NIEse is the standard error for  the indirect effect \cr
#'     NIEci is the 95% confidence interval for the indirect effect \cr
#'     Intref is the Interactive Reference effect \cr
#'     Intrefse is the standard error for Intref \cr
#'     IntrefCI is the CI for Intref \cr
#'     PIE is the pure indirect effect \cr
#'     PIEse is the standard error of PIE \cr
#'     PIECI is the CI for PIE \cr
#'     Intmed is the interactive mediation effect \cr
#'     Intmedse is the error associated with Intmed \cr
#'     IntmedCI is the CI for Intmed \cr
#'     TE is the total effect \cr
#'     TEse is the error of the total effect \cr
#'     TECI is the CI for the total effect \cr
#'     OvInt is the overall additive interaction effect \cr
#'     OvIntse is the standard error for the additive interaction \cr
#'     OvIntCI is the confidence interval for the interaction effect \cr
#'     PropInt is the proportion attributable to the interaction effect \cr
#'     PE is the proportion eliminated
#' @examples
#'     lmoutzimedint(outcome=ContinuousOutcome,mediator=ZICount,exposure=race,n=200,error="Boot")
#'     lmoutzimedint(outcome=data$outcome,mediator=data$mediator,exposure=data$exp,confounder=cbind(data$var1,data$var2),X=10,Xstar=0,C=c(1,3),M=100)
#'     lmoutzimedint(outcome=outcome,mediator=mediator,exposure=race,confounder=sex,C=0,M=0,error='Delta')
#' @export


lmoutzimedint=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,M=NULL,error='Delta'){
  
  interaction=mediator*exposure
  lmout=data.frame(outcome)
  if (is.null(confounder)){
    lmpred=data.frame(exposure,mediator,interaction)
    medreg=mzip(y=mediator,pred=exposure,print=F)
  } else{
    lmpred=data.frame(exposure,mediator,interaction,confounder)
    medreg=mzip(y=mediator,pred=cbind(exposure,confounder),print=F)
  }
  
  lmdata=data.frame(lmout,lmpred)
  f<-stats::as.formula(paste(colnames(lmout),paste(colnames(lmpred),collapse="+"),sep="~"))
  
  k=ncol(lmdata)-2
  r=ncol(lmdata)
  outreg=stats::lm(f,data=lmdata)
  
  if (!is.null(confounder)){
    if (is.null(C)){
      confounder=cbind(confounder)
      C=colMeans(confounder)
    }
  }
  if (is.null(M)){
    M=mean(mediator)
  }
  
  
  #Controlled Direct Effect
  CDE=(outreg$coefficients[[2]]+outreg$coefficients[[4]]*M)*(X-Xstar)
  
  
  #Natural Direct Effect
  NDE=outreg$coefficients[[2]]*(X-Xstar)+outreg$coefficients[[4]]*(X-Xstar)*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg$Aest[c(3:k)]*C))
  
  
  #Indirect Effect
  if (is.null(confounder)){
    IE=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
  } else{
    IE=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
  }
  
  #Interactive Reference Effect
  Intref=NDE-CDE
  
  #Pure Indirect Effect
  if (is.null(confounder)){
    PIE=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
  } else{
    PIE=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
  }
  
  #Interactive Mediation Effect
  Intmed=IE-PIE
  
  #Interactive Effect
  OvInt=Intmed+Intref
  
  #Total Effect
  TE=NDE+IE
  
  #Proportion Mediated
  PM=IE/(IE+NDE)
  
  #Proportion Interacted
  PI=(Intref+Intmed)/TE
  
  #Proportion Eliminated
  PE=(TE-CDE)/TE
  
  if (error=='Delta'){
    #Standard Errors with Delta Method
    if (is.null(confounder)){
      GamCDE=c(0,0,0,1,M,0)
      
      d1=outreg$coefficients[[4]]*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)
      d2=outreg$coefficients[[4]]*Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)
      d4=0
      d5=1
      d6=0
      d7=exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)
      GamNDE=c(d1,d2,d4,d5,d6,d7)
      
      i1=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      i2=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)*(X*exp(medreg$Aest[1]+medreg$Aest[2]*X)-Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      i4=0
      i5=0
      i6=(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      i7=X*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      GamIE=c(i1,i2,i4,i5,i6,i7)
      
      r1=outreg$coefficients[[4]]*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      r2=outreg$coefficients[[4]]*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))*Xstar
      r4=0
      r5=0
      r6=0
      r7=(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))-M
      GamIntref=c(r1,r2,r4,r5,r6,r7)
      
      p1=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      p2=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)*(X*exp(medreg$Aest[1]+medreg$Aest[2]*X)-Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      p4=0
      p5=0
      p6=(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      p7=(Xstar)*(exp(medreg$Aest[1]+medreg$Aest[2]*X)-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      GamPIE=c(p1,p2,p4,p5,p6,p7)
    } else {
      
      GamCDE=c(0,0,0*C,0,1,0,M,0*C)
      
      d1=outreg$coefficients[[4]]*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg$Aest[c(3:k)]*C))
      d2=outreg$coefficients[[4]]*Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg$Aest[c(3:k)]*C))
      d3=outreg$coefficients[[4]]*C*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg$Aest[c(3:k)]*C))
      d4=0
      d5=1
      d6=0
      d7=exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg$Aest[c(3:k)]*C))
      d8=C*0
      GamNDE=c(d1,d2,d3,d4,d5,d6,d7,d8)
      
      i1=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      i2=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)*(X*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      i3=C*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      i4=0
      i5=0
      i6=(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      i7=X*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      i8=C*0
      GamIE=c(i1,i2,i3,i4,i5,i6,i7,i8)
      
      r1=outreg$coefficients[[4]]*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      r2=outreg$coefficients[[4]]*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))*Xstar
      r3=outreg$coefficients[[4]]*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))*C
      r4=0
      r5=0
      r6=0
      r7=(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))-M
      r8=C*0
      GamIntref=c(r1,r2,r3,r4,r5,r6,r7,r8)
      
      p1=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      p2=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)*(X*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-Xstar*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      p3=(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))*C
      p4=0
      p5=0
      p6=(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      p7=(Xstar)*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:k)]*C))-exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:k)]*C)))
      p8=C*0
      GamPIE=c(p1,p2,p3,p4,p5,p6,p7,p8)
      
    }
    GamIntmed=GamIE-GamPIE
    GamTE=GamNDE*(X-Xstar)+GamIE
    GamInt=GamIntmed+GamIntref
    
    #Covariance Matrices
    MZIPCov=medreg$AlphaCov
    lmCov=stats::vcov(outreg) #uses stats package
    
    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nMZI,nlm)
    Botmleft=matrix(0,nlm,nMZI)
    
    Top=cbind(MZIPCov,Topright)
    Btm=cbind(Botmleft,lmCov)
    
    CovM=rbind(Top,Btm)
    
    #Calculate Standard Errors
    CDEse=sqrt(GamCDE %*% CovM %*% GamCDE)*abs(X-Xstar)
    NDEse=sqrt(GamNDE %*% CovM %*% GamNDE)*abs(X-Xstar)
    IEse=sqrt(GamIE %*% CovM %*% GamIE)
    Intrefse=sqrt(GamIntref %*% CovM %*% GamIntref)*abs(X-Xstar)
    PIEse=sqrt(GamPIE %*% CovM %*% GamPIE)
    Intmedse=sqrt(GamIntmed %*% CovM %*% GamIntmed )
    TEse=sqrt(GamTE %*% CovM %*% GamTE)
    OvIntse=sqrt(GamInt %*% CovM %*% GamInt)
    
    #Confidence Intervals
    LCDEci=CDE-1.96*CDEse
    UCDEci=CDE+1.96*CDEse
    CDECI=c(LCDEci,UCDEci)
    
    LNDEci=NDE-1.96*NDEse
    UNDEci=NDE+1.96*NDEse
    NDECI=c(LNDEci,UNDEci)
    
    LIEci=IE-1.96*IEse
    UIEci=IE+1.96*IEse
    IECI=c(LIEci,UIEci)
    
    LIntrefci=Intref-1.96*Intrefse
    UIntrefci=Intref+1.96*Intrefse
    IntrefCI=c(LIntrefci,UIntrefci)
    
    LPIEci=PIE-1.96*PIEse
    UPIEci=PIE+1.96*PIEse
    PIECI=c(LPIEci,UPIEci)
    
    LIntmedci=Intmed-1.96*Intmedse
    UIntmedci=Intmed+1.96*Intmedse
    IntmedCI=c(LIntmedci,UIntmedci)
    
    LTEci=TE-1.96*TEse
    UTEci=TE+1.96*TEse
    TECI=c(LTEci,UTEci)
    
    LOvIntci=OvInt-1.96*OvIntse
    UOvIntci=OvInt+1.96*OvIntse
    OvIntCI=c(LOvIntci,UOvIntci)
  }
  
  if (error=='Boot'){
    datab=list()
    datab2=list()
    outregb=list()
    medregb=list()
    CDEb=list()
    NDEb=list()
    IEb=list()
    confb=list()
    Intrefb=list()
    TEb=list()
    Intmedb=list()
    Intb=list()
    PIEb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(lmdata),replace=T)
      datab2[[i]]=lmdata[datab[[i]],]
      outregb[[i]]=stats::lm(f,data=datab2[[i]])
      if (is.null(confounder)){
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=datab2[[i]][["exposure"]],print=F)
        IEb[[i]]=(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[4]*X)*(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X)-exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar))
        PIEb[[i]]=(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[4]*Xstar)*(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X)-exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar))
      } else {
        confb[[i]]=data.matrix(datab2[[i]][c(5:r)])
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=cbind(datab2[[i]][["exposure"]],confb[[i]]),print=F)
        IEb[[i]]=(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*X)*(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+sum(medregb[[i]][["Aest"]][c(3:k)]*C))-exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:k)]*C)))
        PIEb[[i]]=(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)*(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+sum(medregb[[i]][["Aest"]][c(3:k)]*C))-exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:k)]*C)))
      }
      CDEb[[i]]=(outregb[[i]]$coefficients[[2]]+outregb[[i]]$coefficients[[4]]*M)*(X-Xstar)
      NDEb[[i]]=outregb[[i]]$coefficients[[2]]*(X-Xstar)+outregb[[i]]$coefficients[[4]]*(X-Xstar)*exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]]$Aest[c(3:k)]*C))
      Intrefb[[i]]=NDEb[[i]]-CDEb[[i]]
      TEb[[i]]=NDEb[[i]]+IEb[[i]]
      Intmedb[[i]]=IEb[[i]]-PIEb[[i]]
      
      Intb[[i]]=Intmedb[[i]]+Intrefb[[i]]
    }
    
    #Controlled Direct Effect
    CDEse=matrixStats::colSds(do.call(rbind,CDEb))
    CDECI=stats::quantile(do.call(rbind,CDEb),c(0.025,.975))
    
    #Natural Direct Effect
    NDEse=matrixStats::colSds(do.call(rbind,NDEb))
    NDECI=stats::quantile(do.call(rbind,NDEb),c(0.025,.975))
    
    #Indirect Effect 
    IEse=matrixStats::colSds(do.call(rbind,IEb))
    IECI=stats::quantile(do.call(rbind,IEb),c(0.025,.975))
    
    #Interactive Reference Effect
    Intrefse=matrixStats::colSds(do.call(rbind,Intrefb))
    IntrefCI=stats::quantile(do.call(rbind,Intrefb),c(0.025,.975))
    
    #Pure Indirect Effect
    PIEse=matrixStats::colSds(do.call(rbind,PIEb))
    PIECI=stats::quantile(do.call(rbind,PIEb),c(0.025,.975))
    
    #Interactive Mediation Effect
    Intmedse=matrixStats::colSds(do.call(rbind,Intmedb))
    IntmedCI=stats::quantile(do.call(rbind,Intmedb),c(0.025,.975))
    
    #Overall Interactive Effect
    OvIntse=matrixStats::colSds(do.call(rbind,Intb))
    OvIntCI=stats::quantile(do.call(rbind,Intb),c(0.025,.975))
    
    #Total Effect
    TEse=matrixStats::colSds(do.call(rbind,TEb))
    TECI=stats::quantile(do.call(rbind,TEb),c(0.025,.975))
  }
  
  output=list(lm=outreg,mzip=medreg,CDE=CDE,NDE=NDE,NIE=IE,CDEse=CDEse,CDECI=CDECI,NDEse=NDEse,NDECI=NDECI,NIEse=IEse,NIECI=IECI,
              Intref=Intref,Intrefse=Intrefse,IntrefCI=IntrefCI,Intmed=Intmed,Intmedse=Intmedse,IntmedCI=IntmedCI,
              PIE=PIE,PIEse=PIEse,PIECI=PIECI,TE=TE,TEse=TEse,TECI=TECI,
              OvInt=OvInt,OvIntse=OvIntse,OvIntCI=OvIntCI,PM=PM,PropInt=PI, PE=PE)
}




#' Mediation Analysis for Zero-Inflated Count Mediators using MZIP (Binary or Count Outcome)
#'
#' This function incorporates the MZIP model into the counterfactual approach to mediation analysis
#' as proposed by Vanderweele when the mediator is a Zero-Inflated count variable for cases with
#' binary or count outcome using a Poisson regression with robust standard errors. Standard Errors for
#' direct and indirect effects are computed using delta method or bootstrapping. Note: This function
#' assumes that the outcome is continuous and all exposure, mediator, outcome, and confounder variables
#' have the same sample size. Binary variables must be dummy coded prior.
#' A Poisson regression with robust standard errors were used to obtain direct and indirect
#' effect estimates on a risk ratio scale because odds ratios are a non-collapsible measure which
#' can cause issues in a mediation framework (see Vanderweele 2016).
#' @param outcome is the binary or count outcome variable
#' @param mediator is the zero-inflated mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetition if bootstrapped errors are used. Default is 1000
#' @param C is a vector for theoretical values of each confounder. By default each each value of C will be the mean value of each confounder.
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @return The function will return a list of 12 elements.
#'     GLM is the results of regressing the mediator+exposure+confounder on the outcome using a Poisson model with robust standard errors \cr
#'     MZIP is the results of regressing the exposure and confounders on the mediator using the MZIP model \cr
#'     RRNDE is the risk ratio of the direct effect \cr
#'     RRNIE is the risk ratio of the indirect effect. \cr
#'     logRRNDEse is the standard error for the log risk ratio of NDE \cr
#'     RRNDEci is the 95% confidence interval for the direct effect risk ratio\cr
#'     logRRNIEse is the standard error for  the indirect effect log risk ratio \cr
#'     RRNIEci is the 95% confidence interval for the indirect effect risk ratio \cr
#'     RRTE is the total effect risk ratio \cr
#'     logRRTEse is the standard error for the total effect log risk ratio\cr
#'     RRTECI is the confidence interval for the total effect risk ratio \cr
#'     PM is the proportion mediated
#' @examples
#'     binoutzimed(outcome=BinaryOutcome,mediator=ZICount,exposure=race,n=1000,error='Boot')
#'     binoutzimed(outcome=data$outcome,mediator=data$mediator,exposure=data$exp,confounder=cbind(data$var1,data$var2),X=10,Xstar=0,C=c(1,3))
#'     binoutzimed(outcome=countoutcome,mediator=mzipmediator,exposure=race,confounder=sex,C=0,error='Delta')
#' @export


binoutzimed=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,error='Delta'){
  glmout=data.frame(outcome)
  if (is.null(confounder)){
    glmpred=data.frame(exposure,mediator)
    medreg=mzip(y=mediator,pred=exposure,print=F)
  } else{
    glmpred=data.frame(exposure,mediator,confounder)
    medreg=mzip(y=mediator,pred=cbind(exposure,confounder),print=F)
  }
  #as.formula and lm part of stats package
  glmdata=data.frame(glmout,glmpred)
  f<-stats::as.formula(paste(colnames(glmout),paste(colnames(glmpred),collapse="+"),sep="~"))
  
  m=ncol(glmdata)-1
  r=ncol(glmdata)
  outreg=robust::glmRob(f,data=glmdata,family=poisson())
  
  if (!is.null(confounder)){
    if (is.null(C)){
      confounder=cbind(confounder)
      C=colMeans(confounder)
    }
  }
  
  #Direct Effect
  RRDE=exp(outreg$coefficients[[2]]*(X-Xstar))
  
  #Indirect effect
  if (is.null(confounder)){
    RRIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*X)+
                                                          exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]])-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                   exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]])-1))))
  } else{
    RRIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                          exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                   exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))))
  }
  #Proportion Mediated
  RRTE=RRIE*RRDE
  PM=(RRDE*(RRIE-1))/(RRIE*RRDE-1)
  
  if (error=='Delta'){
    #Standard Errors Delta Method
    
    #Gamma for delta involves gamma coefficients from MZIP
    if (is.null(confounder)){
      GamDE=c(0,0,0,0,0,X-Xstar,0)
      P=exp(medreg$Gest[1]+medreg$Gest[2]*X)
      S=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)
      V=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]])-1))
      W=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]])-1))
      H=P+V
      J=S+W
      
      r1=(1/(1+S))-(1/(1+P))+((P-(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H)-
        ((S-(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J)
      r2=Xstar*((1/(1+S))-((S-(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J))+
        X*(-(1/(1+P))+((P-(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H))
      a1=-(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)+
        (((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      a2=-X*(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)+
        Xstar*(((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      t1=0
      t2=0
      t3=-(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]))/H)
      +(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]))/J)
      
      
      GamIE=c(r1,r2,a1,a2,t1,t2,t3)
      
    } else {
      GamDE=c(0,0,0*C,0,0,0*C,0,X-Xstar,0,0*C)
      
      P=exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))
      S=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))
      V=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))
      W=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))
      H=P+V
      J=S+W
      
      r1=(1/(1+S))-(1/(1+P))+((P-(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H)-
        ((S-(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J)
      r2=Xstar*((1/(1+S))-((S-(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J))+
        X*(-(1/(1+P))+((P-(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H))
      r3=C*r1
      a1=-(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)+
        (((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      a2=-X*(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)+
        Xstar*(((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      a3=C*a1
      t1=0
      t2=0
      t3=-(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]))/H)+
        +(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]))/J)
      t4=C*0
      
      
      GamIE=c(r1,r2,r3,a1,a2,a3,t1,t2,t3,t4)
    }
    
    GamTE=GamIE+GamDE
    #Extract Covariance matrices
    MZIPCov=medreg$Cov
    lmCov=outreg$cov #uses stats package
    
    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nMZI,nlm)
    Botmleft=matrix(0,nlm,nMZI)
    
    Top=cbind(MZIPCov,Topright)
    Btm=cbind(Botmleft,lmCov)
    
    CovM=rbind(Top,Btm)
    
    logRRDEse=sqrt(GamDE %*% CovM %*% GamDE)
    logRRIEse=sqrt(GamIE %*% CovM %*% GamIE)
    logRRTEse=sqrt(GamTE %*% CovM %*% GamTE)
    
    LRRDEci=exp(log(RRDE)-1.96*logRRDEse)
    URRDEci=exp(log(RRDE)+1.96*logRRDEse)
    RRDECI=c(LRRDEci,URRDEci)
    
    LRRIEci=exp(log(RRIE)-1.96*logRRIEse)
    URRIEci=exp(log(RRIE)+1.96*logRRIEse)
    RRIECI=c(LRRIEci,URRIEci)
    
    LRRTEci=exp(log(RRTE)-1.96*logRRTEse)
    URRTEci=exp(log(RRTE)+1.96*logRRTEse)
    RRTECI=c(LRRTEci,URRTEci)
  }
  
  if (error=='Boot'){
    datab=list()
    datab2=list()
    outregb=list()
    medregb=list()
    RRDEb=list()
    RRIEb=list()
    confb=list()
    RRTEb=list()
    logRRDEb=list()
    logRRIEb=list()
    logRRTEb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(glmdata),replace=T)
      datab2[[i]]=glmdata[datab[[i]],]
      outregb[[i]]=robust::glmRob(f,data=datab2[[i]],family=poisson())
      if (is.null(confounder)){
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=datab2[[i]][["exposure"]],print=F)
        RRIEb[[i]]=((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X)+
                                                                                exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X))))*(exp(outregb[[i]]$coefficients[[3]])-1))))/
          ((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar)+
                                                                   exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))))*(exp(outregb[[i]]$coefficients[[3]])-1))))
      } else {
        confb[[i]]=data.matrix(datab2[[i]][c(4:r)])
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=cbind(datab2[[i]][["exposure"]],confb[[i]]),print=F)
        RRIEb[[i]]=((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                                                                                      exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]])-1))))/
          ((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                                                                         exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]])-1))))
      }
      RRDEb[[i]]=exp(outregb[[i]]$coefficients[[2]]*(X-Xstar))
      RRTEb[[i]]=RRIEb[[i]]*RRDEb[[i]]
      
      logRRDEb[[i]]=log(RRDEb[[i]])
      logRRIEb[[i]]=log(RRIEb[[i]])
      logRRTEb[[i]]=log(RRTEb[[i]])
    }
    
    #quantile part of stats package. colSds part of matrixStats package (not base)
    #Risk Ratio Direct Effect
    logRRDEse=matrixStats::colSds(do.call(rbind,logRRDEb))
    RRDECI=stats::quantile(do.call(rbind,RRDEb),c(0.025,.975))
    
    #Risk Ratio Indirect Effect
    logRRIEse=matrixStats::colSds(do.call(rbind,logRRIEb))
    RRIECI=stats::quantile(do.call(rbind,RRIEb),c(0.025,.975))
    
    logRRTEse=matrixStats::colSds(do.call(rbind,logRRTEb))
    RRTECI=stats::quantile(do.call(rbind,RRTEb),c(0.025,.975))
  }
  
  output=list(GLM=outreg,MZIP=medreg,RRNDE=RRDE,RRNIE=RRIE,PM=PM,logRRNDEse=logRRDEse,RRNDEci=RRDECI,logRRNIEse=logRRIEse,RRNIEci=RRIECI,
              RRTE=RRTE,logRRTEse=logRRTEse,RRTEci=RRTECI)
}



#' Mediation Analysis for Zero-Inflated Count Mediators using MZIP with Exposure-Mediator Interactions (Binary/Count Outcome)
#'
#' This function will do the same thing as the binoutzimed function, but includes an exposure-mediator interaction.
#' 4-way decomposition of total effect (Vanderweele) are included in the output.
#' @param outcome is the continuous outcome variable
#' @param mediator is the zero-inflated mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetitions for bootstrapping. Default is 1000. Setting n when using delta method errors will have no effect on output.
#' @param C is a vector for theoretical values of each confounder. If left out the default will be set to the mean of each confounder giving marginal effects
#' @param M is a fixed value for the mediator, M. If M is not specified, M will be set to its mean value 
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @return The function will return a list of 34 elements.
#'     GLM is the results of regressing the mediator+exposure+confounder on the outcome using a Poisson model with robust standard errors. To assess interaction effect individually look in the glm statement at the 4th parameter estimate \cr
#'     MZIP is the results of regressing the exposure and confounders on the mediator using the MZIP model \cr
#'     RRCDE is the controlled direct effect risk ratio \cr
#'     RRNDE is the natural direct effect risk ratio \cr
#'     RRNIE is the indirect effect risk ratio. \cr
#'     PM is the proportion mediated\cr
#'     logRRCDEse is the standard error for the  controlled direct effect log risk ratio \cr
#'     RRCDEci is the 95% confidence interval for the controlled direct effect risk raito\cr
#'     logRRNDEse is the standard error for the  natural direct effect log risk ratio \cr
#'     RRNDEci is the 95% confidence interval for the natural direct effect risk ratio\cr
#'     logRRNIEse is the standard error for  the indirect effect log risk ratio \cr
#'     RRNIEci is the 95% confidence interval for the indirect effect risk ratio\cr
#'     Intref is the Interactive Reference effect (not a risk ratio) \cr
#'     Intrefse is the standard error for Intref \cr
#'     IntrefCI is the CI for Intref \cr
#'     RRPIE is the pure indirect effect risk ratio \cr
#'     logRRPIEse is the standard error of PIE log risk ratio \cr
#'     RRPIECI is the CI for PIE risk ratio \cr
#'     Intmed is the interactive mediation effect (not a risk ratio) \cr
#'     Intmedse is the error associated with Intmed \cr
#'     IntmedCI is the CI for Intmed \cr
#'     RRTE is the total effect risk ratio \cr
#'     logRRTEse is the error of the total effect log risk ratio \cr
#'     RRTECI is the CI for the total effect risk ratio\cr
#'     OvInt is the overall additive interaction effect \cr
#'     OvIntse is the standard error for the additive interaction \cr
#'     OvIntCI is the confidence interval for the interaction effect \cr
#'     PAINT is the proportion attributable to the interaction effect \cr
#'     PE is the proportion eliminated \cr
#'     PACDE is the proportion of the total effect due to neither mediation nor interaction \cr
#'     PAIntref is the proportion of the total effect due to just interaction \cr
#'     PAIntmed is the proportion of the total effect attributable to the joint effect of mediation and interaction \cr
#'     PAPIE is the proportion of the total effect attributable to just mediation \cr
#'     terr is the total excess relative risk 
#' @examples
#'     binoutzimedint(outcome=ContinuousOutcome,mediator=ZICount,exposure=race,n=200,error="Boot")
#'     binoutzimedint(outcome=data$outcome,mediator=data$mediator,exposure=data$exp,confounder=cbind(data$var1,data$var2),X=10,Xstar=0,C=c(1,3),M=100)
#'     binoutzimedint(outcome=countoutcome,mediator=zimediator,exposure=race,confounder=sex,C=0,M=0,error='Delta')
#' @export

binoutzimedint=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,M=NULL,error='Delta'){
  #lm, as.formula, quantile in stats, colSds in matrixStats
  interaction=mediator*exposure
  glmout=data.frame(outcome)
  if (is.null(confounder)){
    glmpred=data.frame(exposure,mediator,interaction)
    medreg=mzip(y=mediator,pred=exposure,print=F)
  } else{
    glmpred=data.frame(exposure,mediator,interaction,confounder)
    medreg=mzip(y=mediator,pred=cbind(exposure,confounder),print=F)
  }
  
  glmdata=data.frame(glmout,glmpred)
  f<-stats::as.formula(paste(colnames(glmout),paste(colnames(glmpred),collapse="+"),sep="~"))
  
  m=ncol(glmdata)-2
  r=ncol(glmdata)
  outreg=robust::glmRob(f,data=glmdata,family=poisson())
  
  if (!is.null(confounder)){
    if (is.null(C)){
      confounder=cbind(confounder)
      C=colMeans(confounder)
    }
  }
  if (is.null(M)){
    M=mean(mediator)
  }
  
  #Controlled Direct Effect
  RRCDE=exp((outreg$coefficients[[2]]+outreg$coefficients[[4]]*M)*(X-Xstar))
  
  
  #Natural Direct Effect
  if (is.null(confounder)){
    RRNDE=((exp(outreg$coefficients[[2]]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((exp(outreg$coefficients[[2]]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  } else{
    RRNDE=((exp(outreg$coefficients[[2]]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((exp(outreg$coefficients[[2]]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  }
  
  #Indirect Effect
  if (is.null(confounder)){
    RRIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*X)+
                                                          exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                   exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))
  } else{
    RRIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                          exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                   exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))
  }
  
  
  
  #Pure Indirect Effect
  if (is.null(confounder)){
    RRPIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*X)+
                                                           exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                   exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  } else{
    RRPIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                           exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                   exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  }
  
  
  #Proportion Mediated
  PM=RRNDE*(RRIE-1)/(RRIE*RRNDE-1)
  
  #kappa
  kappa=(exp(outreg$coefficients[[3]]*M+outreg$coefficients[[4]]*Xstar*M)*(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/
    ((exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))+exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)))
  
  
  #Total Effect
  RRTE=RRNDE*RRIE
  
  
  #Interactive Reference Effect
  RRIntref=((RRNDE-1)/kappa)-RRCDE+1
  
  
  #Interactive Mediation Effect
  RRIntmed=(RRTE-RRNDE-RRPIE+1)/kappa
  
  #Proportion Eliminated
  PE=(RRTE-RRCDE)/(RRTE-1)
  
  terr=kappa*(RRCDE-1)+kappa*RRIntref+kappa*RRIntmed+(RRPIE-1)
  
  PACDE=(kappa*(RRCDE-1))/terr
  PAIntref=(kappa*(RRIntref))/terr
  PAIntmed=(kappa*(RRIntmed))/terr
  PAPIE=(RRPIE-1)/terr
  PAINT=PAIntref+PAIntmed
  
  #Interactive Effect
  RRInt=RRIntmed+RRIntref
  
  if (error=='Delta'){
    #Standard Errors with Delta Method
    if (is.null(confounder)){
      GamCDE=c(0,0,0,0,0,1,M,0)
      
      Q=(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
           exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)))
      R=(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
           exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)))
      B=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)
      D=Q-B
      G=R-B
      
      dr1=((B-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*B*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/Q)-((B-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*B*
                                                               exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/R)
      dr2=dr1*Xstar
      da1=((-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*(B+1)*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/Q)+
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*(B+1)*
            exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/R)
      da2=da1*Xstar
      dt1=0
      dt2=X-Xstar
      dt3=((-(B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)+
        (((B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      dt4=((-X*(B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)+
        ((Xstar*(B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      GamNDE=c(dr1,dr2,da1,da2,dt1,dt2,dt3,dt4)
      
      
      P=exp(medreg$Gest[1]+medreg$Gest[2]*X)
      S=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)
      V=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      W=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      H=P+V
      J=S+W
      
      nr1=(1/(1+S))-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H)-
        ((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J)
      nr2=Xstar*((1/(1+S))-((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J))+
        X*(-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H))
      na1=-(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)+
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      na2=-X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)+
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      nt1=0
      nt2=0
      nt3=-(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/H)
      +(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/J)
      nt4=X*nt3
      
      GamIE=c(nr1,nr2,na1,na2,nt1,nt2,nt3,nt4)
      
      V2=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      W2=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      H2=P+V2
      J2=S+W2
      
      pr1=(1/(1+S))-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H2)-
        ((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J2)
      pr2=Xstar*((1/(1+S))-((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J2))+
        X*(-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H2))
      pa1=-(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H2)+
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J2)
      pa2=-X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H2)+
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J2)
      pt1=0
      pt2=0
      pt3=-(((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/H2)+
        +(((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/J2)
      pt4=Xstar*pt3
      
      GamPIE=c(pr1,pr2,pa1,pa2,pt1,pt2,pt3,pt4)
      
      Z=exp((X-Xstar)*outreg$coefficients[[2]]-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      
      mr1=((-Z*H*P)/((1+P)^2))+(Z*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((Z*J*S)/((1+S)^2))-(Z*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((-exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)
      mr2=((-X*Z*H*P)/((1+P)^2))+(Z*X*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((Xstar*Z*J*S)/((1+S)^2))-(Z*Xstar*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)))/(1+S)+
        ((X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((-Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)
      ma1=(Z/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))-
        (Z/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      ma2= (X*Z/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))-
        (Xstar*Z/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))-
        (X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      mt1=0
      mt2=((X-Xstar)*Z*H/(1+P))-((X-Xstar)*Z*J/(1+S))
      mt3=-M*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+P))*(-(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      mt4=-M*Xstar*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+P))*(-(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z*X/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      
      GamIntmed=c(mr1,mr2,ma1,ma2,mt1,mt2,mt3,mt4)
      
      rr1=(-(Z*J*S)/((1+S)^2))+(Z*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)
      rr2=Xstar*rr1
      ra1=(Z/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      ra2=Xstar*ra1
      rt1=0
      rt2=(Xstar*Z*J/(1+S))-(X-Xstar)*RRCDE
      rt3=-M*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      rt4=-M*Xstar*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      -M*(X-Xstar)*RRCDE
      
      GamIntref=c(rr1,rr2,ra1,ra2,rt1,rt2,rt3,rt4)
      
    } else {
      
      GamCDE=c(0,0,0*C,0,0,0*C,0,1,0,M,0*C)
      
      Q=(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
           exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)))
      R=(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
           exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)))
      B=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))
      D=Q-B
      G=R-B
      
      dr1=((B-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*B*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/Q)-((B-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*B*
                                                                                               exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/R)
      dr2=dr1*Xstar
      dr3=dr1*C
      da1=((-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*(B+1)*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/Q)+
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*(B+1)*
            exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/R)
      da2=da1*Xstar
      da3=da1*C
      dt1=0
      dt2=X-Xstar
      dt3=((-(B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)+
        (((B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      dt4=((-X*(B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)+
        ((Xstar*(B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      dt5=C*0
      GamNDE=c(dr1,dr2,dr3,da1,da2,da3,dt1,dt2,dt3,dt4,dt5)
      
      P=exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))
      S=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))
      V=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      W=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      H=P+V
      J=S+W
      
      nr1=(1/(1+S))-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H)-
        ((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J)
      nr2=Xstar*((1/(1+S))-((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J))+
        X*(-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H))
      nr3=C*nr1
      na1=-(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)+
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      na2=-X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)+
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      na3=C*na1
      nt1=0
      nt2=0
      nt3=-(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/H)
      +(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/J)
      nt4=X*nt3
      nt5=C*0
      
      
      GamIE=c(nr1,nr2,nr3,na1,na2,na3,nt1,nt2,nt3,nt4,nt5)
      
      V2=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      W2=exp(-(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      H2=P+V2
      J2=S+W2
      
      pr1=(1/(1+S))-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H2)-
        ((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J2)
      pr2=Xstar*((1/(1+S))-((S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J2))+
        X*(-(1/(1+P))+((P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H2))
      pr3=C*pr1
      pa1=-(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H2)+
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J2)
      pa2=-X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H2)+
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J2)
      pa3=C*pa1
      pt1=0
      pt2=0
      pt3=-(((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/H2)
      +(((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/J2)
      pt4=Xstar*pt3
      pt5=C*0
      
      GamPIE=c(pr1,pr2,pr3,pa1,pa2,pa3,pt1,pt2,pt3,pt4,pt5)
      
      Z=exp((X-Xstar)*outreg$coefficients[[2]]-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      
      mr1=((-Z*H*P)/((1+P)^2))+(Z*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)+
        ((Z*J*S)/((1+S)^2))-(Z*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)-
        ((-exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)
      mr2=((-X*Z*H*P)/((1+P)^2))+(Z*X*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)+
        ((Xstar*Z*J*S)/((1+S)^2))-(Z*Xstar*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)+
        ((X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)-
        ((-Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)
      mr3=mr1*C
      ma1=(Z/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (Z/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      ma2= (X*Z/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (Xstar*Z/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      ma3=ma1*C
      mt1=0
      mt2=((X-Xstar)*Z*H/(1+P))-((X-Xstar)*Z*J/(1+S))
      mt3=-M*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+P))*(-(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      mt4=-M*Xstar*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+P))*(-(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z*X/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*(-(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      mt5=0*C
      
      GamIntmed=c(mr1,mr2,mr3,ma1,ma2,ma3,mt1,mt2,mt3,mt4,mt5)
      
      rr1=(-(Z*J*S)/((1+S)^2))+(Z*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)
      rr2=Xstar*rr1
      rr3=C*rr1
      ra1=(Z/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      ra2=Xstar*ra1
      ra3=C*ra1
      rt1=0
      rt2=(Xstar*Z*J/(1+S))-(X-Xstar)*RRCDE
      rt3=-M*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      rt4=-M*Xstar*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+S))*(-(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*(-(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      -M*(X-Xstar)*RRCDE
      rt5=C*0
      
      GamIntref=c(rr1,rr2,rr3,ra1,ra2,ra3,rt1,rt2,rt3,rt4,rt5)
      
    }
    GamTE=GamNDE+GamIE
    GamInt=GamIntref+GamIntmed
    
    #Covariance Matrices
    MZIPCov=medreg$Cov
    lmCov=outreg$cov #uses stats package
    
    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nMZI,nlm)
    Botmleft=matrix(0,nlm,nMZI)
    
    Top=cbind(MZIPCov,Topright)
    Btm=cbind(Botmleft,lmCov)
    
    CovM=rbind(Top,Btm)
    
    #Now calculate standard errors
    logRRCDEse=sqrt(GamCDE %*% CovM %*% GamCDE)*(X-Xstar)
    logRRNDEse=sqrt(GamNDE %*% CovM %*% GamNDE)
    logRRIEse=sqrt(GamIE %*% CovM %*% GamIE)
    RRIntrefse=sqrt(GamIntref %*% CovM %*% GamIntref)
    RRIntmedse=sqrt(GamIntmed %*% CovM %*% GamIntmed)
    logRRPIEse=sqrt(GamPIE %*% CovM %*% GamPIE)
    logRRTEse=sqrt(GamTE %*% CovM %*% GamTE)
    RRIntse=sqrt(GamInt %*% CovM %*% GamInt)
    
    LRRCDEci=exp(log(RRCDE)-1.96*logRRCDEse)
    URRCDEci=exp(log(RRCDE)+1.96*logRRCDEse)
    RRCDECI=c(LRRCDEci,URRCDEci)
    
    LRRNDEci=exp(log(RRNDE)-1.96*logRRNDEse)
    URRNDEci=exp(log(RRNDE)+1.96*logRRNDEse)
    RRNDECI=c(LRRNDEci,URRNDEci)
    
    LRRIEci=exp(log(RRIE)-1.96*logRRIEse)
    URRIEci=exp(log(RRIE)+1.96*logRRIEse)
    RRIECI=c(LRRIEci,URRIEci)
    
    LRRIntrefci=(RRIntref)-1.96*RRIntrefse
    URRIntrefci=(RRIntref)+1.96*RRIntrefse
    RRIntrefCI=c(LRRIntrefci,URRIntrefci)
    
    LRRIntmedci=(RRIntmed)-1.96*RRIntmedse
    URRIntmedci=(RRIntmed)+1.96*RRIntmedse
    RRIntmedCI=c(LRRIntmedci,URRIntmedci)
    
    LRRPIEci=exp(log(RRPIE)-1.96*logRRPIEse)
    URRPIEci=exp(log(RRPIE)+1.96*logRRPIEse)
    RRPIECI=c(LRRPIEci,URRPIEci)
    
    LRRTEci=exp(log(RRTE)-1.96*logRRTEse)
    URRTEci=exp(log(RRTE)+1.96*logRRTEse)
    RRTECI=c(LRRTEci,URRTEci)
    
    LRRIntci=(RRInt)-1.96*RRIntse
    URRIntci=(RRInt)+1.96*RRIntse
    RRIntCI=c(LRRIntci,URRIntci)
  }
  
  if (error=='Boot'){
    datab=list()
    datab2=list()
    outregb=list()
    medregb=list()
    RRCDEb=list()
    RRNDEb=list()
    RRIEb=list()
    confb=list()
    RRIntrefb=list()
    RRTEb=list()
    RRIntmedb=list()
    RRIntb=list()
    RRPIEb=list()
    kappab=list()
    logRRNDEb=list()
    logRRIEb=list()
    logRRTEb=list()
    logRRCDEb=list()
    logRRPIEb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(glmdata),replace=T)
      datab2[[i]]=glmdata[datab[[i]],]
      outregb[[i]]=robust::glmRob(f,data=datab2[[i]],family=poisson())
      if (is.null(confounder)){
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=datab2[[i]][["exposure"]],print=F)
        RRNDEb[[i]]=((exp(outregb[[i]]$coefficients[[2]]*X))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar)+
                                                                exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*X)-1))))/
          ((exp(outregb[[i]]$coefficients[[2]]*Xstar))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar)+
                                                          exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1))))
        
        RRIEb[[i]]=((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X)+
                                                                                exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*X)-1))))/
          ((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar)+
                                                                   exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*X)-1))))
        RRPIEb[[i]]=((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X)+
                                                                                 exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1))))/
          ((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar)+
                                                                   exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1))))
        kappab[[i]]=(exp(outregb[[i]]$coefficients[[3]]*M+outregb[[i]]$coefficients[[4]]*Xstar*M)*(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar)))/
          ((exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))+exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1)))
      } else {
        confb[[i]]=data.matrix(datab2[[i]][c(5:r)])
        medregb[[i]]=mzip(y=datab2[[i]][["mediator"]],pred=cbind(datab2[[i]][["exposure"]],confb[[i]]),print=F)
        RRNDEb[[i]]=((exp(outregb[[i]]$coefficients[[2]]*X))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                                exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*X)-1))))/
          ((exp(outregb[[i]]$coefficients[[2]]*Xstar))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                          exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1))))
        
        RRIEb[[i]]=((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                                                                                      exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*X)-1))))/
          ((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                                                                         exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*X)-1))))
        RRPIEb[[i]]=((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                                                                                       exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*X+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1))))/
          ((1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*X+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))*(exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C))+
                                                                                                         exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1))))
        kappab[[i]]=(exp(outregb[[i]]$coefficients[[3]]*M+outregb[[i]]$coefficients[[4]]*Xstar*M)*(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C))))/
          ((exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))+exp(-(exp(medregb[[i]]$Aest[1]+medregb[[i]]$Aest[2]*Xstar+sum(medregb[[i]][["Aest"]][c(3:m)]*C)+log(1+exp(medregb[[i]]$Gest[1]+medregb[[i]]$Gest[2]*Xstar+sum(medregb[[i]][["Gest"]][c(3:m)]*C)))))*(exp(outregb[[i]]$coefficients[[3]]+outregb[[i]]$coefficients[[4]]*Xstar)-1)))
        
      }
      RRCDEb[[i]]=exp((outregb[[i]]$coefficients[[2]]+outregb[[i]]$coefficients[[4]]*M)*(X-Xstar))
      RRTEb[[i]]=RRIEb[[i]]*RRNDEb[[i]]
      RRIntrefb[[i]]=((RRNDEb[[i]]-1)/kappab[[i]])-RRCDEb[[i]]+1
      RRIntmedb[[i]]=(RRTEb[[i]]-RRNDEb[[i]]-RRPIEb[[i]]+1)/kappab[[i]]
      RRIntb[[i]]=RRIntmedb[[i]]+RRIntrefb[[i]]
      
      logRRCDEb[[i]]=log(RRCDEb[[i]])
      logRRNDEb[[i]]=log(RRNDEb[[i]])
      logRRIEb[[i]]=log(RRIEb[[i]])
      logRRTEb[[i]]=log(RRTEb[[i]])
      logRRPIEb[[i]]=log(RRPIEb[[i]])
    }
    
    #RRCDE
    logRRCDEse=matrixStats::colSds(do.call(rbind,logRRCDEb))
    RRCDECI=stats::quantile(do.call(rbind,RRCDEb),c(0.025,.975))
    
    #Risk Ratio Natural Direct Effect
    logRRNDEse=matrixStats::colSds(do.call(rbind,logRRNDEb))
    RRNDECI=stats::quantile(do.call(rbind,RRNDEb),c(0.025,.975))
    
    #Risk Ratio Indirect Effect
    logRRIEse=matrixStats::colSds(do.call(rbind,logRRIEb))
    RRIECI=stats::quantile(do.call(rbind,RRIEb),c(0.025,.975))
    
    #Interactive Reference Effect
    RRIntrefse=matrixStats::colSds(do.call(rbind,RRIntrefb))
    RRIntrefCI=stats::quantile(do.call(rbind,RRIntrefb),c(0.025,.975))
    
    #Pure Indirect Effect
    logRRPIEse=matrixStats::colSds(do.call(rbind,logRRPIEb))
    RRPIECI=stats::quantile(do.call(rbind,RRPIEb),c(0.025,.975))
    
    #Interactive Mediation Effect
    RRIntmedse=matrixStats::colSds(do.call(rbind,RRIntmedb))
    RRIntmedCI=stats::quantile(do.call(rbind,RRIntmedb),c(0.025,.975))
    
    #Total Effect
    logRRTEse=matrixStats::colSds(do.call(rbind,logRRTEb))
    RRTECI=stats::quantile(do.call(rbind,RRTEb),c(0.025,.975))
    
    #Interaction
    RRIntse=matrixStats::colSds(do.call(rbind,RRIntb))
    RRIntCI=stats::quantile(do.call(rbind,RRIntb),c(0.025,.975))
  }
  
  output=list(GLM=outreg,MZIP=medreg,RRCDE=RRCDE,RRNDE=RRNDE,RRNIE=RRIE,logRRCDEse=logRRCDEse,
              logRRNDEse=logRRNDEse,logRRNIEse=logRRIEse,RRCDEci=RRCDECI,RRNDEci=RRNDECI,RRNIEci=RRIECI, PM=PM,
              Intref=RRIntref,Intrefse=RRIntrefse,Intrefci=RRIntrefCI, Intmed=RRIntmed,Intmedse=RRIntmedse,
              Intmedci=RRIntmedCI,RRPIE=RRPIE,logRRPIEse=logRRPIEse,RRPIEci=RRPIECI,
              RRTE=RRTE,logRRTEse=logRRTEse,RRTEci=RRTECI,OvInt=RRInt,
              OvIntse=RRIntse,OvIntci=RRIntCI,PE=PE,PACDE=PACDE,PAIntref=PAIntref,PAIntmed=PAIntmed,PAPIE=PAPIE,PAINT=PAINT,terr=terr)
}
