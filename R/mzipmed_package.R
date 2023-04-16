#MZIP Estimation Function




#' Marginalized Zero-Inflated Poisson Regression Model
#'
#' This function uses the MZIP model to allow you to fit counts variables
#'    with excess zeroes
#'    while allowing for easy interpretations. This function assumes that
#'    the outcome and covariates are all the same sample size without missing
#'    data. Covariates must be numerical, so binary predictors such as
#'    gender or race need to be dummy coded with zeroes and ones. For more
#'    information about this model and interpretations see Long, D Leann et al. "A marginalized
#'    zero-inflated Poisson regression model with overall exposure effects." Statistics in
#'    medicine vol. 33,29 (2014): 5151-65. doi:10.1002/sim.6293.
#'    Note: BFGS likelihood optimization was used for this R package
#' @param y is the outcome variable
#' @param pred is a vector of covariates (use cbind for multiple)
#' @param print if =TRUE will give model parameters estimates and overall mean relative risks. Default =FALSE
#' @return The function will return a list of 22 elements.
#'     In the list G(Gamma) refers to the excess zero/logistic part of the model \cr
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
#'     Cov is the covariance matrix for the MZIP model \cr
#'     RobAlphaCov robust covariance matrix for the Poisson component of MZIP \cr
#'     RobCov robust covariance matrix
#' @examples
#'     test=mzip(y=mzipmed_data$ziY1,pred=cbind(mzipmed_data$X,mzipmed_data$C1,
#'               mzipmed_data$C2),print=FALSE)
#'
#'    \dontrun{
#'    test= mzip(y=mzipmed_data$ziY1,pred=cbind(X=mzipmed_data$X,C1=mzipmed_data$C1,
#'               C2=mzipmed_data$C2),print=TRUE)
#'               }
#' @export



mzip = function(y,pred,print=FALSE){

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

  # diag_gg	= (psi_hat^2*(1-psi_hat)*(psi_hat*psi_hat2*nu_hat+1)*(exp(nu_hat*psi_hat2)-nu_hat*(psi_hat2)-1))/(psi_hat*exp(nu_hat*psi_hat2)+(1-psi_hat))
  #
  # diag_aa = (nu_hat*(psi_hat*(exp(nu_hat*psi_hat2)-nu_hat*psi_hat2-1)+1))/(psi_hat*exp(nu_hat*psi_hat2)+(1-psi_hat))
  #
  # diag_ga = nu_hat*psi_hat*(1-exp(-nu_hat*psi_hat2)+(-(1+nu_hat*psi_hat*psi_hat2+nu_hat*(psi_hat*psi_hat2)^2+exp(-nu_hat*psi_hat2))/((psi_hat*psi_hat2)*exp(nu_hat*psi_hat2)+1)))
  #
  # #test for error
  # dd_gg=diag(as.vector(diag_gg))
  # dd_aa=diag(as.vector(diag_aa))
  #
  # tx_dd_gg=t(Z)%*%diag(as.vector(diag_gg))
  # dd_gg_x=diag(as.vector(diag_gg))%*%(Z)
  #
  # I_gg			= t(Z)%*%diag(as.vector(diag_gg))%*%(Z)
  # I_aa			= t(X)%*%diag(as.vector(diag_aa))%*%(X)
  # I_ga			= t(X)%*%diag(as.vector(diag_ga))%*%(Z)
  # I_ag			= t(I_ga)
  #
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
  RobaCov=robust[w:k,w:k]
  m_se			= sqrt(diag(Inv_inform))

  r_se			= sqrt(diag(robust))

  mupper = c(gamma_hat,alpha_hat) + 1.96*m_se
  mlower = c(gamma_hat,alpha_hat) - 1.96*m_se
  rupper = c(gamma_hat,alpha_hat) + 1.96*r_se
  rlower = c(gamma_hat,alpha_hat) - 1.96*r_se

  GWald=(gamma_hat/m_se[1:dim(Z)[2]])^2
  GPval=ifelse(1-stats::pchisq(q=(gamma_hat/m_se[1:dim(Z)[2]])^2,df=1)<.0001,"<.0001",round(1-stats::pchisq(q=(gamma_hat/m_se[1:dim(Z)[2]])^2,df=1),digits=5))
  GRobWald=(gamma_hat/r_se[1:dim(Z)[2]])^2
  GRobPval=ifelse(1-stats::pchisq(q=(gamma_hat/r_se[1:dim(Z)[2]])^2,df=1)<.0001,"<.0001",round(1-stats::pchisq(q=(gamma_hat/r_se[1:dim(Z)[2]])^2,df=1),digits=5))
  AWald=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2
  APval=ifelse(1-stats::pchisq(q=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1)<.0001,"<.0001",round(1-stats::pchisq(q=(alpha_hat/m_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1),digits=5))
  ARobWald=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2
  ARobPval=ifelse(1-stats::pchisq(q=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1)<.0001,"<.0001",round(1-stats::pchisq(q=(alpha_hat/r_se[(dim(Z)[2]+1):(dim(Z)[2]+dim(X)[2])])^2,df=1),digits=5))


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
                 AlphaCov=aCov,Cov=Inv_inform,RobCov=robust,RobAlphaCov=RobaCov)

  if(print){
    name=data.frame(intercept,pred)
    varname=colnames(name)
    MZIP_Mean=data.frame(Variable=varname,Alpha_Estimate=round(output$Aest,digits=3),SE=round(output$AModelSE,digits=3),
                         P_Value=round(output$AModelZpval,digits=4),Robust_SE=round(output$ARobustSE,digits=3),
                         Robust_P_Value=round(output$ARobustZpval,digits=4))
    print("Overall Mean Estimates from MZIP")
    print(MZIP_Mean)

    MZIP_ExcessZero=data.frame(Variable=varname,Gamma_Estimate=round(output$Gest,digits=3),SE=round(output$GModelSE,digits=3),
                               P_Value=round(output$GModelZpval,digits=4),Robust_SE=round(output$GRobustSE,digits=3),
                               Robust_P_Value=round(output$GRobustZpval,digits=4))
    print("Excess Zero Estimates from MZIP")
    print(MZIP_ExcessZero)

    rel_risk=exp(output$Aest)
    rrlowci=exp(output$AModelLower)
    rruppci=exp(output$AModelUpper)
    rrroblowci=exp(output$ARobustLower)
    rrrobuppci=exp(output$ARobustUpper)
    Relative_Risk=data.frame(Variable=varname,rel_risk=round(rel_risk,digits=3),
                             Lower_CI=round(rrlowci,digits=3),Upper_CI=round(rruppci,digits=3),
                             RobLower_CI=round(rrroblowci,digits=3),RobUpper_CI=round(rrrobuppci,digits=3))
    print("Overall Mean Relative Risk Estimates")
    print(Relative_Risk)
  }
  return(output)
}




#' Mediation Analysis for Zero-Inflated Count Mediators using MZIP (Continuous Outcome)
#'
#' This function incorporates the MZIP model into the counterfactual approach to mediation analysis
#' as proposed by Vanderweele when the mediator is a Zero-Inflated count variable. Errors for
#' direct and indirect effects are computed using delta method or bootstrap. Note: This function
#' assumes that the outcome is continuous and all exposure, mediator, outcome, and covariates
#' have the same sample size. Binary variables must be dummy coded prior.
#' @param outcome is the continuous outcome variable
#' @param mediator is the zero-inflated mediator variable, currently only 1 mediator allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetition if bootstrapped errors are used
#' @param C is a vector for theoretical values of each confounder. By default each each value of C will be the mean value of each confounder.
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
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
#'     #Example with delta method
#'     zimed=lmoutzimed(outcome=mzipmed_data$lmY,mediator=mzipmed_data$ziM,
#'                  exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                  mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0)
#'
#'     #Example using bootstrapping, 10 iterations used for succinctness
#'     zimed2=lmoutzimed(outcome=mzipmed_data$lmY,mediator=mzipmed_data$ziM,
#'                   exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                    mzipmed_data$C2),error="Boot",n=10,C=c(0,0.5))
#' @export



lmoutzimed=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,error='Delta',robust=FALSE){
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
   if (robust){
     MZIPCov=medreg$RobAlphaCov
   } else {
     MZIPCov=medreg$AlphaCov
   }

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

  outprint=round(matrix(c(output$NDE,output$NDEse,output$NDEci[1],output$NDEci[2],
                          output$NIE,output$NIEse,output$NIEci[1],output$NIEci[2],
                          output$TE,output$TEse,output$TEci[1],output$TEci[2],
                          output$PM,NA,NA,NA),nrow=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect","Natural Indirect Effect","Total Effect","Proportion Mediated")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
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
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
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
#'     Int is the overall additive interaction effect \cr
#'     Intse is the standard error for the additive interaction \cr
#'     IntCI is the confidence interval for the interaction effect \cr
#'     PAINT is the proportion attributable to the interaction effect \cr
#'     PE is the proportion eliminated
#' @examples
#'    #Example with exposure-mediator interaction
#'    #This builds upon function without interaction
#'     zimmed=lmoutzimedint(outcome=mzipmed_data$lmY,mediator=mzipmed_data$ziM,
#'                   exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                   mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0,M=NULL,C=NULL)
#' @export


lmoutzimedint=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,M=NULL,error='Delta',robust=FALSE){

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
    if (robust){
      MZIPCov=medreg$RobAlphaCov
    } else {
      MZIPCov=medreg$AlphaCov
    }
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
              Int=OvInt,Intse=OvIntse,IntCI=OvIntCI,PM=PM,PAINT=PI, PE=PE)

  outprint=round(matrix(c(output$NDE,output$NDEse,output$NDECI[1],output$NDECI[2],
                          output$NIE,output$NIEse,output$NIECI[1],output$NIECI[2],
                          output$TE,output$TEse,output$TECI[1],output$TECI[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$CDE,output$CDEse,output$CDECI[1],output$CDECI[2],
                          output$PIE,output$PIEse,output$PIECI[1],output$PIECI[2],
                          output$Intref,output$Intrefse,output$IntrefCI[1],output$IntrefCI[2],
                          output$Intmed,output$Intmedse,output$IntmedCI[1],output$IntmedCI[2],
                          output$Int,output$Intse,output$IntCI[1],output$IntCI[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$PM,NA,NA,NA,
                          output$PE,NA,NA,NA,
                          output$PAINT,NA,NA,NA),ncol=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect","Natural Indirect Effect","Total Effect",
                        "","4-way Decomposition","Controlled Direct Effect","Pure Indirect Effect",
                        "Interactive Reference Effect","Interactive Mediation Effect","Total Additive Interaction",
                        "","Proportions","Proportion Mediated","Proportion Eliminated","Proportion from Interaction")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
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
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
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
#'     #Example with delta method
#'     zimed=binoutzimed(outcome=mzipmed_data$binY,mediator=mzipmed_data$ziM,
#'                      exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                      mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0)
#'
#'     #Example using bootstrapping, 10 iterations are used for succinctness
#'     zimed2=binoutzimed(outcome=mzipmed_data$binY,mediator=mzipmed_data$ziM,
#'                    exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                    mzipmed_data$C2),error="Boot",n=10,C=c(0,0.5))
#' @export


binoutzimed=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,error='Delta',robust=FALSE){
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
  outreg=robust::glmRob(f,data=glmdata,family=stats::poisson())

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
                                                          exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]])-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                   exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]])-1))))
  } else{
    RRIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                          exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                   exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))))
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
      V=exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]])-1))
      W=exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]])-1))
      H=P+V
      J=S+W

      r1=(S/(1+S))-(P/(1+P))+((P+(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H)-
        ((S+(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J)
      r2=Xstar*((S/(1+S))-((S+(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J))+
        X*(-(P/(1+P))+((P+(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H))
      a1=(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)-
        (((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      a2=X*(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)-
        Xstar*(((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      t1=0
      t2=0
      t3=(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]))/H)
      -(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]))/J)


      GamIE=c(r1,r2,a1,a2,t1,t2,t3)

    } else {
      GamDE=c(0,0,0*C,0,0,0*C,0,X-Xstar,0,0*C)

      P=exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))
      S=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))
      V=exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))
      W=exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]])-1))
      H=P+V
      J=S+W

      r1=(S/(1+S))-(P/(1+P))+((P+(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H)-
        ((S+(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J)
      r2=Xstar*((S/(1+S))-((S+(exp(outreg$coefficients[[3]])-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J))+
        X*(-(P/(1+P))+((P+(exp(outreg$coefficients[[3]])-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H))
      r3=C*r1
      a1=(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)-
        (((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      a2=X*(((exp(outreg$coefficients[[3]])-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)-
        Xstar*(((exp(outreg$coefficients[[3]])-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      a3=C*a1
      t1=0
      t2=0
      t3=(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]))/H)+
        -(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]))/J)
      t4=C*0


      GamIE=c(r1,r2,r3,a1,a2,a3,t1,t2,t3,t4)
    }

    GamTE=GamIE+GamDE
    #Extract Covariance matrices
    if (robust){
      MZIPCov=medreg$RobCov
    } else {
      MZIPCov=medreg$Cov
    }
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
      outregb[[i]]=robust::glmRob(f,data=datab2[[i]],family=stats::poisson())
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

  outprint=round(matrix(c(output$RRNDE,output$logRRNDEse,output$RRNDEci[1],output$RRNDEci[2],
                          output$RRNIE,output$logRRNIEse,output$RRNIEci[1],output$RRNIEci[2],
                          output$RRTE,output$logRRTEse,output$RRTEci[1],output$RRTEci[2],
                          output$PM,NA,NA,NA),nrow=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","log SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect (RR)","Natural Indirect Effect (RR)","Total Effect (RR)","Proportion Mediated")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
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
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
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
#'     Int is the overall additive interaction effect \cr
#'     Intse is the standard error for the additive interaction \cr
#'     IntCI is the confidence interval for the interaction effect \cr
#'     PAINT is the proportion attributable to the interaction effect \cr
#'     PE is the proportion eliminated \cr
#'     PACDE is the proportion of the total effect due to neither mediation nor interaction \cr
#'     PAIntref is the proportion of the total effect due to just interaction \cr
#'     PAIntmed is the proportion of the total effect attributable to the joint effect of mediation and interaction \cr
#'     PAPIE is the proportion of the total effect attributable to just mediation \cr
#'     terr is the total excess relative risk
#' @examples
#'    #Example with exposure-mediator interaction
#'    #This builds upon function without interaction
#'     zimmed=binoutzimedint(outcome=mzipmed_data$binY,mediator=mzipmed_data$ziM,
#'                    exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                    mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0,M=NULL,C=NULL)

#' @export

binoutzimedint=function(outcome,mediator,exposure,confounder=NULL,C=NULL,n=1000,X=1,Xstar=0,M=NULL,error='Delta',robust=FALSE){
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
  outreg=robust::glmRob(f,data=glmdata,family=stats::poisson())

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
                                                exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((exp(outreg$coefficients[[2]]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  } else{
    RRNDE=((exp(outreg$coefficients[[2]]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((exp(outreg$coefficients[[2]]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  }

  #Indirect Effect
  if (is.null(confounder)){
    RRIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*X)+
                                                          exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                   exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))
  } else{
    RRIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                          exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                   exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))))
  }



  #Pure Indirect Effect
  if (is.null(confounder)){
    RRPIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))*(exp(medreg$Gest[1]+medreg$Gest[2]*X)+
                                                           exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
                                                   exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  } else{
    RRPIE=((1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                           exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))/
      ((1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))*(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
                                                                                   exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))))
  }


  #Proportion Mediated
  PM=RRNDE*(RRIE-1)/(RRIE*RRNDE-1)

  #kappa
  kappa=(exp(outreg$coefficients[[3]]*M+outreg$coefficients[[4]]*Xstar*M)*(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/
    ((exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))+exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)))


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
           exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)))
      R=(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)+
           exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)))
      B=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)
      D=Q-B
      G=R-B

      dr1=((B+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*B*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/Q)-((B+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*B*
                                                               exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/R)
      dr2=dr1*Xstar
      da1=(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*(B+1)*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/Q)-
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*(B+1)*
            exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/R)
      da2=da1*Xstar
      dt1=0
      dt2=X-Xstar
      dt3=(((B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)-
        (((B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      dt4=((X*(B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)-
        ((Xstar*(B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      GamNDE=c(dr1,dr2,da1,da2,dt1,dt2,dt3,dt4)


      P=exp(medreg$Gest[1]+medreg$Gest[2]*X)
      S=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar)
      V=exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      W=exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      H=P+V
      J=S+W

      nr1=(S/(1+S))-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H)-
        ((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J)
      nr2=Xstar*((S/(1+S))-((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J))+
        X*(-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H))
      na1=(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)-
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      na2=X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H)-
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J)
      nt1=0
      nt2=0
      nt3=(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/H)
      -(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/J)
      nt4=X*nt3

      GamIE=c(nr1,nr2,na1,na2,nt1,nt2,nt3,nt4)

      V2=exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      W2=exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      H2=P+V2
      J2=S+W2

      pr1=(S/(1+S))-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H2)-
        ((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J2)
      pr2=Xstar*((S/(1+S))-((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+medreg$Gest[1]+medreg$Gest[2]*Xstar)))/J2))+
        X*(-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+medreg$Gest[1]+medreg$Gest[2]*X)))/H2))
      pa1=(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H2)-
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J2)
      pa2=X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))/H2)-
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))/J2)
      pt1=0
      pt2=0
      pt3=(((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/H2)+
        -(((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/J2)
      pt4=Xstar*pt3

      GamPIE=c(pr1,pr2,pa1,pa2,pt1,pt2,pt3,pt4)

      Z=exp((X-Xstar)*outreg$coefficients[[2]]-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))

      mr1=((-Z*H*P)/((1+P)^2))+(Z*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((Z*J*S)/((1+S)^2))-(Z*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((-exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+S)
      mr2=((-X*Z*H*P)/((1+P)^2))+(Z*X*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((Xstar*Z*J*S)/((1+S)^2))-(Z*Xstar*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)))/(1+S)+
        ((X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+P)+
        ((-Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+S)
      ma1=(Z/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))-
        (Z/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      ma2= (X*Z/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X))-
        (Xstar*Z/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))-
        (X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      mt1=0
      mt2=((X-Xstar)*Z*H/(1+P))-((X-Xstar)*Z*J/(1+S))
      mt3=-M*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+P))*((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      mt4=-M*Xstar*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+P))*((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z*X/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))

      GamIntmed=c(mr1,mr2,ma1,ma2,mt1,mt2,mt3,mt4)

      rr1=(-(Z*J*S)/((1+S)^2))+(Z*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar)))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X)))/(1+S)
      rr2=Xstar*rr1
      ra1=(Z/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar))
      ra2=Xstar*ra1
      rt1=0
      rt2=(Xstar*Z*J/(1+S))-(X-Xstar)*RRCDE
      rt3=-M*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      rt4=-M*Xstar*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      -M*(X-Xstar)*RRCDE

      GamIntref=c(rr1,rr2,ra1,ra2,rt1,rt2,rt3,rt4)

    } else {

      GamCDE=c(0,0,0*C,0,0,0*C,0,1,0,M,0*C)

      Q=(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
           exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)))
      R=(exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))+
           exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)))
      B=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))
      D=Q-B
      G=R-B

      dr1=((B+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*B*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/Q)-((B+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*B*
                                                                                               exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/R)
      dr2=dr1*Xstar
      dr3=dr1*C
      da1=(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*D*(B+1)*
              exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/Q)-
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*G*(B+1)*
            exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/R)
      da2=da1*Xstar
      da3=da1*C
      dt1=0
      dt2=X-Xstar
      dt3=(((B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)-
        (((B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      dt4=((X*(B+1)*D*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/Q)-
        ((Xstar*(B+1)*G*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/R)
      dt5=C*0
      GamNDE=c(dr1,dr2,dr3,da1,da2,da3,dt1,dt2,dt3,dt4,dt5)

      P=exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))
      S=exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))
      V=exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      W=exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1))
      H=P+V
      J=S+W

      nr1=(S/(1+S))-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H)-
        ((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J)
      nr2=Xstar*((S/(1+S))-((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J))+
        X*(-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H))
      nr3=C*nr1
      na1=(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)-
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      na2=X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H)-
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J)
      na3=C*na1
      nt1=0
      nt2=0
      nt3=(((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/H)
      -(((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))/J)
      nt4=X*nt3
      nt5=C*0


      GamIE=c(nr1,nr2,nr3,na1,na2,na3,nt1,nt2,nt3,nt4,nt5)

      V2=exp((exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      W2=exp((exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+log(1+exp(medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C)))))*(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1))
      H2=P+V2
      J2=S+W2

      pr1=(S/(1+S))-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H2)-
        ((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J2)
      pr2=Xstar*((S/(1+S))-((S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*(exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*Xstar+sum(medreg[["Gest"]][c(3:m)]*C))))/J2))+
        X*(-(P/(1+P))+((P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*(exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+medreg$Gest[1]+medreg$Gest[2]*X+sum(medreg[["Gest"]][c(3:m)]*C))))/H2))
      pr3=C*pr1
      pa1=(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H2)-
        (((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J2)
      pa2=X*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))/H2)-
        Xstar*(((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))/J2)
      pa3=C*pa1
      pt1=0
      pt2=0
      pt3=(((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/H2)
      -(((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/J2)
      pt4=Xstar*pt3
      pt5=C*0

      GamPIE=c(pr1,pr2,pr3,pa1,pa2,pa3,pt1,pt2,pt3,pt4,pt5)

      Z=exp((X-Xstar)*outreg$coefficients[[2]]-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))

      mr1=((-Z*H*P)/((1+P)^2))+(Z*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)+
        ((Z*J*S)/((1+S)^2))-(Z*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)-
        ((-exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)
      mr2=((-X*Z*H*P)/((1+P)^2))+(Z*X*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*V*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)+
        ((Xstar*Z*J*S)/((1+S)^2))-(Z*Xstar*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)+
        ((X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2*P)/((1+P)^2))-(X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*V2*P*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+P)-
        ((-Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))+(Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)
      mr3=mr1*C
      ma1=(Z/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (Z/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      ma2= (X*Z/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (Xstar*Z/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (X*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      ma3=ma1*C
      mt1=0
      mt2=((X-Xstar)*Z*H/(1+P))-((X-Xstar)*Z*J/(1+S))
      mt3=-M*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+P))*((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      mt4=-M*Xstar*((Z*H)/(1+P)-(Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*H2)/(1+P)+(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+P))*((P+1)*V*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Z*X/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+P))*((P+1)*V2*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))+
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      mt5=0*C

      GamIntmed=c(mr1,mr2,mr3,ma1,ma2,ma3,mt1,mt2,mt3,mt4,mt5)

      rr1=(-(Z*J*S)/((1+S)^2))+(Z*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X)-1)*W*S*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)+
        ((exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2*S)/((1+S)^2))-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+(exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar)-1)*W2*S*exp(medreg$Aest[1]+medreg$Aest[2]*X+sum(medreg[["Aest"]][c(3:m)]*C))))/(1+S)
      rr2=Xstar*rr1
      rr3=C*rr1
      ra1=(Z/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))*(S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((exp(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*(S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)))
      ra2=Xstar*ra1
      ra3=C*ra1
      rt1=0
      rt2=(Xstar*Z*J/(1+S))-(X-Xstar)*RRCDE
      rt3=-M*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      rt4=-M*Xstar*((Z*J)/(1+S)-(exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))*J2)/(1+S))+
        (Z*X/(1+S))*((S+1)*W*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*X))-
        (Xstar*exp(-M*(outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))/(1+S))*((S+1)*W2*exp(medreg$Aest[1]+medreg$Aest[2]*Xstar+sum(medreg[["Aest"]][c(3:m)]*C)+outreg$coefficients[[3]]+outreg$coefficients[[4]]*Xstar))
      -M*(X-Xstar)*RRCDE
      rt5=C*0

      GamIntref=c(rr1,rr2,rr3,ra1,ra2,ra3,rt1,rt2,rt3,rt4,rt5)

    }
    GamTE=GamNDE+GamIE
    GamInt=GamIntref+GamIntmed

    #Covariance Matrices
    if (robust){
      MZIPCov=medreg$RobCov
    } else {
      MZIPCov=medreg$Cov
    }
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
      outregb[[i]]=robust::glmRob(f,data=datab2[[i]],family=stats::poisson())
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
              RRTE=RRTE,logRRTEse=logRRTEse,RRTEci=RRTECI,Int=RRInt,
              Intse=RRIntse,Intci=RRIntCI,PE=PE,PACDE=PACDE,PAIntref=PAIntref,PAIntmed=PAIntmed,PAPIE=PAPIE,PAINT=PAINT,terr=terr)

  outprint=round(matrix(c(output$RRNDE,output$logRRNDEse,output$RRNDEci[1],output$RRNDEci[2],
                          output$RRNIE,output$logRRNIEse,output$RRNIEci[1],output$RRNIEci[2],
                          output$RRTE,output$logRRTEse,output$RRTEci[1],output$RRTEci[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$RRCDE,output$logRRCDEse,output$RRCDEci[1],output$RRCDEci[2],
                          output$RRPIE,output$logRRPIEse,output$RRPIEci[1],output$RRPIEci[2],
                          output$Intref,output$Intrefse,output$Intrefci[1],output$Intrefci[2],
                          output$Intmed,output$Intmedse,output$Intmedci[1],output$Intmedci[2],
                          output$Int,output$Intse,output$Intci[1],output$Intci[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$PM,NA,NA,NA,
                          output$PE,NA,NA,NA,
                          output$PAINT,NA,NA,NA),ncol=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","log SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect (RR)","Natural Indirect Effect (RR)","Total Effect (RR)",
                        "","4-way Decomposition","Controlled Direct Effect (RR)","Pure Indirect Effect (RR)",
                        "Interactive Reference Effect","Interactive Mediation Effect","Total Additive Interaction",
                        "","Proportions","Proportion Mediated","Proportion Eliminated","Proportion from Interaction")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
  }







#' Mediation Analysis for Zero-Inflated Count Outcomes using MZIP
#'
#' This function incorporates the MZIP model into the counterfactual approach to mediation analysis
#' as proposed by Vanderweele when the outcome is a Zero-Inflated count variable for cases with
#' continuous mediators. Standard Errors for
#' direct and indirect effects are computed using delta method or bootstrapping. Note: This function
#' assumes that the outcome is continuous and all exposure, mediator, outcome, and confounder variables
#' have the same sample size. Binary variables must be dummy coded prior.
#' @param outcome is the zero-inflated count outcome variable
#' @param mediator is the continuous mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetition if bootstrapped errors are used. Default is 1000
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
#' @return The function will return a list of 12 elements.
#'     LM is the linear model regressing the exposure and covariates on the continuous mediator \cr
#'     MZIP is the results of regressing the exposure, covariates, and mediator on the outcome using the MZIP model \cr
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
#'     #Example using delta method
#'     ziout=zioutlmmed(outcome=mzipmed_data$ziY1,mediator=mzipmed_data$lmM,
#'                  exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                  mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0)
#'
#'    #Example using boostrapping, 10 iterations used for succinctness
#'    ziout2=zioutlmmed(outcome=mzipmed_data$ziY1,mediator=mzipmed_data$lmM,
#'                  exposure=mzipmed_data$X, confounder=cbind(mzipmed_data$C1,
#'                  mzipmed_data$C2),error="Boot",n=10)
#' @export

zioutlmmed=function(outcome,mediator,exposure,confounder=NULL,X=1,Xstar=0,error='Delta',n=1000,robust=FALSE){
  lmout=data.frame(mediator)
  if (is.null(confounder)){
    lmpred=data.frame(exposure)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator),print=F)
    mzipdata=data.frame(outcome,exposure,mediator)
  } else{
    lmpred=data.frame(exposure,confounder)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator,confounder),print=F)
    mzipdata=data.frame(outcome,exposure,mediator,confounder)
  }
  lmdata=data.frame(lmout,lmpred)

  #as.formula part of stats package
  f=stats::as.formula(paste(colnames(lmout),paste(colnames(lmpred),collapse="+"),sep="~"))

  #lm part of stats package
  medreg=stats::lm(f,data=lmdata)
  m=ncol(mzipdata)


  #Risk Ratio Direct Effect
  RRDE=exp(outreg$Aest[2]*(X-Xstar))

  #Risk Ratio Indirect Effect
  RRIE=exp(outreg$Aest[3]*medreg$coefficients[[2]]*(X-Xstar))

  #Proportion Mediated
  PM=RRDE*(RRIE-1)/(RRIE*RRDE-1)

  RRTE=RRIE*RRDE

  if (error=='Delta'){
    #Delta Method SE
    #Gamma
    if (is.null(confounder)){
      GamDE=c(0,0,0,1,0)
      GamIE=c(0,outreg$Aest[3],0,0,medreg$coefficients[[2]])
    } else{
      confounder=cbind(confounder)
      CL=ncol(confounder)
      C=rep(0,CL)
      GamDE=c(0,0,C,0,1,0,C)
      GamIE=c(0,outreg$Aest[3],C,0,0,medreg$coefficients[[2]],C)
    }
    GamTE=GamDE+GamIE

    #Covariance Matrices
    lmCov=stats::vcov(medreg) #uses stats package
    if (robust){
      MZIPCov=outreg$RobAlphaCov
    } else {
      MZIPCov=outreg$AlphaCov
    }

    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nlm,nMZI)
    Botmleft=matrix(0,nMZI,nlm)

    Top=cbind(lmCov,Topright)
    Btm=cbind(Botmleft,MZIPCov)

    CovM=rbind(Top,Btm)


    logRRDEse=sqrt(GamDE %*% CovM %*% GamDE)*abs(X-Xstar)
    logRRIEse=sqrt(GamIE %*% CovM %*% GamIE)*abs(X-Xstar)
    logRRTEse=sqrt(GamTE %*% CovM %*% GamTE)*abs(X-Xstar)

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
    PMb=list()
    confb=list()
    RRTEb=list()
    logRRDEb=list()
    logRRIEb=list()
    logRRTEb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(mzipdata),replace=T)
      datab2[[i]]=mzipdata[datab[[i]],]
      if (is.null(confounder)){
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]]),print=F)
      } else{
        confb[[i]]=data.matrix(datab2[[i]][c(4:m)])
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]],confb[[i]]),print=F)
      }

      medregb[[i]]=stats::lm(f,data=datab2[[i]])
      RRDEb[[i]]=exp(outregb[[i]]$Aest[2]*(X-Xstar))
      RRIEb[[i]]=exp(outregb[[i]]$Aest[3]*medregb[[i]]$coefficients[[2]]*(X-Xstar))
      PMb[[i]]=RRDEb[[i]]*(RRIEb[[i]]-1)/(RRIEb[[i]]*RRDEb[[i]]-1)
      RRTEb[[i]]=RRDEb[[i]]*RRIEb[[i]]

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

  output=list(MZIP=outreg,LM=medreg,RRNDE=RRDE,RRNIE=RRIE,PM=PM,logRRNDEse=logRRDEse,RRNDEci=RRDECI,logRRNIEse=logRRIEse,RRNIEci=RRIECI,
              RRTE=RRTE,logRRTEse=logRRTEse,RRTEci=RRTECI)

  outprint=round(matrix(c(output$RRNDE,output$logRRNDEse,output$RRNDEci[1],output$RRNDEci[2],
                          output$RRNIE,output$logRRNIEse,output$RRNIEci[1],output$RRNIEci[2],
                          output$RRTE,output$logRRTEse,output$RRTEci[1],output$RRTEci[2],
                          output$PM,NA,NA,NA),nrow=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","log SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect (RR)","Natural Indirect Effect (RR)","Total Effect (RR)","Proportion Mediated")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
}




#' Mediation Analysis for Zero-Inflated Count Outcomes using MZIP with Exposure-Mediator Interactions
#'
#' This function will do the same thing as the zioutlmmed function, but includes an exposure-mediator interaction.
#' 4-way decomposition of total effect (Vanderweele) are included in the output.
#' @param outcome is the zero-inflated count outcome variable
#' @param mediator is the continuous mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetitions for bootstrapping. Default is 1000. Setting n when using delta method errors will have no effect on output.
#' @param C is a vector for theoretical values of each confounder. If left out the default will be set to the mean of each confounder giving marginal effects
#' @param M is a fixed value for the mediator, M. If M is not specified, M will be set to its mean value
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
#' @return The function will return a list of 34 elements.
#'     MZIP is the results of regressing the mediator+exposure+confounder on the outcome using MZIP. To assess interaction effect individually look in the glm statement at the 4th parameter estimate \cr
#'     LM is the results of regressing the exposure and confounders on the mediator using linear regression \cr
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
#'     Int is the overall additive interaction effect \cr
#'     Intse is the standard error for the additive interaction \cr
#'     IntCI is the confidence interval for the interaction effect \cr
#'     PAINT is the proportion attributable to the interaction effect \cr
#'     PE is the proportion eliminated \cr
#'     PACDE is the proportion of the total effect due to neither mediation nor interaction \cr
#'     PAIntref is the proportion of the total effect due to just interaction \cr
#'     PAIntmed is the proportion of the total effect attributable to the joint effect of mediation and interaction \cr
#'     PAPIE is the proportion of the total effect attributable to just mediation \cr
#'     terr is the total excess relative risk
#' @examples
#' zimout=zioutlmmedint(outcome=mzipmed_data$ziY1,mediator=mzipmed_data$lmM,
#'              exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'              mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0,M=NULL,C=NULL)
#' @export

zioutlmmedint=function(outcome,mediator,exposure,confounder=NULL,n=1000,M=NULL,X=1,Xstar=0,C=NULL,error='Delta',robust=FALSE){
  interaction=mediator*exposure
  lmout=data.frame(mediator)
  if (is.null(confounder)){
    lmpred=data.frame(exposure)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator,interaction),print=F)
    mzipdata=data.frame(outcome,exposure,mediator,interaction)
  } else{
    lmpred=data.frame(exposure,confounder)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator,interaction,confounder),print=F)
    mzipdata=data.frame(outcome,exposure,mediator,interaction,confounder)
  }
  lmdata=data.frame(lmout,lmpred)


  f<-stats::as.formula(paste(colnames(lmout),paste(colnames(lmpred),collapse="+"),sep="~"))

  medreg=stats::lm(f,data=lmdata)
  r=ncol(lmdata)

  if (!is.null(confounder)){
    if (is.null(C)){
      confounder=cbind(confounder)
      C=colMeans(confounder)
    }
  }
  if (is.null(M)){
    M=mean(mediator)
  }

  RSS=sum(medreg[["residuals"]]^2)
  sigma2=RSS/(nrow(lmdata)-r-1)

  #RRCDE
  RRCDE=exp((outreg$Aest[2]+outreg$Aest[4]*M)*(X-Xstar))
  logRRCDE=log(RRCDE)


  #Risk Ratio Natural Direct Effect
  RRNDE=exp((outreg$Aest[2]+outreg$Aest[4]*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2))
            *(X-Xstar)+0.5*(outreg$Aest[4]^2)*sigma2*((X^2)-(Xstar^2)))
  logRRNDE=log(RRNDE)

  #Risk Ratio Indirect Effect
  RRIE=exp((outreg$Aest[3]*medreg$coefficients[[2]]+outreg$Aest[4]*medreg$coefficients[[2]]*X)*(X-Xstar))
  logRRIE=log(RRIE)


  #Some definitions
  Q=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M+(outreg$Aest[3]+outreg$Aest[4]*X)*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))
        +0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*sigma2)
  R=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M+(outreg$Aest[3]+outreg$Aest[4]*Xstar)*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))
        +0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*sigma2)
  S=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M+(outreg$Aest[3]+outreg$Aest[4]*X)*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))
        +0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*sigma2)
  U=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M+(outreg$Aest[3]+outreg$Aest[4]*Xstar)*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))
        +0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*sigma2)

  #Interactive Reference Effect
  RRIntref= Q-R-exp((outreg$Aest[2]+outreg$Aest[4]*M)*(X-Xstar))+1

  #Pure Indirect Effect
  RRPIE=exp((outreg$Aest[3]*medreg$coefficients[[2]]+outreg$Aest[4]*medreg$coefficients[[2]]*Xstar)*(X-Xstar))
  logRRPIE=log(RRPIE)

  #Interactive Mediation Effect
  RRIntmed=S-U-Q+R

  #Total Effect
  logRRTE=logRRNDE+logRRIE
  RRTE=exp(logRRTE)

  #Interaction Effect
  RRInt=RRIntmed+RRIntref

  #Proportion Mediated
  PM=RRNDE*(RRIE-1)/(RRIE*RRNDE-1)

  #Proportion Eliminated
  PE=(RRTE-RRCDE)/(RRTE-1)

  #Other Proportions
  kappa=exp(outreg$Aest[3]*M+outreg$Aest[4]*Xstar*M-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))-
              0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*sigma2)

  terr=kappa*(RRCDE-1)+kappa*RRIntref+kappa*RRIntmed+(RRPIE-1)

  PACDE=(kappa*(RRCDE-1))/terr
  PAIntref=(kappa*(RRIntref))/terr
  PAIntmed=(kappa*(RRIntmed))/terr
  PAPIE=(RRPIE-1)/terr
  PAINT=PAIntref+PAIntmed

  if (error=='Delta'){
    #Delta Method Standard Errors
    if (is.null(confounder)){
      GamCDE=c(0,0,0,1,0,M,0)
      GamNDE=c(outreg$Aest[4],outreg$Aest[4]*Xstar,0,1,outreg$Aest[4]*sigma2,medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]*sigma2+outreg$Aest[3]*sigma2*(X+Xstar),
               0.5*(outreg$Aest[4]^2)*(X+Xstar))
      GamIE=c(0,outreg$Aest[3]+outreg$Aest[4]*X,0,0,medreg$coefficients[[2]],medreg$coefficients[[2]]*X,0)

      GamPIE=c(0,outreg$Aest[3]+outreg$Aest[4]*Xstar,0,0,medreg$coefficients[[2]],medreg$coefficients[[2]]*Xstar,0)

      r1=(outreg$Aest[3]+outreg$Aest[4]*X)*Q-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R
      r2=(outreg$Aest[3]+outreg$Aest[4]*X)*Xstar*Q-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*Xstar*R
      r4=0
      r5=(X-Xstar)*Q-(X-Xstar)*exp((outreg$Aest[2]+outreg$Aest[4]*M)*(X-Xstar))
      r6=(M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q-
        (-M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R
      r7=(M*Xstar+X*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q-
        (-M*Xstar+Xstar*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R-
        M*(X-Xstar)*exp((outreg$Aest[2]+outreg$Aest[4]*M)*(X-Xstar))
      r9=0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*Q-0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*R

      GamIntref=c(r1,r2,r4,r5,r6,r7,r9)

      w1=(outreg$Aest[3]+outreg$Aest[4]*X)*S-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*U-(outreg$Aest[3]+outreg$Aest[4]*X)*Q+(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R
      w2=(outreg$Aest[3]+outreg$Aest[4]*X)*S*X-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*U*X-(outreg$Aest[3]+outreg$Aest[4]*X)*Q*Xstar+(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R*Xstar
      w4=0
      w5=(X-Xstar)*(S-Q)
      w6=(M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*S-
        (-M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*U-
        (M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q+
        (-M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R
      w7=(M*Xstar+X*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*S-
        (-M*Xstar+Xstar*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*U-
        (M*Xstar+X*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q+
        (-M*Xstar+Xstar*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R
      w9=0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*S-0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*U-
        0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*Q+0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*R

      GamIntmed=c(w1,w2,w4,w5,w6,w7,w9)
    } else {

      GamCDE=c(0,0,0*C,0,1,0,M,0*C,0)
      GamNDE=c(outreg$Aest[4],outreg$Aest[4]*Xstar,outreg$Aest[4]*C,0,1,outreg$Aest[4]*sigma2,medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2+outreg$Aest[3]*sigma2*(X+Xstar),
               0*C,0.5*(outreg$Aest[4]^2)*(X+Xstar))
      GamIE=c(0,outreg$Aest[3]+outreg$Aest[4]*X,0*C,0,0,medreg$coefficients[[2]],medreg$coefficients[[2]]*X,0*C,0)

      GamPIE=c(0,outreg$Aest[3]+outreg$Aest[4]*Xstar,0*C,0,0,medreg$coefficients[[2]],medreg$coefficients[[2]]*Xstar,0*C,0)

      r1=(outreg$Aest[3]+outreg$Aest[4]*X)*Q-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R
      r2=(outreg$Aest[3]+outreg$Aest[4]*X)*Xstar*Q-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*Xstar*R
      r3=(outreg$Aest[3]+outreg$Aest[4]*X)*C*Q-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R*C
      r4=0
      r5=(X-Xstar)*Q-(X-Xstar)*exp((outreg$Aest[2]+outreg$Aest[4]*M)*(X-Xstar))
      r6=(M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q-
        (-M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R
      r7=(M*Xstar+X*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q-
        (-M*Xstar+Xstar*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R-
        M*(X-Xstar)*exp((outreg$Aest[2]+outreg$Aest[4]*M)*(X-Xstar))
      r8=0
      r9=0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*Q-0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*R

      GamIntref=c(r1,r2,r3,r4,r5,r6,r7,0*C,r9)

      w1=(outreg$Aest[3]+outreg$Aest[4]*X)*S-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*U-(outreg$Aest[3]+outreg$Aest[4]*X)*Q+(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R
      w2=(outreg$Aest[3]+outreg$Aest[4]*X)*S*X-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*U*X-(outreg$Aest[3]+outreg$Aest[4]*X)*Q*Xstar+(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R*Xstar
      w3=C*((outreg$Aest[3]+outreg$Aest[4]*X)*S-(outreg$Aest[3]+outreg$Aest[4]*Xstar)*U-(outreg$Aest[3]+outreg$Aest[4]*X)*Q+(outreg$Aest[3]+outreg$Aest[4]*Xstar)*R)
      w4=0
      w5=(X-Xstar)*(S-Q)
      w6=(M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*S-
        (-M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*U-
        (M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q+
        (-M+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R
      w7=(M*Xstar+X*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*S-
        (-M*Xstar+Xstar*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*U-
        (M*Xstar+X*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))+outreg$Aest[3]*sigma2+outreg$Aest[4]*X*sigma2)*Q+
        (-M*Xstar+Xstar*(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))+outreg$Aest[3]*sigma2+outreg$Aest[4]*Xstar*sigma2)*R
      w8=0
      w9=0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*S-0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*U-
        0.5*(outreg$Aest[3]+outreg$Aest[4]*X)^(2)*Q+0.5*(outreg$Aest[3]+outreg$Aest[4]*Xstar)^(2)*R

      GamIntmed=c(w1,w2,w3,w4,w5,w6,w7,0*C,w9)

    }

    GamTE=GamNDE+GamIE
    GamInt=GamIntref+GamIntmed
    #Covariance Matrices
    lmCov=stats::vcov(medreg) #uses stats package
    if (robust){
      MZIPCov=outreg$RobAlphaCov
    } else {
    MZIPCov=outreg$AlphaCov
    }

    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nlm,nMZI)
    Botmleft=matrix(0,nMZI,nlm)

    Top=cbind(lmCov,Topright)
    Btm=cbind(Botmleft,MZIPCov)
    z1=matrix(0,nlm+nMZI,1)

    Comb=rbind(Top,Btm)
    Comb2=cbind(Comb,z1)


    sigma2var=2*(sigma2^2)/((nrow(lmdata)-r-1))
    z2=matrix(0,1,nlm+nMZI)
    sigmarow=cbind(z2,sigma2var)

    CovM=rbind(Comb2,sigmarow)

    #Now calculate standard errors
    logRRCDEse=sqrt(GamCDE %*% CovM %*% GamCDE)*abs(X-Xstar)
    logRRNDEse=sqrt(GamNDE %*% CovM %*% GamNDE)*abs(X-Xstar)
    logRRIEse=sqrt(GamIE %*% CovM %*% GamIE)*abs(X-Xstar)
    RRIntrefse=sqrt(GamIntref %*% CovM %*% GamIntref)
    RRIntmedse=sqrt(GamIntmed %*% CovM %*% GamIntmed)
    logRRPIEse=sqrt(GamPIE %*% CovM %*% GamPIE)*abs(X-Xstar)
    logRRTEse=sqrt(GamTE %*% CovM %*% GamTE)*abs(X-Xstar)
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

    LRRPIEci=exp(logRRPIE-1.96*logRRPIEse)
    URRPIEci=exp(logRRPIE+1.96*logRRPIEse)
    RRPIECI=c(LRRPIEci,URRPIEci)

    LRRTEci=exp(logRRTE-1.96*logRRTEse)
    URRTEci=exp(logRRTE+1.96*logRRTEse)
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
    RRPIEb=list()
    PMb=list()
    confb=list()
    RSSb=list()
    sigma2b=list()
    RRIntrefb=list()
    RRIntmedb=list()
    PIb=list()
    RRTEb=list()
    Qb=list()
    Rb=list()
    Sb=list()
    Ub=list()
    RRIntb=list()
    logRRNDEb=list()
    logRRIEb=list()
    logRRTEb=list()
    logRRCDEb=list()
    logRRPIEb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(mzipdata),replace=T)
      datab2[[i]]=mzipdata[datab[[i]],]
      if (is.null(confounder)){
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]],datab2[[i]][["interaction"]]),print=F)
      } else{
        confb[[i]]=data.matrix(datab2[[i]][c(5:r)])
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]],datab2[[i]][["interaction"]],confb[[i]]),print=F)
      }

      medregb[[i]]=stats::lm(f,data=datab2[[i]])
      RSSb[[i]]=sum(medregb[[i]][["residuals"]]^2)
      sigma2b[[i]]=RSSb[[i]]/(nrow(lmdata)-r-1)
      RRCDEb[[i]]=exp((outregb[[i]]$Aest[2]+outregb[[i]]$Aest[4]*M)*(X-Xstar))
      RRNDEb[[i]]=exp((outregb[[i]]$Aest[2]+outregb[[i]]$Aest[4]*(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3]*sigma2b[[i]]))
                      *(X-Xstar)+0.5*(outregb[[i]]$Aest[4]^2)*sigma2b[[i]]*((X^2)-(Xstar^2)))
      RRIEb[[i]]=exp((outregb[[i]]$Aest[3]*medregb[[i]]$coefficients[[2]]+outregb[[i]]$Aest[4]*medregb[[i]]$coefficients[[2]]*X)*(X-Xstar))
      RRPIEb[[i]]=exp((outregb[[i]]$Aest[3]*medregb[[i]]$coefficients[[2]]+outregb[[i]]$Aest[4]*medregb[[i]]$coefficients[[2]]*Xstar)*(X-Xstar))

      PMb[[i]]=RRNDEb[[i]]*(RRIEb[[i]]-1)/(RRIEb[[i]]*RRNDEb[[i]]-1)
      RRTEb[[i]]=exp(log(RRIEb[[i]])+log(RRNDEb[[i]]))
      Qb[[i]]=exp(outregb[[i]]$Aest[2]*(X-Xstar)-outregb[[i]]$Aest[3]*M-outregb[[i]]$Aest[4]*Xstar*M+(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*X)*(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C))
                  +0.5*(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*X)^(2)*sigma2b[[i]])
      Rb[[i]]=exp(-outregb[[i]]$Aest[3]*M-outregb[[i]]$Aest[4]*Xstar*M+(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar)*(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C))
                  +0.5*(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar)^(2)*sigma2b[[i]])
      Sb[[i]]=exp(outregb[[i]]$Aest[2]*(X-Xstar)-outregb[[i]]$Aest[3]*M-outregb[[i]]$Aest[4]*Xstar*M+(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*X)*(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C))
                  +0.5*(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*X)^(2)*sigma2b[[i]])
      Ub[[i]]=exp(-outregb[[i]]$Aest[3]*M-outregb[[i]]$Aest[4]*Xstar*M+(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar)*(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C))
                  +0.5*(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar)^(2)*sigma2b[[i]])
      RRIntrefb[[i]]=Qb[[i]]-Rb[[i]]-exp((outregb[[i]]$Aest[2]+outregb[[i]]$Aest[4]*M)*(X-Xstar))+1
      RRIntmedb[[i]]=Sb[[i]]-Ub[[i]]-Qb[[i]]+Rb[[i]]

      RRIntb[[i]]=RRIntrefb[[i]]+RRIntmedb[[i]]

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

  output=list(MZIP=outreg,LM=medreg,RRCDE=RRCDE,RRNDE=RRNDE,RRNIE=RRIE,logRRCDEse=logRRCDEse,
              logRRNDEse=logRRNDEse,logRRNIEse=logRRIEse,RRCDEci=RRCDECI,RRNDEci=RRNDECI,RRNIEci=RRIECI, PM=PM,
              Intref=RRIntref,Intrefse=RRIntrefse,Intrefci=RRIntrefCI,Intmed=RRIntmed,Intmedse=RRIntmedse,
              Intmedci=RRIntmedCI,RRPIE=RRPIE,logRRPIEse=logRRPIEse,RRPIEci=RRPIECI,
              RRTE=RRTE,logRRTEse=logRRTEse,RRTEci=RRTECI,Int=RRInt,
              Intse=RRIntse,Intci=RRIntCI,PE=PE,PACDE=PACDE,PAIntref=PAIntref,PAIntmed=PAIntmed,PAPIE=PAPIE,PAINT=PAINT,terr=terr)

  outprint=round(matrix(c(output$RRNDE,output$logRRNDEse,output$RRNDEci[1],output$RRNDEci[2],
                          output$RRNIE,output$logRRNIEse,output$RRNIEci[1],output$RRNIEci[2],
                          output$RRTE,output$logRRTEse,output$RRTEci[1],output$RRTEci[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$RRCDE,output$logRRCDEse,output$RRCDEci[1],output$RRCDEci[2],
                          output$RRPIE,output$logRRPIEse,output$RRPIEci[1],output$RRPIEci[2],
                          output$Intref,output$Intrefse,output$Intrefci[1],output$Intrefci[2],
                          output$Intmed,output$Intmedse,output$Intmedci[1],output$Intmedci[2],
                          output$Int,output$Intse,output$Intci[1],output$Intci[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$PM,NA,NA,NA,
                          output$PE,NA,NA,NA,
                          output$PAINT,NA,NA,NA),ncol=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","log SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect (RR)","Natural Indirect Effect (RR)","Total Effect (RR)",
                        "","4-way Decomposition","Controlled Direct Effect (RR)","Pure Indirect Effect (RR)",
                        "Interactive Reference Effect","Interactive Mediation Effect","Total Additive Interaction",
                        "","Proportions","Proportion Mediated","Proportion Eliminated","Proportion from Interaction")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
  }






#' Mediation Analysis for Zero-Inflated Count Outcomes using MZIP with binary mediators
#'
#' This function incorporates the MZIP model into the counterfactual approach to mediation analysis
#' as proposed by Vanderweele when the outcome is a Zero-Inflated count variable for cases with
#' binary mediators using a logistic regression mediator model. Standard Errors for
#' direct and indirect effects are computed using delta method or bootstrapping. Note: This function
#' assumes that the outcome is continuous and all exposure, mediator, outcome, and confounder variables
#' have the same sample size. Binary variables must be dummy coded prior.
#' @param outcome is the zero-inflated count outcome variable
#' @param mediator is the binary mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param C is a vector for theoretical values of each confounder. If left out the default will be set to the mean of each confounder giving marginal effects
#' @param n is the number of repetition if bootstrapped errors are used. Default is 1000
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
#' @return The function will return a list of 12 elements.
#'     GLM is the logistic model regressing the exposure and covariates on the continuous mediator \cr
#'     MZIP is the results of regressing the exposure, covariates, and mediator on the outcome using the MZIP model \cr
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
#'     #Example using delta method
#'     ziout=zioutbinmed(outcome=mzipmed_data$ziY2,mediator=mzipmed_data$binM,
#'                    exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                    mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0)
#' \dontrun{
#'     #Example using bootstrapping with 10 iterations
#'     ziout2=zioutbinmed(outcome=mzipmed_data$ziY2,mediator=mzipmed_data$binM,
#'                    exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                    mzipmed_data$C2),error="Boot",n=10,C=c(0,0.5))
#'    }
#' @export

zioutbinmed=function(outcome,mediator,exposure,confounder=NULL,n=1000,X=1,Xstar=0,C=NULL,error='Delta',robust=FALSE){
  glmout=data.frame(mediator)
  if (is.null(confounder)){
    glmpred=data.frame(exposure)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator),print=F)
    mzipdata=data.frame(outcome,exposure,mediator)
  } else{
    glmpred=data.frame(exposure,confounder)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator,confounder),print=F)
    mzipdata=data.frame(outcome,exposure,mediator,confounder)
  }
  glmdata=data.frame(glmout,glmpred)
  r=ncol(glmdata)

  #as.formula part of stats package
  f=stats::as.formula(paste(colnames(glmout),paste(colnames(glmpred),collapse="+"),sep="~"))

  #lm part of stats package
  medreg=stats::glm(f,data=glmdata,family=stats::binomial)
  m=ncol(mzipdata)

  if (!is.null(confounder)){
    if (is.null(C)){
      confounder=cbind(confounder)
      C=colMeans(confounder)
    }
  }

  #Risk Ratio Direct Effect
  RRDE=exp(outreg$Aest[2]*(X-Xstar))

  #Risk Ratio Indirect Effect
  RRIE=((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))*
          (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3])))/
    ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))*
       (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3])))


  #Proportion Mediated
  PM=RRDE*(RRIE-1)/(RRIE*RRDE-1)

  RRTE=RRIE*RRDE

  if (error=='Delta'){
    #Delta Method SE
    #Gamma
    if (is.null(confounder)){

      D=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))
      R=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X++outreg$Aest[3]))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]))
      K=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))
      Fe=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]))

      GamDE=c(0,0,0,X-Xstar,0)
      GamIE=c(D+R-Fe-K,Xstar*(D-Fe)+X*(R-K),0,0,R-Fe)


    } else {

      D=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))
      R=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]))
      K=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))
      Fe=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]))

      GamDE=c(0,0,C*0,0,X-Xstar,0,C*0)
      GamIE=c(D+R-Fe-K,Xstar*(D-Fe)+X*(R-K),C*(D+R-Fe-K),0,0,R-Fe,0*C)
    }
    GamTE=GamDE+GamIE

    #Covariance Matrices
    lmCov=stats::vcov(medreg) #uses stats package
    if (robust){
      MZIPCov=outreg$RobAlphaCov
    } else {
      MZIPCov=outreg$AlphaCov
    }

    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nlm,nMZI)
    Botmleft=matrix(0,nMZI,nlm)

    Top=cbind(lmCov,Topright)
    Btm=cbind(Botmleft,MZIPCov)

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
    PMb=list()
    confb=list()
    RRTEb=list()
    logRRDEb=list()
    logRRIEb=list()
    logRRTEb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(mzipdata),replace=T)
      datab2[[i]]=mzipdata[datab[[i]],]
      if (is.null(confounder)){
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]]),print=F)
      } else{
        confb[[i]]=data.matrix(datab2[[i]][c(4:m)])
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]],confb[[i]]),print=F)
      }

      medregb[[i]]=stats::glm(f,data=datab2[[i]],family=stats::binomial)
      RRDEb[[i]]=exp(outregb[[i]]$Aest[2]*(X-Xstar))
      if (is.null(confounder)){
        RRIEb[[i]]=((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar))*
                      (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+outregb[[i]]$Aest[3])))/
          ((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X))*
             (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+outregb[[i]]$Aest[3])))
      } else{
        RRIEb[[i]]=((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))*
                      (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3])))/
          ((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))*
             (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3])))

      }

      PMb[[i]]=RRDEb[[i]]*(RRIEb[[i]]-1)/(RRIEb[[i]]*RRDEb[[i]]-1)
      RRTEb[[i]]=RRDEb[[i]]*RRIEb[[i]]

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

  output=list(MZIP=outreg,GLM=medreg,RRNDE=RRDE,RRNIE=RRIE,PM=PM,logRRNDEse=logRRDEse,RRNDEci=RRDECI,logRRNIEse=logRRIEse,RRNIEci=RRIECI,
              RRTE=RRTE,logRRTEse=logRRTEse,RRTEci=RRTECI)

  outprint=round(matrix(c(output$RRNDE,output$logRRNDEse,output$RRNDEci[1],output$RRNDEci[2],
                          output$RRNIE,output$logRRNIEse,output$RRNIEci[1],output$RRNIEci[2],
                          output$RRTE,output$logRRTEse,output$RRTEci[1],output$RRTEci[2],
                          output$PM,NA,NA,NA),nrow=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","log SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect (RR)","Natural Indirect Effect (RR)","Total Effect (RR)","Proportion Mediated")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
}






#' Mediation Analysis for Zero-Inflated Count Outcomes using MZIP with Exposure-Mediator Interactions (Binary Outcome)
#'
#' This function will do the same thing as the zioutbinmed function, but includes an exposure-mediator interaction.
#' 4-way decomposition of total effect (Vanderweele) are included in the output.
#' @param outcome is the zero-inflated count outcome variable
#' @param mediator is the binary mediator variable, currently only 1 mediator variable is allowed
#' @param exposure is the primary exposure being considered, only 1 is allowed
#' @param confounder is a vector of confounder variables. If no confounder variables are needed then confounder is set to NULL. If more than 1 confounder is being considered then use the cbind function, e.g. cbind(var1,var2)
#' @param X is the theoretical value for the exposure variable to be set at. The default is to 1
#' @param Xstar is the theoretical value for the exposure variable to be compared to X. The default is 0, so direct, indirect, and proportion mediated values will be for a 1 unit increase in the exposure variable.
#' @param n is the number of repetitions for bootstrapping. Default is 1000. Setting n when using delta method errors will have no effect on output.
#' @param C is a vector for theoretical values of each confounder. If left out the default will be set to the mean of each confounder giving marginal effects
#' @param M is a fixed value for the mediator, M. If M is not specified, M will be set to its mean value
#' @param error ='Delta' for delta method standard errors and ='Boot' for bootstrap. Default is delta method
#' @param robust indicates if a robust covariance matrix should be used for MZIP in delta method derivations. Default is FALSE.
#' @return The function will return a list of 34 elements.
#'     MZIP is the results of regressing the mediator+exposure+confounder on the outcome using MZIP. To assess interaction effect individually look in the glm statement at the 4th parameter estimate \cr
#'     GLM is the results of regressing the exposure and confounders on the mediator using logistic regression \cr
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
#'     Int is the overall additive interaction effect \cr
#'     Intse is the standard error for the additive interaction \cr
#'     IntCI is the confidence interval for the interaction effect \cr
#'     PAINT is the proportion attributable to the interaction effect \cr
#'     PE is the proportion eliminated \cr
#'     PACDE is the proportion of the total effect due to neither mediation nor interaction \cr
#'     PAIntref is the proportion of the total effect due to just interaction \cr
#'     PAIntmed is the proportion of the total effect attributable to the joint effect of mediation and interaction \cr
#'     PAPIE is the proportion of the total effect attributable to just mediation \cr
#'     terr is the total excess relative risk
#' @examples
#'     zimout=zioutbinmedint(outcome=mzipmed_data$ziY2,mediator=mzipmed_data$binM,
#'                    exposure=mzipmed_data$X,confounder=cbind(mzipmed_data$C1,
#'                    mzipmed_data$C2),error="Delta",robust=FALSE,X=1,Xstar=0,M=NULL,C=NULL)
#' @export

zioutbinmedint=function(outcome,mediator,exposure,confounder=NULL,n=1000,M=NULL,X=1,Xstar=0,C=NULL,error='Delta',robust=FALSE){
  #lm,quantile,as.formula in Stats
  #colSds in matrixStats
  interaction=mediator*exposure
  glmout=data.frame(mediator)
  if (is.null(confounder)){
    glmpred=data.frame(exposure)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator,interaction),print=F)
    mzipdata=data.frame(outcome,exposure,mediator,interaction)
  } else{
    glmpred=data.frame(exposure,confounder)
    outreg=mzip(y=outcome,pred=cbind(exposure,mediator,interaction,confounder),print=F)
    mzipdata=data.frame(outcome,exposure,mediator,interaction,confounder)
  }
  glmdata=data.frame(glmout,glmpred)


  f<-stats::as.formula(paste(colnames(glmout),paste(colnames(glmpred),collapse="+"),sep="~"))

  medreg=stats::glm(f,data=glmdata,family=stats::binomial)
  r=ncol(glmdata)

  if (!is.null(confounder)){
    if (is.null(C)){
      confounder=cbind(confounder)
      C=colMeans(confounder)
    }
  }
  if (is.null(M)){
    M=mean(mediator)
  }


  #RRCDE
  RRCDE=exp((outreg$Aest[2]+outreg$Aest[4]*M)*(X-Xstar))
  logRRCDE=log(RRCDE)


  #Risk Ratio Natural Direct Effect
  RRNDE=exp(outreg$Aest[2]*(X-Xstar))*(1+exp(outreg$Aest[3]+outreg$Aest[4]*X+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))/
    (1+exp(outreg$Aest[3]+outreg$Aest[4]*Xstar+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))
  logRRNDE=log(RRNDE)

  #Risk Ratio Indirect Effect
  RRIE=((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))*
          (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X)))/
    ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))*
       (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X)))
  logRRIE=log(RRIE)


  #Pure Indirect Effect
  RRPIE=((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar)))/
    ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))*
       (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar)))
  logRRPIE=log(RRPIE)

  #Total Effect
  logRRTE=logRRNDE+logRRIE
  RRTE=exp(logRRTE)

  #kappa
  kappa=exp(outreg$Aest[3]*M+outreg$Aest[4]*Xstar*M)/
    ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar))/
       (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))))


  #Interactive Reference Effect
  RRIntref=((RRNDE-1)/kappa)-RRCDE+1


  #Interactive Mediation Effect
  RRIntmed=(RRTE-RRNDE-RRPIE+1)/kappa

  #Interaction Effect
  RRInt=RRIntmed+RRIntref

  #Proportion Mediated
  PM=RRNDE*(RRIE-1)/(RRIE*RRNDE-1)

  #Proportion Eliminated
  PE=(RRTE-RRCDE)/(RRTE-1)

  #Other Proportions
  terr=kappa*(RRCDE-1)+kappa*RRIntref+kappa*RRIntmed+(RRPIE-1)

  PACDE=(kappa*(RRCDE-1))/terr
  PAIntref=(kappa*(RRIntref))/terr
  PAIntmed=(kappa*(RRIntmed))/terr
  PAPIE=(RRPIE-1)/terr
  PAINT=PAIntref+PAIntmed

  if (error=='Delta'){
    #Delta Method Standard Errors
    if (is.null(confounder)){
      GamCDE=c(0,0,0,1,0,M)
      Q=exp(outreg$Aest[3]+outreg$Aest[4]*X+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)/
        (1+exp(outreg$Aest[3]+outreg$Aest[4]*X+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))
      B=exp(outreg$Aest[3]+outreg$Aest[4]*Xstar+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)/
        (1+exp(outreg$Aest[3]+outreg$Aest[4]*Xstar+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))

      GamNDE=c(Q-B,Xstar*(Q-B),0,X-Xstar,Q-B,X*Q-Xstar*B)

      D=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))
      R=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*X))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*X))
      K=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))
      Fe=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*X))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*X))

      GamIE=c(D+R-K-Fe,Xstar*(D-Fe)+X*(R-K),0,0,R-Fe,X*(R-Fe))

      R2=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*Xstar))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*Xstar))
      Fe2=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*Xstar))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*Xstar))


      GamPIE=c(D+R2-K-Fe2,Xstar*(D-Fe2)+X*(R2-K),0,0,R2-Fe2,Xstar*(R2-Fe2))

      E=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*X)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))
      G=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X++outreg$Aest[3]+outreg$Aest[4]*X)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))
      H=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*Xstar)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))
      J=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*Xstar)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))

      E1=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*X))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*X))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))^2)
      G1=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*X))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*X))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))^2)
      H1=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X))^2)
      J1=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar))^2)


      E2=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*X))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)))
      G2=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*X))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)))
      H2=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+outreg$Aest[3]+outreg$Aest[4]*Xstar))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X)))
      J2=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+outreg$Aest[3]+outreg$Aest[4]*Xstar))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar)))

      GamIntmed=c(E1-G1-H1+J1,Xstar*(E1+J1)-X*(G1+H1),0,(X-Xstar)*(E-G),
                  -M*(E-G-H+J)+E2-G2-H2-J2,-Xstar*M*(E-G-H+J)+X*(E-G)+Xstar*(J-H))



      GamIntref=c(H1-J1,X*H1-Xstar*J1,0,-(X-Xstar)*RRCDE,
                  -M*(H-J)+H2-J2,-Xstar*M*(H-J)+Xstar*(H-J)-M*(X-Xstar)*RRCDE)
    } else {

      GamCDE=c(0,0,0*C,0,1,0,M,0*C)
      Q=exp(outreg$Aest[3]+outreg$Aest[4]*X+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))/
        (1+exp(outreg$Aest[3]+outreg$Aest[4]*X+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))
      B=exp(outreg$Aest[3]+outreg$Aest[4]*Xstar+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))/
        (1+exp(outreg$Aest[3]+outreg$Aest[4]*Xstar+medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))

      GamNDE=c(Q-B,Xstar*(Q-B),C*(Q-B),0,X-Xstar,Q-B,X*Q-Xstar*B,C*0)

      D=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))
      R=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*X))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*X))
      K=exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))
      Fe=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*X))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*X))

      GamIE=c(D+R-K-Fe,Xstar*(D-Fe)+X*(R-K),C*(D+R-K-Fe),0,0,R-Fe,X*(R-Fe),C*0)

      R2=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*Xstar))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*Xstar))
      Fe2=(exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*Xstar))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[3]*Xstar))

      GamPIE=c(D+R2-K-Fe2,Xstar*(D-Fe2)+X*(R2-K),C*(D+R2-K-Fe2),0,0,R2-Fe2,Xstar*(R2-Fe2),C*0)


      E=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))
      G=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))
      H=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))
      J=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar)))/
        (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))

      E1=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))^2)
      G1=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))^2)
      H1=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)))^2)
      J1=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        ((exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))-
           (1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar))*
           (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)))^2)


      E2=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))))
      G2=exp(outreg$Aest[2]*(X-Xstar)-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*X))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))))
      H2=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*X+sum(medreg[["coefficients"]][c(3:r)]*C))))
      J2=exp(-outreg$Aest[3]*M-outreg$Aest[4]*Xstar*M)*
        (exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C)+outreg$Aest[3]+outreg$Aest[4]*Xstar))/
        ((1+exp(medreg$coefficients[[1]]+medreg$coefficients[[2]]*Xstar+sum(medreg[["coefficients"]][c(3:r)]*C))))

      GamIntmed=c(E1-G1-H1+J1,Xstar*(E1+J1)-X*(G1+H1),C*(E1-G1-H1+J1),0,(X-Xstar)*(E-G),
                  -M*(E-G-H+J)+E2-G2-H2-J2,-Xstar*M*(E-G-H+J)+X*(E-G)+Xstar*(J-H),C*0)



      GamIntref=c(H1-J1,X*H1-Xstar*J1,C*(H1-J1),0,-(X-Xstar)*RRCDE,
                  -M*(H-J)+H2-J2,-Xstar*M*(H-J)+Xstar*(H-J)-M*(X-Xstar)*RRCDE,C*0)
    }

    GamTE=GamNDE+GamIE
    GamInt=GamIntref+GamIntmed
    #Covariance Matrices
    lmCov=stats::vcov(medreg) #uses stats package
    if (robust){
      MZIPCov=outreg$RobAlphaCov
    } else {
      MZIPCov=outreg$AlphaCov
    }

    nlm=nrow(lmCov)
    nMZI=nrow(MZIPCov)
    Topright=matrix(0,nlm,nMZI)
    Botmleft=matrix(0,nMZI,nlm)

    Top=cbind(lmCov,Topright)
    Btm=cbind(Botmleft,MZIPCov)




    CovM=rbind(Top,Btm)

    #Now calculate standard errors
    logRRCDEse=sqrt(GamCDE %*% CovM %*% GamCDE)
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

    LRRPIEci=exp(logRRPIE-1.96*logRRPIEse)
    URRPIEci=exp(logRRPIE+1.96*logRRPIEse)
    RRPIECI=c(LRRPIEci,URRPIEci)

    LRRTEci=exp(logRRTE-1.96*logRRTEse)
    URRTEci=exp(logRRTE+1.96*logRRTEse)
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
    PMb=list()
    confb=list()
    RRIntrefb=list()
    RRIntmedb=list()
    RRPIEb=list()
    kappab=list()
    RRTEb=list()
    logRRNDEb=list()
    logRRIEb=list()
    logRRTEb=list()
    logRRCDEb=list()
    logRRPIEb=list()
    RRIntb=list()
    for (i in 1:n){
      datab[[i]]=sample(1:nrow(mzipdata),replace=T)
      datab2[[i]]=mzipdata[datab[[i]],]
      if (is.null(confounder)){
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]],datab2[[i]][["interaction"]]),print=F)
      } else{
        confb[[i]]=data.matrix(datab2[[i]][c(5:r)])
        outregb[[i]]=mzip(y=datab2[[i]][["outcome"]],pred=cbind(datab2[[i]][["exposure"]],datab2[[i]][["mediator"]],datab2[[i]][["interaction"]],confb[[i]]),print=F)
      }

      medregb[[i]]=stats::glm(f,data=datab2[[i]],family=stats::binomial)

      RRCDEb[[i]]=exp((outregb[[i]]$Aest[2]+outregb[[i]]$Aest[4]*M)*(X-Xstar))
      RRNDEb[[i]]=exp(outregb[[i]]$Aest[2]*(X-Xstar))*(1+exp(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*X+medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))/
        (1+exp(outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar+medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))
      RRIEb[[i]]=((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))*
                    (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*X)))/
        ((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))*
           (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*X)))
      PMb[[i]]=RRNDEb[[i]]*(RRIEb[[i]]-1)/(RRIEb[[i]]*RRNDEb[[i]]-1)
      RRPIEb[[i]]=((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))*
                     (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar)))/
        ((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*X+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)))*
           (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar)))
      RRTEb[[i]]=RRIEb[[i]]*RRNDEb[[i]]

      kappab[[i]]=exp(outregb[[i]]$Aest[3]*M+outregb[[i]]$Aest[4]*Xstar*M)/
        ((1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C)+outregb[[i]]$Aest[3]+outregb[[i]]$Aest[4]*Xstar))/
           (1+exp(medregb[[i]]$coefficients[[1]]+medregb[[i]]$coefficients[[2]]*Xstar+sum(medregb[[i]][["coefficients"]][c(3:r)]*C))))

      RRIntrefb[[i]]=((RRNDEb[[i]]-1)/kappab[[i]])-RRCDEb[[i]]+1
      RRIntmedb[[i]]=(RRTEb[[i]]-RRNDEb[[i]]-RRPIEb[[i]]+1)/kappab[[i]]
      RRIntb[[i]]=RRIntrefb[[i]]+RRIntmedb[[i]]
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

  output=list(MZIP=outreg,LM=medreg,RRCDE=RRCDE,RRNDE=RRNDE,RRNIE=RRIE,logRRCDEse=logRRCDEse,
              logRRNDEse=logRRNDEse,logRRNIEse=logRRIEse,RRCDEci=RRCDECI,RRNDEci=RRNDECI,RRNIEci=RRIECI, PM=PM,
              Intref=RRIntref,Intrefse=RRIntrefse,Intrefci=RRIntrefCI,Intmed=RRIntmed,Intmedse=RRIntmedse,
              Intmedci=RRIntmedCI,RRPIE=RRPIE,logRRPIEse=logRRPIEse,RRPIEci=RRPIECI,
              RRTE=RRTE,logRRTEse=logRRTEse,RRTEci=RRTECI,Int=RRInt,
              Intse=RRIntse,Intci=RRIntCI,PE=PE,PACDE=PACDE,PAIntref=PAIntref,PAIntmed=PAIntmed,PAPIE=PAPIE,PAINT=PAINT,terr=terr)

  outprint=round(matrix(c(output$RRNDE,output$logRRNDEse,output$RRNDEci[1],output$RRNDEci[2],
                          output$RRNIE,output$logRRNIEse,output$RRNIEci[1],output$RRNIEci[2],
                          output$RRTE,output$logRRTEse,output$RRTEci[1],output$RRTEci[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$RRCDE,output$logRRCDEse,output$RRCDEci[1],output$RRCDEci[2],
                          output$RRPIE,output$logRRPIEse,output$RRPIEci[1],output$RRPIEci[2],
                          output$Intref,output$Intrefse,output$Intrefci[1],output$Intrefci[2],
                          output$Intmed,output$Intmedse,output$Intmedci[1],output$Intmedci[2],
                          output$Int,output$Intse,output$Intci[1],output$Intci[2],
                          NA,NA,NA,NA,
                          NA,NA,NA,NA,
                          output$PM,NA,NA,NA,
                          output$PE,NA,NA,NA,
                          output$PAINT,NA,NA,NA),ncol=4,byrow=TRUE),digits=3)
  colnames(outprint)<-c("Estimate","log SE","Lower CI","Upper CI")
  rownames(outprint)<-c("Natural Direct Effect (RR)","Natural Indirect Effect (RR)","Total Effect (RR)",
                        "","4-way Decomposition","Controlled Direct Effect (RR)","Pure Indirect Effect (RR)",
                        "Interactive Reference Effect","Interactive Mediation Effect","Total Additive Interaction",
                        "","Proportions","Proportion Mediated","Proportion Eliminated","Proportion from Interaction")
  outprint<-as.table(outprint)
  print(outprint)

  return(output)
}
