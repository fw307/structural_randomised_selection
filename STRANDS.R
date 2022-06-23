library(glmnet)
srs=function(X,Y,times=300,sel=NA,rho=0.5,alpha_glmnet=1,nfolds)        

  {
L=list()
n=nrow(X);p=ncol(X) 
resample=function(x, ...) x[sample.int(length(x), ...)]       
if(rho>1 | rho<0) {stop("rho should be between 0 and 1!")}
if(!is.na(sel)){if(sel>=1 | sel<=0) {stop("sel should be between 0 and 1!")}}   
if(is.na(sel)) {sel=0.5}

Y=Y-mean(Y); cv0=cv.glmnet(x=X,y=Y,alpha=alpha_glmnet,nfolds=nfolds,maxit=1000)
lasso0=cv0$glmnet.fit
Beta=predict(lasso0,type="coef",s=cv0$lambda.min)
beta.estimate=Beta[-1]
intercept=Beta[1]


sds=apply(X,2,sd)
X=scale(X)


COR=abs(cor(X));  diag(COR)=0 
no_cor_group=0; group0=1:p; used.idx=c()
relevant.idx=order(abs(beta.estimate),decreasing=TRUE)[0:sum(beta.estimate!=0)]



     for(idx in relevant.idx)
  {
if(length(setdiff(which(COR[idx,]>=rho),used.idx))==0 | idx %in% used.idx) {next}
new.group=idx; no_cor_group=no_cor_group+1; flag=0
         while(flag==0) 
     {
       rest.idx=setdiff(group0,new.group)
       if(length(rest.idx)==0) {break}
       if(length(rest.idx)==1) {m.cor=mean(COR[rest.idx,new.group])}
       if(length(rest.idx)>1) {SUB.COR=as.matrix(COR[rest.idx,new.group]); m.cor=apply(SUB.COR,1,mean)}

       if(max(m.cor)<rho){flag=1}
       if(max(m.cor)>=rho) {new.idx=rest.idx[which.max(m.cor)]; new.group=union(new.group,new.idx)}
       assign(sprintf("group%s",no_cor_group),new.group); group0=setdiff(group0,new.group); used.idx=union(used.idx,new.group)
     }

  }


##################################Step 1##################################
beta_coef=rep(0,p)                         
beta_freq=rep(0,p)                         
beta_times=rep(0,p)           


lambdas=cv0$lambda.min


iter=1     


                                 while(iter<=times)
               {
         idx=c()     

         for(j in 0:no_cor_group)    
     {
    groupj=get(sprintf("group%s",j)); group.size=length(groupj)
    if(group.size>0){idx=c(idx,resample(groupj,size=resample(0:group.size,size=1),replace=FALSE))}
     } 
    if(length(idx)<=1) {next}


    X_iter=X[,idx]; Y_iter=Y
    cv_iter=cv.glmnet(x=X_iter,y=Y_iter,alpha=alpha_glmnet,nfolds=nfolds,maxit=1000000) 
    lasso_iter=cv_iter$glmnet.fit
    Beta=predict(lasso_iter,type="coef",s=cv_iter$lambda.min)
    beta_estimate=Beta[-1]
   
    beta_freq[idx]=beta_freq[idx]+abs(sign(beta_estimate))
    beta_coef[idx]=beta_coef[idx]+abs(beta_estimate)
    beta_times[idx]=beta_times[idx]+1


    lambdas=c(lambdas,cv_iter$lambda.min)
    iter=iter+1
               }     

if(prod(beta_times)==0) {stop("Need to run more times!")}
weights=abs(beta_freq/beta_times)
coef.weights=beta_coef/beta_times

if(sum(abs(weights))==0) {print("Null model!")}
L$"times"=beta_times
############################################################################


lambdas=unique(lambdas[lambdas<=max(cv0$lambda) & lambdas>=min(cv0$lambda)])
lambdas=sort(lambdas,decreasing=TRUE)


#######################################Step 2############################################
beta_times=rep(0,p);beta_freq=rep(0,p);beta_coef=rep(0,p)
iter=1
model.size=ceiling(sum(weights)) 
print(sum(weights))


                                            while(iter<=times)
               { 
   core=c()
   if(model.size>0) {core=resample(1:p,size=model.size,prob=abs(coef.weights*weights),replace=FALSE)}

   idx=core
   if(length(idx)==0) {beta_times=beta_times+1;iter=iter+1;next}
   if(length(idx)==1) 
{
X_iter=X[,idx]; Y_iter=Y 
linearmodel=lm(Y_iter~X_iter)
p.val=summary(linearmodel)$coefficients[2,4]
if(p.val>=0.05) {beta_times=beta_times+1;iter=iter+1;next}
    if(p.val<=0.05) 
  {
     beta_freq[idx]=beta_freq[idx]+1
     beta_coef[idx]=beta_coef[idx]+summary(linearmodel)$coefficients[2,1]
     beta_times=beta_times+1    
     iter=iter+1
  }
next
}


X_iter=X[,idx]; Y_iter=Y 
cv_iter=cv.glmnet(x=X_iter,y=Y_iter,alpha=alpha_glmnet,nfolds=nfolds,maxit=1000000,lambda=lambdas)  
lasso_iter=cv_iter$glmnet.fit
Beta=predict(lasso_iter,type="coef",s=cv_iter$lambda.min)



beta_estimate=Beta[-1] 
if(sum(is.na(beta_estimate))!=0) {next}

beta_freq[idx]=beta_freq[idx]+abs(sign(beta_estimate))
beta_coef[idx]=beta_coef[idx]+beta_estimate 
beta_times=beta_times+1    
iter=iter+1
               }          

####################################################################################################################


result=coef=ifelse(beta_times!=0,beta_coef/beta_times,0)
beta_freq=ifelse(beta_times!=0,beta_freq/beta_times,0)
L$"intercept"=intercept

select.idx1=which(beta_freq>=sel); s.hat=length(select.idx1)
select.idx2=order(abs(result),decreasing=TRUE)[0:s.hat]
select.idx=union(select.idx1,select.idx2)
if(length(select.idx)>0) {result[-select.idx]=0}
if(length(select.idx)==0) {result[1:p]=0}

L$"cv0"=cv0
L$"beta0"=beta.estimate
L$"lambdas"=lambdas
L$"coef"=coef/sds
L$"beta"=result/sds
L$"freq"=beta_freq
L$"weights"=weights 
L$"coef.weights"=coef.weights



L
  }
