
#R functions for app
library(lavaan)
# library(ggplot2)
# library(plotly)


#--------------------------------------------------------------------------------------#
#this function creates a lavaan model with k=2 factors and p/k indicators per factor 

model.2f<-function(p){
k=2
if(p%%k != 0) { stop("p is not a multiple of k")}  
v1=paste("x",1:(p/k-1),"+",sep="",collapse="")
v2=paste("x",(p/k+1):(p-1),"+",sep="",collapse="")
model=paste("f1=~",v1,"x",p/k," \n ", "f2=~",v2,"x",p,sep="",collapse="")
}
#---------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------#
#this function creates a lavaan model with k = 3 factors and p/k indicators per factor 

model.3f<-function(p){
   k=3
   if(p%%k != 0) { stop("p is not a multiple of k")} 
   v1=paste("x",1:(p/k-1),"+",sep="",collapse="")
   v2=paste("x",(p/k+1):(2*p/k-1),"+",sep="",collapse="")
   v3=paste("x",(2*p/k+1):(p-1),"+",sep="",collapse="")
   model=paste("f1=~",v1,"x",p/k," \n ", "f2=~",v2,"x",2*p/k," \n ","f3=~",v3,"x",p,sep="",collapse="")
}
#---------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------#
#this function fits a lavaan model to a p x p population covariance matrix and computes pop fit indices
#the matrix must have 1s on diagonal or else everything defaults to NA

get.fits<-function(Sigmastar,lavmodel){
      p<-nrow(Sigmastar)
      colnames(Sigmastar)<-rownames(Sigmastar)<-paste0("x",1:p) 
      
      fit<-try(cfa(lavmodel,sample.nobs=100,sample.cov=Sigmastar,check.post=FALSE,warn=FALSE),
               silent=TRUE) #n.obs is arbitrary
      
      if (class(fit)=="try-error" | inspect(fit,"converged")==FALSE){
        fits=rep(NA,4)
        } else {
          conv=inspect(fit,"converged")
          df=fitmeasures(fit)["df"]
          fmin<-fitmeasures(fit)["fmin"]*2 #lavaan saves 1/2 of Fml
          rmsea<-sqrt(fmin/df) 
      
          srmr<-fitmeasures(fit)["srmr"] #can be taken directly from lavaan
      
          #below follows the traditional computation of the CFI, although we could have just used 
          #determinant of Sigmastar, as shown in the appendix
          fit_null<- lavaan:::lav_object_independence(fit) #independence model
          fmin_null<-fitmeasures(fit_null)["fmin"]*2 #lavaan saves 1/2 of Fml
          df_null<-fitmeasures(fit_null)["df"]
          cfi<-1-fmin/fmin_null
          fits<-c(conv,rmsea,srmr, cfi)
      }
      names(fits)<-c("conv","rmsea","srmr","cfi")
      return(fits) 
}  
    
#---------------------------------------------------------------------------------------------------#  
#the main function below creates population covariance matrices for a 2- or a 3-factor model 
#with an increasing number of cross-loadings, then use helper function model.2f and model.3f to create
#lavaan model files and use helper function get.fits to fit the model without crossloadings to 
#all the population covaraince matrices and compute rmsea, cfi, and srmr
#k is number of factors, p is the number of variables, fcor is factor correlation, 
#l is the vector of main loadings for all the variables, cld is the vector of cross-loadings  

#rmsea.2f has been renamed into main.2f

main.2f <- function(p=30,fcor=.1,l,cld){ #length of l and cld depends on p

k=2
if(p%%k != 0) { print(p); stop("p is not a multiple of k")}  #is quit better?   
model<-model.2f(p)

#factor loadings matrix      
llist<-split(l, rep(1:k,each=p/k)) #will generate error if p is not a multiple of k
l1<-llist[[1]] 
l2<-llist[[2]]
L<-matrix(c(l1,rep(0,p),l2),nrow=p,ncol=k)
Lsame<-Ldif<-L #create two alternating orders
Phi<-matrix(fcor,nrow=k,ncol=k)  #factor cor matrix
diag(Phi)<-1
diagPsi <- 1-diag(L%*%Phi%*%t(L)) #residual variances 

tosave<-c("number_crossloadings","conv_same","rmsea_same_f","srmr_same_f","cfi_same_f",
           "conv_dif","rmsea_dif_f","srmr_dif_f","cfi_dif_f",paste0("x",1:p,"same"),paste0("x",1:p,"dif"))
results <- data.frame(matrix(vector(), 0, length(tosave), dimnames=list(c(), tosave)))
results[1,]<-c(0,1,0,0,1,1,0,0,1,diagPsi,diagPsi) #for zero cross-loadings

for (t in (1:p)) { 
    #add one cross-loading at a time, to the same factor
    if (t<(p/k+1)) {Lsame[t,2]=cld[t]} else {Lsame[t,1]=cld[t]}
    Sigmastar_same<-Lsame%*%Phi%*%t(Lsame) #true cov matrix
    diagPsi_same<-1-diag(Sigmastar_same)
    
   if (max(diag(Sigmastar_same))>.99) {fits_same<-rep(NA,4)} else {
      diag(Sigmastar_same)<-1; fits_same<-get.fits(Sigmastar_same,model)
      }  
      
   #add one cross-loading at a time, to alternating factors
   if (t %% 2 ==0) {j=1; a<-(.5*p+.5*t)} else {j=2; a<-(.5*(t-1)+1)}
   Ldif[a,j]<-cld[t]
   Sigmastar_dif<-Ldif%*%Phi%*%t(Ldif) #true cov matrix
   diagPsi_dif<-1-diag(Sigmastar_dif)
   
   if (max(diag(Sigmastar_dif))>.99){fits_dif<-rep(NA,4)} else {
      diag(Sigmastar_dif)<-1 #adding error variances in the remaining amounts
      fits_dif<-get.fits(Sigmastar_dif,model)  
       }
  
   results[t+1,]<-c(t,fits_same,fits_dif,diagPsi_same,diagPsi_dif) 
  
}  #end of t loop
#print(dim(results))

#print(dim(results))
results.full<-list(results) #now includes orders
names(results.full)<-c("results")
return(results.full) #now returns a list with two elements
}  


#k=3
main.3f <- function(p=30,fcor=.1,l,cld){ #length of l must be p, length of cld must be 2*p
   
   k=3
   if(p%%k != 0) { stop("p is not a multiple of k")}  #is quit better?   
   model<-model.3f(p)
   
   #factor loadings matrix      
   llist<-split(l, rep(1:k,each=p/k)) #will generate error if p is not a multiple of k
   l1<-llist[[1]] 
   l2<-llist[[2]]
   l3<-llist[[3]]
   L<-matrix(c(l1,rep(0,p),l2,rep(0,p),l3),nrow=p,ncol=k)
   Lsame1<-Lsame2<-Ldif1<-Ldif2<-L
   Phi<-matrix(fcor,nrow=k,ncol=k)  #factor cor matrix
   diag(Phi)<-1
   diagPsi <- 1-diag(L%*%Phi%*%t(L)) #residual variances 
   
   tosave<-c("number_crossloadings","conv_same1","rmsea_same1_f","srmr_same1_f","cfi_same1_f",
             "conv_same2","rmsea_same2_f","srmr_same2_f","cfi_same2_f",
              "conv_dif1","rmsea_dif1_f","srmr_dif1_f","cfi_dif1_f",
             "conv_dif2","rmsea_dif2_f","srmr_dif2_f","cfi_dif2_f",paste0("x",1:p,"same1"),
             paste0("x",1:p,"same2"),paste0("x",1:p,"dif1"),paste0("x",1:p,"dif2"))
  
   results <- data.frame(matrix(vector(), 0, length(tosave), dimnames=list(c(), tosave)))
   results[1,]<-c(0,rep(c(1,0,0,1),4),rep(diagPsi,4)) #for zero cross-loadings
   
   #here we explicitly calculate the column/row index for each of the four orders
   #Partition the matrix L into 9 nine cells where each has p/k rows and 1 column:
   # L= * 1 2
   #    3 * 4
   #    5 6 *
   # asterisks indicate main loadings. The numbered groups of loadings, 1-6, are cycled through
   # in different orders in the four configurations below
   
   #first "same" order is columnwise: all loadings in 3, then in 5, then in 1, then in 6, 
   #then in 2, then in 4
   row.same1<-c((p/k+1):p,1:(p/k),(2*p/k+1):p,1:(2*p/k)) 
   col.same1<-c(rep(1,2*p/k),rep(2,2*p/k),rep(3,2*p/k))
      
   #second "same" order is to columnwise until indicators of just 1 factor have been covered:
   #all loadings in 3, then in 6, then in 2, then in 5, then in 1, then in 4
   row.same2<-row.same1 #row order is the same
   col.same2<-c(rep(1,p/k),rep(2,p/k),rep(3,p/k),rep(1,p/k),rep(2,p/k),rep(3,p/k))
   
   #the first "different" order is filling out all the zeros going horizontally by row
   #that is add by row all loadings to 1, 2, then to 3, 4; then to 5, 6 
   row.dif1<-c(rep(1:(p/k),each=2),rep((p/k+1):(2*p/k),each=2),rep((2*p/k+1):p,each=2)) 
   col.dif1<-c(rep(c(2,3),p/k),rep(c(1,3),p/k),rep(c(1,2),p/k))
   
   #the second "different" order is as follows: 
   #one loading from each of: 3, 6, 2, 5, 1, 4; then repeat
  
   ni<-p/k #number of items per factor
   #te<-c(11,21,1,12,22,2) #the first six elements, specific to p=30
   te<-c(ni+1,ni*2+1,1,ni+2,ni*2+2,2) #the first six elements
   
   t2<-c(rep(0:(ni-1),each=6))
   ind<-te+t2
   ind[(2*p-2):(2*p)]<-ind[1:3] #awkward but ok
   
   row.dif2<-ind
   col.dif2<-rep(c(1,2,3,3,1,2),p/k) 
   
   #added 8/15/2022
   orders<-cbind(cld,row.same1,col.same1,row.same2,col.same2,row.dif1,col.dif1,row.dif2,col.dif2)
   
   for (t in (1:length(row.same1))) {  #cycling through all (p-1)*k crossloadings to add
      
      #print(t)
      #add one cross-loading at a time, in all four orders
      Lsame1[row.same1[t],col.same1[t]] <- cld[t]
      Lsame2[row.same2[t],col.same2[t]] <- cld[t]
      Ldif1[row.dif1[t],col.dif1[t]] <- cld[t]
      Ldif2[row.dif2[t],col.dif2[t]] <- cld[t]
       # print(Lsame1)
       # print(Lsame2)
       # print(Ldif1)
       # print(Ldif2)
      
      #true score cov matrices
      Sigmastar_same1<-Lsame1%*%Phi%*%t(Lsame1) 
      Sigmastar_same2<-Lsame2%*%Phi%*%t(Lsame2) 
      Sigmastar_dif1<-Ldif1%*%Phi%*%t(Ldif1) 
      Sigmastar_dif2<-Ldif2%*%Phi%*%t(Ldif2)
      
      #debug (delete later)
      if (t==15){save(Sigmastar_same1,Sigmastar_same2,file="debug.RData")}
      
      #residual variances (computing for all t)
      diagPsi_same1<-1-diag(Sigmastar_same1)
      diagPsi_same2<-1-diag(Sigmastar_same2)
      diagPsi_dif1<-1-diag(Sigmastar_dif1)
      diagPsi_dif2<-1-diag(Sigmastar_dif2)
      
      # rather than checking the eigenvalues, lets check if true score variances are too close to 1
      #is this always the same for these types of models? 
      if (max(diag(Sigmastar_same1))>.99) {fits_same1<-rep(NA,4)} else {
         diag(Sigmastar_same1)<-1; fits_same1<-get.fits(Sigmastar_same1,model)
      }  
      
      if (max(diag(Sigmastar_same2))>.99) {fits_same2<-rep(NA,4)} else {
         diag(Sigmastar_same2)<-1; fits_same2<-get.fits(Sigmastar_same2,model) 
      } 
      
      if (max(diag(Sigmastar_dif1))>.99) {fits_dif1<-rep(NA,4)} else {
         diag(Sigmastar_dif1)<-1; fits_dif1<-get.fits(Sigmastar_dif1,model)
      } 
      
      if (max(diag(Sigmastar_dif2))>.99) {fits_dif2<-rep(NA,4)} else {
         diag(Sigmastar_dif2)<-1; fits_dif2<-get.fits(Sigmastar_dif2,model)
      } 
      
      results[t+1,]<-c(t,fits_same1,fits_same2,fits_dif1,fits_dif2,diagPsi_same1,diagPsi_same2,
                       diagPsi_dif1,diagPsi_dif2) 
      
    }  #end of t loop
   
   #print(dim(results))
   results.full<-list(results,orders) #now includes orders
   names(results.full)<-c("results","orders")
   return(results.full) #now returns a list with two elements
}  

#constraint: l^2+cld^2+2*phi*ld*cld < 1 

#test of main.2f
#l=rep(.7,10)
#cld=rep(.3,10)
#out<-main.2f(p=10,l=l,cld=cld,fcor=0) 
#out

#test of main.3f
# p=15
# fcor=0
# l=rep(.7,p)
# cld<-seq(.3,.89,by=.01)
# cld=rep(.7,2*p)
# out1<-main.3f(fcor=fcor,p=p,l=l,cld=cld)
# out1$results  
# out1$orders

#test 
# p=8
# l=c(0.608, 0.915, 0.568, 0.605, 0.753, 0.902, 0.590, 0.458)
# cld=c(.529, .366, .538, .464, .393,.405,.539, .521)
# fcor=.2
# results<-main.2f(p=p,fcor=fcor,l=l,cld=cld)
# results

#p=20
#l=runif(p,min=.7-.3/2,max=.7+.3/2) 
#l=rep(.7,p)
#cld=rep(.3,p)
#rmseas<-rmsea.2f(p=p,fcor=.5,l,cld=rep(.3,p))

# p=20
# fcor=.2
# l=rep(.7,p)
# cld=rep(.4,p)
# results<-main.2f(p=p,fcor=fcor,l=l,cld=cld)
# results
# round(results$rmsea_dif_f,3)
# #max(results$rmsea_dif_f)
# 
# plot(results$number_crossloadings,results$rmsea_dif_f,type="l")
# lines(results$number_crossloadings,results$rmsea_same_f,col="red")
# abline(h=.05)
# 
# plot(results$number_crossloadings,results$srmr_dif_f,type="l")
# lines(results$number_crossloadings,results$srmr_same_f,col="red")
# abline(h=.05)
# 
# plot(results$number_crossloadings,results$cfi_dif_f,type="l")
# lines(results$number_crossloadings,results$cfi_same_f,col="red")
# abline(h=.95)
# 
# 
# #  
# #ggplot
# ggplot(data=results) + 
# geom_line(mapping = aes(x=number_crossloadings, y=rmsea_same_f,color="To 1st factor, then 2d"),size=1)+ 
# geom_line(mapping = aes(x=number_crossloadings, y=rmsea_dif_f,color="To alternating factors"),size=1)+ 
# geom_abline(color="grey",slope=0, intercept=0.05) + labs(color = "How crossloadings are added") +
# xlab("Number of crossloadings in the true model")+
# ylab("RMSEA for the model with no crossloadings")
# 
#  
