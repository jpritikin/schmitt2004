################################
# concatenatedata.R
#
# Assembles parameter estimates into a single vector
#
#

#variance components
a<-mxEval(Al,fit)
c<-mxEval(Cl,fit)#excluded
e<-mxEval(El,fit)
t<-mxEval(Tl,fit)#excluded

# x * x' variance components
A<-mxEval(A,fit)
C<-mxEval(C,fit)#excluded
E<-mxEval(E,fit)
T<-mxEval(T,fit)#excluded

#total varaince

V<-(A+C+E+T)
er<-mxEval(err,fit) #error parameter
b <- mxEval(betas,fit) #regression parameter estimates

LL<-mxEval(objective,fit) #log likelihood

#Proportional (co)variance components
a2<-A/V
c2<-C/V
t2<-T/V
e2<-E/V

#decompose matrices into vectors

aV<-vech(a)
cV<-vech(c)
eV<-vech(e)
tV<-vech(t)
vV<-vech(V)

AV<-vech(A)
CV<-vech(C)
EV<-vech(E)
TV<-vech(T)
VV<-vech(V)

a2V<-vech(a2)
c2V<-vech(c2)
e2V<-vech(e2)
t2V<-vech(t2)


LL<-as.numeric(LL)

#concatenate vectors
ROIout<-c(i,"voxR",aV,cV,eV,tV,vV,AV,CV,EV,TV,VV,a2V,c2V,e2V,t2V,b,er,LL)

#clean up
rm(a,c,e,t,A,C,E,T,V,a2,c2,t2,e2,aV,cV,eV,tV,vV,AV,CV,EV,TV,VV,a2V,c2V,e2V,t2V)



