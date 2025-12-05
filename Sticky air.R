#Demi-angle de l'entonnoir
a=c(1e-2,2e-2,3e-2) #en m
inca=1e-2 #en m
h=c(sqrt(5)*1e-2,sqrt(13)*1e-2,sqrt(153)*1e-2) #en m
inch=1e-2 #en m
o=c(2e-2,3e-2,12e-2) #en m
inco=1e-2 #en m
alpha=NULL
for (i in 1:length(a)){alpha[i]=mean(c(acos(a[i]/h[i]),asin(o[i]/h[i])))}
r=c(1e-2,2e-2,3e-2) #en m
incr=1e-2 #en m
L=c(10e-2,20e-2,30e-2) #en m
incL=1e-2 #en m
incalpha=alpha*sqrt((inca/a)^2+(inch/h)^2+(inco/o)^2)

#Mesures du temps
Z01=c(0,1,2) #en m
Z02=c(0,1,2) #en m
Z03=c(0,1,2) #en m
incZ0=1e-2 #en m
tf1=c(0,1,2)
inctf1=0.01 #en s
tf2=c(3,4,5)
inctf2=inctf1
tf3=c(6,7,8)
inctf3=inctf1
data1=data.frame(tf1,r)
data2=data.frame(tf2,L)
data3=data.frame(tf3,alpha)

par(mfrow=c(2,2))
#Ajustements linéaires
lmtf1=lm(tf1~r,data=data1)
summary(lmtf1)
a1=lmtf1[["coefficients"]][["r"]]
b1=lmtf1[["coefficients"]][["(Intercept)"]]
plot(r,tf1,type="p",
     xlab="r",ylab="tf",main="tf en fonction de r",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement linéaire: tf=",sprintf("%.3f",a1),"r+",sprintf("%.3f",b1),"s"))
abline(lmtf1)
arrows(r, tf1-inctf1, r, tf1+inctf1,
       angle=90, code=3, length=0.05)

lmtf2=lm(tf2~L,data=data2)
summary(lmtf2)
a2=lmtf2[["coefficients"]][["L"]]
b2=lmtf2[["coefficients"]][["(Intercept)"]]
plot(L,tf2,type="p",
     xlab="L",ylab="tf",main="tf en fonction de L",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement linéaire: tf=",sprintf("%.3f",a2),"L+-",sprintf("%.3f",b2),"s"))
abline(lmtf2)
arrows(L, tf2-inctf2, L, tf2+inctf2,
       angle=90, code=3, length=0.05)

lmtf3=lm(tf3~alpha,data=data3)
summary(lmtf3)
a3=lmtf3[["coefficients"]][["alpha"]]
b3=lmtf3[["coefficients"]][["(Intercept)"]]
plot(alpha,tf3,type="p",
     xlab="alpha",ylab="tf",main="tf en fonction de alpha",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement linéaire: tf=",sprintf("%.3f",a3),"alpha+-",sprintf("%.3f",b3),"s"))
abline(lmtf3)
arrows(alpha, tf3-inctf3, alpha, tf3+inctf3,
       angle=90, code=3, length=0.05)

#Viscosité
gamma=1e-3 #en N.m^-1
mu1=(tf1*gamma*r^4)/(L[1]*Z01^4*(1-cos(alpha[1])))
incmu1=mu1*sqrt((incr/r)^2+(incL/L[1])^2+(incalpha/alpha[1]))
mu2=(tf2*gamma*r[1]^4)/(L*Z01^4*(1-cos(alpha[1])))
incmu2=mu2*sqrt((incr/r[1])^2+(incL/L)^2+(incalpha/alpha[1]))
mu3=(tf3*gamma*r[1]^4)/(L[1]*Z01^4*(1-cos(alpha)))
incmu3=mu3*sqrt((incr/r[1])^2+(incL/L[1])^2+(incalpha/alpha))

par(mfrow=c(2,2))
plot(r,mu1,type="o",
     xlab="r",ylab="mu(r)",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue")
arrows(r, mu1-incmu1, r, mu1+incmu1,
       angle=90, code=3, length=0.05)
plot(L,mu2,type="o",
     xlab="L",ylab="mu(L)",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue")
arrows(L, mu2-incmu2, L, mu2+incmu2,
       angle=90, code=3, length=0.05)
plot(alpha,mu3,type="o",
     xlab="alpha",ylab="mu(alpha)",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue")
arrows(alpha, mu3-incmu3, alpha, mu3+incmu3,
       angle=90, code=3, length=0.05)


matplot(r,cbind(mu1,mu2,mu3),type="o", lwd=1.5,
     xlab="r",ylab="mu",main="mu en fonction de r, alpha et L",
     col=c("red","lightblue","brown"),
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue")
legend("topright", legend=c("mu(r)","mu(L)","mu(alpha)"),
       col=c("red","lightblue","brown"),lty=1, lwd=1.5)
arrows(r, mu1-incmu1, r, mu1+incmu1,
       angle=90, code=3, length=0.05)
arrows(r, mu2-incmu2, r, mu2+incmu2,
       angle=90, code=3, length=0.05)
arrows(r, mu3-incmu3, r, mu3+incmu3,
       angle=90, code=3, length=0.05)

