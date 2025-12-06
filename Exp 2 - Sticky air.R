#Demi-angle de l'entonnoir
alpha=c(40,30,20,5)*(pi/180)
incalpha=pi/360
r=c(6,8,10,12)*((1e-3)/2) #en m
incr=2e-3 #en m
L=c(21,50,10.8,30)*1e-2 #en m
incL=incr #en m

#Mesures du temps
r=sort(r) #Variation de r
Z01=c(95-8,95-6,95-9,105-6)*1e-3 #en m
L=sort(L) #Variation de L
Z02=c(95-12,95-6,95-6,95-7)*1e-3 #en m
alpha=sort(alpha) #Variation de alpha
Z03=c(95-5,95-6,75-7,87-16)*1e-3 #en m
incZ0=1e-2 #en m
tf1=c(43.49,27.5,11.42,11.29)
inctf1=250e-3 #en s
tf2=c(29.77,38.56,40.07,46.53)
inctf2=inctf1
tf3=c(16.62,15.51,10.77,35.57)
inctf3=inctf1
data2=data.frame(tf2,L)

par(mfrow=c(2,2))
#Ajustements linéaires
r2=1/(r^4)
data1=data.frame(tf1,r2)
lmtf1=lm(tf1~r2,data=data1)
summary(lmtf1)
a1=lmtf1[["coefficients"]][["r2"]]
b1=lmtf1[["coefficients"]][["(Intercept)"]]
plot(1/(r^4),tf1,type="p",
     xlab="r",ylab="tf",main="tf en fonction de r",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     xlim=c(0,max(1/(r^4))),
     ylim=c(0,max(tf1)),
     sub=paste("Ajustement linéaire: tf=",sprintf("%.3f",a1),"r+-",sprintf("%.3f",b1),"s"))
abline(lmtf1)
arrows(r2, tf1-inctf1, r2, tf1+inctf1,
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
     xlim=c(0,max(L)),
     ylim=c(0,max(tf2)),
     sub=paste("Ajustement linéaire: tf=",sprintf("%.3f",a2),"L+-",sprintf("%.3f",b2),"s"))
abline(lmtf2)
arrows(L, tf2-inctf2, L, tf2+inctf2,
       angle=90, code=3, length=0.05)

tf32=tf3/(Z03^4)
inctf32=tf32*sqrt((inctf3/tf3)^2+16*(incZ0/Z03)^2)
alpha2=1-cos(alpha)
data3=data.frame(tf32,alpha2)

lmtf3=lm(tf32~alpha2,data=data3)
summary(lmtf3)
a3=lmtf3[["coefficients"]][["alpha2"]]
b3=lmtf3[["coefficients"]][["(Intercept)"]]
plot(alpha2,tf32,type="p",
     xlab="alpha",ylab="tf",main="tf en fonction de alpha",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     xlim=c(0,max(1-cos(alpha))),
     ylim=c(0,max(tf3/(Z03^4))),
     sub=paste("Ajustement linéaire: tf=",sprintf("%.3f",a3),"alpha+-",sprintf("%.3f",b3),"s"))
abline(lmtf3)
arrows(alpha2, tf32-inctf32, alpha2, tf32+inctf32,
       angle=90, code=3, length=0.05)

#Viscosité
gamma=25e-3 #en N.m^-1
mu1=(tf1*gamma*r^4)/(L[4]*Z01^4*(1-cos(alpha[2])))
incmu1=mu1*sqrt((incr/r)^2+(incL/L[1])^2+(incalpha/alpha[1]))
mu2=(tf2*gamma*r[1]^4)/(L*Z01^4*(1-cos(alpha[2])))
incmu2=mu2*sqrt((incr/r[1])^2+(incL/L)^2+(incalpha/alpha[1]))
mu3=(tf3*gamma*r[2]^4)/(L[4]*Z01^4*(1-cos(alpha)))
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

