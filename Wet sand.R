#Saisir les masses mesurées
tq=as.integer(seq(1,120,5)) #en s
m=1:23 #en kg
rho=999.19 #en kg.m^-3
tq=tq[1:length(m)]
Q=numeric(length(tq))
for (i in 2:length(m)) {Q[i]=(m[i]-m[i-1])/((tq[i]-tq[i-1])*rho)}
print(Q)
Q=Q[-1]
tq=tq[-1]
data1=data.frame(Q,tq)

#Incertitudes
incm=0.01 #en kg
inctq=0.01 #en s (incertitude du lecteur vidéo)
incQ=Q*sqrt((incm/m[-1])^2+(inctq/tq)^2)

#Ajustement linéaire
lmQ=lm(Q~tq,data=data1)
summary(lmQ)
a1=lmQ[["coefficients"]][["tq"]]
b1=lmQ[["coefficients"]][["(Intercept)"]]
plot(tq,Q,,type="p",
     xlab="t",ylab="Q(t)",main="Débit volumique en fonction du temps",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement: Q(t)=",a1,"t+",b1))
abline(lmQ)
arrows(tq, Q-incQ, tq, Q+incQ,
       angle=90, code=3, length=0.05)
Qadj=a1*tq+b1
incQadj=mean(incQ)

#Saisir les hauteurs mesurées
th=as.integer(seq(1,120,5)) #en s
h=rev(1:23) #en m
g=9.81 #en m.s^-2
th=th[1:length(h)]
dP=rev(-rho*g*h)
print(dP)
data2=data.frame(dP,th)

#Incertitudes
incth=inctq
inch=1e-3 #en m
incdP=dP*sqrt((incth/th)^2+(inch/h)^2)

#Ajustement linéaire
lmdP=lm(dP~th,data=data2)
summary(lmdP)
a2=lmdP[["coefficients"]][["th"]]
b2=lmdP[["coefficients"]][["(Intercept)"]]
plot(th,dP,,type="p",
     xlab="t",ylab="dP(t)",main="Différence de pression près du poreux en fonction du temps",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement: dP(t)=",a2,"t+",b2))
abline(lmdP)
arrows(th, dP-incdP, th, dP+incdP,
       angle=90, code=3, length=0.05)
dPadj=a2*th+b2
incdPadj=mean(incdP)

#Proportionnalité entre Q et dP
prop=function(t){(a1*t+b1)/a2*t}
x=seq(0,120,1)
val=prop(x)
data3=data.frame(val,x)

#Incertitudes
incprop=val*sqrt((incQ/Q)^2+(incdP/dP)^2)

#Ajustement linéaire
lmx=lm(val~x,data=data3)
summary(lmx)
a3=lmx[["coefficients"]][["x"]]
b3=lmx[["coefficients"]][["(Intercept)"]]
plot(x,val,type="p",
     xlab="t",ylab="Q/dP",main="Rapport entre débit volumique et différence de pression",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement: Q/dP=",a3,"t+",b3))
abline(lmx)
arrows(x, val-incprop, x, val+incprop,
       angle=90, code=3, length=0.05)
Rh=abs(1/a3)
incRh=Rh*sqrt((incQadj/Qadj)^2+(incdPadj/dPadj)^2)
incRh=mean(incRh)
print(paste("Rh=",Rh,"+-",incRh,"m.s"))

#Reynolds
D=10e-2 #en m
incD=1e-2 #en m
dp=1e-2 #en m
incdp=1e-2 #en m
mu=1.1375e-3 #en Pa.s
Re=(4*rho*D*Qadj)/(pi*mu*(dp^2))
incRe=Re*sqrt((incD/D)^2+4*(incdp/dp)^2+(incQadj/Qadj)^2)
Re=mean(Re)
incRe=mean(incRe)
print(paste("Re=",Re,"+-",incRe))


#Perméabilité k
k=(4*mu)/(pi*(dp^2)*Rh)
inck=k*sqrt(4*(incdp/dp)^2+(incRh/Rh)^2)
print(paste("k=",k,"+-",inck,"m^2"))

#Taille des grains
phi=c(0.65,0.74)
d=(phi/(dp*(1-phi)))*sqrt(720*mu/(pi*Rh*(1-phi)))
incd=d*sqrt((incdp/dp)^2+(incRh/Rh)^2)
print(paste("d=",d,"+-",incd,"m"))
ht0=(1-exp(-1))*h[1] #On détermine tau expérimentalement
tau=10 #en s
inctau=0.01 #en s
d=(phi*D/(dp*(1-phi)))*sqrt(180*tau*mu/(rho*g*(1-phi)))
incd=d*sqrt((incD/d)^2+(incdp/dp)^2)
print(paste("d=",d,"+-",incd,"m"))
