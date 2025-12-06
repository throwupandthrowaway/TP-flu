#Saisir les masses mesurées
tq=c(3.44,6.70,10.30,13.7,17.35,20.68,24.51,28.22,31.70,35.74,39.75,43.55,46.90,51.13,55.13,60.65,64.92,69.72,74.71,80.65,87.42,93.56,100.16,107.90,117.12,126.32,142.22,153.92) #en s
m=c(9,24,45,65,85,100,125,140,160,180,198,223,238,255,273,295,312,330,350,368,385,408,423,439,455,467,485,500)*1e-3 #en kg
rho=1000 #en kg.m^-3
Q=numeric(length(tq))
for (i in 2:length(m)) {Q[i]=(m[i]-m[i-1])/((tq[i]-tq[i-1])*rho)}
print(Q)
Q=Q[-1]
tq=tq[-1]
data1=data.frame(Q,tq)

#Incertitudes
incm=2.5e-3 #en kg
inctq=250e-3 #en s (incertitude du lecteur vidéo)
incQ=Q*sqrt((incm/m)^2+(inctq/tq)^2)

par(mfrow=c(2,2))
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
th=tq #en s
h=rev(seq(0,27,1))*1e-3 #en m
h=h[1:length(th)]
g=9.81 #en m.s^-2
dP=rho*g*h
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

data3=data.frame(Q,dP)
lmx=lm(Q~0+dP,data=data3)
summary(lmx)
a3=lmx[["coefficients"]][["dP"]]

plot(dP,Q,type="p",
     xlab="dP",ylab="Q(dP)",main="Rapport entre débit volumique et différence de pression",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     xlim=c(0,max(dP)),
     ylim=c(0,max(Q)),
     sub=paste("Ajustement: Q(dP)=",a3,"dP"))

arrows(dP, Q-incQ, dP, Q+incQ,
       angle=90, code=3, length=0.05)
abline(lmx)

Rh=abs(1/a3)
incRh=Rh*sqrt((incQadj/Qadj)^2+(incdPadj/dPadj)^2)
incRh=mean(incRh)
print(paste("Rh=",Rh,"+-",incRh,"m.s"))

#Reynolds
D=4e-3 #en m
incD=0.5e-3 #en m
d=2.7e-2 #en m
incd=0.5e-3 #en m
mu=1e-3 #en Pa.s
Re=(4*rho*D*max(Q))/(pi*mu*(d^2))
incRe=Re*sqrt((incD/D)^2+4*(incd/d)^2+(incQadj/max(Q))^2)
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
tau=39.75 #en s
inctau=0.01 #en s
d=(phi*D/(dp*(1-phi)))*sqrt(180*tau*mu/(rho*g*(1-phi)))
incd=d*sqrt((incD/d)^2+(incdp/dp)^2)
print(paste("d=",d,"+-",incd,"m"))
