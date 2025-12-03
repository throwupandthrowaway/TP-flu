#Saisir les masses mesurées
tq=as.integer(seq(0,120,5)) #en s
m=0:23 #en kg
rho=999.19 #en kg.m^-3
tq=tq[1:length(m)]
Q=numeric(length(tq))
for (i in 2:length(m)) {Q[i]=(m[i]-m[i-1])/((tq[i]-tq[i-1])*rho)}
print(Q)
Q=Q[-1]
tq=tq[-1]
data1=data.frame(Q,tq)
#Ajustement linéaire
lmQ=lm(Q~tq,data=data1)
summary(lmQ)
a1=lmQ[["coefficients"]][["tq"]]
b1=lmQ[["coefficients"]][["(Intercept)"]]
plot(tq,Q,,type="o",
     xlab="t",ylab="Q(t)",main="Débit volumique en fonction du temps",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement: Q(t)=",a1,"t+",b1))
abline(lmQ)
Qadj=a1*tq+b1


#Saisir les hauteurs mesurées
th=as.integer(seq(0,120,5)) #en s
h=0:23 #en m
g=9.80665 #en m.s^-2
th=th[1:length(h)]
th=th[-1]
dP=-rho*g*th
print(dP)
data2=data.frame(dP,th)
#Ajustement linéaire
lmdP=lm(dP~th,data=data2)
summary(lmdP)
a2=lmdP[["coefficients"]][["th"]]
b2=lmdP[["coefficients"]][["(Intercept)"]]
plot(th,dP,,type="o",
     xlab="t",ylab="dP(t)",main="Différence de pression près du poreux en fonction du temps",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement: dP(t)=",a2,"t+",b2))
abline(lmdP)
dPadj=a2*th+b2

#Proportionnalité entre Q et dP
prop=function(t){(a1*t+b1)/a2*t}
x=seq(0,120,1)
val=prop(x)
data3=data.frame(val,x)
lmx=lm(val~x,data=data3)
summary(lmx)
a3=lmx[["coefficients"]][["x"]]
b3=lmx[["coefficients"]][["(Intercept)"]]
plot(x,val,type="o",
     xlab="t",ylab="Q/dP",main="Rapport entre débit volumique et différence de pression",
     col="red",
     fg="darkgreen",
     col.axis="blue",
     col.lab="blue",
     sub=paste("Ajustement: Q/dP=",a3,"t+",b3))
abline(lmx)
Rh=abs(a3)

#Reynolds
D=10e-2 #en m
dp=1e-2 #en m
mu=1.1375e-3 #en Pa.s
Re=(4*rho*D*Qadj)/(pi*mu*(dp^2))
sdRe=sd(Re)
Re=mean(Re)
print(Re)


#Perméabilité k
k=(4*mu)/(pi*(dp^2)*Rh)

#Taille des grains
phi=c(0.65,0.74)
d=(phi/(dp*(1-phi)))*sqrt(720*mu/(pi*Rh*(1-phi)))
ht0=(1-exp(-1))*h[1] #On détermine tau expérimentalement
tau=10 #en s
d=(phi*D/(dp*(1-phi)))*sqrt(180*tau*mu/(rho*g*(1-phi)))

#Incertitudes
