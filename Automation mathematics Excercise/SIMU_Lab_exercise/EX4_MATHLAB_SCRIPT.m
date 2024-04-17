%%
close all
P=50
I=10
D=5
b=1
n=10
K=4
Km=0.1
Ra=20
J=0.1
Kb=1
La=1

num1=[K*Km]
den1=[J*Ra,Kb*Km,K*Km]
poles1=roots(den1)

num2=[Km,0]
den2=[J*Ra,Kb*Km,K*Km]


plot(real(poles1),imag(poles1),'ro','linewidth',1.5)
grid
xlabel('real')
ylabel('imag')
title('poles of G')
