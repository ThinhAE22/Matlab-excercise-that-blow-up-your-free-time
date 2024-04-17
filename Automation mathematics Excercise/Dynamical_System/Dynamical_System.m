%% ex 1
clear
close all
R=1
L=3
U=5
%L*i'(t) + R*i(t) =U
%q(0)=0

%first degree equation
%a*y'(t)+b*y(t)=A,y(0)=y0
%solution_formulas.pdf, 1b i) p.2
a=L
b=R
A=U
y0=0
la=-b/a
%solution
%i(t)
i=@(t) A/b+(y0-A/b)*exp(la*t)
%i'(t)
di=@(t) (y0-A/b)*exp(la*t)*la

tau=-1/la %time constant
tmax=5*tau
t=0:tau/100:tmax;
plot(t,R*i(t),'linewidth',1.5)
hold on
plot(t,L*di(t),'linewidth',1.5)
hold off
grid
xlabel('time t')
xlim([0,tmax])
legend({'Ur','UL'},'fontsize',12)

%% ex 2
clear
close all
%x'(t) + k*x(t)/V = k*(C*sin(wt)+D)
%x(0)=0


%first degree equation
%a*y'(t)+b*y(t)=A,y(0)=y0
%Initiate
k=50
V=1000
C=2
D=2 
w=0.2

a=1
b=k/V
p=@(t) C*sin(w*t)+D
y0=0
la=-b/a
%amplification
K= 1/sqrt(b^2 + (a*w)^2)
%phase shift
theta = atan2(-a*w, b)
%transient
Ctr =@(t) y0 - K*k*C*sin(theta) - k*D/b
yh =@(t) Ctr(t)*exp(la*t)
%solution
x=@(t) K*k*C*sin(w*t+theta) + k*D/b + yh(t)

tau=V/k %time constant
tmax=10*tau
t=0:tau/100:tmax;

%omega
dw = 10
wmax = 10*dw
w_range = 0:dw/1000:wmax;

%amplification
K_w =@(w) 1./sqrt(b^2 + (a*w).^2)
theta_w = @(w) atan2(-a*w, b)

figure(1)
plot(t,p(t),'linewidth',1.5)
hold
plot(t,x(t)/V,'linewidth',1.5)
hold off
grid
xlabel('time (min)')
ylabel('suolapitoisuudet (kg/litra)')
xlim([0,tmax])

figure(2)
plot(t,K*k/V*C*sin(w*t)+D,'linewidth',1.5)
hold
plot(t,x(t)/V,'linewidth',1.5)
hold off
grid
xlabel('time (min)')
ylabel('suolapitoisuudet (kg/litra)')
title('the x(t)/V = K*C*sin(w*t +theta)+D as t increase')
xlim([0,tmax])

figure(3)
subplot(2,1,1)
plot(w_range,K_w(w_range)*k/V,'linewidth',1.5)
hold off
grid
title('vahvitus K')
xlim([0,wmax]/10)
ylim([0,1])
subplot(2,1,2)
plot(w_range,theta_w(w_range),'linewidth',1.5)
hold off
grid
xlabel('kulmataajuus w')
title('vaihesiirto theta')
xlim([0,wmax]/10)
ylim([-2,0])

%ex3
clear
close all

%L*q''(t) + R*q'(t) + 1/C * q(t) = U
%q'(0) = 0
%q(0) = 0

%Initiate
R = 1
L = 0.2
L1 = 2
C = 0.5
U = 1
q0 = 0
qd0 = 0
%function second_degree_homogeneous.m
[x,dx,d2x,tau]=second_degree_homogeneous(L,R,1/C,q0-U/(1/C),qd0)
[x1,dx1,d2x1,tau1]=second_degree_homogeneous(L1,R,1/C,q0-U/(1/C),qd0)

q  =  @(t) U/(1/C)+ x(t)
qd =  @(t) dx(t)
q2d = @(t) d2x(t)

q1  =  @(t) U/(1/C)+ x1(t)
qd1 =  @(t) dx1(t)
q2d1 = @(t) d2x1(t)

tmax=10*tau
t=0:tau/100:tmax;

tmax1=10*tau1
t1=0:tau1/100:tmax1;

figure(1)
plot(t,q(t)/C,'linewidth',1.5)
grid
hold
plot(t,qd(t)*R,'linewidth',1.5)
plot(t,q2d(t)*L,'linewidth',1.5,'g')
hold off
xlabel('time t')
legend({'U_C','U_R','U_L'},'fontsize',12)
title('R=1, L=2, C = 0.5, U = 1')
xlim([0,tmax])

figure(2)
plot(t1,q1(t1)/C,'linewidth',1.5)
grid
hold
plot(t1,qd1(t1)*R,'linewidth',1.5)
plot(t1,q2d1(t1)*L1,'linewidth',1.5,'g')
hold off
xlabel('time t')
legend({'U_C','U_R','U_L'},'fontsize',12)
title('R=1, L=0.2, C = 0.5, U = 1')
xlim([0,tmax1])

%ex 4
clear
close all
%Formular
%M*y''(t) + b*y'(t) + k*y(t) = m*R*w^2*sin(w*t)
%y(0) = 0
%y'(0) = 0

%Initial
M = 100
m= 10
k = 500000
b = 100
R = 0.1
w = 100
y0 = 0
yd0 = 0

%aplification
K = 1/sqrt((b*w)^2 + (k-M*w^2).^2)

%phase shift
theta = atan2(-b*w,k-M*w^2)

%intial condition
yh0 = y0 -K*(m*R*w^2)*sin(theta)
ydh0 = yd0 -K*(m*R*w^2)*w*cos(theta)

[x,dx,d2x,tau]=second_degree_homogeneous(M,b,k,yh0,ydh0)
y =@(t) K*(m*R*w^2)*sin(w*t + theta) + x(t)

tmax=20*tau
t=0:tau/1000:tmax;

%Amplitude
w_range = 0:(10*sqrt(k/M))/700:10*sqrt(k/M);

K_w = @(w_range) 1./(sqrt((b*w_range).^2 + (k-M*w_range.^2).^2))

figure(1)
plot(t,y(t),'linewidth',1.5)
grid
hold
plot(t,K*(m*R*w^2)*sin(w*t),'linewidth',1.5)
hold off
xlabel('time t')
ylabel('paikka y(t)')
title('M =100, m = 10,k = 500000, b= 100, R = 0.1, w= 100')
xlim([0,1])
ylim([-0.05,0.05])

figure(2)
plot(t,y(t),'linewidth',1.5)
grid
hold
plot(t,K*(m*R*w^2)*sin(w*t),'linewidth',1.5)
hold off
xlabel('time t')
ylabel('paikka y(t)')
title('M =100, m = 10,k = 500000, b= 100, R = 0.1, w= 100')
xlim([20,21])
ylim([-0.02,0.02])

figure(3)
plot(w_range,(m*R*w^2)*K_w(w_range),'linewidth',1.5)
grid
hold
plot(100,(m*R*w^2)*K_w(100),'markersize',20)
hold off
xlabel('kulmanopeus omega')
ylabel('y(t):n amplitudi')
title('M =100, m = 10,k = 500000, b= 100, R = 0.1, w= 100')
ylim([0,0.8])
xlim([0,800])

%ex 5
%case 1
clear
close all
L=2
R=1
C=0.1
U0=10
T=5 %period
q0=0 %x(0)
qd0=0 %x'(0)=v(0)


dt=0.001 %time step
%increasing part
t1=0:dt:T;
U1=(U0/T)*t1;
%decreasing part
t2=T+dt:dt:T+dt;
U2=-U0/T*(t2-T)+U0;
%one sawtooth
t=[t1,t2];
U=[U1,U2];

%5 sawtooth
t=[t,t+T,t+2*T,t+3*T,t+4*T];
U=[U,U,U,U,U];

figure(1)
subplot(2,1,1)
plot(t,U,'linewidth',1.5)
grid
title('jannite F(t)')

%Formular
% L*q''(t) + R*q'(t) + 1/C*q(t) = u(t)

%% values of x(t),x'(t),x''(t) numerically
N=length(t)
q=zeros(1,N);
qd=zeros(1,N);
q2d=zeros(1,N-1);
q(1)=q0;
qd(1)=qd0;


for n=1:N-1   
   q2d(n)=1/L*(U(n)-R*qd(n)-(1/C)*q(n)); %x''(t)=a(t)
   q(n+1)=q(n)+qd(n)*dt+1/2*q2d(n)*dt^2;
   qd(n+1)=qd(n)+q2d(n)*dt;%x'(t)=v(t)
end



%% 
figure(2) 
subplot(3,1,1)
plot(t,(1/C)*q,'linewidth',1.5,'r')
grid
title('L = 2, R = 1, C= 0.1, U = 10, T=5')
ylabel('jannite U_c')
subplot(3,1,2)
plot(t,qd*R,'linewidth',1.5,'g')
grid
ylabel('jannite U_R')
subplot(3,1,3)
plot(t(1:N-1),q2d*L,'linewidth',1.5,'b')
grid
xlabel('time t')
ylabel('jannite U_L')

%case 2
clear
close all
L=0.5
R=1
C=2
U0=10
T=2 %period
q0=0 %x(0)
qd0=0 %x'(0)=v(0)


dt=0.001 %time step
%increasing part
t1=0:dt:T;
U1=(U0/T)*t1;
%decreasing part
t2=T+dt:dt:T+dt;
U2=-U0/T*(t2-T)+U0;
%one sawtooth
t=[t1,t2];
U=[U1,U2];

%5 sawtooth
t=[t,t+T,t+2*T,t+3*T,t+4*T];
U=[U,U,U,U,U];

figure(1)
subplot(2,1,1)
plot(t,U,'linewidth',1.5)
grid
title('jannite F(t)')

%Formular
% L*q''(t) + R*q'(t) + 1/C*q(t) = u(t)

%% values of x(t),x'(t),x''(t) numerically
N=length(t)
q=zeros(1,N);
qd=zeros(1,N);
q2d=zeros(1,N-1);
q(1)=q0;
qd(1)=qd0;


for n=1:N-1   
   q2d(n)=1/L*(U(n)-R*qd(n)-(1/C)*q(n)); %x''(t)=a(t)
   q(n+1)=q(n)+qd(n)*dt+1/2*q2d(n)*dt^2;
   qd(n+1)=qd(n)+q2d(n)*dt;%x'(t)=v(t)
end



%% 
figure(2) 
subplot(3,1,1)
plot(t,(1/C)*q,'linewidth',1.5,'r')
grid
title('L = 0.5, R = 1, C= 2, U = 10, T=2')
ylabel('jannite U_c')
subplot(3,1,2)
plot(t,qd*R,'linewidth',1.5,'g')
grid
ylabel('jannite U_R')
subplot(3,1,3)
plot(t(1:N-1),q2d*L,'linewidth',1.5,'b')
grid
xlabel('time t')
ylabel('jannite U_L')



%ex 6
%J*L*w''(t) + (R*J +b*L)*w'(t) + (b*R +K^2)w(t) = K*u(t)

%w't = (K*u(t) - J*L*w''(t)-(b*R +K^2)w(t))/(R*J +b*L)

%w''(t) = (K*u(t) - (R*J +b*L)*w'(t)-(b*R +K^2)w(t))/(J*L)
clear
close all
%DC motor
J=0.1
L=0.5
R=0.3
bm= 0.08
K=0.1
%PID
%fuction u(k)=Kp*error+Ki*integral+Kd*derivative
a= J*L
b= (R*J +bm*L)
c =(b*R +K^2)
Kp=5
Ki=2
Kd=0.001
dt=0.01 %time step
N=1000
t=0:dt:N*dt;
r=ones(1,N+1);%setpoint
y=zeros(1,N+1);
dy=zeros(1,N);
u=zeros(1,N);
e=zeros(1,N);
y(1)=0 
dy(1) = 0
integral=0 
previous_error=0

for k=1:N
    error=r(k)-y(k);
    e(k)=error;
    integral=integral+error*dt;
    derivative=(error-previous_error)/dt;
    u(k)=Kp*error+Ki*integral+Kd*derivative;%output of PID
    d2y = (K*u(k) - (R*J +bm*L)*dy(k)-(b*R +K^2)*y(k))./(J*L);
    dy(k+1)=dy(k)+d2y*dt;
    y(k+1)=y(k)+dy(k)*dt;
    previous_error=error;
end


figure(1)
subplot(3,1,1)
plot(t(1:end-1),e,'linewidth',1.5)
grid
title('error e(t) = r(t)-y(t)')

subplot(3,1,2)
plot(t(1:end-1),u,'linewidth',1.5)
grid
title('u(t) = ouput of PID ')

subplot(3,1,3)
plot(t,y,'linewidth',1.5)
grid
xlabel('time t')
title('w(t) = output of the system')
