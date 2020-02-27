%Implicit NR solver for spreading-vanishing dichotomy on an n-sphere
%Mat Simpson February 2020
clear all
aa=0.0; %aa=(d-1) where d is dimension.
dxi=0.001; %dxi
dt=0.0001;
tf=100.00;
kappa=1.00;
tol=1e-10; 
dldt=0.0;
l=1.0;
pl=1.0;
lcrit=1.570796327; %1.570796327; %2.404825558 %3.141592654;  
maxsteps=tf/dt;

Lrecord=zeros(maxsteps+1,1);
Lrecord(1,1)=1.0;
wsrecord=zeros(maxsteps+1,1);
wsrecord(1,1)=0;

N=1.0/dxi+1;

u=zeros(1,N);
pu=zeros(1,N);
delu=ones(1,N);
xi=zeros(1,N);
for i=1:N
    xi(1,i)=0.0+(i-1)*dxi;
end

for i=1:round(1.0/dxi)
    u(1,i)=1.0; %Initial condition
    pu(1,i)=u(1,i);
end

a=zeros(1,N);
b=zeros(1,N);
c=zeros(1,N);
d=zeros(1,N);

t=0.0;

for i=1:maxsteps
    t=t+dt;
    kk=0;
    delu=ones(1,N);
    while norm(delu,Inf) > tol
        kk=kk+1;
        
a(1,1)=0.0;
b(1,1)=-1.0;
c(1,1)=1.0;
d(1,1)=-1*(u(1,1)-u(1,2)); %LHS Boundary Condition

a(1,N)=0.0;
b(1,N)=1.0;
c(1,N)=0.0;
d(1,N)=0.0; %RHS boundary condition

dldt=kappa*u(1,N-1)/(l*dxi);
for j=2:N-1 %Discretized equations for NR iterations written in Tridiagonal form.
a(1,j)=1.0/(l*dxi)^2-(xi(1,j)*dldt/l + aa/(xi(1,j)*l^2))/(2*dxi);
b(1,j)=-1.0/dt-2.0/(l*dxi)^2+1.0*(1.0-2.0*u(1,j));
c(1,j)=1.0/(l*dxi)^2+(xi(1,j)*dldt/l + aa/(xi(1,j)*l^2))/(2*dxi);
d(1,j)=(u(1,j)-pu(1,j))/dt-(u(1,j-1)-2.0*u(1,j)+u(1,j+1))/(l*dxi)^2 ...
-((u(1,j+1)-u(1,j-1))/(2.0*dxi))*(aa/(xi(1,j)*l^2)+xi(1,j)*dldt/l) ...
-1.0*(u(1,j)*(1-u(1,j)));
end

delu = thomas(N,a,b,c,d);

u(1,:)=u(1,:)+delu(1,:); 
l=pl+dt*dldt;

    end
    
    if mod(i,1000)==0
    fprintf('Time %d\n',t); 
    fprintf('Iteration %d\n',kk);
    fprintf('Wave speed is %d\n',(l-pl)/dt);
    end
Lrecord(i+1,1)=l; 
wsrecord(i+1,1)=(l-pl)/dt;
pu(1,:)=u(1,:);
pl=l;
Lrecord(i+1,1)=l;
end

for i=1:N
    xx(1,i)=l*xi(1,i);
end
plot(xx,u)
xlabel('x')
ylabel('u(x,t)')
%axis([0 70 0 1.2])

figure
plot(dt*(1:maxsteps+1),Lrecord)
hold on
plot([0, tf], [lcrit, lcrit],'-r')
axis([0 tf 0 1.2*max(max(Lrecord),lcrit*1.2)])
xlabel('time')
ylabel('L(t)')

figure
plot(1:(maxsteps+1),wsrecord)
xlabel('time')
ylabel('c')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Algorithm
function x = thomas(N,a,b,c,d)
x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end
    
    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end



