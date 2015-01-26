%This code is to study how to improve POD in model reduction
function pod()
initialx=rand(6,1);
deltat=0.01;
tspan=0:deltat:5;
odeset('MaxStep',10^-6);
[t,x] = ode45(@rhs,tspan,initialx);
[eigv,vec]=eig(corr(x));
rho=eigv(:,1:3)';
%rho=[eye(3) zeros(3)];
P=rho'*rho;


A1=[-0.1 0 0;0 -0.1732 2;0 -2 -0.1732];
A12=[0.3893 0.5179 -1.543;1.39 1.3 0.8841;0.0629 -0.9078 -1.184];
A2=[-1 0 0;0 -1.226 -0.7080;0 0.7080 -1.226];
Jacobian=[A1 A12;zeros(size(A1,1),size(A2,2)) A2];
A=blkdiag(rho(1,:),rho(2,:),rho(3,:))
B1=rho*Jacobian;
B=blkdiag(B1(1,:),B1(2,:),B1(3,:));
%rho*rho'=I; rho*Jacobian*rho' min.
% rho*X=I;
% rho*Jacobian*X min
[t,xhat]=ode45(@(t,xhat)reduced(t,xhat,P),tspan,initialx);

rhoc=[zeros(3) eye(3)];
 plot(t,xhat(:,1:3),'--');
 xtheuta=P*(x-repmat(mean(x),size(x,1),1))'+repmat(mean(x),size(x,1),1)';
[V,D]= eig(cov(x))
 hold on;plot(t,xtheuta(1:3,:)')
 hold off;
 Hoang.standardize_fig(get(0,'CurrentFigure'))
 Hoang.save_standized_fig(get(0,'CurrentFigure'),'ex1')
 
 %close all;
end
function rhs=rhs(t,x)
A1=[-0.1 0 0;0 -0.1732 2;0 -2 -0.1732];
A12=[0.3893 0.5179 -1.543;1.39 1.3 0.8841;0.0629 -0.9078 -1.184];
A2=[-1 0 0;0 -1.226 -0.7080;0 0.7080 -1.226];
A=[A1 A12;zeros(size(A1,1),size(A2,2)) A2];
rhs=A*x;
end
function rhsreduced=reduced(t,xhat,P)
rhsreduced=P*rhs(t,xhat);
end