%
dbstop at 72
clc
clearvars -except output tu tl
tl1=100;tf = 700;
tu1=400;
close all
C0_data = dlmread('initial_concentration.txt');

C0 = zeros(34,1);

C0(1,1) =  25E-12; % TF
% alpha=0.5;
% C0(2,1) = alpha*C0_data(1,1); % VII
% C0(8,1) = alpha*C0_data(1,2); % X
% C0(11,1) = alpha*C0_data(1,3); % IX
% C0(14,1) = alpha*C0_data(1,4); % II
% C0(15,1) = alpha*C0_data(1,5); % VIII
% C0(21,1) = alpha*C0_data(1,6); % V
% C0(26,1) = alpha*C0_data(1,7); % TFPI
% C0(29,1) = alpha*C0_data(1,8); % AT
% % % 
C0(2,1) = 0.87*C0_data(1,1); % VII
C0(8,1) = 1.1*C0_data(1,2); % X
C0(11,1) = 1.07*C0_data(1,3); % IX
C0(14,1) = 1.15*C0_data(1,4); % II
C0(15,1) = 0.8*C0_data(1,5); % VIII
C0(21,1) = 1.05*C0_data(1,6); % V
C0(26,1) = 0.63*C0_data(1,7); % TFPI
C0(29,1) = 1.13*C0_data(1,8); % AT

% C0(2,1) = 0.87*C0_data(1,1); % VII
% C0(8,1) = 1.1*C0_data(1,2); % X
% C0(11,1) = 1.07*C0_data(1,3); % IX
% C0(14,1) = 1.15*C0_data(1,4); % II
% C0(15,1) = 0.04*C0_data(1,5); % VIII
% C0(21,1) = 1.05*C0_data(1,6); % V
% C0(26,1) = 0.63*C0_data(1,7); % TFPI
% C0(29,1) = 1.13*C0_data(1,8); % AT

% C0(4,1) = C0(2,1)*0; % VIIa
C0(4,1) = C0(2,1)*0.01; % VIIa

t0 = 0; 
deltat=0.5;
tspan = [t0:deltat:tf];

options = odeset('RelTol', 1e-10,'AbsTol',1e-12*ones(1,34));
[T,C] = ode23t(@reaction_rates,tspan, C0, options);

% [T,C] = ode23s(@reaction_rates,tspan, C0);
B=C(tl1/deltat:tu1/deltat,1:34);
t0=tl1/deltat;
tem=JacobianFunc(t0-1,C(t0-1,:));
Jtotal=tem';

for t=tl1/deltat:tu1/deltat
    tem=JacobianFunc(t,C(t,:));
Jtotal=[Jtotal;tem'];
end
%tl=[tl tl1];
%tu=[tu tu1];
%B=C(50:600,1:34); -> pretty good
%B=C(50:700,1:34);->very good
% [coeff,score,latent,tsquared,explained,mu] = pca(C(50:250,1:34));
% norm(B-B*coeff(:,1)*coeff(:,1)'-B*coeff(:,2)*coeff(:,2)'-B*coeff(:,3)*coeff(:,3)','fro')/norm(B,'fro')
%figure; stem(explained);
Jtotal=Jtotal';

nullJtotal=null(Jtotal)';
Ainv=nullJtotal(1:6,:);
h=figure;plot(T,C);

xlabel('t');ylabel('y_i')
%set(gca, 'YScale', 'log')
%Hoang.output_standardize_figure(h,'H:\RESEARCH\thrombin\figures\timeseries','.eps');

sim=corr(B);
h=figure;plot(sim);
ylabel('Correlation');xlabel('i - index for variable')
%Hoang.output_standardize_figure(h,'H:\RESEARCH\thrombin\figures\correlation','.eps');
% 14->25 are highly correlated, ~1
ind=[9:10,14:25];
%% svd to do dimension reduction
 X=B;
 [U,S,V] = svd(X);
for t1=6:6
clear('St','Vt','Ut','z','Xreconstruct');
Ut=U(:,1:t1);
St=S(1:t1,1:t1);
Vt=V(:,1:t1);
Xsvdt=Ut*St*Vt';
A=(St*Vt')';
for i1=1:length(tspan)
fAz(:,i1)=A'*reaction_rates(i1,C(i1,:));
end
errorate_svd(t1)=1-norm(X-Xsvdt,'fro')/norm(X,'fro');
optionst = odeset('RelTol', 1e-9,'AbsTol',1e-10*ones(1,t1));
[T1,z] = ode23s(@(t,z)reaction_rates_new(t,z,A,Ainv),tspan, Ainv*C0);
Xreconstruct=A*z';
sizedata(t1)=size(z,1);
h=figure;
hold on;
plot(Xreconstruct(7,:),'r')
plot(C(:,7));
plot(Xreconstruct(7,:)'-C(1:length(Xreconstruct(7,:)),7),'g')
legend('Reconstructed','Original','Error')

errorate_reconstruct(t1)=norm(C(1:size(Xreconstruct,2),:)-Xreconstruct','fro')/norm(C(1:size(Xreconstruct,2),:),'fro');
end
output1=[sizedata;errorate_reconstruct]
output=[output;output1];
%save('output.mat','output','tu','tl')
h=figure;stem(1-errorate_svd);
t='Error rate';
%='$\frac{||S-\hat{S}||_{F}}{||S||_{F}}$';
ylabel(t);xlabel('Number of eigenvectors');
Hoang.output_standardize_figure(h,'H:\RESEARCH\thrombin\figures\errorrate','.eps');

h=figure;
subplot(2,1,1);plot(Xreconstruct');title('Reconstructed time series');ylabel('y_i')
subplot(2,1,2);plot(C);title('Original time series');
xlabel('t');ylabel('y_i')
Hoang.output_standardize_figure(h,'H:\RESEARCH\thrombin\figures\timeseries_reconstruct','.eps');
h=figure;plot(Xreconstruct'-C);Hoang.output_standardize_figure(h,'H:\RESEARCH\thrombin\figures\timeseries_reconstruct_error','.eps');

h=figure;
hold on;
plot(Xreconstruct(14,:),'r')
plot(C(:,14));
plot(Xreconstruct(14,:)'-C(:,14),'g')
legend('Reconstructed','Original','Error')

Hoang.output_standardize_figure(h,'H:\RESEARCH\thrombin\figures\timeseries_reconstruct_thrombin_4components_use3rdinital_original','.eps');

% save('SolutionData.mat','T', 'C')
% filename = 'SolutionData.mat';
% m = matfile(filename);