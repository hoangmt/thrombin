% Model reduction demo for the transient evolution of the temperature field
% around a cylinder immersed into fluid.

%  Antonio Cosmin Ionita
%  Rice University
%  2013/05/30

% clean up
clear all, close all, format compact, format short e, clc
% setup plots
set(0,'DefaultFigurePosition', [100 100 940 500]);
set(0,'defaultlinelinewidth',2)
set(0,'defaultlinemarkersize',20)
set(0,'defaultaxesfontsize',16)
set(0,'defaulttextinterpreter','latex')

% path to plotting functions adapted from the rbMIT software package, see
% http://augustine.mit.edu/ for details and other interesting examples.
addpath('./conv_diff_plotting/')
load conv_diff_VIS.mat

echo on;
% Transient evolution of the temperature field around a cylinder:
%
% state-space matrices of the full-order finite element model
%
%     E dx/dt = A x + B u
%           y = C x + D u
%
% where A = A1+1/mu*A2, u is the input, x is the state, and y is the output
%
echo off;
load FEM_system_878.mat E A1 A2 B C D
whos('E','A1','A2','B','C','D')


% the finite element discretization has n = 878 states
n = length(B);

% the Peclet number (set in the interval [0.1, 10])
mu = 0.1;

% the corresponding state-space matrix
A = A1+1/mu*A2;

% setup time grid
dt = .01;         % time step
t = [0:dt:1].';   % time grid

% the input excitation
u = 10*t;

%
% evolve the finite element (FE) model in time using Backward Euler
%

% initialization
xFE(:,1) = zeros(n,1);           % zero initial conditions
yFE(:,1) = C*xFE(:,1) + u(1)*D;
% evolution
[L,U] = lu(E-dt*A);
for i = 2:length(t)
   % the state trajectory x
   xFE(:,i) = U\(L\( E*xFE(:,i-1) + dt*B*u(i) )); 
   % the outputs y
   yFE(:,i) = C*xFE(:,i) + u(i)*D;
end

%
% from the state trajectories, extract a low order approximation,
% i.e., Proper Orthogonal Decomposition (POD) approach
%
[Y,S,~] = svd(xFE);

figure, clf
   semilogy(diag(S)./S(1,1),'.')
   axis tight
   title('Singular values of state snapshots')

%
% since the singular values exhibit a fast decay, we can select a low value
% for the reduced order k
%
k = 10;

%
% project to get the reduced-order state-space matrices 
%
Yk = Y(:,1:k);

Ek = Yk'*E*Yk;
Ak = Yk'*A*Yk;
Bk = Yk'*B;
Ck = C*Yk;

%
% evolve the POD reduced-order model in time
%

% initialization
xPOD(:,1) = zeros(k,1);           % zero initial conditions
yPOD(:,1) = Ck*xPOD(:,1) + u(1)*D;
% evolution
[L,U] = lu(Ek-dt*Ak);
for i = 2:length(t)
   % the state trajectory x
   xPOD(:,i) = U\(L\( Ek*xPOD(:,i-1) + dt*Bk*u(i) )); 
   % the outputs y
   yPOD(:,i) = Ck*xPOD(:,i) + u(i)*D;
end

%
% compare the norm of the error between the full-order outputs and
% reduced-order outputs
%
figure,
   semilogy(t,sum(abs(yFE-yPOD)))
   xlabel('time $t$')
   title('Error between full-order outputs and reduced-order outputs, $\| y_n- y_k\|_2$')
   xlim([t(1),t(end)])

%
% extract a .gif animation from the state snapshots, showing the
% temperature evolution around the cylinder as the system evolves in time
%

% use the plotting functions from the rbMIT software package
tt = conv_diff_VIS.t;
pmu=geoplot('conv_diff',mu);   
pp=pmnode(pmu',tt(1:3,:)');


% filename of the gif animation
filename = 'convection_diffusion_POD.gif';

% the animation shows only 11 snapshots
idx = [2; [11:10:length(t)]'];
t = t(idx);          % truncate to 11 time stamps
yn = yFE(:,idx);     % truncate to 11 snapshots
yk = yPOD(:,idx);    % truncate to 11 snapshots

for i = 1:length(t)

   % plot the FE full-order model
   figure(3), subplot(1,3,1)
      plotsol(pp,tt(1:3,:)',yn(:,i)); 
      axis equal; axis tight; axis on; 
      hcb = colorbar('WestOutside');
      set(hcb, 'Visible', 'off')
      %set(hcb,'YTick',floor([min(yn(:,i)),max(yn(:,i))]*10000)/10000)
      set(gca,'XTick',[])
      set(gca,'YTick',[])
      xlabel(['time t = ', num2str(floor(t(i)*100)/100)])
      title(['Full-order model $y_n$', sprintf('\n'), ...
         '($n$ = ', num2str(n), ' states)'])
      set(gcf,'renderer','painters')

   % plot the POD reduced-order model
   figure(3), subplot(1,3,2)
      plotsol(pp,tt(1:3,:)',yk(:,i)); 
      axis equal; axis tight; axis on; 
      hcb = colorbar('EastOutside');
      set(hcb, 'Visible', 'off')
      set(gca,'XTick',[])
      set(gca,'YTick',[])
      xlabel(['time t = ', num2str(floor(t(i)*100)/100)])
      title(['Reduced-order model $y_k$ ', sprintf('\n'), ...
         '($k$ = ', num2str(k), ' states)'])
      set(gcf,'renderer','painters')

   % also plot the error
   figure(3), subplot(1,3,3)      
      plotsol(pp,tt(1:3,:)',abs(yn(:,i)-yk(:,i))); 
      axis equal; axis tight; axis on; 
      set(gca,'CLim',[0, max(max(abs(yn-yk)))])
      colorbar('EastOutside');
      set(gca,'XTick',[])
      set(gca,'YTick',[])
      xlabel(['time t = ', num2str(floor(t(i)*100)/100)])
      title(['$| y_n-y_k |$', sprintf('\n'),' '])
      set(gcf,'renderer','painters')
   
   % add titles
   annotation(gcf,'textbox',[0.28 0.905 0.246808510638298 0.076],...
      'String',{'Temperature field'},'FontWeight','bold',...
      'FontSize',18,'LineStyle','none');
   annotation(gcf,'textbox',[0.69 0.905 0.278723404255319 0.072],...
      'String',{'Approximation error'},'FontWeight','bold',...
      'FontSize',18,'LineStyle','none');

   drawnow

   % write to .gif file
   frame = getframe(gcf);
   [im,map] = frame2im(frame);
   [imind,cm] = rgb2ind(im,256);
   if i == 1;
      imwrite(imind,cm,filename,'gif','DelayTime',0.75,'Loopcount',inf);
   else
      imwrite(imind,cm,filename,'gif','DelayTime',0.75,'WriteMode','append');
   end

end
