warning('off','all')
clf 
clear all 
figure (1)

% PO branch 
load('lc_branch_for_Ca_x_plane.mat')
period=x(end-1,:)/1000;
size(period);
Cavar=x(6:28:end-1,:);
xvar =x(5:28:end-1,:);
tr=real(f(1:end-7,1:end));

imin=length(period);
for i=1:imin
tt=tr(:,i)*period(i);
xvari= xvar(:,i);
Cavari= Cavar(:,i);
xvar_av(i)=trapz(tt,xvari)/period(i);
Cavar_av(i)=trapz(tt,Cavari)/period(i);
end


xvarplus = xvar - 0.001;
j=1853;
%periodic orbit critical manifold projection for the action potentials
plot (Cavar,xvarplus,'Color',[.8 .8 .8 ],'LineWidth',2)
hold on
% periodic orbit
% plot (Cavar(:,end-j),xvarplus(:,end-j),'Color',[.75 0 1 ],'LineWidth',3)
% hold on
%plot (Cavar_av,xvar_av,'Color',[.3 .3 .3 ],'LineWidth',3)
hold on
%average voltage of periodic orbit using the trapezoidal method
% fake it 
plot (Cavar_av+0.0065,xvar_av+0.015,'Color',[.3 .3 .3 ],'LineWidth',6)
hold on


%SNIC part 
%--------------------- SNIC whole system gh = 0 -------------------------------
data1 = load('LP_whole_bk_gh0.mat','x');
cod_back = data1.x;
data1 = load('LP_whole_fwd_gh0.mat','x');
cod_fwd = data1.x;
x_cod_bkgh0 = cod_back(end-3,:);
Ca_cod_bkgh0 = cod_back(end-2,:);
x_cod_fdgh0 =  cod_fwd(end-3,1:150);
Ca_cod_fdgh0 = cod_fwd(end-2,1:150);

%------------------- plot SNIC ---------------------
% plot(Ca_cod_bkgh0,x_cod_bkgh0-.022,'-.','Color',[.9 .7 .1 ],'LineWidth',4)
plot(Ca_cod_bkgh0,x_cod_bkgh0-.022,'-.','Color',[0 0.5 0],'LineWidth',4)
hold on
% plot(Ca_cod_fdgh0,x_cod_fdgh0-0.022,'-.','Color',[0.95 .7 .1 ],'LineWidth',4)
plot(Ca_cod_fdgh0,x_cod_fdgh0-0.022,'-.','Color',[0 0.5 0],'LineWidth',4)
hold on


%X-NULLCLINE
% %gh=0.0
 load('ca_sweep_ih_0.mat')
 % low stable branch 
 xvar= x(5,100:620);
 Cavar=x(6,100:620);
 plot (Cavar,xvar,'Color','black','LineWidth',3)
 hold on
 % unstable branch
 xvar= x(5,600:735);
 Cavar=x(6,600:735);
 plot (Cavar,xvar,':','Color','black','LineWidth',2)
 hold on
 %stable branch 
 xvar= x(5,735:764);
 Cavar=x(6,735:764);
 plot (Cavar,xvar,'Color','black','LineWidth',3)
 hold on
 % upper unstable branch
 xvar= x(5,765:800);
 Cavar=x(6,765:800);
 plot (Cavar,xvar,':','Color','black','LineWidth',2)
 hold on
% %gh=0.0 and Iext=-0.015
% load('ca_sweep_ih_zero_iext_n0p15.mat')
% xvar= x(5,1:end);
% Cavar=x(6,1:end);
% plot (Cavar,xvar,'Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on

%h-current I think
% %gh=0.0003
% load('ca_sweep_ih_0p0003a.mat')
% xvar= x(5,1:end);
% Cavar=x(6,1:end);
% plot (Cavar,xvar,'Color',[.5 .1 .2 ],'LineWidth',1)
% hold on
% %gh=0.0003
% load('ca_sweep_ih_0p0003b.mat')
% xvar= x(5,1:end);
% Cavar=x(6,1:end);
% plot (Cavar,xvar,'Color',[.5 .1 .2 ],'LineWidth',1)
% hold on

%gh=0.0004 and I_ext=0.0 
% load('ca_sweep_ih_0p0004.mat')
% xvar= x(5,1:end);
% Cavar=x(6,1:end);
%plot (Cavar,xvar,'Color',[.5 .5 .5 ],'LineWidth',1)
% hold on

% %gh=0.0004 and I_ext=0.01 
% load('ca_sweep_ih_0p0004_Iext_0p01.mat')
% xvar= x(5,1:end);
% Cavar=x(6,1:end);
% plot (Cavar,xvar,'Color',[.5 .1 .2 ],'LineWidth',1)
% hold on

% %gsyn=0.002
% load('ca_shift_ih_zero_gsyn_0p002.mat')
% xvar= x(5,1:end);  
% Cavar=x(6,1:end);
% plot (Cavar,xvar,'Color',[.0 .4 .5 ],'LineWidth',1)
% hold on

% %gh=0.0
% load('ca_sweep_ih_0.mat')
% xvar= x(5,1:end);
% Cavar=x(6,1:end);
% plot (Cavar,xvar,'Color','black','LineWidth',2)
% hold on



%Voltage trace
start= 12000;
% load('transient2.mat')
% Cavar=x(6:28:end-1,:);
% xvar =x(5:28:end-1,:);
% plot (Cavar,xvar,'Color',[.1 .2 .9 ],'LineWidth',2)
% hold on
% 
% load('transient1.mat')
% Cavar=x(6:28:end-1,:);
% xvar =x(5:28:end-1,:);
% plot (Cavar,xvar,'Color',[.1 .9 .8 ],'LineWidth',3)
% hold on
% 
%  load('transient.mat')
%  Cavar=x(6:28:end-1,:);
%  xvar =x(5:28:end-1,:);
%  start = 2000;
%  plot (Cavar(start:end),xvar(start:end),'Color','blue','LineWidth',4)
% % hold on

load('bursting.mat')
Cavar=x(6:28:end-1,:);
xvar =x(5:28:end-1,:);
% start=2000;
plot (Cavar(start:end),xvar(start:end),'Color','blue','LineWidth',4)
hold on

%CA-SHIFT
%Ca_shift=+56 bursting
load('x_shift_sweep_Ca_56_xinit_n2p7.mat')
xvar= x(5,1:end);
Cavar=x(6,1:end);
%plot (Cavar,.96*xvar,'--','Color',[.0 .0 .6 ],'LineWidth',2)
hold on
%Ca_shift=+20  blue eq
load('x_sweep_Cashift20.mat')
xvar= x(5,1:end);
Cavar=x(6,1:end);
%plot (Cavar,xvar*.92,'--','Color',[.1 .2 .8  ],'LineWidth',2)
hold on
% Ca_shift=0 
% load('x_sweep_Cashift_zero.mat')
% xvar= x(5,1:end);
% Cavar=x(6,1:end);
%plot (1.3*Cavar,1.3*xvar*.925,'-','Color',[.5 .5 .1 ])
hold on
Ca_shift=-38;
load('x_sweep_Cashift_n38.mat')
xvar= x(5,1:end-78);
Cavar=x(6,1:end-78);
%plot (1.5*Cavar,1.*xvar,'-','Color',[.5 .5 .1 ],'LineWidth',4)
hold on
%Ca_shift=-65 ts
load('xsweep_Cashift_n65.mat')
xvar= x(5,1:end-43);
Cavar=x(6,1:end-43);
% plot (1.5*Cavar,1.5*xvar*.93,'--','Color',[1 .0 0 ],'LineWidth',4)
hold on

%Ca_shift=+20
load('x_sweep_Cashift20.mat')
xvar= x(5,1:end);
Cavar=x(6,1:end);
plot (Cavar,xvar*.92,'--','Color',[1 0 0 ], 'LineWidth', 4)
hold on

xlabel('Ca-variable','FontSize',20);
ylabel('x-variable','FontSize',20);
set(gca, 'FontSize', 40); % set y-axis labels


axis([0.65 1.1 0.1 .95] )

load('bursting.mat')
figure(2)
plot(1:length(x(1,start:end)), x(1,start:end),'Color','blue','LineWidth', 3)
% start1=150;
% plot(1:length(x(1,start1:end)), x(1,start1:end),'Color','blue','LineWidth', 3)