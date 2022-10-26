clf
clear all
figure(1)
% load('lc_period_4000_gh_zero.mat')
% xshift=x(end-1,:);
% Cashift=x(end,:);
% hold on
% plot (Cashift,xshift,'Color',[.5 .5 .5 ],'LineWidth',2)
% hold on
% load('lc_period_4000_gh_0p0001.mat')
% xshift=x(end-1,:);
% Cashift=x(end,:);
% hold on
% plot (Cashift,xshift, 'Color',[.5 .5 .5 ],'LineWidth',2)
% hold on



% period=x(end-1,:);
% for i=1:length(param)
%     tt=t(:,i)*period(i);
% end
% load('period2000b.mat')
% xshift= x(end-1,:);
% Cashift=x(end,:);
% plot (Cashift,xshift,'Color',[.6 .6 .6 ],'LineWidth',2)
% hold on
% load('period2000a.mat')
% xshift=x(end-1,:);
% Cashift=x(end,:);
% plot (Cashift,xshift,'Color',[.5 .5 .5 ],'LineWidth',2)
% hold on

% figure(2)
% v=x(1:28:end-1,:);
% t=real(f(1:end-7,:));
% plot (t,v,'Color','r')
% hold on

%load('hopf_0p003.mat')
%xshift=x(end-2,:);
%Cashift=x(end-1,:);
%plot (Cashift,xshift,'Color',[.1 .1 .9 ],'LineWidth',1)
hold on


%SNIC part 
%--------------------- SNIC whole system gh = 0 -------------------------------
data1 = load('LP_whole_bk_gh0.mat','x');
cod_back = data1.x;
data1 = load('LP_whole_fwd_gh0.mat','x');
cod_fwd = data1.x;
x_cod_bkgh0 = cod_back(end-1,:);
Ca_cod_bkgh0 = cod_back(end,:);
x_cod_fdgh0 = cod_fwd(end-1,:);
Ca_cod_fdgh0 = cod_fwd(end,:);

%------------------- plot SNIC ---------------------
plot(Ca_cod_fdgh0(1:end),x_cod_fdgh0(1:end),'-','Color',[.1 .1 .1 ],'LineWidth',2)
hold on
plot(Ca_cod_bkgh0+1.2,x_cod_bkgh0,'-','Color',[.75 .75 .75 ],'LineWidth',12)
hold on
plot(Ca_cod_fdgh0(1:end-227)+1.2,x_cod_fdgh0(1:end-227),'-','Color',[.7 .7 .7 ],'LineWidth',12)
hold on
plot(Ca_cod_bkgh0,x_cod_bkgh0,'-','Color',[.0 .0 .0 ],'LineWidth',2)
hold on
plot(Ca_cod_fdgh0(1:end-225),x_cod_fdgh0(1:end-225),'-','Color',[.0 .0 .0 ],'LineWidth',2)
hold on


load('hopf_bif_dagram_ih_zero.mat')
xshift= x(end-2,1:end);
Cashift=x(end-1,1:end);

stop1=560;
plot (Cashift(stop1:end-10)+1.2,xshift(stop1:end-10)-.07,'Color',[.75 .75 .75 ],'LineWidth',12)
hold on

stop1=558;
plot (Cashift(1:stop1),xshift(1:stop1)+.06,'Color',[.75 .75 .99 ],'LineWidth',12)
hold on


plot (Cashift,xshift,'Color',[.1 .1 .8 ],'LineWidth',4)
hold on
plot (Cashift(stop1:end),xshift(stop1:end),'.-','MarkerSize',6,'Color',[.9 .9 .9 ],'LineWidth',2)
hold on
plot(Ca_cod_bkgh0-.2,x_cod_bkgh0,'-','Color',[.0 .0 .0 ],'LineWidth',2)
hold on
plot(Ca_cod_fdgh0(1:end)-.2,x_cod_fdgh0(1:end),'-','Color',[.1 .1 .1 ],'LineWidth',2)
hold on

plot(30.007875,-2.77,'.','MarkerSize',80,'Color',[.0 .0 .5 ])
hold on

% -2.779353 

ax = gca;
ax.FontSize = 16; 

xlabel('\Delta_{Ca}','FontSize',26);
ylabel('\Delta_{x}','FontSize',26);
axis([-80 370 -8 8])
axis([-95 50 -5 2])
box on
%axis tight 
