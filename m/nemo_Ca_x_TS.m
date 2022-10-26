clear all 

figure(11)
set(gcf,'position',[150,5,800,800])
clf
colorbackground = [1 1 1]*0.8;
alphabackground = 0.10;
colornemo1 = [220 220 0]/255;
% colornemo2 = [255 128 0]/255;
colornemo2 = [255 150 0]/255;
alphafish = 0.3;

load('LC.mat')
%-------------- 6 variables, last one is s ----------------
period=x(end-1,:)/1000;
%size(period)
Vvar=x (1:6*4:end-1,:);
Cavar=x(6:6*4:end-1,:);
xvar =x(5:24:end-1,:);
nvar =x(2:24:end-1,:);
tr=real(f(1:end-6,1:end));
%=================

imin=length(period);
for i=1:imin
tt=tr(:,i)*period(i);
xvari= xvar(:,i);
Cavari= Cavar(:,i);
xvar_av(i)=trapz(tt,xvari)/period(i);
Cavar_av(i)=trapz(tt,Cavari)/period(i);
end

plot3(Cavar_av*1.01,0*xvar_av,xvar_av+0.015,'Color',[.0 .0 .0 ],'LineWidth',2)
hold on

intv = 20;
 st=1;
% kend=st+150;
%plot3(Cavar(:,st:intv:end), Vvar(:,st:intv:end), xvar(:,st:intv:end),'Color',[0.5 0.5 0.5])
hold on


%  
 intv = 2;
 st=850;

 coorV= Vvar   (:,1:intv:end);
 coorCa=Cavar(:,1:intv:end);
 coorx =xvar   (:,1:intv:end);
 coorn =nvar   (:,1:intv:end);
zstep =3;     
strip1 = floor(size(coorV,2)/zstep);
%% number of strip
colo = colornemo2;

for ii = 1:6
    Castrip = coorCa(:,(ii-1)*zstep+1:ii*zstep+1);
    nstrip = coorn(:,(ii-1)*zstep+1:ii*zstep+1);
    xstrip = coorx(:,(ii-1)*zstep+1:ii*zstep+1);
    Vstrip = coorV(:,(ii-1)*zstep+1:ii*zstep+1);
    surf0 = surf(Castrip,Vstrip,xstrip);
    set(surf0,'EdgeColor','n','FaceColor',[.8 .3 0],'FaceAlpha',1);
    hold on
end
% 
zstep = 5; 
strip1 = floor(size(coorV,2)/zstep);
for ii = 1:strip1-1
    Castrip = coorCa(:,(ii-1)*zstep+1:ii*zstep+1);
    nstrip = coorn(:,(ii-1)*zstep+1:ii*zstep+1);
    xstrip = coorx(:,(ii-1)*zstep+1:ii*zstep+1);
    Vstrip = coorV(:,(ii-1)*zstep+1:ii*zstep+1);
    surf0 = surf(Castrip,Vstrip,xstrip);
    set(surf0,'EdgeColor','none','FaceColor',colo,'FaceAlpha',alphafish);
    if (colo == colornemo2)
        colo = colornemo1;
    else
        colo = colornemo2;
    end;
    hold on
end
hold on
ii = strip1;
Castrip = coorCa(:,(ii-1)*zstep+1:end);
nstrip = coorn(:,(ii-1)*zstep+1:end);
xstrip = coorx(:,(ii-1)*zstep+1:end);
Vstrip = coorV(:,(ii-1)*zstep+1:end);
sstrip = surf(Castrip,Vstrip,xstrip);
set(sstrip,'EdgeColor','none','FaceColor',colo,'FaceAlpha',alphafish);
    if (colo == colornemo2)
        colo = colornemo1;
    else
        colo = colornemo2;
    end;
hold on

view (29,5)
view (0,0)
% 
clear all 
% % 
load('eq1.mat')
% %  % low stable branch 
Vvar= x(1,1:end);
 nvar=x(2,1:end); 
 xvar= x(5,1:end);
 Cavar=x(6,1:end);
 end1=1200;
 plot3 (Cavar(1:end1),Vvar(1:end1),xvar(1:end1),'Color','black','LineWidth',2)
 hold on
 end2=1460;
 plot3 (Cavar(end1:end2),Vvar(end1:end2),xvar(end1:end2),':','Color','black','LineWidth',2)
 hold on
 end3=1582;
 plot3 (Cavar(end2:end3),Vvar(end2:end3),xvar(end2:end3),'Color','black','LineWidth',2)
 hold on
 end4=1700;
 plot3 (Cavar(end3:end4),Vvar(end3:end4),xvar(end3:end4),':','Color','black','LineWidth',2)
 hold on

 clear all 
% % 
load('eq2.mat')
% %  % low stable branch 
Vvar= x(1,1:end);
 nvar= x(2,1:end); 
 xvar= x(5,1:end);
 Cavar=x(6,1:end);
 end1=30;
 plot3 (Cavar(1:end1),Vvar(1:end1),xvar(1:end1),'Color','black','LineWidth',2)
 hold on
 end2=270;
 plot3 (Cavar(end1:end2),Vvar(end1:end2),xvar(end1:end2),':','Color','black','LineWidth',2)
 hold on

 
clear all 
% % 
load('EQgsyn0p01.mat')
% %  % low stable branch 
Vvar= x(1,1:end);
 nvar= x(2,1:end); 
 xvar= x(5,1:end);
 Cavar=x(6,1:end);
 end1=306;
 plot3 (Cavar(1:end1),Vvar(1:end1),xvar(1:end1),'Color',[0.6 0.6 0.6],'LineWidth',1.2)
 hold on
 end2=399;
  plot3 (Cavar(end1:end2),Vvar(end1:end2),xvar(end1:end2),'.','Color',[0.5 0.5 0.5],'LineWidth',1.9)
 hold on
 
 clear all 
% % 
load('EQgsyn0p02.mat')
% %  % low stable branch 
Vvar= x(1,1:end);
 nvar= x(2,1:end); 
 xvar= x(5,1:end);
 Cavar=x(6,1:end);
 end1=460;
 plot3 (Cavar(1:end1),Vvar(1:end1),xvar(1:end1),'Color',[0.6 0.6 0.6],'LineWidth',1.9)
 hold on
 
%  clear all 
% load('quescence1.mat')
% Vvar=  x(1,260:end);
%  nvar= x(2,260:end); 
%  xvar= x(5,260:end);
%  Cavar=x(6,260:end);
% %plot3(Cavar,Vvar,xvar,'Color',[.0 .1 .9],'LineWidth',2)
% hold on

%  clear all 
% load('TS1.mat')
% Vvar=  x(1,260:end);
%  nvar= x(2,260:end); 
%  xvar= x(5,260:end);
%  Cavar=x(6,260:end);
% %plot3(Cavar,Vvar,xvar,'Color',[.0 .2 .8],'LineWidth',2)
% hold on


clear all 
load('caprime_n35.mat')
% %  % low stable branch 
Vvar= x(1,1:end);
 nvar= x(2,1:end); 
 xvar= x(5,1:end);
 Cavar=x(6,1:end);
 end1=460;
plot3 (Cavar(1:end1)*1.45,Vvar(1:end1),1.7*xvar(1:end1)-0.083,'--','Color',[0.99 0.5 0.9],'LineWidth',2)
 hold on




%SNIC part 
%--------------------- SNIC whole system gh = 0 -------------------------------
data1 = load('LP_whole_bk_gh0.mat','x');
cod_back = data1.x;
data1 = load('LP_whole_fwd_gh0.mat','x');
cod_fwd = data1.x;
x_cod_bkgh0 = cod_back(end-3,:)-0.04;
Ca_cod_bkgh0 = cod_back(end-2,:);
x_cod_fdgh0 = cod_fwd(end-3,:)-0.04;
Ca_cod_fdgh0 = cod_fwd(end-2,:);
%------------------- plot SNIC ---------------------
plot3(Ca_cod_bkgh0, 0*Ca_cod_bkgh0,x_cod_bkgh0-0.007,'-.','Color',[.85 .85 .85 ],'LineWidth',4)
hold on
plot3(Ca_cod_fdgh0, 0*Ca_cod_fdgh0,x_cod_fdgh0-0.007,'-.','Color',[.85 .85 .85 ],'LineWidth',4)
hold on

axis([0.45 1.42 -80 40 0.2 .97])
axis on 
view(0,0)
%axis tight 


clear all 
load('TS_Q_after_pulse.mat')
st=  30150;
end1=47000;
plot3(Caa2(st:end1),vv2(st:end1),xx2(st:end1),'Color',[.9 .1 .1],'LineWidth',2)
hold on

st=  47000;
plot3(Caa2(st:end),vv2(st:end),xx2(st:end),'Color',[.0 .2 .8],'LineWidth',2)
hold on

x3=Caa2(st);
y3=vv2(st);
z3=xx2(st);
[X,Y,Z]=sphere(20);
mag0=0.012;

s3=surf (mag0*X+x3, mag0*Y+y3, mag0*Z/1.6+z3,'EdgeColor','none','FaceColor',...
    [200./255  51./255  55./255],'FaceAlpha',0.99) ;
hold on



%camlight headlight
%h=camlight;
%set(h,'position',[10 1 1],'style','local'); 
material shiny;




fftt = 16;
txpo = text(0.64,0.5,0.925,'M_{po}','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(1.01,0.5,0.76,'TS','Fontsize',18,'Color','blue','FontName','Arial')
txpo = text(0.54,0.5,0.87,'\langle x/V \rangle ','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(1.31,0.5,0.87,'SNIC','Fontsize',fftt,'Color',[0.6 0.6 0.6],'FontName','Arial')
txpo = text(0.55,0.5,0.37,'[Ca]^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')

txpo = text(1.3,0.5,0.73,'x^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(1.12,0.5,0.34,'x_{syn1}^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(0.9,0.5,0.305,'x_{syn2}^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')



xlabel('[Ca]-variable','FontSize',18);
zlabel('x-variable','FontSize',18);

%print(gcf,'-djpeg','-r600' ,'nemoCaX_TS_after_pulse.jpeg');


figure(10)
load('TS_Q_after_pulse1.mat')

subplot(2,1,1)
 plot(time-25, vv2,'Color',[0 0.2 0.8],'LineWidth',2);
 hold on
 st=30120;
end1=35040;
 plot(time(st:end1)-25, vv2(st:end1),'Color',[0.9 0.1 0.1],'LineWidth',2);
 hold on
 xlim([25-25 55-25]) 
 ylim([-65 45])
 
xlabel('time [sec]','FontSize',18);
ylabel('Voltage [mV]','FontSize',18);

alpha = 1;
beta =  0.035;

s2=zeros(length(vv2),1);
s=0.0001;
for i=2:1:length(vv1)
s=s+0.1*(alpha*s*(1-s)/(1+exp(-10*(vv2(i)+20)))-beta*s);
s2(i)=s;
end

subplot(2,1,2)
plot(time-25, s2,'Color',[0.1 0.2 0.8],'LineWidth',2);
 hold on

 
 plot(time-25,smoothdata(s2),'Color',[0.9 0.9 0.99],'LineWidth',2);
 hold on
 % plot(time-25,smoothdata(smoothdata(s2)),'Color',[0.9 0.2 0.1],'LineWidth',2);
 hold on
 xlim([25-25 55-25]) 
 ylim([-0.01 1.1])
 
xlabel('time [sec]','FontSize',18);
ylabel('rate s(t)','FontSize',18);




%print(gcf,'-djpeg','-r600' ,'V_TS_after_pulse.jpeg');