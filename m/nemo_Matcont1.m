clear all 

figure(21)
clf 
colorbackground = [1 1 1]*0.8;
alphabackground = 0.10;
colornemo1 = [220 220 0]/255;
% colornemo2 = [255 128 0]/255;
colornemo2 = [255 150 0]/255;
alphafish = 0.5;

load('LC.mat')
%-------------- 6 variables, last one is s ----------------
period=x(end-1,:)/1000;
%size(period)
Vvar=x (1:6*4:end-1,:);
Cavar=x(6:6*4:end-1,:);
xvar =x(5:24:end-1,:);
nvar =x(2:24:end-1,:);
tr=real(f(1:end-7,1:end));

%=================

intv = 20;
 st=1;
% kend=st+150;
plot3(Cavar(:,st:intv:end), nvar(:,st:intv:end), Vvar(:,st:intv:end),'Color',[0.5 0.5 0.5])
hold on

%  
 intv = 2;
 st=850;

 coorV= Vvar   (:,1:intv:end);
 coorCa=Cavar  (:,1:intv:end);
 coorx =xvar   (:,1:intv:end);
 coorn =nvar   (:,1:intv:end);
zstep =3;     
strip1 = floor(size(coorV,2)/zstep)
%% number of strip
colo = colornemo2;

for ii = 1:6
    Castrip = coorCa(:,(ii-1)*zstep+1:ii*zstep+1);
    nstrip = coorn(:,(ii-1)*zstep+1:ii*zstep+1);
    xstrip = coorx(:,(ii-1)*zstep+1:ii*zstep+1);
    Vstrip = coorV(:,(ii-1)*zstep+1:ii*zstep+1);
    surf0 = surf(Castrip,nstrip,Vstrip);
    set(surf0,'EdgeColor','n','FaceColor',[.8 .3 0],'FaceAlpha',1);
    hold on
end
% 
zstep = 30; 
strip1 = floor(size(coorV,2)/zstep)
for ii = 1:strip1-1
    Castrip = coorCa(:,(ii-1)*zstep+1:ii*zstep+1);
    nstrip = coorn(:,(ii-1)*zstep+1:ii*zstep+1);
    xstrip = coorx(:,(ii-1)*zstep+1:ii*zstep+1);
    Vstrip = coorV(:,(ii-1)*zstep+1:ii*zstep+1);
    surf0 = surf(Castrip,nstrip,Vstrip);
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
sstrip = surf(Castrip,nstrip,Vstrip);
set(sstrip,'EdgeColor','none','FaceColor',colo,'FaceAlpha',alphafish);
    if (colo == colornemo2)
        colo = colornemo1;
    else
        colo = colornemo2;
    end;
hold on

view (12,3)
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
 plot3 (Cavar(1:end1),nvar(1:end1),Vvar(1:end1),'Color',[0.8 0.8 0.8],'LineWidth',2)
 hold on
 end2=1460;
 plot3 (Cavar(end1:end2),nvar(end1:end2),Vvar(end1:end2),':','Color',[0.8 0.8 0.8],'LineWidth',2)
 hold on
 end3=1589;
 plot3 (Cavar(end2:end3),nvar(end2:end3),Vvar(end2:end3),'Color',[0.8 0.8 0.8],'LineWidth',2)
 hold on
 end4=1700;
 plot3 (Cavar(end3:end4),nvar(end3:end4),Vvar(end3:end4),':','Color',[0.8 0.8 0.8],'LineWidth',2)
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
 plot3 (Cavar(1:end1),nvar(1:end1),Vvar(1:end1),'Color',[0.8 0.8 0.8],'LineWidth',2)
 hold on
 end2=270;
 plot3 (Cavar(end1:end2),nvar(end1:end2),Vvar(end1:end2),':','Color',[0.8 0.8 0.8],'LineWidth',2)
 hold on

 clear all 
 load('TS1.mat')
 Vvar= x(1, 1000:end);
 nvar= x(2,1000:end); 
 xvar= x(5,1000:end);
 Cavar=x(6,1000:end);
 time=t(1000:end)'; 
 size(time)
 size (Vvar)
 
plot3(Cavar,nvar,Vvar,'Color',[.1 .2 .9],'LineWidth',2)
hold on


fftt = 15;
txpo = text(0.4,0.5,-58,'M_{po}','Fontsize',fftt,'Color','black','FontName','Arial')
%txpo = text(0.7,0.5, 0.85,'\langle x/V \rangle ','Fontsize',fftt,'Color','black','FontName','Arial')
%txpo = text(0.7,0.5, 0.68,'SNIC','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(1.5,0.7,-48,'M_{eq}','Fontsize',fftt,'Color','black','FontName','Arial')


xlabel('[Ca]-variable','FontSize',16)
ylabel('n','FontSize',16)
zlabel('Voltage [mV]','FontSize',16)
axis([-.2 1.9 0 1 -70 30])

print(gcf,'-djpeg','-r600' ,'nemo2.jpeg');
% 
% figure(3)
%  clf
%  plot(time,Vvar,'Color',[.1 .2 .9 ],'LineWidth',1.5)
% xlabel('time [mS]','FontSize',16);
% ylabel('Voltage [mV]','FontSize',16);
% print(gcf,'-djpeg','-r600' ,'burst.jpeg');
