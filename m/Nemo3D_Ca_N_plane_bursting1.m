
 
%clear all 
figure (24)
clf

colorbackground = [1 1 1]*0.8;
alphabackground = 0.10;
colornemo1 = [220 220 0]/255;
% colornemo2 = [255 128 0]/255;
colornemo2 = [255 150 0]/255;
alphafish = 0.3;

% PO branch 
load('lc_branch_for_Ca_x_plane.mat')
period=x(end-1,:)/1000;
Vvar= x(1:28:end-1,:);
%nvar= x(2:28:end-1,:);
Cavar=x(6:28:end-1,:);
xvar =x(5:28:end-1,:);
nvar =x(2:28:end-1,:);
tr=real(f(1:end-7,1:end));

imin=length(period);
for i=1:imin
tt=tr(:,i)*period(i);
xvari= xvar(:,i);
nvari= nvar(:,i);
Cavari= Cavar(:,i);
Vvari=  Vvar(:,i);
xvar_av(i)=trapz(tt,xvari)/period(i);
nvar_av(i)=trapz(tt,nvari)/period(i);
Cavar_av(i)=trapz(tt,Cavari)/period(i);
Vvar_av(i)=trapz(tt,Vvari)/period(i);
end

plot3(Cavar_av+0.0065,nvar_av+0.015,Vvar_av,'Color',[.3 .3 .3 ],'LineWidth',2)
hold on

%st=  30300;
%plot3(Cavar(st:end) ,xvar(st:end),Vvar(st:end),'Color',[.3 .3 .3 ],'LineWidth',0.1)
hold on

 st=  3000;
 intv = 1;

 coorV= Vvar   (:,1:intv:end);
 coorCa=Cavar  (:,1:intv:end);
 coorx =xvar   (:,1:intv:end);
 coorn =nvar   (:,1:intv:end);
zstep =3;
 size(coorV)
 
zstep =10;     
strip1 = floor(size(coorV,2)/zstep)
%% number of strip
colo = colornemo2;

for ii = 72:85
    Castrip = coorCa(:,(ii-1)*zstep+1:ii*zstep+1);
    nstrip = coorn(:,(ii-1)*zstep+1:ii*zstep+1);
    xstrip = coorx(:,(ii-1)*zstep+1:ii*zstep+1);
    Vstrip = coorV(:,(ii-1)*zstep+1:ii*zstep+1);
    surf0 = surf(Castrip,nstrip,Vstrip);
    set(surf0,'EdgeColor','n','FaceColor',[.9 .6 0],'FaceAlpha',.99);
    hold on
end
% 
% 
zstep = 90; 
strip1 = floor(size(coorV,2)/zstep);
%% number of strip
colo = colornemo2;

for ii = 10:strip1-10
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

% 
% 
% % %gh=0.0
  load('ca_sweep_ih_0.mat')
  % low stable branch 
  xvar= x(5,100:620);
  nvar= x(2,100:620);
  Cavar=x(6,100:620);
  Vvar=x(1,100:620);
  plot3 (Cavar,nvar,Vvar,'Color','black','LineWidth',2)
  hold on
  % unstable branch
  xvar= x(5,600:735);
  nvar= x(2,600:735);
  Cavar=x(6,600:735);
  Vvar= x(1,600:735);
  plot3(Cavar,nvar,Vvar,':','Color','black','LineWidth',2)
  hold on
%  %stable branch 
  xvar= x(5,735:764);
  nvar= x(2,735:764);
  Cavar=x(6,735:764);
  Vvar= x(1,735:764);
  plot3 (Cavar,nvar,Vvar,'Color','black','LineWidth',2)
  hold on
%  % upper unstable branch
  xvar= x(5,765:800);
  nvar= x(2,765:800);
  Cavar=x(6,765:800);
  Vvar=x(1, 765:800);
  plot3 (Cavar,nvar,Vvar,':','Color','black','LineWidth',2)
  hold on
  
  load('eq2.mat')
% %  % low stable branch 
Vvar= x(1,1:end);
 nvar= x(2,1:end); 
 xvar= x(5,1:end);
 Cavar=x(6,1:end);
 end1=30;
 plot3 (Cavar(1:end1),nvar(1:end1),Vvar(1:end1),'Color',[0.2 0.2 0.2],'LineWidth',2)
 hold on
 end2=270;
 plot3 (Cavar(end1:end2),nvar(end1:end2),Vvar(end1:end2),':','Color',[0.2 0.2 0.2],'LineWidth',2)
 hold on
  

load('bursting.mat')
Cavar=x(6:28:end-1,:);
xvar =x(5:28:end-1,:);
nvar =x(2:28:end-1,:);
Vvar =x(1:28:end-1,:);
plot3 (Cavar(700:end),nvar(700:end),Vvar(700:end),'Color',[.0 .2 .85 ],'LineWidth',2)
hold on

fftt = 16;
txpo = text(0.4,0,20,'M_{po}','Fontsize',fftt,'Color','black','FontName','Arial');
txpo = text(0.4,0,-23,'\langle V \rangle ','Fontsize',fftt,'Color','black','FontName','Arial');
%txpo = text(0.7,0.68,'SNIC','Fontsize',fftt,'Color','black','FontName','Arial')
%txpo = text(0.68,0.39,'[Ca]^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')

txpo = text(1.25,0,-28,'M_{eq}=0','Fontsize',fftt,'Color','black','FontName','Arial');


xlabel('[Ca]-variable','FontSize',16);
ylabel('n-variable','FontSize',16);
zlabel('Voltage [mv]','FontSize',16)

view (25,5)
%view (0,90)

axis([-.2 1.2 0.0 1.5 -72 35] )
%axis tight 

 print(gcf,'-djpeg','-r600' ,'nemo3D_Ca_n1.jpeg');
 
  figure(25)
   clf 
  load('bursting.mat')
Vvar =x(1:28:end-1,:);
time1=1:1:length(Vvar(700:end));

plot (time1*0.0001-1,Vvar(700:end),'Color',[.0 .2 .85 ],'LineWidth',2)
hold on
xlim([0 6])
xlabel('time [sec]','FontSize',16);
ylabel('Voltage [mV]','FontSize',16);

print(gcf,'-djpeg','-r600' ,'burst_time.jpeg');

set(gcf, 'color', 'none');    
set(gca, 'color', 'none'); 
% exportgraphics(gca,'burst_time.png','BackgroundColor','none')
print(gcf,'-dpng','-r600' ,'burst_time.png');

exportgraphics(gcf,'transparent.png',...   % since R2020a
    'ContentType','vector',...
    'BackgroundColor','none')


% 
%  clf
%  time=t(1:end)'/10000; 
%  plot(time-1, Vvar,'Color',[0 0.2 0.8],'LineWidth',2);
%  hold on
%  xlim([0 7]) 
%  ylim([-70 40])
% %axis tight  
% xlabel('time [sec]','FontSize',18);
% ylabel('Voltage [mV]','FontSize',18);
% 
% %print(gcf,'-djpeg','-r600' ,'bursting.jpeg');
