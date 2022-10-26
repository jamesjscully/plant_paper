
 
clear all 
figure (27)
clf

colorbackground = [1 1 1]*0.8;
alphabackground = 0.10;
colornemo1 = [220 220 0]/255;
% colornemo2 = [255 128 0]/255;
colornemo2 = [255 150 0]/255;
alphafish = 0.5;

% PO branch 
load('lc_branch_for_Ca_x_plane.mat')
period=x(end-1,:)/1000;

Vvar= x(1:28:end-1,:);
%nvar= x(2:28:end-1,:);
Cavar=x(6:28:end-1,:);
xvar =x(5:28:end-1,:);
tr=real(f(1:end-7,1:end));

imin=length(period);
for i=1:imin
tt=tr(:,i)*period(i);
xvari= xvar(:,i);
Cavari= Cavar(:,i);
Vvari=  Vvar(:,i);
xvar_av(i)=trapz(tt,xvari)/period(i);
Cavar_av(i)=trapz(tt,Cavari)/period(i);
Vavar_av(i)=trapz(tt,Vvari)/period(i);
end

hold on
% fake it 
plot (Cavar_av+0.0065,xvar_av+0.015,'Color',[.3 .3 .3 ],'LineWidth',3)
hold on


[x y] = meshgrid(-.5:0.1:1.5); % Generate x and y data
z = -.2+zeros(size(x, 1)); % Generate z data
splane=surf(x, y, z);
set(splane,'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',.3);
hold on


intv = 10;
st=800;
 coorV= Vvar (:,st:1:end);
 coorCa=Cavar(:,st:1:end);
 coorx =xvar (:,st:1:end);
 

 zstep = 60; 
strip1 = floor(size(coorV,2)/zstep);
%% number of strip
colo = colornemo2;

for ii = 1:strip1-1
    Castrip = coorCa(:,(ii-1)*zstep+1:ii*zstep+1);
    %nstrip = coorn(:,(ii-1)*zstep+1:ii*zstep+1);
    xstrip = coorx(:,(ii-1)*zstep+1:ii*zstep+1);
    Vstrip = coorV(:,(ii-1)*zstep+1:ii*zstep+1);
    surf0 = surf(Castrip,xstrip,0*Vstrip);
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
%nstrip = coorn(:,(ii-1)*zstep+1:end);
xstrip = coorx(:,(ii-1)*zstep+1:end);
Vstrip = coorV(:,(ii-1)*zstep+1:end);
sstrip = surf(Castrip,xstrip,0*Vstrip);
set(sstrip,'EdgeColor','none','FaceColor',colo,'FaceAlpha',alphafish);
    if (colo == colornemo2)
        colo = colornemo1;
    else
        colo = colornemo2;
    end;
hold on

view (29,5)
view (20,90)




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
plot(Ca_cod_bkgh0, x_cod_bkgh0,'-.','Color',[.8 .8 .8 ],'LineWidth',4)
hold on
plot(Ca_cod_fdgh0, x_cod_fdgh0,'-.','Color',[.8 .8 .8 ],'LineWidth',4)
hold on


% %gh=0.0
 load('ca_sweep_ih_0.mat')
 % low stable branch 
 xvar= x(5,100:620);
 Cavar=x(6,100:620);
 plot (Cavar,xvar,'Color','black','LineWidth',2)
 hold on
 % unstable branch
 xvar= x(5,600:735);
 Cavar=x(6,600:735);
 plot (Cavar,xvar,':','Color','black','LineWidth',2)
 hold on
 %stable branch 
 xvar= x(5,735:764);
 Cavar=x(6,735:764);
 plot (Cavar,xvar,'Color','black','LineWidth',2)
 hold on
 % upper unstable branch
 xvar= x(5,765:800);
 Cavar=x(6,765:800);
 plot (Cavar,xvar,':','Color','black','LineWidth',2)
 hold on



%Ca_shift=+56 bursting
load('x_shift_sweep_Ca_56_xinit_n2p7.mat')
xvar= x(5,1:end);
Cavar=x(6,1:end);
plot (Cavar,xvar,'--','Color',[.9 .0 .6 ],'LineWidth',2)
hold on

load('bursting.mat')
Cavar=x(6:28:end-1,:);
xvar =x(5:28:end-1,:);
plot (Cavar,xvar,'Color',[.0 .2 .85 ],'LineWidth',2)
hold on

fftt = 15;
txpo = text(0.77,0.925,'M_{po}','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(0.7,0.85,'\langle x/V \rangle ','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(0.7,0.65,'SNIC','Fontsize',fftt,'Color','black','FontName','Arial')
txpo = text(0.68,0.39,'[Ca]^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')

txpo = text(1.04,0.63,'x^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')


xlabel('[Ca]-variable','FontSize',16);
ylabel('x-variable','FontSize',16);


axis([0.65 1.12 0.1 .95] )
view(0,90)

daspect([1 1 20])
view(20,30)
camproj perspective
axis on

%axis off

%print(gcf,'-djpeg','-r600' ,'bursting_xCa_skewed1.jpeg');
