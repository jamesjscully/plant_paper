clear all 
clf 
colorbackground = [1 1 1]*0.8;
alphabackground = 0.10;

colornemo1 = [220 220 0]/255;
% colornemo2 = [255 128 0]/255;
colornemo2 = [255 150 0]/255;
alphafish = 0.5;



figure(11)
load('lc_branch_for_Ca_x_plane.mat')
%-------------- 7 variables, last one is s ----------------
period=x(end-1,:)/1000;
%size(period)
Vvar=x (1:28:end-1,:);
Cavar=x(6:28:end-1,:);
xvar =x(5:28:end-1,:);
nvar =x(2:28:end-1,:);
%tr=real(f(1:end-7,1:end));
%=================
size(Cavar)

% intv = 20;
% st=700;
% kend=st+150;
% plot3(Cavar(:,st:intv:kend), nvar(:,st:intv:kend), Vvar(:,st:intv:kend),'Color',[0.5 0.5 0.5])
% hold on

%  start from AH to SNIC
st=850;
intv = 70;
plot3(Cavar(:,st:intv:end-750), nvar(:,st:intv:end-750), Vvar(:,st:intv:end-750),'Color',[0.5 0.5 0.5])
hold on


xlabel('Ca')
ylabel('x')
zlabel('V')'

intv = 2;
st=850;
%  coorV= x(1:kk:end-1,1:10:end-10);
%  coorCa=x(6:kk:end-1,1:10:end-10);
%  coorx =x(5:kk:end-1,1:10:end-10);
%  coorn =x(2:kk:end-1,1:10:end-10);
 coorV= Vvar (:,700:intv:end-1000);
 coorCa=Cavar(:,700:intv:end-1000);
 coorx =xvar (:,700:intv:end-1000);
 coorn =nvar (:,700:intv:end-1000);
zstep = 30;     
strip1 = floor(size(coorV,2)/zstep)
%% number of strip
colo = colornemo2;

for ii = 1:3
    Castrip = coorCa(:,(ii-1)*zstep+1:ii*zstep+1);
    nstrip = coorn(:,(ii-1)*zstep+1:ii*zstep+1);
    xstrip = coorx(:,(ii-1)*zstep+1:ii*zstep+1);
    Vstrip = coorV(:,(ii-1)*zstep+1:ii*zstep+1);
    surf0 = surf(Castrip,nstrip,Vstrip);
    set(surf0,'EdgeColor','n','FaceColor',[.6 .3 0],'FaceAlpha',1);
    hold on
end

zstep = 50; 
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

view (29,5)

% clear all 
% 
% load('ca_sweep_ih_0.mat')
%  % low stable branch 
%  Vvar= x(1,100:620);
%  xvar= x(5,100:620);
%  Cavar=x(6,100:620);
%  plot3 (Cavar,xvar,Vvar,'Color','black','LineWidth',3)
%  hold on
%  % unstable branch
%  Vvar= x(1,600:735);
%  xvar= x(5,600:735);
%  Cavar=x(6,600:735);
%  plot3 (Cavar,xvar, Vvar,':','Color','black','LineWidth',2)
%  hold on
%  %stable branch 
%  Vvar= x(1,735:764);
%  xvar= x(5,735:764);
%  Cavar=x(6,735:764);
%  plot3(Cavar,xvar,Vvar,'Color','black','LineWidth',3)
%  hold on
%  % upper unstable branch
%  Vvar= x(1,765:800);
%  xvar= x(5,765:800);
%  Cavar=x(6,765:800);
%  plot3 (Cavar,xvar,Vvar,':','Color','black','LineWidth',2)
%  hold on