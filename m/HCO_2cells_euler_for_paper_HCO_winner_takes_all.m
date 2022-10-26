
clear all 

tmax=50; % max time in seconds


%--------- green, blue: ts ---------------
alpha1 = 0.021
beta1 =  0.0013;
scale1=alpha1/(alpha1+beta1)

alpha2 = 0.03;
beta2 =  0.0011;
scale2=alpha2/(alpha2+beta2)

g12=.1;
g21=1.;


gelec=0.0002;


%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.000;
%----------------------

%Parameters3
Ca_shift1 = -50;
Ca_shift2 = -40;

x_shift =  -4;

Iapp=0.0;

t1=00*1000; % this time is in msec
t2=00*1000; % this time is in msec
      
% H-current      
    gh    = 0.0004;
    Vhh   = -55;

% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.1;       % time step in msec
tsamp=2; % sampling interval to save data in msec



cut1=cutoff/1000; % cutoff frequency converted to msec
% Initial values
V1= -44; V2= -44; Ca1=.2; Ca2=1.3; h1 =0; h2 =0; n1 =0; n2 =0; 
x1 =0.8; x2 =.82; y1 =0; y2 =0; s1 =0; s2 =0; m1 =0; m2 =0;


tic
[rnd1,nt] = bandlimnoise (cut1,step,tf);
[rnd2,nt] = bandlimnoise (cut1,step,tf);
rnd1=rnd1*noise; rnd2=rnd2*noise;

nt1=round(tmax/tsamp)+1;
is=round(tsamp/step);

time=zeros(nt1,1);vv1=time;vv2=time;ss1=time;ss2=time;mm1=time;mm2=time;Caa1=time;
Caa2=time;xx1=time;xx2=time;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CALCULATE heaviside functions in advance! This speeds up calculation!

heav1=zeros(nt,1); heav2=heav1;
j1=round(t1/step)+1; j2=round(t2/step)+1;
heav1(j1:nt)=1; heav2(j2:nt)=1; heav=heav1-heav2; 
clear heav1 heav2;

j=0;
for i=1:nt
 tt=(i-1)*step;
V1 =V1 +step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
    +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1) +gh*((1/(1+exp(-(V1+63)/7.8)))^3)*y1*(+120-V1)...
-g21*(V1+80)*s2/scale1+gelec*(V2-V1)+rnd1(i));
V2= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
    +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) +gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)...
-g12*(V2-40)*s1/scale1+gelec*(V1-V2)+rnd2(i));
Ca1=Ca1+step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
Ca2=Ca2+step*(0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/235);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/235);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
y1 =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
y2 =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));
s1 =s1+step*(alpha1*s1*(1-s1)/(1+exp(-20*(V1+20)))-beta1*(s1-0.00001));
s2 =s2+step*(alpha2*s2*(1-s2)/(1+exp(-20*(V2+20)))-beta2*(s2-0.00001));
% m1 =m1+step*((1/(1+exp(-(V1+20)))-m1)/taum);
%m2 =m2+step*((1/(1+exp(-(V2+20)))-m2)/taum);


% Do not need to save every point!
if mod(i,is) ==0
    j=j+1;
    time(j)=tt; vv1(j)=V1; vv2(j)=V2; ss1(j)=s1;
    ss2(j)=s2; mm1(j)=m1; mm2(j)=m2; Caa1(j)=Ca1;
    Caa2(j)=Ca2; xx1(j)=x1; xx2(j)=x2; heavy(j)=heav(i);
end
end  
toc    

%  
%  [vp1,tp1] = findpeaks(vv1,time,'MinPeakHeight',-20);
%  fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%  
%  % convert back to sec
%  tp1=tp1/1000;
 time=time/1000;
%%
 %     
 figure(4)
 clf
 subplot(4,1,1)

 plot(time,vv1,'Color',[0 0  ,.6]','LineWidth',1.5)
 hold on
%  plot(time,heavy*10,'Color',[1 0  1]','LineWidth',1.5)
%  hold on
 xlim([0 tmax]) 
 ylim([-85 40])
  fftt = 28;
txpo = text(-4.2,25,'B','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')
   
   xlabel('Time [sec]','Fontsize', 14),ylabel('Voltage','Fontsize', 16)

 subplot(4,1,2)
 plot(time,vv2,'Color',[0 .5 0]','LineWidth',1.5)
 hold on
  xlabel('Time [sec]','Fontsize', 14),ylabel('Voltage','Fontsize', 16)
 xlim([0 tmax]) 
 ylim([-65 40])
 
  %xlim([350 time(end)])

  subplot(4,1,3)
 % plot(time,10*mm1.*ss1,'blue')
  hold on
  plot(time,ss1/scale1,'Color',[0.1 0.1 0.6],'LineWidth',1.5)
  hold on
  plot(time,ss2/scale2,'Color',[0.1 0.6 0.1],'LineWidth',1.5)
  hold on
  ylim([-0.02  .1])
  xlim([0 tmax])
  box on
   xlabel('Time [sec]','Fontsize', 16),ylabel('Release rate','Fontsize', 16)
  
subplot(4,1,4)
  plot(time,Caa1,'Color',[0. 0.0 0.8],'Linewidth',1.5)
  hold on
  plot(time,Caa2,'Color',[0.0 0.7 0.0],'Linewidth',1.5)
  hold on
  %plot(time,xx1,'Color',[0 .1 1])
  hold on
  %plot(time,xx2,'Color',[0 1 .1])
  hold on
  ylim([0.0 1.4])
  xlim([0 tmax])
 xlabel('Time [sec]','Fontsize', 14),ylabel('[Ca]','Fontsize', 16)
  %title('Ca and Ca-activaterd possasion gating varaible in time', 'Fontsize', 11);
  
   print(gcf,'-djpeg','-r600' ,'HCO_winner1_B.jpeg')

  toc 
 
  figure(5)
 clf
 st=1;
 plot(Caa1(st:end),xx1(st:end),'Color',[0. 0. 0.6],'LineWidth',1.5)
  hold on
  plot(Caa2(st:end),xx2(st:end),'Color',[0.1 0.6 0.1],'LineWidth',1.5)
  hold on
  %title('Ca vs x', 'Fontsize', 11);
  

 xlabel('[Ca]-variable','Fontsize', 16),
  ylabel('x-variable','Fontsize', 16)  
  %title('Ca vs Ca-activaterd possasion gating varaible', 'Fontsize', 11);
  
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
plot(Ca_cod_bkgh0,x_cod_bkgh0,'-.','Color',[.85 .85 .85 ],'LineWidth',4)
hold on
plot(Ca_cod_fdgh0,x_cod_fdgh0,'-.','Color',[.85 .85 .85 ],'LineWidth',4)
hold on

axis([0.01 1.4 0.01 1.0])

fftt = 28;
txpo = text(-.1,.98,'А','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')

fftt = 15;
txpo = text(0.18,.1,'[Ca]^\prime=0','Fontsize',fftt,'Color','red','FontName','Arial')

fftt = 15;
%txpo = text(1.2,.51,'[Ca]^\prime=0','Fontsize',fftt,'Color','red','FontName','Arial')


fftt = 16;
txpo = text(1.26,.15,'x^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')

fftt = 16;
txpo = text(.14,.38,'SNIC','Fontsize',fftt,'Color',[0.6 0.6 0.6],'FontName','Arial')




pos=length(Caa1);
  plot(Caa1(pos),xx1(pos),'.','Color',[0. 0.0 1],'MarkerSize',50)
  hold on
  plot(Caa2(pos),xx2(pos),'.','Color',[0. 0.99 0.],'MarkerSize',50)
  hold on
plot(Caa1(pos),xx1(pos),'.','Color',[0.5 0.7 1],'MarkerSize',40)
  hold on
  plot(Caa2(pos),xx2(pos),'.','Color',[0. 0.7 0.],'MarkerSize',40)
  hold on
   plot(Caa1(pos),xx1(pos),'.','Color',[0.1 0.1 .5],'MarkerSize',20)
  hold on
  plot(Caa2(end),xx2(end),'.','Color',[0. 0.5 0.],'MarkerSize',20)
  hold on 
  

dca=-43;
dx=-4;
% 
%x1=0.1:0.01:.999;
%V1=dx-50-log(1./x1-1)/0.15;
 V1=-68:.1:-30;
x1=1./(1+exp(-0.15*(V1+50-dx)));
% 
gsyn=-0.000;
Iapp=+0.00;

VI=30;
VK=-75;
Vs=127*V1/(VI-VK)-(115*VK+VI*12)/(VI-VK); 
amm=0.1*(50-Vs)./(exp((50-Vs)/10)-1);
bmm=4*exp((25-Vs)/18);                                                                                                                                         
ah=0.07*exp((25-Vs)/20);   
bh=1.0./(1+exp((55-Vs)/10));          
an=0.01*(55-Vs)./(exp((55-Vs)/10)-1);
bn=0.125*exp((45-Vs)/80);  
minf=amm./(amm+bmm);  
hinf=ah./(ah+bh); 
ninf=an./(an+bn);  
tauh=12.5./(ah+bh);  
taun=12.5./(an+bn);                                                                                                                                             
taux=100;                                                                                                                                                                                                                                                                                                                                                                                                     
h1=hinf;                                                                                                                                             
n1=ninf;  
Vhh=-50;
y1=1./(1+exp(10*(V1-Vhh)));   

%Jack 
% r = 127*V1/105 + 8265/105;
% amm=0.1*(50-r)./(exp((50-r)./10)-1);
% bmm=4*exp((25-r)./18);
% minfty = amm./(amm+bmm);

% TTX current
gT=0.01;
% little change from -0.15 as above
xinf=1./(1+exp(-0.15*(V1+50-dx)));  
Ittx=gT*xinf.*(V1-VI);

% leak current 
gL=0.003; VL=-40; 
Ileak=gL*(V1-VL);

%Calcium current 
rho=0.0003; gKCa=0.03; VCa=140;
Ca1=0.0085*x1.*(VCa-V1+dca);
ICa=gKCa*Ca1.*(V1-VK)./(0.5+Ca1);

%very small conntributions 
% Potassium currebt
gK=0.3; 
IK=gK*n1.^4.*(V1-VK);
%Sodium current
gI=4;
INa=-gI*minf.^3.*h1.*(30-V1);

% h-current
gh=0.0000;
Ihh=gh./((1.+exp(-(V1+63.)/7.8))).^3.*y1.*(V1-70);


%all currents  
Q1 =-gsyn*(V1+70) -Iapp -INa  -IK -Ittx -Ileak -Ihh;
 
% solve for Q1 = -0.03*Ca1./(Ca1+0.5).*(V1+75)
Ca2=0.5*Q1./(0.03*(V1-VK)-Q1);

%x22=Ca1./(0.0085*(140-V1+dca));
plot (Ca2,x1,'black','LineWidth',1)
hold on

% nullcline Ca'=0
Q11=-gsyn*(V1+70)-Iapp+gI*minf.^3.*h1.*(VI-V1) + gK*n1.^4.*(VK-V1) +0.03*Ca1.*(V1+75)./(0.5+Ca1) -gL*(V1-VL) -gh*(((1./(1.+exp(-(V1+63.)/7.8))).^3).*y1.*(V1-70));
x11=Q11./(gT*(VI-V1));
Ca3=0.0085*x11.*(140-V1+dca);
plot (Ca3,x11,'r')
hold on



dca=-50;
dx=-4;
% 
%x1=0.1:0.01:.999;
%V1=dx-50-log(1./x1-1)/0.15;
 V1=-71.1:1:-17;
x1=1./(1+exp(-0.15*(V1+50-dx)));
% 
gsyn=0.0083;
Iapp=+0.00;

VI=30;
VK=-75;
Vs=127*V1/(VI-VK)-(115*VK+VI*12)/(VI-VK); 
amm=0.1*(50-Vs)./(exp((50-Vs)/10)-1);
bmm=4*exp((25-Vs)/18);                                                                                                                                         
ah=0.07*exp((25-Vs)/20);   
bh=1.0./(1+exp((55-Vs)/10));          
an=0.01*(55-Vs)./(exp((55-Vs)/10)-1);
bn=0.125*exp((45-Vs)/80);  
minf=amm./(amm+bmm);  
hinf=ah./(ah+bh); 
ninf=an./(an+bn);  
tauh=12.5./(ah+bh);  
taun=12.5./(an+bn);                                                                                                                                             
taux=100;                                                                                                                                                                                                                                                                                                                                                                                                     
h1=hinf;                                                                                                                                             
n1=ninf;  
Vhh=-50;
y1=1./(1+exp(10*(V1-Vhh)));   

%Jack 
% r = 127*V1/105 + 8265/105;
% amm=0.1*(50-r)./(exp((50-r)./10)-1);
% bmm=4*exp((25-r)./18);
% minfty = amm./(amm+bmm);

% TTX current
gT=0.01;
% little change from -0.15 as above
xinf=1./(1+exp(-0.15*(V1+50-dx)));  
Ittx=gT*xinf.*(V1-VI);

% leak current 
gL=0.003; VL=-40; 
Ileak=gL*(V1-VL);

%Calcium current 
rho=0.0003; gKCa=0.03; VCa=140;
Ca1=0.0085*x1.*(VCa-V1+dca);
ICa=gKCa*Ca1.*(V1-VK)./(0.5+Ca1);

%very small conntributions 
% Potassium currebt
gK=0.3; 
IK=gK*n1.^4.*(V1-VK);
%Sodium current
gI=4;
INa=-gI*minf.^3.*h1.*(30-V1);

% h-current
gh=0.0000;
Ihh=gh./((1.+exp(-(V1+63.)/7.8))).^3.*y1.*(V1-70);


%all currents  
Q1 =-gsyn*(V1+80) -Iapp -INa  -IK -Ittx -Ileak -Ihh;
 
% solve for Q1 = -0.03*Ca1./(Ca1+0.5).*(V1+75)
Ca2=0.5*Q1./(0.03*(V1-VK)-Q1);

%x22=Ca1./(0.0085*(140-V1+dca));
plot (Ca2,x1,'Color',[0 0.2 .99],'LineWidth',1)
hold on

% nullcline Ca'=0
Q11=-gsyn*(V1+80)-Iapp+gI*minf.^3.*h1.*(VI-V1) + gK*n1.^4.*(VK-V1) +0.03*Ca1.*(V1+75)./(0.5+Ca1) -gL*(V1-VL) -gh*(((1./(1.+exp(-(V1+63.)/7.8))).^3).*y1.*(V1-70));
x11=Q11./(gT*(VI-V1));
Ca3=0.0085*x11.*(140-V1+dca);
plot (Ca3,x11,'r')
hold on


 print(gcf,'-djpeg','-r600' ,'HCO_winner1_A.jpeg');


