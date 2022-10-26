clear all 

% Vs=127/(VI-VK)*V1-(115*VK+VI*12)/(VI-VK);                                                                                                                       
% amm=0.1*(50-Vs)/(exp((50-Vs)/10)-1);                                                                                                                            
% bmm=4*exp((25-Vs)/18);                                                                                                                                          
% ah=0.07*exp((25-Vs)/20);                                                                                                                                        
% bh=1.0/(1+exp((55-Vs)/10));                                                                                                                                     
% an=0.01*(55-Vs)/(exp((55-Vs)/10)-1);                                                                                                                            
% bn=0.125*exp((45-Vs)/80);                                                                                                                                       
% minf=amm/(amm+bmm);                                                                                                                                             
% hinf=ah/(ah+bh);                                                                                                                                                
% ninf=an/(an+bn);                                                                                                                                                
% xinf=1/(exp(0.2*(-V1-50+xshift))+1);                                                                                                                            
% tauh=12.5/(ah+bh);                                                                                                                                              
% taun=12.5/(an+bn);                                                                                                                                              
% taux=100;                                                                                                                                                       
% V1'=-gsyn*(V1+60)+Iapp+gI*minf^3*h1*(VI-V1)+gK*n1^4*(VK-V1)+gT*x1*(VI-V1)+gKCa*Ca1/(.5+Ca1)*(VK-V1)+gL*(VL-V1)+gh*(((1./(1.+exp(-(V1+63.)/7.8)))^3)*y1*(70-V1));
% h1'=(hinf-h1)/tauh;                                                                                                                                             
% n1'=(ninf-n1)/taun;                                                                                                                                             
% y1'=0.5*((1./(1+exp(10*(V1+50.))))-y1)/(7.1+10.4/(1+exp((V1+68.)/2.2)));                                                                                        
% x1'=(xinf-x1)/taux;                                                                                                                                             
% Ca1'=rho*( Kc*x1*(VCa-V1+Cashift)-Ca1 );                                                                                                                          
dca=-35;
dx=-4;
s=0;
% 
%x1=0.1:0.01:.999;
%V1=dx-50-log(1./x1-1)/0.15;
 V1=-68:1:-15;
x1=1./(1+exp(-0.15*(V1+50-dx)));
% 
 figure(5)
 clf 
gsyn=0.000;
Iapp=+0.00;

VI=30;
VK=-75;
Vs=127*V1/(VI-VK)-(115*VK+VI*12)/(VI-VK); 
amm=0.1*(50-Vs)./(exp((50-Vs)/10)-1);
bmm=4*exp((25-Vs)/18);                                                                                                                                         
ah=0.07*exp((25-Vs)/20);   
size(ah)
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
plot (Ca2,x1,'blue','LineWidth',2)
hold on

% nullcline Ca'=0
Q11=-gsyn*(V1+70)-Iapp+gI*minf.^3.*h1.*(VI-V1) + gK*n1.^4.*(VK-V1) +0.03*Ca1.*(V1+75)./(0.5+Ca1) -gL*(V1-VL) -gh*(((1./(1.+exp(-(V1+63.)/7.8))).^3).*y1.*(V1-70));

x11=Q11./(gT*(VI-V1));
Ca3=0.0085*x11.*(140-V1+dca);
plot (Ca3,x11,'r')
hold on

load('ca_sweep_ih_0.mat')
% %  % low stable branch 
Vvar1= x(1,1:end);
nvar= x(2,1:end); 
xvar1= x(5,1:end);
Cavar1=x(6,1:end);
plot(Cavar1, xvar1,'Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
load('eq1.mat')
Vvar2= x(1,1:end);
nvar= x(2,1:end); 
xvar2= x(5,1:end);
Cavar2=x(6,1:end);
plot(Cavar2,xvar2,'Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
load('eq2.mat')
Vvar3= x(1,1:end);
nvar= x(2,1:end); 
xvar3= x(5,1:end);
Cavar3=x(6,1:end);
plot(Cavar3,xvar3,'Color',[0.8 0.8 0.8],'LineWidth',1)
hold on

axis([-0.2 1.6 -0.1 1.2])


 figure(6)
 clf 
 plot(Cavar2,Vvar2,'Color',[0.8 0.8 0.8],'LineWidth',2)
 hold on
plot(Cavar3,Vvar3,'Color',[0.8 0.8 0.8],'LineWidth',2)
 hold on
 plot(Ca2,V1,'Color',[0.8 0. 0.],'LineWidth',2)
% 
% 
% 
 axis([-0.2 1.6 -70 0 ])



%figure(8)
%clf 
%plot (INa,V1,'--')
hold on
%plot (IK,V1,'.')
hold on
%plot (Ittx,V1)
% hold on
% plot (INa,V1)
% hold on
% plot (IK,V1)
% hold on
