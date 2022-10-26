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
% Ca1'=rho*(Kc*x1*(VCa-V1+Cashift)-Ca1);                                                                                                                          

gI=4;
gK=0.3;
gT=0.01;
gKCa=0.03;
gL=0.003;
gh=0.0005;
VI=30;
VK=-75;
VL=-40;
VCa=140;
Vca=140;
Kc=0.0085;
rho=0.0003;
gsyn=0.00;
Iapp=-0.0;

dca=-40;
dx=-4;
g21=0;
s2=0;
gh=0.000;
%dca_shift=-100:1:100;

V1=-68:0.1:-40;
figure(5)
clf 
h1=(0.07*exp((25-(127*V1/105+8265/105))/20))./(0.07*exp((25-(127*V1/105 + 8265/105))/20) + (1.0./(1 + exp((55 - (127*V1/105 + 8265/105))/10)))) ;
n1= (0.01*(55 -(127*V1/105+8265/105))/(exp((55-(127*V1/105+8265/105))/10)- 1))./((0.01*(55 - (127*V1/105 + 8265/105))./(exp((55-(127*V1/105+8265/105))/10)-1))+(0.125*exp((45-(127*V1/105+8265/105))/80)) )  ;
Vhh   = -50;
y1=1./(1+exp(10*(V1-Vhh)));
x1=1./(1+exp(-0.15*(V1+50-dx)));  
Q=4*((0.1*(50-(127*V1./105+8265/105))/(exp((50 - (127*V1./105 ...
    +8265/105))/10) - 1))./((0.1*(50 - (127*V1./105 + 8265/105))./(exp((50 - (127*V1./105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18)))).^3.*h1.*(30 - V1) +gh*1./((1+exp(-(V1+63)/7.8))).^3.*y1.*(-V1+70) + 0.01*x1.*(30-V1) ...
    + 0.003*(-40-V1) + 0.3*n1.^4.*(-75 - V1) -g21*(V1+80);
Ca1=0.5*Q./(-Q+0.03*(75+V1));
plot (Ca1,x1,'b','LineWidth',5)
hold on
plot (Ca1,x1,'b','LineWidth',5)
hold on

clear x1 Ca1 V1
V1=-69:0.1:40;
Vs=127/(VI-VK)*V1-(115*VK+VI*12)/(VI-VK);                                                                                                                      
amm=0.1*(50-Vs)/(exp((50-Vs)/10)-1);                                                                                                                           
bmm=4*exp((25-Vs)/18);                                                                                                                                         
ah=0.07*exp((25-Vs)/20);                                                                                                                                       
bh=1.0./(1+exp((55-Vs)/10));                                                                                                                                    
an=0.01*(55-Vs)/(exp((55-Vs)/10)-1);                                                                                                                           
bn=0.125*exp((45-Vs)/80);                                                                                                                                      
minf=amm./(amm+bmm);                                                                                                                                            
hinf=ah./(ah+bh);                                                                                                                                               
ninf=an./(an+bn);                                                                                                                                               
tauh=12.5./(ah+bh);                                                                                                                                             
taun=12.5./(an+bn);                                                                                                                                             
taux=100;                                                                                                                                                                                                                                                                                                                                                                                                     
h1=hinf;                                                                                                                                             
n1=ninf;                                                                                                                                             
y1=1./(1+exp(10*(V1+50.)));    
xinf=1./(1+exp(-0.2*(V1+50-dx)));                                                                                                                            
Q1=-gsyn*(V1+80)+Iapp+gI*minf.^3.*h1.*(VI-V1)+gK*n1.^4.*(VK-V1) + gT*xinf.*(VI-V1)+gL*(VL-V1)+gh*(((1./(1.+exp(-(V1+63.)/7.8))).^3).*y1.*(70-V1));
% solve for Q1 = -0.03*Ca1./(Ca1+0.5).*(V1+75)
Ca1=0.0085*xinf.*(140-V1+dca);
Ca2=0.5*Q1./(0.03*(V1+75)-Q1);
x22=Ca2./(0.0085*(140-V1+dca));
plot (Ca2,xinf,'green')

%Ca2=16.6667*Q1./(V1-33.333.*Q1+75);

% nullcline Ca'=0
Ca1=0.0085*xinf.*(140-V1+dca);
Q11=-gsyn*(V1+80)+Iapp+gI*minf.^3.*h1.*(VI-V1)+gK*n1.^4.*(VK-V1)+0.03*Ca1./(0.5+Ca1).*(V1+75) +gL*(VL-V1)+gh*(((1./(1.+exp(-(V1+63.)/7.8))).^3).*y1.*(70-V1));
x11=Q11./(gT*(VI-V1));
Ca3=0.0085*x11.*(140-V1+dca);
plot (Ca3,x11,'r')
hold on

load('ca_sweep_ih_0.mat')
% %  % low stable branch 
Vvar= x(1,1:end);
nvar= x(2,1:end); 
xvar= x(5,1:end);
Cavar=x(6,1:end);
plot(Cavar, xvar,'Color',[0.8 0.8 0.8],'LineWidth',2)
hold on
load('eq1.mat')
Vvar= x(1,1:end);
nvar= x(2,1:end); 
xvar= x(5,1:end);
Cavar=x(6,1:end);
plot(Cavar,xvar,'Color',[0.8 0.8 0.8],'LineWidth',2)
hold on
load('eq2.mat')
Vvar= x(1,1:end);
nvar= x(2,1:end); 
xvar= x(5,1:end);
Cavar=x(6,1:end);
plot(Cavar,xvar,'Color',[0.8 0.8 0.8],'LineWidth',2)
hold on

axis([0.65 1.1 0.1 .95] )


figure(6)
clf 
plot (Ca1, V1)
hold on 
plot(Cavar,Vvar,'Color',[0.8 0.8 0.8],'LineWidth',2)
hold on

axis([-0.2 2.6 -80 40])