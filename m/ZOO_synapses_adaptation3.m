% Model 
% V : Membrane Potential(mV)
% t : Time (ms)33
% C_m : Membrane capacitance(pico farad/cm^2)
% g_i : Maximum conductance 
% V_i : Nerst reversal potential(mV) 
% x, n : dimensionless ativation variables
% h : Dimensionless inactivation variable
% m_inf : Instantenous activation function for the fast inward current
% Ca : The dimensionless calcium ion concentration at the innner cell
clear all 

global alpha beta alpha1 beta1 g12 g21 gelec Iapp Iapp1 Iapp2 t1 t2  t3 t4  t5 t6 Ca_shift x_shift gh Vhh Ca_shift2

alpha = 0.005;
beta =  0.002;

g12=.02;
g21=0;
gelec=0.0002;


Ca_shift = -8;
x_shift =  -4;
Ca_shift2 =-10;

% Positive Constant stimuli
Iapp=-.2; %1Hz
t1=5*1000;
t2=15*1000;

% H-current      
    gh    = 0.0001;
    Vhh   = -53;

% Initial Conditions
%      V Ca  h  n  x y s r m   
s0 = [-43.85 1.1 0 0 1 0 0 0 0 0 0 -44 1.25 0 0 0.82 0 0.0];
opts = odeset('Stats','off','RelTol',1e-4,'AbsTol',1e-5);
[t,s] = ode15s(@Plant_Model,[0:70*1000], s0,opts);
t=t/1000;
%  
% [vp1,tp1] = findpeaks(s(:,1),t,'MinPeakHeight',-20);
% fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%     

 
%load('post_alpha1.mat')
%load('post_dyn1.mat')

tmax=65;
tmin=2;

figure(4)
clf 
ax = gca;
ax.FontSize = 14;
 [ha, pos] = tight_subplot(5,1,[.01 .01],[.05 .02],[.06 .01]); 

axes(ha(1));
plot(t,s(:,1),'b','LineWidth',1.5)
hold on  
 xlim([tmin tmax])
 ylim([-80 40])
 ylabel('pre-Voltage','Fontsize',18)
ax = gca;
ax.FontSize = 14;
 
 axes(ha(2));
 plot(t,s(:,12),'r','LineWidth',1.5)
 hold on
 %plot(t,post_dyn,'Color',[.1 .7 .7],'LineWidth',1.5)
 hold on
 %plot(t,post_alpha,'Color',[.1 .1 .1],'LineWidth',1.5)
 hold on
 ax = gca;
ax.FontSize = 14;
 
 xlim([tmin tmax])  

% post_alpha=s(:,12);
% post_dyn=s(:,12);
% save ('post_dyn1.mat','t','post_dyn') 
%save ('post_alpha1.mat','t','post_alpha') 
 
 ylabel('post-voltage','Fontsize',18)
%  xlim([0 42]) 
  ylim([-55 -42]) ; xlim([tmin tmax])  
  axis on
 
  clear vp vp1 tp tp1 avs tav
  
 [vp,tp] =  findpeaks(-s(:,7)-.002,t,'MinPeakHeight',-2);
 [vp1,tp1] = findpeaks(s(:,7),t,'MinPeakHeight', 0.002);

 for i =1:min(length(vp),length(vp1))
    avs(i)=(-vp(i)+vp1(i))/2;
    tav(i)=tp(i);
 end

[vp3,tp3] = findpeaks( -s(:,8)-0.002,t,'MinPeakHeight', -2);
[vp4,tp4] = findpeaks(  s(:,8)+0.001,t,'MinPeakHeight',0.002);
size(vp3)
size(vp4)
% 
for i =1:min(length(vp3),length(vp4))
    avs3(i)=(-vp3(i)+vp4(i))/2;
    tav3(i)=tp3(i);
end

 axes(ha(3));
 plot(t,s(:,7),'Color',[0.7 0.7 0.7],'LineWidth',2.)
 hold on
 %plot(tp1,vp1,'.','MarkerSize',10','Color',[1 0 0],'LineWidth',1.5)
 hold on
 %plot(tp,-vp,'.','MarkerSize',10','Color',[0 0 1],'LineWidth',1.5)
 hold on
 plot(tav,avs,'-','Color',[1 1 1],'LineWidth',1.5)
 hold on
 plot(t,s(:,8),'Color',[.1 .1 .1],'LineWidth',2.)
 hold on
% plot(tp3,vp3,'.','MarkerSize',10','Color',[1 0 0],'LineWidth',1.5)
% hold on
% plot(tp4,-vp4,'.','MarkerSize',10','Color',[0 0 1],'LineWidth',1.5)
% hold on
 
 plot(tav3,avs3,'-','Color',[.99 .99 .99],'LineWidth',1.5)
 hold on
 ylabel('\alpha-synapse','Fontsize',18)
 ylim([0 1]); ; xlim([tmin tmax])   
 set(ha(1:5),'XTickLabel','');
 ax = gca;
ax.FontSize = 16;

axis tight 
 
 [vp4,tp4] = findpeaks( -s(:,7).*s(:,10),t,'MinPeakHeight', -2);
 [vp5,tp5] = findpeaks  (s(:,7).*s(:,10),t,'MinPeakHeight',0.0);
 
     
for i =1:min(length(vp),length(vp1))
    avs1(i)=(-vp4(i)+vp5(i))/2;
    tav1(i)=tp5(i);
end

 axes(ha(4));
 plot(t,7*s(:,7).*s(:,10),'Color',[.1 .7 .7],'LineWidth',2.)
 hold on
 hold on
 plot(tav1,7*avs1,'-','Color',[1 1 1],'LineWidth',1.)
  hold on
 ylabel('dyn-synapse','Fontsize',18)
 ylim([0 1]) ; xlim([tmin tmax])  
 ax = gca;
ax.FontSize = 14;
axis tight 


 [vp,tp] = findpeaks(-s(:,11),t,'MinPeakHeight', -2);
 [vp1,tp1] = findpeaks(s(:,11),t,'MinPeakHeight',0);
% fr1=length(tp1)/(tp1(end)-tp1(1))*1000


for i =1:min(length(vp),length(vp1))
    avs2(i)=(-vp(i)+vp1(i))/2;
    tav2(i)=tp1(i);
end
    
 axes(ha(5));
 plot(t,5.4*s(:,11),'red','LineWidth',2.)
   hold on
% %plot(tp1,1.4*vp1)
%   hold on
% %plot(tp,-1.4*vp)
%   hold on
 plot(tav2,5.4*avs2,'-','Color',[1 1 1],'LineWidth',1.)
   hold on  
%   
   ax = gca;
 ax.FontSize = 14;
% 
  xlabel('time','Fontsize',18),ylabel('sigm-synapse','Fontsize',18)
  ylim([-.0 1]); xlim([tmin tmax])  
  axis tight 

set(ha(1:5),'XTickLabel',''); 
%set(ha,'YTickLabel','')
set(gca,'XTickLabel',cellstr(num2str(get(gca,'XTick')'))) 
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')'))) 
ax = gca;
% 


%print(gcf,'-djpeg','-r600' ,'PIR_adaptation.jpeg');

% 
% figure(2)
% clf
% plot(s(:,13),s(:,16),'r','LineWidth',1.5)
%  hold on

function [qd] = Plant_Model(t,q)
global alpha beta alpha1 beta1 g12 g21 gelec Iapp Iapp1 Iapp2 t1 t2 t3 t4 t5 t6 Ca_shift Ca_shift2 x_shift gh Vhh
  qd = zeros(18,1);
          
      V1 = q(1);                       % Fast
     Ca1 = q(2);                       % Slow
      h1 = q(3);                       % Fast 
      n1 = q(4);                       % Fast
      x1 = q(5);                       % Slow
      y1 = q(6);         
      s1 = q(7);
      r1 = q(8);
      r2 = q(9);
      m1 = q(10); 
      s3 = q(11);                       % Fast
     
      V2 = q(12);
      Ca2 =q(13);                       % Slow
      h2 = q(14);                       % Fast 
      n2 = q(15);                       % Fast
      x2 = q(16);                       % Slow
      y2 = q(17);         
      Iap=q(18);
    
qd(1) =Iapp*heaviside(t-t1)*heaviside(t2-t)+ 4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - ...
    (127*V1/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 -...
    (127*V1/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30- V1)...
    + 0.3*n1^4*(-75 - V1) + 0.01*x1*(30-V1) +0.03*Ca1/(.5 + Ca1)*(-75 - V1) +...
    0.003*(-40 - V1) +gh*((1/(1+exp(-(V1+63)/7.8)))^3)*y1*(+120-V1);  
qd(2)  =0.00003*(0.0085*x1*(140 - V1+Ca_shift)-Ca1);
qd(3)  =((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5;
qd(4)  =((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5;
qd(5)  =((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100;
qd(6)  =.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2)));
qd(7)  =alpha*(1-s1)/(1+exp(-100*(V1+40)))-beta*s1;
qd(8)  =alpha*(1-r1)*s1-beta*r1;
qd(9) = alpha*(1-r2)*r1-beta*r2;
qd(10)  =(.0+1.0/(1+exp(-1000*(V1+10)))-0.8*m1)/1200;
qd(11) =1.3*alpha*s3*(1-s3)/(1+exp(-1000*(V1+40)))-.0016*(s3-0.004);

qd(12)= 4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))...
    /((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2) +...
    0.01*x2*(30-V2) +0.03*Ca2/(.5 + Ca2)*(-75 - V2) + 0.008*(-40 - V2) +...
    gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)  -g12*(V2+55)*s1;
qd(13) =0.00003*(0.0085*x2*(140 - V2+Ca_shift2)-Ca2);
qd(14) =((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5;
qd(15) =((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5;
qd(16) =((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100;
qd(17) =.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2)));

qd(18)=0.0000;

end


