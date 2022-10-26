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
function Melibe_Fig


global alpha beta alpha1 beta1 g12 g21 gelec Iapp Iapp1 Iapp2 t1 t2  t3 t4  t5 t6 Ca_shift x_shift gh Vhh

alpha = 0.01;
beta =  0.008;

g12=.2;
g21=0;
gelec=0.0002;


Ca_shift = -20;
x_shift =  -4;

% Positive Constant stimuli
Iapp=.013; %1Hz
Iapp2=.04;
Iapp1=.15;

%Iapp=.03; %  2Hz
%Iapp=.115; % 5Hz
%Iapp=.655; % 10Hz


% Negative constant stimuli
%Iapp=-.025; %  




t1=26*1000;
t2=28*1000;

t5=32*1000;
t6=34*1000;

t3=38*1000;
t4=40*1000;




alpha1 = 0.002;
beta1 = 0.002;

% H-current      
    gh    = 0.0001;
    Vhh   = -53;

% Initial Conditions
%      V Ca  h  n  x y s r m   
s0 = [-43.85 1.1 0 0 1 0   0 0 0 0 0 43.85 1.2 0. 0 1 0 ];
opts = odeset('Stats','off','RelTol',1e-7,'AbsTol',1e-7);
[t,s] = ode15s(@Plant_Model,[0:50*1000], s0,opts);
t=t/1000;
%  
% [vp1,tp1] = findpeaks(s(:,1),t,'MinPeakHeight',-20);
% fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%     

load('post_alpha.mat')
load('post_dyn.mat')


clf 
 figure(2)
 ax = gca;
ax.FontSize = 14;
 [ha, pos] = tight_subplot(5,1,[.01 .01],[.05 .02],[.05 .01]) 
axes(ha(1));
 plot(t,s(:,1),'b','LineWidth',2)
hold on  
 xlim([25 42])
 ylabel('pre-Voltage','Fontsize',16)

 
 axes(ha(2));
  ax = gca;
ax.FontSize = 14;
 plot(t,s(:,12),'r','LineWidth',2)
 hold on
 plot(t,post_dyn,'Color',[.1 .7 .7],'LineWidth',2)
 hold on
 plot(t,post_alpha,'Color',[.1 .1 .1],'LineWidth',2)
 hold on

 
 %post_alpha=s(:,12);
 %save ('post_dyn.mat','t','post_dyn') 
 %save ('post_alpha.mat','t','post_alpha') 
 
 ylabel('post-Voltage','Fontsize',16)
  xlim([25 42]) 
  ylim([-55 -42])
  axis on
  
  axes(ha(3));
   ax = gca;
ax.FontSize = 14
   plot(t,1.6*s(:,7),'Color',[0.1 0.1 0.1],'LineWidth',1.5)
 hold on
 plot(t,1.6*s(:,8),'Color',[.5 .5 .5],'LineWidth',2.5)
 hold on
 ylabel('\alpha-synapse','Fontsize',16)
 xlim([25 42]) 
  ylim([-0.01 1])
 set(ha(1:5),'XTickLabel','');
;


 axes(ha(4));
  ax = gca;
ax.FontSize = 14
 plot(t,20*s(:,7).*s(:,10),'Color',[.1 .7 .7],'LineWidth',2)
 hold on
 ylabel('dyn-synapse','Fontsize',16)
 xlim([25 42]) 
 ylim([-0.01 1]) 


axes(ha(5));
plot(t,1.4*s(:,11),'red','LineWidth',2)
hold on
xlim([25 42])   
 ylim([-0.01 1])
set(ha(1:5),'XTickLabel',''); 
ax = gca;
ax.FontSize = 12

xlabel('time [sec]','Fontsize',16),ylabel('sigm-synapse','Fontsize',16)

%set(ha,'YTickLabel','')
set(gca,'XTickLabel',cellstr(num2str(get(gca,'XTick')'))) 
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')'))) 
ax = gca;


print(gcf,'-djpeg','-r600' ,' ZOO_synapses_squre_tight.jpeg');

end

function [qd] = Plant_Model(t,q)
global alpha beta alpha1 beta1 g12 g21 gelec Iapp Iapp1 Iapp2 t1 t2 t3 t4 t5 t6 Ca_shift x_shift gh Vhh
  qd = zeros(17,1);
          
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

    
qd(1) =Iapp*heaviside(t-t1)*heaviside(t2-t)+Iapp1*heaviside(t-t3)*heaviside(t4-t)+...
      Iapp2*heaviside(t-t5)*heaviside(t6-t)+ 4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - ...
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
qd(10)  =(.01+0.9/(1+exp(-1000*(V1+40)))-2*m1)/8000;
qd(11) =2.8*alpha*s3*(1-s3)/(1+exp(-1000*(V1+40)))-.006*(s3-0.001);

qd(12)= 4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))...
    /((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2) +...
    0.01*x2*(30-V2) +0.03*Ca2/(.5 + Ca2)*(-75 - V2) + 0.008*(-40 - V2) +...
    gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)  -g12*(V2+55)*s3;


qd(13) =0.00003*(0.0085*x2*(140 - V2+Ca_shift)-Ca2);
qd(14) =((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5;
qd(15) =((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5;
qd(16) =((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100;
qd(17) =.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2)));
end


     %clf
%  subplot(5,1,1)
%  plot(t,s(:,1),'b','LineWidth',1.5)
%  hold on
%  % plot(t,s(:,12),'r','LineWidth',1)
%  hold on
%  xlim([25 42]) 
% 
%  subplot(5,1,2)
%  %plot(t,s(:,1),'b','LineWidth',1.5)
%  hold on
%   plot(t,s(:,12),'r','LineWidth',1)
%  hold on
%  xlim([25 42]) 
%  
%  subplot(5,1,3)
%  plot(t,s(:,7),'Color',[0 0 0],'LineWidth',1.5)
%  hold on
%  plot(t,s(:,8),'Color',[.1 .6 .9],'LineWidth',1.5)
%  hold on
% % plot(t,s(:,9),'Color',[.6 .0 .0],'LineWidth',1.5)
%  hold on
%  ylabel('\alpha-synapse')
%  xlim([25 42])   
%  subplot(5,1,4)
%  plot(t, s(:,7).*s(:,10),'Color',[.1 .7 .7],'LineWidth',1.5)
%  hold on
%  % plot(t,s(:,10),'Color',[.8 .8 .8],'LineWidth',1)
%  hold on
%  xlabel('Time'),ylabel('Voltage')
%  xlim([25 42])   
%  subplot(5,1,5)
%   plot(t,s(:,11),'red','LineWidth',1.5)
%  hold on
%   plot(t,movavg(s(:,11)),'red','LineWidth',1.5)
%  hold on
 
%  xlabel('Time'),ylabel('Voltage')
%  xlim([25 42])   
%  

