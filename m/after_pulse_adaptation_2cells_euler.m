
%clear all 

alpha = 0.02;
%alpha synapse 
%slow m*s
%fast
beta =  0.01;
taum=100;
%slow
beta =  0.005;

g12=0.1;
%
%g21=0.26;
g21=0.1;
gelec=0.0000;

noise=.0001;

%Parameters3
Ca_shift1= -30.;
Ca_shift2 = -34;
x_shift =  -4;

% Constant stimuli if amy
Iapp=.008; %1Hz
Iapp=.022; %  2Hz
%Iapp=.115; % 5Hz
%Iapp=.655; % 10Hz
Iapp=0.00;

gsyn=0.0215;

t1=35*1000;
t2=45*1000;
      
% H-current      
    gh    = 0.0001;
    Vhh   = -53;

% Intergration    
tf=200;   
step=1.;

% Initial values
V1= -40;
V2= -40; 
Ca1=1.1;
Ca2=1.;
h1 =0;
h2 =0;
n1 =0;
n2 =0;
x1 =0.85;
x2 =0.85;
y1 =0;
y2 =0;
s1 =0;
s2 =0;
m1 =0;
m2 =0;
i=0;
tt=0;

tic
while (tt < tf*1000)
 
V1 =V1 +step*(-gsyn*(V1+60)*heaviside(tt-t1)*heaviside(t2-tt)+4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1) +gh*((1/(1+exp(-(V1+63)/7.8)))^3)*y1*(+120-V1)-g21*(V1+55)*s2*m2+gelec*(V2-V1)+noise*(rand-0.5));
V2= V2 +step*(-gsyn*(V2+60)*heaviside(tt-t1)*heaviside(t2-tt)+4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) +gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)-g12*(V2-55)*s1*m1+gelec*(V1-V2)+noise*(rand-0.5));
Ca1=Ca1+step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
Ca2=Ca2+step*(0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
y1 =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
y2 =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));
s1 =s1+step*(alpha*(1-s1)/(1+exp(-10*(V1+20)))-beta*s1);
s2 =s2+step*(alpha*(1-s2)/(1+exp(-10*(V2+20)))-beta*s2);
tt=tt+step;
i=i+1;
time(i)=tt;
vv1(i)=V1;
vv2(i)=V2;
ss1(i)=s1;
ss2(i)=s2;

Caa1(i)=Ca1;
Caa2(i)=Ca2;
xx1(i)=x1;
xx2(i)=x2;

end  
toc    
    
%  
%  [vp1,tp1] = findpeaks(vv1,time,'MinPeakHeight',-20);
%  fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%  tp1=tp1/1000;
 time=time/1000;

%     
 figure(4)
 clf
 subplot(4,1,1)
 % plot(tp1, vp1,'.','Color',[0 0 0]);
  hold on
 plot(time,vv1,'Color',[0 0  1]','LineWidth',1.5)
 hold on
 %plot(time,vv2-5,'Color',[0 1 0]','LineWidth',1.5)
 hold on
 xlim([0 tf]) 
 ylim([-65 45])
 
 subplot(4,1,2)
 plot(time,vv2,'Color',[0 1 0]','LineWidth',1.)
 hold on
  xlabel('Time'),ylabel('Voltage')
 xlim([0 tf]) 
 ylim([-65 45])

  subplot(4,1,3)
  plot(time,10*mm1.*ss1,'blue')
  hold on
  plot(time,ss1,'Color',[0.5 0.5 0.5])
  hold on
  plot(time,10*mm1,'Color',[0 1 1])
  hold on
  plot(time,10*ss2.*mm2,'Color',[0.0 0.5 0.0])
  hold on
  %ylim([0 .02])
  xlim([0 tf])
  title('s(t) probability of neuro-transmitter release', 'Fontsize', 11);
  
  
  
  subplot(4,1,4)
  plot(time,Caa1,'blue')
  hold on
  plot(time,Caa2,'Color',[0.0 1 0])
  hold on
  plot(time,xx1,'Color',[0 .1 1])
  hold on
  plot(time,xx2,'Color',[0 1 .1])
  hold on
  ylim([0.5 1.05])
  xlim([0 tf])
  
  title('Ca and Ca-activaterd possasion gating varaible in time', 'Fontsize', 11);
  
  
  toc 
 
  figure(5)
  clf
  plot(Caa1,xx1,'b','LineWidth',1.5)
  hold on
  plot(Caa2,xx2,'green','LineWidth',1.5)
  hold on
  xlabel('Ca'),ylabel('x-variable')  
  title('Ca vs Ca-activaterd possasion gating varaible', 'Fontsize', 11);
 
  axis([0.45 1.42 0.2 .97])
 
  
%  save('TS_Q_after_pulse1.mat','time','vv1','vv2','Caa1','Caa2','xx1','xx2')

