
clear all 

tmax=40; % max time in seconds


%--------- green, left small canard; blue: ts ---------------
alpha1 = 0.02;
beta1 =  0.002;

alpha2 = 0.01;
beta2 =  0.001;

g12=.02;
g21=.01;


gelec=0.0001;


%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.004;
%----------------------

%Parameters3
Ca_shift1 = -60;
Ca_shift2 = -20;

x_shift =  -4;

% Constant stimuli if amy
Iapp=.008; %1Hz
Iapp=.022; %  2Hz
%Iapp=.115; % 5Hz
%Iapp=.655; % 10Hz
Iapp=0.0;

t1=00*1000; % this time is in msec
t2=00*1000; % this time is in msec
      
% H-current      
    gh    = 0.0004;
    Vhh   = -52;

% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.2;       % time step in msec
tsamp=2; % sampling interval to save data in msec



cut1=cutoff/1000; % cutoff frequency converted to msec
% Initial values
V1= -44; V2= -44; Ca1=.4; Ca2=.9; h1 =0; h2 =0; n1 =0; n2 =0; x1 =0.8;
x2 =.0; y1 =0; y2 =0; s1 =0; s2 =0; m1 =0; m2 =0;


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
-g21*(V1+80)*s2+gelec*(V2-V1)+rnd1(i));
V2= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
    +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) +gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)...
-g12*(V2-40)*s1+gelec*(V1-V2)+rnd2(i));
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
s1 =s1+step*(alpha1*(1-s1)/(1+exp(-50*(V1+30)))-beta1*(s1-0.000));
s2 =s2+step*(alpha2*s2*(1-s2)/(1+exp(-50*(V2+20)))-beta2*(s2-0.00002));
% m1 =m1+step*((1/(1+exp(-(V1+20)))-m1)/taum);
% m2 =m2+step*((1/(1+exp(-(V2+20)))-m2)/taum);

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
 subplot(3,1,1)
%  plot(tp1, vp1,'.','Color',[0 0 0]);
  hold on
 plot(time,vv1,'Color',[0 0  1]','LineWidth',1.5)
 hold on
%  plot(time,heavy*10,'Color',[1 0  1]','LineWidth',1.5)
%  hold on
 xlim([0 time(end)]) 
 ylim([-85 40])
 
 %xlim([350 time(end)])
 
 subplot(3,1,2)
 plot(time,vv2,'Color',[0 1 0]','LineWidth',1.5)
 hold on
 xlabel('Time'),ylabel('Voltage')
 xlim([0 time(end)]) 
 ylim([-85 40])
 
  %xlim([350 time(end)])

  subplot(3,1,3)
 % plot(time,10*mm1.*ss1,'blue')
  hold on
  plot(time,ss1,'Color',[0.1 0.1 0.8],'LineWidth',1.5)
  hold on
  hold on
  plot(time,ss2,'Color',[0.1 0.8 0.1],'LineWidth',1.5)
  hold on
  %ylim([0 .02])
  xlim([0 time(end)])
  
  title('s(t) probability of neuro-transmitter release', 'Fontsize', 11);
  
  
  
%   subplot(4,1,4)
%   plot(time,Caa1,'blue')
%   hold on
%   plot(time,Caa2,'Color',[0.0 1 0])
%   hold on
%   plot(time,xx1,'Color',[0 .1 1])
%   hold on
%   plot(time,xx2,'Color',[0 1 .1])
%   hold on
%   ylim([0.5 1.05])
%   xlim([0 time(end)])
%   
%   xlim([0 time(end)])
%   title('Ca and Ca-activaterd possasion gating varaible in time', 'Fontsize', 11);
  
  
  toc 
 
  figure(5)
  clf
  plot(Caa1,xx1,'b','LineWidth',1.5)
  hold on
  plot(Caa2,xx2,'green','LineWidth',1.5)
  hold on
  xlabel('Ca'),ylabel('x-variable')  
  title('Ca vs Ca-activaterd possasion gating varaible', 'Fontsize', 11);
  
%   figure(6)
% 
%    plot(ss1,ss2,'Color',[0.2 0.2 0.2],'LineWidth',1.5)