
clear all 

tmax=1800; % max time in seconds


%--------- green, left small canard; blue: ts ---------------
alpha = 0.03;
beta =  0.003;
g12=.0;
g21=.0;
gelec=0;

%----------------------
cutoff=20.0; % cutoff frequency in Hz
noise=.000;
%----------------------

%Parameters3
Ca_shift = -8;
x_shift =  -1.86;
x_shift =  -1.8;

x_shift =  -2;
Ca_shift = -4.6;

% Constant stimuli if amy
Iapp=.008; %1Hz
Iapp=.022; %  2Hz
%Iapp=.115; % 5Hz
%Iapp=.655; % 10Hz
Iapp=0.0;

t1=00*1000; % this time is in msec
t2=00*1000; % this time is in msec
      
% H-current      
    gh    = 0.00;
    Vhh   = -55;

% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.25;       % time step in msec
tsamp=1; % sampling interval to save data in msec



cut1=cutoff/1000; % cutoff frequency converted to msec
% Initial values
V1= -48; V2= 20; Ca1=0.91; Ca2=0.92; h1 =0; h2 =0; n1 =0; n2 =0; 
x1 =0.6;x2 =.65; y1 =0; y2 =0; s1 =0; s2 =0; m1 =0; m2 =0;


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
    +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1) +gh*((1/(1+exp(-(V1+63)/7.8)))^3)*y1*(+70-V1)...
-g21*(V1+80)*s2+gelec*(V2-V1)+rnd1(i));
V2= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
    +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) +gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+70)...
-g12*(V2+80)*s1+gelec*(V1-V2)+rnd2(i));
Ca1=Ca1+step*(0.0003*(0.0085*x1*(140-V1+Ca_shift)-Ca1));
Ca2=Ca2+step*(0.0003*(0.0085*x2*(140-V2+Ca_shift)-Ca2));
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
%y1 =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
%y2 =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));
%s1 =s1+step*(alpha*s1*(1-s1)/(1+exp(-10*(V1+20)))-beta*(s1-0.001));
%s2 =s2+step*(alpha*s2*(1-s2)/(1+exp(-10*(V2+20)))-beta*(s2-0.001));
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
 time=time/1000;
[v1min,tp1min] = findpeaks(-vv1,time,'MinPeakProminence',1,'MinPeakHeight',52);
[v2min,tp2min] = findpeaks(-vv2,time,'MinPeakProminence',1,'MinPeakHeight',52);
%  fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%  
%  % convert back to sec



%%
 %     
 figure(4)
clf
 subplot(4,1,1)
 plot(time,vv1,'Color',[0 0  .8]','LineWidth',1.5)
 hold on
 plot(tp1min,-v1min,'.','Color',[1 0   0]','MarkerSize',10)
 hold on
 
 xlabel('Time [sec]','Fontsize', 14),ylabel('Voltage_1','Fontsize', 16)
 xlim([0 tmax]) 
 ylim([-75 40])
  fftt = 28;
txpo = text(-82.2,25,'B','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')
 
 
 %xlim([350 time(end)])
 
 subplot(4,1,2)
 plot(time,vv2,'Color',[0 .7 0]','LineWidth',1.5)
 hold on
 xlabel('Time [sec]','Fontsize', 14),ylabel('Voltage_2','Fontsize', 16)
 xlim([0 tmax]) 
 ylim([-75 40])
 box on 
  
 subplot(4,1,3)
 plot(time,Caa1,'Color',[0 .0 .8]','LineWidth',1.5)
 hold on
[Ca1min,tp1min] = findpeaks(-Caa1,time,'MinPeakProminence',.01);
[Ca2min,tp2min] = findpeaks(-Caa2,time,'MinPeakProminence',.01);
 plot(tp1min,-Ca1min,'.','Color',[1 0   0]','MarkerSize',10)
 hold on

 xlabel('Time [sec]','Fontsize', 14),ylabel('Ca_1','Fontsize', 16)
 xlim([0 tmax]) 
 ylim([0.4 1.2])
 box on 
 
  subplot(4,1,4)
 plot(time,xx1,'Color',[0 .0 .8]','LineWidth',1.5)
 hold on
[x1min,tp1min] = findpeaks(-xx1,time,'MinPeakHeight',-0.5);
[x2min,tp2min] = findpeaks(-xx2,time,'MinPeakHeight',-0.5);
 plot(tp1min,-x1min,'.','Color',[1 0   0]','MarkerSize',10)
 hold on

 xlabel('Time [sec]','Fontsize', 14),ylabel('x_1','Fontsize', 16)
 xlim([0 tmax]) 
 ylim([0.0 1.])
 box on 
 


   %print(gcf,'-djpeg','-r600' ,'HCO1a.jpeg');
  
  toc 
 
  figure(5)
  clf
  plot(Caa1,xx1,'Color',[0. 0.0 0.8],'Linewidth',1.0)
  hold on
 % plot(Caa2,xx2,'Color',[0. 0.7 0.],'Linewidth',1.)
  hold on
  %plot(Caa2(end-1000:end),xx2(end-1000:end),'.','Color',[1. 0.0 0.],'Linewidth',3.)
  hold on
  
   
  xlabel('[Ca]-variable','Fontsize', 16),
  ylabel('x-variable','Fontsize', 16)  
  %title('Ca vs Ca-activaterd possasion gating varaible', 'Fontsize', 11);
  
  %SNIC part 
%--------------------- SNIC whole system gh = 0 -------------------------------
% data1 = load('LP_whole_bk_gh0.mat','x');
% cod_back = data1.x;
% data1 = load('LP_whole_fwd_gh0.mat','x');
% cod_fwd = data1.x;
% x_cod_bkgh0 = cod_back(end-3,:);
% Ca_cod_bkgh0 = cod_back(end-2,:);
% x_cod_fdgh0 =  cod_fwd(end-3,1:150);
% Ca_cod_fdgh0 = cod_fwd(end-2,1:150);
% 
% %------------------- plot SNIC ---------------------
% plot(Ca_cod_bkgh0,x_cod_bkgh0+0.015,'-.','Color',[.85 .85 .85 ],'LineWidth',2)
% hold on
% plot(Ca_cod_fdgh0,x_cod_fdgh0+0.015,'-.','Color',[.85 .85 .85 ],'LineWidth',2)
% hold on

axis([0.5 1.1 0.1 1.0])

fftt = 28;
txpo = text(0.,.98,'А','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')

pos=1;
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
  plot(Caa2(pos),xx2(pos),'.','Color',[0. 0.5 0.],'MarkerSize',20)
  hold on 
  

  v1min1=circshift(v1min,1);
  v2min1=circshift(v2min,1);
figure(7)
clf
%plot (-v1min1,-v1min,'Color',[0.7 0.7 0.7])
hold on

for i =1:length(v1min)-2
plot ([- v1min1(i),-v1min(i)],[-v1min(i),-v1min(i)],'Color',[0.2 0.2 0.2],'Linewidth',1)
hold on
plot ([-v1min(i),-v1min(i)],[-v1min1(i+1),-v1min(i+1)],'Color',[0.2 0.2 0.2],'Linewidth',1)
hold on
end

plot (-v1min1(1),-v1min(1),'.','Color','r','Markersize',30)
hold on
plot (-v1min1,-v1min,'.','Color','b','Markersize',8)
hold on
plot (-v2min1,-v2min,'.','Color','g','Markersize',8)
hold on

% for i=1:length(v1min)-1
% plot (-v1min(i),-v1min(i+1),'x','Color','r','Markersize',30)
% hold on 
% end

plot (-v1min1(200:end),-v1min(200:end),'.','Color',[.6 0 .6],'Markersize',10)
hold on
plot (-v2min1(200:end),-v2min(200:end),'.','Color',[.6 0 .6],'Markersize',10)
hold on
plot ([1.02*min(-v1min1),.98*max(-v1min1)],[1.02*min(-v1min1),.98*max(-v1min1)],'r')
hold on
axis tight 

 xlabel('V(n)','Fontsize', 16),
  ylabel('V(n+1)','Fontsize', 16)


  Ca1min1=circshift(Ca1min,1);
  Ca2min1=circshift(Ca2min,1);
figure(8)
clf

for i =1:length(Ca1min)-2
plot ([-Ca1min1(i),-Ca1min(i)],[-Ca1min(i),-Ca1min(i)],'Color',[0.2 0.2 0.2],'Linewidth',1)
hold on
plot ([-Ca1min(i),-Ca1min(i)],[-Ca1min1(i+1),-Ca1min(i+1)],'Color',[0.2 0.2 0.2],'Linewidth',1)
hold on
end

%plot (-Ca1min1,-Ca1min,'Color',[0.7 0.7 0.7])
hold on
plot (-Ca1min1,-Ca1min,'.','Color','b','Markersize',10)
hold on


plot (-Ca1min1(1),-Ca1min(1),'.','Color','r','Markersize',30)
hold on
plot (-Ca2min1,-Ca2min,'.','Color','g','Markersize',10)
hold on
plot (-Ca1min1(200:end),-Ca1min(200:end),'.','Color',[.6 0 .6],'Markersize',10)
hold on
plot (-Ca2min1(200:end),-Ca2min(200:end),'.','Color',[.6 0 .6],'Markersize',10)
hold on
plot ([.95*min(-Ca1min1),1.05*max(-Ca1min1)],[.95*min(-Ca1min1),1.05*max(-Ca1min1)],'r')
hold on
axis tight

xlabel('Ca(n)','Fontsize', 16),
ylabel('Ca(n+1)','Fontsize', 16)
 
%%%%%%%%%%%%%%%%%%%%%%%%
x1min1=circshift(x1min,1);
x2min1=circshift(x2min,1);
figure(9)
clf

%plot (-x1min1,-x1min,'Color',[0.7 0.7 0.7])

for i =1:length(x1min1)-2
plot ([-x1min1(i),-x1min(i)],[-x1min(i),-x1min(i)],'Color',[0.2 0.2 0.2],'Linewidth',1)
hold on
plot ([-x1min(i),-x1min(i)],[-x1min1(i+1),-x1min(i+1)],'Color',[0.2 0.2 0.2],'Linewidth',1)
hold on
end

plot (-x1min1(1),-x1min(1),'.','Color','r','Markersize',30)
hold on
plot (-x2min1,-x2min,'.','Color','g','Markersize',10)
hold on
plot (-x1min1(200:end),-x1min(200:end),'.','Color',[.6 0 .6],'Markersize',10)
hold on
plot (-x2min1(200:end),-x2min(200:end),'.','Color',[.6 0 .6],'Markersize',10)
hold on
plot ([.95*min(-x1min1),max(-x1min1)*1.05],[.95*min(-x1min1),max(-x1min1)*1.05],'r')

  xlabel('x(n)','Fontsize', 16),
  ylabel('x(n+1)','Fontsize', 16) 
  axis tight 

 
  figure(11)
  clf
  plot3(Caa1,xx1,vv1,'Color',[0. 0.0 0.8],'Linewidth',1.0)
  hold on
  tran=length(xx1)-622600;
  % plot3(Caa1(tran:end),xx1(tran:end),vv1(tran:end),'Color',[1 0.0 0],'Linewidth',2.0)
  hold on
  view(42,40)
  