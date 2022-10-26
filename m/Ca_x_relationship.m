
clear all
v = -68.5:.1:-27;
r = 127*v./105 + 8265/105;

%sodium current
gna = 4;
% fit with sigmoid --jack's simple approximation
%minfty = 1./(exp((-22-v)./8)+1);
% version that dr shilnikov was using.
amm=0.1*(50-r)./(exp((50-r)./10)-1);
bmm=4*exp((25-r)./18);
minfty = amm./(amm+bmm);

ex_h1 = .07*exp((25-r)./20);
ex_h2 = 1./(1+exp((55-r)./10));
h = ex_h1 ./ (ex_h1 +ex_h2);
INa = gna.*minfty.^3.*h.*(30-v);

%potassium
gn = .3;
ex_n1 = .01*(55-r)./(exp((55-r)./10)-1);
ex_n2 = 0.125*exp((45 - r)./80);
n= ex_n1 ./ (ex_n1 +ex_n2);
Ik = gn.*n.^4.0.*(-75 -v);

%leak
gleak = .003;
Il = gleak.*(-40-v);

%hcurrent
gh=0.0000;
Vh=-50;
y1=1./(1+exp(10*(v-Vh)));
Ih=gh./((1.+exp(-(v+63.)/7.8))).^3.*y1.*(v-70);

% find x nullcline

xs = -4;
x = 1./(exp(.15.*(-50 + xs -v))+1);
Ix = .01.*x.*(30-v);

% v' without ca
q = INa+Ik+Ix+Il+Ih;
c = (-.015.*v-1.125)./(-.03.*v+q-2.25)-1./2;

figure(900)
clf
plot(c,x,'LineWidth',4)
hold on

% find C nullcline
cs = -30;

% previous Ca nullcline
Ca = .0085.*x.*(140-v+cs); % assumes x is constant
%ICa = 0.03*Ca./(0.5+Ca).*(v+75);
%qq = INa+Ik+Il+Ih+ICa;
%x11=-qq./(.01*(30-v));
%Ca1=0.0085*x11.*(140-v+cs);
%plot(Ca1, x11, 'r')
%hold on,

y = INa+ Ik +Il+ Ih;

ex_c1 = sqrt( ...
  (51.*cs.*v + 425.*cs.*(4.*y + 9) + 51.*v.^2 + 5.*v.*(340.*y - 863) - ...
  500.*(476.*y + 1011)).^2 + 6800000.*(v - 30).*y.*(cs + v - 140));
ex_c2 = 51.*cs.*v + 1700.*cs.*y +3825.*cs+ 51.*v.^2 + 1700.*v.*y - 4315.*v - ...
  238000.*y - 505500;


c1 = real(-ex_c1+ex_c2)./(4000.*(v - 30));
c2 = real(ex_c1+ex_c2)./(4000.*(v - 30));

x1 = c1./.0085./(140 - v + cs);
plot(c1, x1)
hold on

x2 = c2./.0085./(140 - v + cs);
plot(c2, x2)

axis([-0.2 2.6 -0.1 1.])
hold on


load('ca_sweep_ih_0.mat')
% %  % low stable branch 
Vvar1= x(1,1:end);
nvar= x(2,1:end); 
xvar1= x(5,1:end);
Cavar1=x(6,1:end);
plot(Cavar1, xvar1,'Color',[0.8 0.8 0.8, .5],'LineWidth',2)
hold on
load('eq1.mat')
Vvar2= x(1,1:end);
nvar= x(2,1:end); 
xvar2= x(5,1:end);
Cavar2=x(6,1:end);
plot(Cavar2,xvar2,'Color',[0.8 0.8 0.8, .5],'LineWidth',2)
hold on
load('eq2.mat')
Vvar3= x(1,1:end);
nvar= x(2,1:end); 
xvar3= x(5,1:end);
Cavar3=x(6,1:end);
plot(Cavar3,xvar3,'Color',[0.8 0.8 0.8, .5],'LineWidth',2)
hold on


  
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
plot(Ca_cod_bkgh0,x_cod_bkgh0+0.015,'-.','Color',[.85 .85 .85 ],'LineWidth',2)
hold on
plot(Ca_cod_fdgh0,x_cod_fdgh0+0.015,'-.','Color',[.85 .85 .85 ],'LineWidth',2)
hold on
  
axis([-0.2 1.6 -0.1 1.2])