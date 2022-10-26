%--------------------- SNIC whole system gh = 0 -------------------------------
data1 = load('LP_whole_bk_gh0.mat','x');
cod_back = data1.x;
data1 = load('LP_whole_fwd_gh0.mat','x');
cod_fwd = data1.x;
x_cod_bkgh0 = cod_back(end-3,:);
Ca_cod_bkgh0 = cod_back(end-2,:);
x_cod_fdgh0 = cod_fwd(end-3,:);
Ca_cod_fdgh0 = cod_fwd(end-2,:);

%------------------- plot SNIC ---------------------
figure
plot(Ca_cod_bkgh0,x_cod_bkgh0,'lineStyle','-')
hold on

plot(Ca_cod_fdgh0,x_cod_fdgh0,'lineStyle','-')
hold on