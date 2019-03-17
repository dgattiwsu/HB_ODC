%% HB ODC KINETIC METHOD

% Copyright (c) 2018, Domenico L. Gatti, Nazzareno Capitanio
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
clear, clc, close all

% Factor to convert O2 concentration to partial pressure under the current
% conditions (37 degC, ~23g/Kg salinity). This factor will be different 
% for solutions with different buffers or at different temperatures.
cO2_to_pO2 = 0.8157;

% Amount of the enzyme (E) cytochrome oxidase present in the added
% mitochondria.
E_init = 0.05;

%% CSV filename
csv_filename = 'O2_pc_hb';

%% Declare global variables
global reox_start_time k1 km1 k2 ko1 kom1 ko2 kom2 ...
       ko3 kom3 ko4 kom4 dk_ref ks1 ks2 ks3 ks4 sorc AIR

%% Load data and parameters from the graphic method run 
data_struct = load([csv_filename '_GRAPHIC_ODC.mat'],...
    'O2_pc','pc_vector_deox_reox','time_vector_deox_reox',...
    'O2_equil','total_Hb','K','deox_start_time', 'linear_start_ind', ...
    'linear_end_ind','deox_start_ind','reox_start_time','reox_start_ind', ...
    'deox_end_time','deox_end_ind','dk_in','dk_out');

O2_pc = data_struct.O2_pc;
O2_equil = data_struct.O2_equil
total_Hb = data_struct.total_Hb
K_graphic = data_struct.K

% Here we plot the original data.
% Kinetic_Method_fig_1 = figure;
% set(gcf,'Unit','Normalized','Position',[0.1 0.4 0.8 0.6]);
% plot(O2_pc(:,1),O2_pc(:,2),'-b','LineWidth',2); hold on
% xlim([-50 1900]);ylim([-20 200]); 
% xlabel('Time (s)');ylabel('[O2] (\muM)'); grid on

% Here we select only the part of the data beyond the initial equilibration
% by choosing a starting index, and we reset the time vector at 0 starting
% from this index. We also read in the start time of the reoxygenation
% phase from the graphic method run.
deox_start_ind = data_struct.deox_start_ind; 
deox_start_time = data_struct.deox_start_time; 
deox_end_ind = data_struct.deox_end_ind;
deox_end_time = data_struct.deox_end_time; 
linear_start_ind = data_struct.linear_start_ind;
linear_end_ind = data_struct.linear_end_ind;
reox_start_ind = data_struct.reox_start_ind; 
reox_start_time = data_struct.reox_start_time;
Product_vector = data_struct.pc_vector_deox_reox;
Time_vector = data_struct.time_vector_deox_reox;

% Kinetic_Method_fig_2 = figure;
% set(gcf,'Unit','Normalized','Position',[0.1 0.4 0.8 0.6]);
% plot(Time_vector,Product_vector,'-b','LineWidth',2);
% vline(reox_start_time,'g'); 
% xlim([-50 1600]);ylim([-20 200]); 
% xlabel('Time (s)');ylabel('[O2] (\muM)'); grid on

%% Set initial parameters

% All 2nd order rate constants
sorc = 100;

% Reference human Hb (Imai)
%     K_ref(1) = 0.0188;
%     K_ref(2) = 0.0566;
%     K_ref(3) = 0.407;
%     K_ref(4) = 4.28;
    
% Recommended reference human Hb for this method.
    K_ref(1) = 0.0273;
    K_ref(2) = 0.0036;
    K_ref(3) = 0.3343;
    K_ref(4) = 0.0667;

% Initial guess for U variables    
    U(1) = 4.8;         % k2
    U(2) = total_Hb;
    
% Here we use the reference values from above 
    U(3) = sorc./K_ref(4);
    U(4) = sorc./K_ref(3);
    U(5) = sorc./K_ref(2);
    U(6) = sorc./K_ref(1);

% Here we use the values from the graphic method
%     U(3) = sorc./data_struct.K(4);
%     U(4) = sorc./data_struct.K(3);
%     U(5) = sorc./data_struct.K(2);
%     U(6) = sorc./data_struct.K(1);

% Recommended initialization for this method: harmonic mean on K's
% between the reference values and the graphic method values.
%     U(3) = (sorc/data_struct.K(4) + sorc/K_ref(4))/2;
%     U(4) = (sorc/data_struct.K(3) + sorc/K_ref(3))/2;
%     U(5) = (sorc/data_struct.K(2) + sorc/K_ref(2))/2;
%     U(6) = (sorc/data_struct.K(1) + sorc/K_ref(1))/2;

%     U(7) = 0.012; % dk_ref
%     U(8) = U(7);  % dk_out

    U(7) = data_struct.dk_in; % dk_ref
    
% Here we define additional four variables that allow shifts in rate
% constants between deoxygenation and reoxygenation phase
    U(8) = 0;
    U(9) = 0;
    U(10) = 0;
    U(11) = 0;
    
% Here we include in the variables also km1
    U(12) = 0.1*sorc;    %   1/second

% Optionally, shake thing up a little
U_std = std(U);
scale = 10^(-8);
U = abs(scale*U_std.*randn(1,length(U)) + U);

% U_struct = load('HB_ODC_KINETIC_METHOD_1.mat','U');
% U = U_struct.U';

U_orig = U;

%% Kinetic model

%    Index:    Reaction:
%    1         E + O <-> EO         % k1 , km1
%    2         EO -> W + E          % k2
%    3         HO4 <-> O + HO3      % ko1 , kom1
%    4         HO3 <-> O + HO2      % ko2 , kom2
%    5         HO2 <-> O + HO       % ko3 , kom3
%    6         HO <-> O + H         % ko4 , kom4
%    7         AIR <-> O
   
k1 = sorc;              %   1/(micromolarity*second)
% km1     =  0.1*sorc;    %   1/second
kom1    =  sorc;        %   1/(micromolarity*second)
kom2    =  sorc;        %   1/(micromolarity*second)
kom3    =  sorc;        %   1/(micromolarity*second)
kom4    =  sorc;        %   1/(micromolarity*second)

k2 = U(1);
ko1 = U(3);    
ko2 = U(4);
ko3 = U(5);
ko4 = U(6);
dk_ref = U(7);
ks1 = U(8);
ks2 = U(9);
ks3 = U(10);
ks4 = U(11);
km1 = U(12);

%

 E = E_init;
 O = Product_vector(1);
 EO = 0;
 HO4 = U(2);
 HO3 = 0;
 HO2 = 0;
 HO = 0;
 H = 0;
 AIR = Product_vector(end);
 
Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H  AIR]';

%% Find most likely Hb concentration
% % EU = HO4;
% % options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',80);
% % fun = @(EU) simul_4(EU,U,Time_vector,Product_vector,E_init);
% % [EU2] = fminsearch(fun,EU,options);
% % HO4 = EU2
% % 
% % Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H  AIR]';

%% Find most likely enzyme concentration
EU = E_init;
options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',60);
fun = @(EU) simul_3(EU,U,Time_vector(linear_start_ind:linear_end_ind), ...
    Product_vector(linear_start_ind:linear_end_ind));
[EU2] = fminsearch(fun,EU,options);
E_init = EU2;
E = EU2

Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H  AIR]';

%% Find most likely enzyme and Hb concentration
% % EU(1) = E_init;
% % EU(2) = HO4;
% % options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',80);
% % fun = @(EU) simul_4(EU,U,Time_vector,Product_vector);
% % [EU2] = fminsearch(fun,EU,options);
% % E = EU2(1)
% % HO4 = EU2(2)
% % 
% % Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H  AIR]';

%% Progress curve
options = odeset('AbsTol',1e-8);
[~,Product_matrix] = ode15s(@reaction,Time_vector,Reactant_vector,options);

% SSE
Trial_vector = Product_matrix(:,2);
g = Trial_vector - Product_vector;
sos1 = g'*g;
sse = sos1

% Rsquared 
c_Product_vector = Product_vector - mean(Product_vector);
sst = c_Product_vector'*c_Product_vector;
rsquare = 1-sse/sst

%% Here we plot the original and the calculated O2 traces.
O2_Trace = figure;
set(gcf,'Unit','Normalized','Position',[0.1 0.4 0.8 0.6]);
plot(Time_vector,Product_vector,'-b','LineWidth',2);hold on
plot(Time_vector,Product_matrix(:,2),'-r','LineWidth',2);hold on

xlim([-50 1900]);ylim([-20 200]); 
xlabel('Time (s)');ylabel('[O2] (\muM)'); grid on

%% Saturation plot with current parameters
Hb_saturation_1 = figure;

% noO2_ind = find(Product_matrix(:,2)<1e-3);
% 
% if isempty(noO2_ind)
%     phi_deox_last = deox_end_ind;
%     phi_reox_first = reox_start_ind;
% else
%     phi_deox_last = noO2_ind(1)-1;
%     phi_reox_first = noO2_ind(end)+1;
% end

phi_deox_last = deox_end_ind;
phi_reox_first = reox_start_ind;

phi_deox_ind = 2:phi_deox_last;
phi_reox_ind = phi_reox_first:size(Product_matrix,1);

phi_deox = (Product_matrix(phi_deox_ind,4)*4+Product_matrix(phi_deox_ind,5)*3+...
    Product_matrix(phi_deox_ind,6)*2+Product_matrix(phi_deox_ind,7))/(Product_matrix(1,4)*4);
plot(Product_matrix(phi_deox_ind,2)*(cO2_to_pO2),phi_deox,'.b','LineWidth',1.5);hold on

phi_reox = (Product_matrix(phi_reox_ind,4)*4+Product_matrix(phi_reox_ind,5)*3+...
    Product_matrix(phi_reox_ind,6)*2+Product_matrix(phi_reox_ind,7))/(Product_matrix(1,4)*4);
plot(Product_matrix(phi_reox_ind,2)*(cO2_to_pO2),phi_reox,'.r','LineWidth',1.5);hold on
legend('Deoxygenation ODC','Reoxygenation ODC','Location','Best');grid on
xlim([0 160])
ylim([-0.02 1.05])
xlabel('P_O_2 (mmHg)')
ylabel('Fractional Saturation')

%% Parameter refinement
nepoch = 3;
nvars = length(U);
U_store = zeros(3,nvars);

for epoch = 1:nepoch
    
% FMINSEARCH
weights = 0.1*ones(4,1);
options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',80);
fun = @(U) simul_2(U,Time_vector,Product_vector,E_init,phi_deox_ind,phi_reox_ind,weights,'sse','do_not_use_K','keep_central');
[U2] = fminsearch(fun,U,options);
U = U2;
U_store(1,:) = U;

% FMINSEARCH
weights = 0.5*ones(4,1);
options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',80);
fun = @(U) simul_2(U,Time_vector,Product_vector,E_init,phi_deox_ind,phi_reox_ind,weights,'sse','do_not_use_K','keep_central');
[U2] = fminsearch(fun,U,options);
U = U2;
U_store(2,:) = U;

% FMINSEARCH
weights = 1*ones(4,1);
options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',80);
fun = @(U) simul_2(U,Time_vector,Product_vector,E_init,phi_deox_ind,phi_reox_ind,weights,'sse','do_not_use_K','keep_central');
[U2,fval] = fminsearch(fun,U,options);
U = U2;
U_store(3,:) = U;

fval_old = fval;
scale = 10^(-8);
golden_ratio = (1 + sqrt(5))/2;

for i = 1:20
weights = 1*ones(4,1);
options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',80);

U_std = std(U_store);
U = abs(scale*U_std.*randn(1,nvars) + U);

% FMINSEARCH
fun = @(U) simul_2(U,Time_vector,Product_vector,E_init,phi_deox_ind,phi_reox_ind,weights,'sse','do_not_use_K','keep_central');
[U2,fval] = fminsearch(fun,U,options);
U = U2;

if fval<=fval_old
scale = scale/(2*golden_ratio);
fval_old = fval;
% U = U2;
end

U_store(1:2,:) = U_store(2:3,:);
U_store(3,:) = U;
end

% LSQNONLIN
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
    'Display','iter','FiniteDifferenceType','central','FunctionTolerance',1e-9);
fun = @(U) simul_2(U,Time_vector,Product_vector,E,phi_deox_ind,phi_reox_ind,weights,'res','do_not_use_K','keep_central');
[U2,resnorm,~,~,~,~,J] = ...
    lsqnonlin(fun,U,[],[],options);
U = U2;

end

%% Mean Squared Error and the Standard Deviation 
nobs = length(Time_vector);
sse = resnorm;
nvars = length(U);
MSE = sse/(nobs-nvars);
RMSE = sqrt(MSE);

% Covariance matrix of the variables from the last value
% of the Jacobian:  J2'*J2 is the 'precision' matrix, its inverse is thus
% the 'covariance' matrix!
J = full(J);
cov_U = inv(J'*J).*MSE;

% Here we fully symmetrize to avoid numerical errors.
cov_U = (cov_U + cov_U')/2

% The covariance matrix can be converted into a correlation matrix; also
% the diagonal of the covariance matrix gives the variances and therefore
% the sigmas.
[corr_U,sigma] = corrcov(cov_U);

% sigma = sqrt(diag(cov_U));

% The 95% confidence intervals is obtained by subtracting or adding
% 1.96*sigma to the 'expected' values of U.
Conf_95 = [U'-1.96*sigma U'+1.96*sigma]

% Rsquared: 
c_Product_vector = Product_vector - mean(Product_vector);
sst = c_Product_vector'*c_Product_vector;
rsquare = 1-sse/sst

%% Progress curve fit and ODC

k2 = U(1);
ko1 = U(3);    
ko2 = U(4);
ko3 = U(5);
ko4 = U(6);
dk_ref = U(7);
ks1 = U(8);
ks2 = U(9);
ks3 = U(10);
ks4 = U(11);
km1 = U(12);

%

 E = E_init;
 O = Product_vector(1);
 EO = 0;
 HO4 = U(2);
 HO3 = 0;
 HO2 = 0;
 HO = 0;
 H = 0;
 % AIR = Product_vector(1);
 
Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H  AIR]';

% Progress curves
options = odeset('AbsTol',1e-8);
[~,Product_matrix] = ode15s(@reaction,Time_vector,Reactant_vector,options);

% SSE without central region
phi_deox_reox_ind = [phi_deox_ind phi_reox_ind];
Trial_vector = Product_matrix(:,2);
g = Trial_vector(phi_deox_reox_ind) - Product_vector(phi_deox_reox_ind);
sos1 = g'*g;
sse = sos1

% Rsquared without central region
c_Product_vector = Product_vector(phi_deox_reox_ind) - mean(Product_vector(phi_deox_reox_ind));
sst = c_Product_vector'*c_Product_vector;
rsquare = 1-sse/sst

% SSE with central region
Trial_vector = Product_matrix(:,2);
g = Trial_vector - Product_vector;
sos1 = g'*g;
sse = sos1

% Rsquared with central region
c_Product_vector = Product_vector - mean(Product_vector);
sst = c_Product_vector'*c_Product_vector;
rsquare = 1-sse/sst

%% Here we plot the original and the calculated O2 traces.
O2_Trace = figure;
set(gcf,'Unit','Normalized','Position',[0.1 0.4 0.8 0.6]);
plot(Time_vector,Product_vector,'-b','LineWidth',2);hold on
plot(Time_vector,Product_matrix(:,2),'-r','LineWidth',2);hold on

xlim([-50 1900]);ylim([-20 200]); 
xlabel('Time (s)');ylabel('[O2] (\muM)'); grid on

%% Saturation plot with refined parameters
Hb_saturation_1 = figure;
% noO2_ind = find(Product_matrix(:,2)<1e-3);
% phi_deox_last = noO2_ind(1)-1;
% phi_reox_first = noO2_ind(end)+1;
% Start at 2 the phi_deox_ind to avoid a point at 100 % saturation that is
% practically not achievable.
phi_deox_ind = 2:phi_deox_last;
phi_reox_ind = phi_reox_first:size(Product_matrix,1);

phi_deox = (Product_matrix(phi_deox_ind,4)*4+Product_matrix(phi_deox_ind,5)*3+...
    Product_matrix(phi_deox_ind,6)*2+Product_matrix(phi_deox_ind,7))/(Product_matrix(1,4)*4);
plot(Product_matrix(phi_deox_ind,2)*(cO2_to_pO2),phi_deox,'.b','LineWidth',1.5);hold on

phi_reox = (Product_matrix(phi_reox_ind,4)*4+Product_matrix(phi_reox_ind,5)*3+...
    Product_matrix(phi_reox_ind,6)*2+Product_matrix(phi_reox_ind,7))/(Product_matrix(1,4)*4);
plot(Product_matrix(phi_reox_ind,2)*(cO2_to_pO2),phi_reox,'.r','LineWidth',1.5);hold on
legend('Deoxygenation ODC','Reoxygenation ODC','Location','Best');grid on
xlim([0 160])
ylim([-0.02 1.05])
xlabel('P_O_2 (mmHg)')
ylabel('Fractional Saturation')

%% Saturation Plot versus partial pressure (mmHg)
Hb_saturation_2 = figure;
set(gcf,'Unit','Normalized','Position',[0 0.4 0.4 0.55])
plot(Product_matrix(phi_deox_ind,2)*(cO2_to_pO2),phi_deox,'.b','MarkerSize',0.2);hold on
plot(Product_matrix(phi_reox_ind,2)*(cO2_to_pO2),phi_reox,'or',...
             'MarkerEdgeColor','r',...
             'MarkerFaceColor','c',...
             'MarkerSize',5); hold on

legend('Deoxygenation ODC','Reoxygenation ODC','Location','Best');grid on
xlim([0 160])
ylim([-0.02 1.05])
xlabel('P_O_2 (mmHg)')
ylabel('Fractional Saturation')

%% Simulation plot versus [O2] micromolar
% Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H  AIR]';

Simulation_plot_1 = figure;
set(gcf,'Unit','Normalized','Position',[0 0.4 0.5 0.6]);
h(1) = plot(Time_vector,Product_vector,'-b','LineWidth',1);hold on
h(2) = plot(Time_vector,Product_matrix(:,2),'--r','LineWidth',2.5);
h(3) = plot(Time_vector,Product_matrix(:,4),'LineStyle','-','Color',[210/255 105/255 30/255],'LineWidth',1.5);
h(4) = plot(Time_vector,Product_matrix(:,5),'LineStyle','-','Color',[1 165/255 0],'LineWidth',1.5);
h(5) = plot(Time_vector,Product_matrix(:,6),'-m','LineWidth',1.5);
h(6) = plot(Time_vector,Product_matrix(:,7),'-g','LineWidth',1.5);
h(7) = plot(Time_vector,Product_matrix(:,8),'LineStyle','-','Color',[30 154 245]/255,'LineWidth',1.5);
grid on
legend('[O2] observed','[O2] fit','[Hb(O2)_4]','[Hb(O2)_3]','[Hb(O2)_2]','[Hb(O2)_1]','[Hb]','Location','Best');
xlabel('Time (s)');
ylabel('[O2] (\muM)');
ylim([-10 190]);

%% Results
% Here we write the results to external files
FIT = [Time_vector Product_vector Product_matrix(:,[2 4:8])];
FIT_table = array2table(FIT,'VariableNames',{'Time','O2_obs','O2_fit','Hb_O2_4','Hb_O2_3','Hb_O2_2','Hb_O2_1','Hb'});
writetable(FIT_table,[csv_filename '_KINETIC_fit.xlsx'],'sheet',1);
PHI = [Product_matrix([phi_deox_ind phi_reox_ind],2) [phi_deox; phi_reox]];
PHI_table = array2table(PHI,'VariableNames',{'O2_conc','Hb_saturation'});
writetable(PHI_table,[csv_filename '_KINETIC_fit.xlsx'],'sheet',2);
PAR = [ U' Conf_95 ];
RowNames = {'kcat_cox';'HbO4';'koff_HbO4';'koff_HbO3';'koff_HbO2';'koff_HbO1';'kdiff_O2';...
    'koff_HbO4_shift';'koff_HbO3_shift';'koff_HbO2_shift';'koff_HbO1_shift';'koff_O2_cox'}
PAR_table = array2table(PAR,'VariableNames',...
    {'refined_parameters','lower_bound','upper_bound'},'RowNames',RowNames);
writetable(PAR_table,[csv_filename '_KINETIC_fit.xlsx'],'sheet',3,'WriteRowNames',true);

%% Equilibrium constants from deox ODC using Adair equation
K_deox_k = sorc./U(6:-1:3);

O2_conc = Product_matrix(phi_deox_ind,2);
fun2 = @(K) simul_1(K,phi_deox,O2_conc);

options = optimoptions('lsqnonlin','FiniteDifferenceType','central',...
    'MaxFunctionEvaluations',40000,'MaxIterations',40000)
[K_deox_adair] = ...
    lsqnonlin(fun2,K_deox_k,[0 0 0 0],[],options);

K_deox_comp = [K_deox_k;K_deox_adair];

% Adair derivation of P50
K = K_deox_adair;

B_num = @(x) K(1)*x + 2*K(1)*K(2)*x.^2 + 3*K(1)*K(2)*K(3)*x.^3 ...
    + 4*K(1)*K(2)*K(3)*K(4)*x.^4;
B_den = @(x) 1 + K(1)*x + K(1)*K(2)*x.^2 + K(1)*K(2)*K(3)*x.^3 ...
    + K(1)*K(2)*K(3)*K(4)*x.^4;
B_ratio = @(x) B_num(x)./(4*B_den(x));

% Saturation at different O2 concentrations
B = B_ratio(O2_conc);

figure;
plot(O2_conc,phi_deox,'-b',O2_conc,B,'-r')

% O2 concentration at 50% saturation
fun3 = @(x) 0.5 - B_ratio(x);

[c50_deox] = lsqnonlin(fun3,20,0,[],options)
p50_deox = cO2_to_pO2*c50_deox

% P50 corrected for CO2
p50_CO2_40_deox = p50_deox + 0.1*(40-0.3) +1.083*10^(-4)*(40-0.3)^2

% Total Hb concentration
Tot_Hb_kin_model = U(2)*4

% Upper limit of the correction for 2,3-DPG binding
% d_23_dpg = 795.63*(3.21*10^(-3)) - 19660*((3.21*10^(-3))^2)

% Equilibrium constants
K_deox_adair

%% Hill plot deox 
% Reconstructed data using Adair equation with refined equilibrium
% constants
HILL_DEOX = figure;
set(gcf,'Unit','Normalized','Position',[0 0.1 0.4 0.6])
% Here we reconstruct the saturation curve using Adair equation
Unbound = 0:0.2:180;
FR1 = B_ratio(Unbound);
X = log(Unbound);
Y = log(FR1./(1 - FR1));
plot(X,Y,'o',...
             'MarkerEdgeColor','b',...
             'MarkerFaceColor','c',...
             'MarkerSize',5);
         
ylim([-8 8])
xlim([-2 5.5])

hold on

f = fittype('a*x + b');

X1 = X(760:end)'; Y1 = Y(760:end)'; 
[Hill_1,GOF_1] = fit(X1,Y1,f,'StartPoint',[2 -8],'Robust','on');
Hill_1
hh1 = plot(Hill_1,'-r');
hh1.LineWidth = 1.5;

X2 = X(135:235)'; Y2 = Y(135:235)';
[Hill_2,GOF_2] = fit(X2,Y2,f,'StartPoint',[1 1],'Robust','on');
Hill_2
hh2 =plot(Hill_2,'-g');
hh2.LineWidth = 1.5;

X3 = X(2:30)'; Y3 = Y(2:30)';
[Hill_3,GOF_3] = fit(X3,Y3,f,'StartPoint',[1 1],'Robust','on');
Hill_3
hh3 = plot(Hill_3,'-b');
hh3.LineWidth = 1.5;
         
xlabel('log[Free O_2]');
ylabel('log(Frac. Sat./(1-Frac. Sat.)');
title('Hill plot (deox phase)');

legend('Hill Plot','Fit1: Kd\_1','Fit2: Kd\_2',...
    'Fit3: Kd\_3','Location','NorthWest');

% Retrieve the fit parameters.
Hill_1_params = coeffvalues(Hill_1);
Hill_2_params = coeffvalues(Hill_2);
Hill_3_params = coeffvalues(Hill_3);
nH_1 = Hill_1_params(1);
nH_2 = Hill_2_params(1);
nH_3 = Hill_3_params(1);

Kd_1 = exp(-Hill_1_params(2));
Kd_2 = exp(-Hill_2_params(2));
Kd_3 = exp(-Hill_3_params(2));

string1 = 'nH\_1 = ';
string2 = num2str(nH_1,'%6.3e\n');
string3 = 'nH\_2 = ';
string4 = num2str(nH_2,'%6.3e\n');
string5 = 'nH\_3 = ';
string6 = num2str(nH_3,'%6.3e\n');

% Create textbox
annotation(HILL_DEOX,'textbox',...
    [0.65 0.4 0.18 0.1],...
    'String',{[string1 string2 ; string3 string4 ; string5 string6]},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 0.800000011920929],...
    'Color',[1 0 0]);

string0 = ' microM';
string1 = 'Kd\_1 = ';
string2 = num2str(Kd_1,'%6.2e\n');
string3 = 'Kd\_2 = ';
string4 = num2str(Kd_2,'%6.2e\n');
string5 = 'Kd\_3 = ';
string6 = num2str(Kd_3,'%6.2e\n');

% Create textbox
annotation(HILL_DEOX,'textbox',...
    [0.60 0.2 0.23 0.1],...
    'String',{[string1 string2 string0; string3 string4 string0; string5 string6 string0]},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 0.800000011920929],...
    'Color',[1 0 0]);

xlabel('log(Free[O_2] \muM)');
ylabel('log(Frac. Sat./(1-Frac. Sat.)');
grid on

%% Equilibrium constants from reox ODC using Adair equation
K_reox_k = sorc./(U(6:-1:3) + U(11:-1:8));

O2_conc = Product_matrix(phi_reox_ind,2);
fun2 = @(K) simul_1(K,phi_reox,O2_conc);

options = optimoptions('lsqnonlin','FiniteDifferenceType','central',...
    'MaxFunctionEvaluations',40000,'MaxIterations',40000)
[K_reox_adair] = ...
    lsqnonlin(fun2,K_reox_k,[0 0 0 0],[],options);

K_reox_comp = [K_reox_k;K_reox_adair];

% Adair derivation of P50
K = K_reox_adair;

B_num = @(x) K(1)*x + 2*K(1)*K(2)*x.^2 + 3*K(1)*K(2)*K(3)*x.^3 ...
    + 4*K(1)*K(2)*K(3)*K(4)*x.^4;
B_den = @(x) 1 + K(1)*x + K(1)*K(2)*x.^2 + K(1)*K(2)*K(3)*x.^3 ...
    + K(1)*K(2)*K(3)*K(4)*x.^4;;
B_ratio = @(x) B_num(x)./(4*B_den(x));

% Saturation at different O2 concentrations
B = B_ratio(O2_conc);

figure;
plot(O2_conc,phi_reox,'-b',O2_conc,B,'-r')

% O2 concentration at 50% saturation
fun3 = @(x) 0.5 - B_ratio(x);

[c50_reox] = lsqnonlin(fun3,20,0,[],options)
p50_reox = cO2_to_pO2*c50_reox

% P50 corrected for CO2
p50_CO2_40_reox = p50_reox + 0.1*(40-0.3) +1.083*10^(-4)*(40-0.3)^2

% Total Hb concentration
Tot_Hb_kin_model = U(2)*4

% Upper limit of the correction for 2,3-DPG binding
% d_23_dpg = 795.63*(3.21*10^(-3)) - 19660*((3.21*10^(-3))^2)

% Equilibrium constants
K_reox_adair

%% Hill plot reox 
% Reconstructed data using Adair equation with refined equilibrium
% constants
HILL_REOX = figure;
set(gcf,'Unit','Normalized','Position',[0 0.1 0.4 0.6])
% Here we reconstruct the saturation curve using Adair equation
Unbound = 0:0.2:180;
FR1 = B_ratio(Unbound);
X = log(Unbound);
Y = log(FR1./(1 - FR1));
plot(X,Y,'o',...
             'MarkerEdgeColor','b',...
             'MarkerFaceColor','c',...
             'MarkerSize',5);
         
ylim([-8 8])
xlim([-2 5.5])

hold on

f = fittype('a*x + b');

X1 = X(760:end)'; Y1 = Y(760:end)'; 
[Hill_1,GOF_1] = fit(X1,Y1,f,'StartPoint',[2 -8],'Robust','on');
Hill_1
hh1 = plot(Hill_1,'-r');
hh1.LineWidth = 1.5;

X2 = X(135:235)'; Y2 = Y(135:235)';
[Hill_2,GOF_2] = fit(X2,Y2,f,'StartPoint',[1 1],'Robust','on');
Hill_2
hh2 =plot(Hill_2,'-g');
hh2.LineWidth = 1.5;

X3 = X(2:40)'; Y3 = Y(2:40)';
[Hill_3,GOF_3] = fit(X3,Y3,f,'StartPoint',[1 1],'Robust','on');
Hill_3
hh3 = plot(Hill_3,'-b');
hh3.LineWidth = 1.5;
         
xlabel('log[Free O_2]');
ylabel('log(Frac. Sat./(1-Frac. Sat.)');
title('Hill plot (reox phase)');

legend('Hill Plot','Fit1: Kd\_1','Fit2: Kd\_2',...
    'Fit3: Kd\_3','Location','NorthWest');

% Retrieve the fit parameters.
Hill_1_params = coeffvalues(Hill_1);
Hill_2_params = coeffvalues(Hill_2);
Hill_3_params = coeffvalues(Hill_3);
nH_1 = Hill_1_params(1);
nH_2 = Hill_2_params(1);
nH_3 = Hill_3_params(1);
Kd_1 = exp(-Hill_1_params(2));
Kd_2 = exp(-Hill_2_params(2));
Kd_3 = exp(-Hill_3_params(2));

string1 = 'nH\_1 = ';
string2 = num2str(nH_1,'%6.3e\n');
string3 = 'nH\_2 = ';
string4 = num2str(nH_2,'%6.3e\n');
string5 = 'nH\_3 = ';
string6 = num2str(nH_3,'%6.3e\n');

% Create textbox
annotation(HILL_REOX,'textbox',...
    [0.65 0.4 0.18 0.1],...
    'String',{[string1 string2 ; string3 string4 ; string5 string6]},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 0.800000011920929],...
    'Color',[1 0 0]);

string0 = ' microM';
string1 = 'Kd\_1 = ';
string2 = num2str(Kd_1,'%6.2e\n');
string3 = 'Kd\_2 = ';
string4 = num2str(Kd_2,'%6.2e\n');
string5 = 'Kd\_3 = ';
string6 = num2str(Kd_3,'%6.2e\n');

% Create textbox
annotation(HILL_REOX,'textbox',...
    [0.60 0.2 0.23 0.1],...
    'String',{[string1 string2 string0; string3 string4 string0; string5 string6 string0]},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 0.800000011920929],...
    'Color',[1 0 0]);

xlabel('log(Free[O_2] \muM)');
ylabel('log(Frac. Sat./(1-Frac. Sat.)');
grid on

%% Equilibrium constants from combined deox-reox ODC using Adair equation
% K_from_rate_constants = sorc./U(6:-1:3);
K_from_rate_constants = mean([K_deox_k;K_reox_k]);

O2_conc = [Product_matrix(phi_deox_ind,2) ; Product_matrix(phi_reox_ind,2)];

fun2 = @(K) simul_1(K,[phi_deox ; phi_reox],O2_conc);
options = optimoptions('lsqnonlin','FiniteDifferenceType','central',...
    'MaxFunctionEvaluations',40000,'MaxIterations',40000)
[K_from_adair,resnorm,residual,exitflag,output,lambda,J] = ...
    lsqnonlin(fun2,K_from_rate_constants,[0 0 0 0],[],options);

% Adair derivation of P50
K = K_from_adair;

B_num = @(x) K(1)*x + 2*K(1)*K(2)*x.^2 + 3*K(1)*K(2)*K(3)*x.^3 ...
    + 4*K(1)*K(2)*K(3)*K(4)*x.^4;
B_den = @(x) 1 + K(1)*x + K(1)*K(2)*x.^2 + K(1)*K(2)*K(3)*x.^3 ...
    + K(1)*K(2)*K(3)*K(4)*x.^4;;
B_ratio = @(x) B_num(x)./(4*B_den(x));

% Saturation at different O2 concentrations
B = B_ratio(O2_conc);

% O2 concentration at 50% saturation
fun3 = @(x) 0.5 - B_ratio(x);

[c50,resnorm_c50,residual_c50,exitflag_c50,output_c50,lambda,J_c50] = ...
    lsqnonlin(fun3,20,[0],[],options);
p50 = cO2_to_pO2*c50

% P50 corrected for CO2
p50_CO2_40 = p50 + 0.1*(40-0.3) +1.083*10^(-4)*(40-0.3)^2

% Total Hb concentration
Tot_Hb_kin_model = U(2)*4

% Upper limit of the correction for 2,3-DPG binding
% d_23_dpg = 795.63*(3.21*10^(-3)) - 19660*((3.21*10^(-3))^2)

% Equilibrium constants
K_adair

%% Hill plot 
% Reconstructed data using Adair equation with refined equilibrium
% constants
HILL = figure;
set(gcf,'Unit','Normalized','Position',[0 0.1 0.4 0.6])
% Here we reconstruct the saturation curve using Adair equation
Unbound = 0:0.2:180;
FR1 = B_ratio(Unbound);
% X = log(Unbound*cO2_to_pO2);
X = log(Unbound);
Y = log(FR1./(1 - FR1));
plot(X,Y,'o',...
             'MarkerEdgeColor','b',...
             'MarkerFaceColor','c',...
             'MarkerSize',5);
         
ylim([-8 8])
xlim([-2 5.5])

hold on

f = fittype('a*x + b');

X1 = X(760:end)'; Y1 = Y(760:end)'; 
[Hill_1,GOF_1] = fit(X1,Y1,f,'StartPoint',[2 -8],'Robust','on');
Hill_1
hh1 = plot(Hill_1,'-r');
hh1.LineWidth = 1.5;

X2 = X(135:235)'; Y2 = Y(135:235)';
[Hill_2,GOF_2] = fit(X2,Y2,f,'StartPoint',[1 1],'Robust','on');
Hill_2
hh2 =plot(Hill_2,'-g');
hh2.LineWidth = 1.5;

X3 = X(2:30)'; Y3 = Y(2:30)';
[Hill_3,GOF_3] = fit(X3,Y3,f,'StartPoint',[1 1],'Robust','on');
Hill_3
hh3 = plot(Hill_3,'-b');
hh3.LineWidth = 1.5;
         
xlabel('log[Free O_2]');
ylabel('log(Frac. Sat./(1-Frac. Sat.)');
title('Hill plot (combined phases)');

% Retrieve the fit parameters.
Hill_1_params = coeffvalues(Hill_1);
Hill_2_params = coeffvalues(Hill_2);
Hill_3_params = coeffvalues(Hill_3);
nH_1 = Hill_1_params(1);
nH_2 = Hill_2_params(1);
nH_3 = Hill_3_params(1);
Kd_1 = exp(-Hill_1_params(2));
Kd_2 = exp(-Hill_2_params(2)/Hill_2_params(1));
Kd_3 = exp(-Hill_3_params(2));

cO2_hill = exp(-Hill_2_params(2)/Hill_2_params(1))
pO2_hill = cO2_to_pO2*cO2_hill

string1 = 'n_H_1 = ';
string2 = num2str(nH_1,'%6.2f\n');
string3 = 'n_H_2 = ';
string4 = num2str(nH_2,'%6.2f\n');
string5 = 'n_H_3 = ';
string6 = num2str(nH_3,'%6.2f\n');

% Create textbox
annotation(HILL,'textbox',...
    [0.65 0.4 0.18 0.1],...
    'String',{[string1 string2 ; string3 string4 ; string5 string6]},...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 0.800000011920929],...
    'Color',[1 0 0]);

string0 = ' \muM';
string00 = ' \muM ';
string1 = 'K_d_\__1 =  ';
string2 = num2str(Kd_1,'%6.2f\n');
string3 = 'K_D     =  ';
string4 = num2str(Kd_2,'%6.2f\n');
string5 = 'K_d_\__3 = ';
string6 = num2str(Kd_3,'%6.2f\n');


% Create textbox
annotation(HILL,'textbox',...
    [0.60 0.2 0.23 0.1],...
    'String',{[string1 string2 string0; string3 string4 string00; string5 string6 string0]},...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 0.800000011920929],...
    'Color',[1 0 0]);


xlabel('log(Free[O_2] \muM)');
ylabel('log(Frac. Sat./(1-Frac. Sat.)');
xlim([-1 5.5])
grid on

line([0 0],[-8 8],'Color','black','LineStyle',':','LineWidth',2)
line([-1 5.485],[0 0],'Color','black','LineStyle',':','LineWidth',2)

legend('Hill Plot','Fit1: K_d_\__1','Fit2: K_D',...
    'Fit3: K_d_\__3','Location','NorthWest');

%%
close all
save([csv_filename '_KINETIC_ODC.mat'])



