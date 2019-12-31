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

%% CSV filename

csv_filename = 'O2_pc_hb';

%% Read in and plot experimental O2 progress curve (O2_pc).
O2_pc = csvread([csv_filename '.csv']);
O2_pc(:,2) = O2_pc(:,2) - min(O2_pc(:,2))

Graphic_method_fig_1 = figure;
set(gcf,'Unit','Normalized','Position',[0.1 0.4 0.8 0.6]);
plot(O2_pc(:,1),O2_pc(:,2),'-b','LineWidth',2);
xlim([-50 1900]);ylim([-20 200]); 
xlabel('Time (s)');ylabel('[O2] (\muM)'); grid on

%% Progress curve: selection of the range used in fitting 

% Index of the equilibration end point and O2 concentration at equilibrium.
equil_end_time = 200;
equil_end_ind = find(O2_pc(:,1) == equil_end_time);
O2_equil = mean(O2_pc(1:equil_end_ind,2));

% Index of deoxygenation curve start point.
deox_start_time = 254;
deox_start_ind = find(O2_pc(:,1) == deox_start_time);

% Index of the beginning and end point of the linear part of the
% deoxygenation curve (as derived from Graphic_Method_fig_1.
linear_start_time = 290;
linear_start_ind = find(O2_pc(:,1) == linear_start_time);

% Use first time at which the concentration of oxygen falls at or below 100
% uM
O2_threshold_1 = 100;
linear_end_ind = find(O2_pc(:,2) <= O2_threshold_1,1)
linear_end_time = O2_pc(linear_end_ind,1);

% Index of the end point of the Hb deoxygenation
deox_end_time = 1208; % 618
deox_end_ind = find(O2_pc(:,1) == deox_end_time);

% Index of the reoxygenation start time
reox_start_time = 1344;
reox_start_ind = find(O2_pc(:,1) == reox_start_time);

% Indices of the start and end points of the progress curve, past the
% reoxygenation of Hb, that will be used to derive a fit.
O2_threshold_2 = 100;
ind1 = find(O2_pc(:,2) < O2_threshold_2,1,'last');
ind2 = find(ind1 > reox_start_ind,1);
reox_fit_start_ind = ind1(ind2);
reox_fit_start_time = O2_pc(reox_fit_start_ind,1);

% Index of the reoxygenation end time (also fit end time)
reox_end_time = 1812; 
reox_end_ind = find(O2_pc(:,1) == reox_end_time);
reox_fit_end_time = reox_end_time;
reox_fit_end_ind = reox_end_ind;

% Plot selected boundaries
figure(Graphic_method_fig_1);
vline(O2_pc(equil_end_ind,1),'k'); 
vline(O2_pc(deox_start_ind,1),'c'); 
vline(O2_pc(linear_start_ind,1),'r'); 
vline(O2_pc(linear_end_ind,1),'r'); 
vline(O2_pc(deox_end_ind,1),'c'); 
vline(O2_pc(reox_start_ind,1),'g'); 
vline(O2_pc(reox_fit_start_ind,1),'m'); 
vline(O2_pc(reox_fit_end_ind,1),'m'); 
vline(O2_pc(reox_end_ind,1),'g'); 

%% Progress curve used for fitting 
% Select the boundaries of the progress curve that will be used for fitting
% and reset times and indices.
time_shift = deox_start_time;
ind_shift = deox_start_ind -1;

pc_vector_deox_reox = O2_pc(deox_start_ind:reox_end_ind,2);
time_vector_deox_reox = O2_pc(deox_start_ind:reox_end_ind,1) - time_shift;

Graphic_method_fig_2 = figure;
set(gcf,'Unit','Normalized','Position',[0.1 0.4 0.8 0.6]);
plot(time_vector_deox_reox,pc_vector_deox_reox,'-b','LineWidth',2);
hold on
xlim([-50 1600]);ylim([-20 200]); 
xlabel('Time (s)');ylabel('[O2] (\muM)'); grid on

% Index of the end point of the linear part of the deoxygenation curve (as 
% derived from Graphic_Method_fig_2.

linear_start_time = linear_start_time - time_shift;
linear_start_ind = linear_start_ind - ind_shift;
linear_end_time = linear_end_time - time_shift;
linear_end_ind = linear_end_ind - ind_shift;
linear_time_vec = time_vector_deox_reox(linear_start_ind:linear_end_ind);
linear_pc_vec = pc_vector_deox_reox(linear_start_ind:linear_end_ind);

% Index of the end point of the Hb deoxygenation
deox_end_time = deox_end_time - time_shift;
deox_end_ind = deox_end_ind - ind_shift;

% Index of the reoxygenation start and end time
reox_start_time = reox_start_time - time_shift;
reox_start_ind = reox_start_ind - ind_shift;
reox_end_time = reox_end_time - time_shift;
reox_end_ind = reox_end_ind - ind_shift;

% Indices of the start and end points of the progress curve, past the
% reoxygenation of Hb, that will be used to derive an exponential fit.
reox_fit_start_time = reox_fit_start_time - time_shift;
reox_fit_start_ind = reox_fit_start_ind - ind_shift;
reox_fit_end_time = reox_fit_end_time - time_shift;
reox_fit_end_ind = reox_fit_end_ind - ind_shift;

% Plot boundaries of the linear and exponential fit
vline(time_vector_deox_reox(1),'c'); 
vline(time_vector_deox_reox(linear_start_ind),'r'); 
vline(time_vector_deox_reox(linear_end_ind),'r'); 
vline(time_vector_deox_reox(deox_end_ind),'c'); 
vline(time_vector_deox_reox(reox_start_ind),'g'); 
vline(time_vector_deox_reox(reox_fit_start_ind),'m'); 
vline(time_vector_deox_reox(reox_fit_end_ind),'m'); 
vline(time_vector_deox_reox(reox_end_ind),'g');

%% Deoxygenation phase fit

% Deoxygenation curve. 
time_vector_deox = time_vector_deox_reox(1:deox_end_ind);
pc_vector_deox = pc_vector_deox_reox(1:deox_end_ind);


% Design matrix for linear regression
lsq_mat = [ones(length(linear_time_vec),1) linear_time_vec];

% Least squares parameters
lsq_param = (lsq_mat'*lsq_mat)\(lsq_mat'*linear_pc_vec);

% Extrapolate O2 consumption until the end point of Hb deoxygantion
% lsq_mat_large = [ones(deox_end_ind,1) time_vector_deox_reox(1:deox_end_ind)];
lsq_mat_large = [ones(deox_end_ind,1) time_vector_deox];
lin_pc_vec_lsq = lsq_mat_large(1:deox_end_ind,:)*lsq_param;

% Plot the deoxygenation phase fit
Graphic_method_fig_3 = figure;
set(gcf,'Unit','Normalized','Position',[0.1 0.4 1 0.6]);
subplot(1,3,1)
plot(time_vector_deox,pc_vector_deox,...
    '-b','LineWidth',3);hold on
plot(time_vector_deox(linear_end_ind:end),lin_pc_vec_lsq(linear_end_ind:end),...
    '--r','LineWidth',1.5);

xlim([-50 deox_end_time + 100]);ylim([-80 200]);
legend('[O2] vs. time','Linear fit', 'Location', 'NorthEast' );

xlabel('Time (s)');ylabel('[O2] (\muM)'); grid on
vline(time_vector_deox(end),'g'); 
hline(0,':k');
hline(O2_threshold_1,':k');

% Take the negative of the difference to derive an upside Hb curve.
Hb_deox_curve =  - (pc_vector_deox_reox(1:deox_end_ind) - lin_pc_vec_lsq);
Hb_deox_curve = Hb_deox_curve - min(Hb_deox_curve);
[Hb_deox_curve_min,Hb_deox_curve_min_ind] = min(Hb_deox_curve);

% Hb concentration based on the deoxygenation curve.
Hb_deox_conc_scale = 1.0 % Alternatively use some fraction of 1, i.e., .999)
Hb_deox_conc = max(Hb_deox_curve)*Hb_deox_conc_scale; 

% Here we divide by the total amount of Hb to obtain the saturation
Hb_deox_sat = Hb_deox_curve/Hb_deox_conc;

%% Reoxygenation phase fit

% Reoxygenation curve. 
time_vector_reox = time_vector_deox_reox(reox_start_ind:reox_end_ind);
pc_vector_reox = pc_vector_deox_reox(reox_start_ind:reox_end_ind);

% Simulate cell:AIR requilibration
global dk_in dk_out  AIR

U(1) = 0.012; 
U(2) = U(1);  

dk_in = U(1);
dk_out = U(2);

time_vector_reox_fit = (reox_fit_start_time:2:reox_end_time)' ;
pc_vector_reox_fit = pc_vector_deox_reox(reox_fit_start_ind:reox_end_ind);

AIR = pc_vector_reox_fit(end);

options = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxIter',200)
fun = @(U) simul_0(U,time_vector_reox_fit,pc_vector_reox_fit);
[U2] = fminsearch(fun,U,options);
U = U2;

%% Progress curve

dk_ref = U(1);
dk_out2 = U(2);
O = pc_vector_reox_fit(1);
 
reactant_vector = [O  AIR]';

% Progress curve
options = odeset('AbsTol',1e-8);
[~,product_matrix] = ode15s(@reaction_0,time_vector_reox_fit,reactant_vector,options);

% SSE
trial_vector = product_matrix(:,1);
g = trial_vector - pc_vector_reox_fit;
sos1 = g'*g;
sse = sos1

% Rsquared 
c_pc_vector_reox_fit = pc_vector_reox_fit - mean(pc_vector_reox_fit);
sst = c_pc_vector_reox_fit'*c_pc_vector_reox_fit;
rsquare = 1-sse/sst

% Extension back in time
[~,product_matrix_rev] = ode15s(@reaction_0,-time_vector_reox_fit,reactant_vector,options);
time_vector_reox_fit_rev = 2*time_vector_reox_fit(1)-time_vector_reox_fit;

[time_vector_reox_fit_rev,I] = sort(time_vector_reox_fit_rev)
product_matrix_rev = product_matrix_rev(I,1);

time_vector_fit = [time_vector_reox_fit_rev; time_vector_reox_fit(2:end)];
product_vector_fit = [product_matrix_rev;product_matrix(2:end,1)];

reox_fit_start_ind = find(time_vector_fit == reox_start_time)

time_vector_fit = time_vector_fit(reox_fit_start_ind:end);
product_vector_fit = product_vector_fit(reox_fit_start_ind:end);

%%

% Plot fit with data.
figure(Graphic_method_fig_3)

subplot(1,3,2)
plot(time_vector_reox,pc_vector_reox,'-b','LineWidth',3);hold on
plot(time_vector_fit,product_vector_fit,'--r','LineWidth',1.5);
xlim([(reox_start_time -50) (reox_end_time + 50)]);ylim([-80 200]);
vline(time_vector_fit(1),'g') 
hline(0,':k')
hline(O2_threshold_2,':k')

legend('[O2] vs. time','Kinetic fit', 'Location', 'NorthEast' );
xlabel('Time (s)')
ylabel('[O2] (\muM)')
grid on

% Take the negative of the difference to derive an upside Hb curve
Hb_reox_curve =  (pc_vector_reox - product_vector_fit);
Hb_reox_curve = -(Hb_reox_curve - Hb_reox_curve(1));
[Hb_reox_curve_max,Hb_reox_curve_max_ind] = max(Hb_reox_curve);

% Divide by the total amount of Hb to obtain the saturation. The scale
% factor can be modified to improve the correspondence between the
% deoxygenation and reoxygention ODC's.

Hb_reox_conc_scale = 1.0; % recommended 1.05
Hb_reox_conc = Hb_reox_curve_max*Hb_reox_conc_scale;
Hb_reox_sat = Hb_reox_curve/Hb_reox_conc;

%% ODC

% Deoxygenation ODC
subplot(1,3,3)
h(1) = plot(pc_vector_deox(1:Hb_deox_curve_min_ind)*(cO2_to_pO2),...
    Hb_deox_sat,'-r','LineWidth',2.0);hold on
xlim([0 100]);ylim([0 1.05]); grid on
xlabel('pO_2 mmHg'); ylabel('Hb saturation');

% Reoxygenation ODC
subplot(1,3,3)
h(2) = plot(pc_vector_reox(1:Hb_reox_curve_max_ind)*(cO2_to_pO2),...
    Hb_reox_sat(1:Hb_reox_curve_max_ind),'--b','LineWidth',3)
xlim([0 100]);ylim([0 1.05]); grid on
xlabel('pO_2 mmHg'); ylabel('Hb saturation');

legend(h(1:2),'Deoxygenation ODC','Reoxygenation ODC', 'Location', 'Best' );

%% Equilibrium constants fit

% Reference human Hb
K_ref(1) = 0.0188;
K_ref(2) = 0.0566;
K_ref(3) = 0.407;
K_ref(4) = 4.28;

% Here we combine the points of the deoxygenation and reoxygenation ODC's.
O2_conc = [pc_vector_deox(1:Hb_deox_curve_min_ind) ; pc_vector_reox(1:Hb_reox_curve_max_ind)];
phi_deox_reox = [Hb_deox_sat(1:Hb_deox_curve_min_ind) ; Hb_reox_sat(1:Hb_reox_curve_max_ind)];

% Here we can cut values below a certain O2 concentration
O2_max = 200;
ind = O2_conc <= O2_max;
fun1 = @(K) simul_1(K,phi_deox_reox(ind),O2_conc(ind));
options = optimoptions('lsqnonlin','FiniteDifferenceType','central',...
    'MaxFunctionEvaluations',40000,'MaxIterations',40000);

% fun2 contains the Adair equation
[K,resnorm_K,residual_K,exitflag_K,output_K,lambda,J_K] = ...
    lsqnonlin(fun1,K_ref,[0 0 0 0],[],options);

% Adair equation
B_num = @(x) K(1)*x + 2*K(1)*K(2)*x.^2 + 3*K(1)*K(2)*K(3)*x.^3 ...
    + 4*K(1)*K(2)*K(3)*K(4)*x.^4;
B_den = @(x) 1 + K(1)*x + K(1)*K(2)*x.^2 + K(1)*K(2)*K(3)*x.^3 ...
    + K(1)*K(2)*K(3)*K(4)*x.^4;
B_ratio = @(x) B_num(x)./(4*B_den(x));
B = B_ratio(O2_conc);

Graphic_method_2 = figure;
set(gcf,'Unit','Normalized','Position',[0.1 0.4 0.35 0.5]);
h3(1) = plot(O2_conc,phi_deox_reox,'.r','MarkerSize',3);hold on
h3(2) = plot(O2_conc,B,'--b','LineWidth',1.5);grid on
xlim([0 200]);ylim([0 1.05]); grid on
legend(h3(1:2),'Deox/reox ODC','Adair fit', 'Location', 'Best' );
xlabel('[O2] (\muM)'); ylabel('Hb saturation');
xticks([0 25 50 75 100 125 150 175 200])

% c50, p50, p50 corrected for CO2 or 2,3-dpg
fun2 = @(x) 0.5 - B_ratio(x)
[c50,resnorm_c50,residual_c50,exitflag_c50,output_c50,lambda,J_c50] = ...
    lsqnonlin(fun2,20,[0],[],options);
p50 = cO2_to_pO2*c50
p50_CO2_40 = p50 + 0.1*(40-0.3) +1.083*10^(-4)*(40-0.3)^2
d_23_dpg = 795.63*(3.21*10^(-3)) - 19660*((3.21*10^(-3))^2)
K
total_Hb = (Hb_deox_conc + Hb_reox_conc)/8

%%
close all
save([csv_filename '_GRAPHIC_ODC.mat'])

