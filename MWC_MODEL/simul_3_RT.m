 
 function loss = simul_3_RT(EU,U,Time_vector,Product_vector)
 
global k_m1 k_2 T_k_off R_k_off ...
       k_r0 k_r1 k_r2 k_r3 k_r4 dk_ref AIR 
         
k_2 = U(1);
T_k_off = U(3); % k_off T state
R_k_off = U(4); % k_off R state

% K_T = sorc/U(3)
% K_R = sorc/U(4)
% c = K_T/K_R
c = U(4)/U(3);

k_r0 = U(5); % k for R_0 --> T_0 transition

k_r1 = k_r0*c; % k for R_1 --> T_1 transition
k_r2 = k_r0*c^2; % k for R_2 --> T_2 transition
k_r3 = k_r0*c^3; % k for R_3 --> T_3 transition
k_r4 = k_r0*c^4; % k for R_4 --> T_4 transition

dk_ref = U(6);
k_m1 = U(7);

E = EU;
O = Product_vector(1);
EO = 0;
T_O4 = 0;
T_O3 = 0;
T_O2 = 0;
T_O1 = 0;
T = 0;
R_O4 = U(2);
R_O3 = 0;
R_O2 = 0;
R_O1 = 0;
R = 0;
% AIR = Product_vector(end);

Reactant_vector = [E  O  EO  T_O4  T_O3  T_O2  T_O1  T  R  R_O1  R_O2  R_O3  R_O4  AIR]';


options = odeset('AbsTol',1e-8);
[~,Product_matrix] = ode15s(@reaction_RT,Time_vector,Reactant_vector,options);


loss = sum( ( (Product_matrix(:,2) - ...
    Product_vector ).^2) );
        
 end
 

