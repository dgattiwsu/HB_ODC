function du = reaction_RT(t,s) 

global reox_start_time k_1 k_m1 k_2 T_k_off R_k_off k_on ...
       k_f k_r0 k_r1 k_r2 k_r3 k_r4 dk_ref sorc AIR 
    
     if	t < reox_start_time                  	
         k_1 = sorc;                      
         k_in = 0;
         k_out = 0;
     else                   	
         k_1 = 0;                        
         k_in = dk_ref;
         k_out = dk_ref;
     end

% Reactant_vector = [E  O  EO  T_O4  T_O3  T_O2  T_O1  T  R  R_O1  R_O2  R_O3  R_O4  AIR]';


E = s(1);
O = s(2);
EO = s(3);
T_O4 = s(4);
T_O3 = s(5);
T_O2 = s(6);
T_O1 = s(7);
T = s(8);
R = s(9);
R_O1 = s(10);
R_O2 = s(11);
R_O3 = s(12);
R_O4 = s(13);

% AIR = s(14);

%% DIRECT 
% % Fluxes:
% reaction_1 = (k1*E*O)-(km1*EO);
% reaction_2 = (k2*EO);
% reaction_3 = (kt_off*HtO4)-(kt_on*O*HtO3);
% reaction_4 = (kt_off*HtO3)-(kt_on*O*HtO2);
% reaction_5 = (kt_off*HtO2)-(kt_on*O*HtO1);
% reaction_6 = (kt_off*HtO1)-(kt_on*O*Ht);
% reaction_7 = (kr_off*HrO1)-(kr_on*O*Hr);
% reaction_8 = (kr_off*HrO2)-(kr_on*O*HrO1);
% reaction_9 = (kr_off*HrO3)-(kr_on*O*HrO2);
% reaction_10 = (kr_off*HrO4)-(kr_on*O*HrO3);
% reaction_11 = (k_tr_f*Ht)-(k_tr_r*Hr);
% reaction_12 = (k_tr_f*HtO1)-(k_tr_r*HrO1);
% reaction_13 = (k_tr_f*HtO2)-(k_tr_r*HrO2);
% reaction_14 = (k_tr_f*HtO3)-(k_tr_r*HrO3);
% reaction_15 = (k_tr_f*HtO4)-(k_tr_r*HrO4);
% O2_diffusion = (dk_in*AIR)-(dk_out*O);

reaction_1 = (k_1*O*E)-(k_m1*EO);
reaction_2 = (k_2*EO);
reaction_3 = (T_k_off*T_O4)-(k_on*T_O3*O);
reaction_4 = (T_k_off*T_O3)-(k_on*T_O2*O);
reaction_5 = (T_k_off*T_O2)-(k_on*T_O1*O);
reaction_6 = (T_k_off*T_O1)-(k_on*T*O);
reaction_7 = (k_in*AIR)-(k_out*O);
reaction_8 = (R_k_off*R_O1)-(k_on*R*O);
reaction_9 = (R_k_off*R_O2)-(k_on*R_O1*O);
reaction_10 = (R_k_off*R_O3)-(k_on*R_O2*O);
reaction_11 = (R_k_off*R_O4)-(k_on*R_O3*O);
reaction_12 = (k_f*T)-(k_r0*R);
reaction_13 = (k_f*T_O1)-(k_r1*R_O1);
reaction_14 = (k_f*T_O2)-(k_r2*R_O2);
reaction_15 = (k_f*T_O3)-(k_r3*R_O3);
reaction_16 = (k_f*T_O4)-(k_r4*R_O4);
             
%% ODEs:

% d(E)/dt = 1/cell*(reaction_2 - reaction_1)
% d(O)/dt = 1/cell*(reaction_3 + reaction_4 + reaction_5 + reaction_6 + reaction_7 - reaction_1 + reaction_8 + reaction_9 + reaction_10 + reaction_11)
% d(EO)/dt = 1/cell*(-reaction_2 + reaction_1)
% d(W)/dt = 1/cell*(reaction_2)
% d(T_O4)/dt = 1/cell*(-reaction_3 - reaction_16)
% d(T_O3)/dt = 1/cell*(reaction_3 - reaction_4 - reaction_15)
% d(T_O2)/dt = 1/cell*(reaction_4 - reaction_5 - reaction_14)
% d(T_O1)/dt = 1/cell*(reaction_5 - reaction_6 - reaction_13)
% d(T)/dt = 1/cell*(reaction_6 - reaction_12)
% d(R)/dt = 1/cell*(reaction_8 + reaction_12)
% d(R_O1)/dt = 1/cell*(-reaction_8 + reaction_9 + reaction_13)
% d(R_O2)/dt = 1/cell*(-reaction_9 + reaction_10 + reaction_14)
% d(R_O3)/dt = 1/cell*(-reaction_10 + reaction_11 + reaction_15)
% d(R_O4)/dt = 1/cell*(-reaction_11 + reaction_16)
% d(AIR)/dt = 0

% Reactant_vector = [E  O  EO  T_O4  T_O3  T_O2  T_O1  T  R  R_O1  R_O2  R_O3  R_O4  AIR]';


du = [     
        (reaction_2 - reaction_1)
        (reaction_3 + reaction_4 + reaction_5 + reaction_6 + reaction_7 - reaction_1 + reaction_8 + reaction_9 + reaction_10 + reaction_11)
        (-reaction_2 + reaction_1)
        (-reaction_3 - reaction_16)
        (reaction_3 - reaction_4 - reaction_15)
        (reaction_4 - reaction_5 - reaction_14)
        (reaction_5 - reaction_6 - reaction_13)
        (reaction_6 - reaction_12)
        (reaction_8 + reaction_12)
        (-reaction_8 + reaction_9 + reaction_13)
        (-reaction_9 + reaction_10 + reaction_14)
        (-reaction_10 + reaction_11 + reaction_15)
        (-reaction_11 + reaction_16)
        0
            ];     
     
end
  
  

 
      

     
