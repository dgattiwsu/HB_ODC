
function Product_matrix_two_koff = two_koff_ODC(U,Time_vector,Product_vector, ...
     E_init,reox_start_time,sorc)

k2 = U(1);
kt_on = sorc;
kt_off = U(3); % k_off T state
dk_ref = U(4);
km1 = U(5);
kt_off_2 = U(6); % k_off_2 T state

E = E_init;
O = Product_vector(1);
EO = 0;
HtO4 = U(2);
HtO3 = 0;
HtO2 = 0;
HtO1 = 0;
Ht = 0;
AIR = Product_vector(end);

Reactant_vector = [E  O  EO  HtO4  HtO3  HtO2  HtO1  Ht AIR]';


options = odeset('AbsTol',1e-8);
[~,Product_matrix_two_koff] = ode15s(@(t,s) reaction_local(t,s,reox_start_time, km1, k2, kt_on, kt_off, ...
    kt_off_2, dk_ref, sorc, AIR),Time_vector,Reactant_vector,options);

end

  
function du = reaction_local(t,s,reox_start_time, km1, k2, kt_on, kt_off, ...
       kt_off_2, dk_ref, sorc, AIR) 


     if	t < reox_start_time                  	
         k1 = sorc;                      
         dk_in = 0;
         dk_out = 0;
     else                   	
         k1 = 0;                        
         dk_in = dk_ref;
         dk_out = dk_ref;
     end

% Reactant_vector = [E  O  EO  HtO4  HtO3  HtO2  HtO1  Ht AIR]';

E = s(1);
O = s(2);
EO = s(3);
HtO4 = s(4);
HtO3 = s(5);
HtO2 = s(6);
HtO1 = s(7);
Ht = s(8);
% AIR = s(9);
        
%% DIRECT 
% Fluxes:
reaction_1 = (k1*E*O)-(km1*EO);
reaction_2 = (k2*EO);
reaction_3 = (kt_off_2*HtO4)-(kt_on*O*HtO3);
reaction_4 = (kt_off_2*HtO3)-(kt_on*O*HtO2);
reaction_5 = (kt_off*HtO2)-(kt_on*O*HtO1);
reaction_6 = (kt_off*HtO1)-(kt_on*O*Ht);
O2_diffusion = (dk_in*AIR)-(dk_out*O); 

%% ODEs:
% d(E)/dt = (-reaction_1 + reaction_2)
% d(O)/dt = (-reaction_1 + reaction_3 + reaction_4 + reaction_5 + reaction_6 + O2_diffusion)
% d(EO)/dt = (reaction_1 - reaction_2)
% d(W)/dt = (reaction_2)
% d(HtO4)/dt = (-reaction_3)
% d(HtO3)/dt = (reaction_3 - reaction_4)
% d(HtO2)/dt = (reaction_4 - reaction_5)
% d(HtO1)/dt = (reaction_5 - reaction_6)
% d(Ht)/dt = (reaction_6)
% dAIR = 0;

du = [     
        (-reaction_1 + reaction_2)
        (-reaction_1 + reaction_3 + reaction_4 + reaction_5 + reaction_6 + O2_diffusion)
        (reaction_1 - reaction_2)
        (-reaction_3)
        (reaction_3 - reaction_4)
        (reaction_4 - reaction_5)
        (reaction_5 - reaction_6)
        (reaction_6)
        0
            ];     
     
end
  

 
      

     
