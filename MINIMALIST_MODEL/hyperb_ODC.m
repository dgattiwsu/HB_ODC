

function Product_matrix_ind = hyperb_ODC(U,Time_vector,Product_vector,E_init,...
     reox_start_time,sorc)
          
k_2 = U(1);
k_off_1 = U(3);    
k_off_2 = U(3);
k_off_3 = U(3);
k_off_4 = U(3);
dk_ref = U(4);
k_m1 = U(5);
k_on = sorc;

%

 E = E_init;
 O = Product_vector(1);
 EO = 0;
 Hb4_O = U(2);
 Hb3_O = U(2);
 Hb2_O = U(2);
 Hb1_O = U(2);
 Hb4 = 0;
 Hb3 = 0;
 Hb2 = 0;
 Hb1 = 0;
 AIR = Product_vector(end);
 
Reactant_vector = [E  O  EO  Hb4_O  Hb3_O  Hb2_O  Hb1_O  Hb4  Hb3  Hb2  Hb1  AIR]';


options = odeset('AbsTol',1e-8);
[~,Product_matrix_ind] = ode15s(@(t,s) reaction(t,s,reox_start_time, k_m1, k_2, k_on, k_off_1, k_off_2, ...
       k_off_3, k_off_4, dk_ref, sorc, AIR),Time_vector,Reactant_vector,options);

end

function du = reaction(t,s,reox_start_time, k_m1, k_2, k_on, k_off_1, k_off_2, ...
       k_off_3, k_off_4, dk_ref, sorc, AIR) 
    
     if	t < reox_start_time                  	
         k_1 = sorc;                      
         k_in = 0;
         k_out = 0;
     else                   	
         k_1 = 0;                        
         k_in = dk_ref;
         k_out = dk_ref;
     end
     
%      k_on = sorc;

 E = s(1);
 O = s(2);
 EO = s(3);
 Hb4_O = s(4);
 Hb3_O = s(5);
 Hb2_O = s(6);
 Hb1_O = s(7);
 Hb4 = s(8);
 Hb3 = s(9);
 Hb2 = s(10);
 Hb1 = s(11);
 %AIR = s(12);
 
% Reactant_vector = [E  O  EO  Hb4_O  Hb3_O  Hb2_O  Hb1_O  Hb4  Hb3  Hb2  Hb1  AIR]';


%% DIRECT 
% Fluxes:

reaction_1 = (k_1*O*E)-(k_m1*EO);
reaction_2 = (k_2*EO);
reaction_3 = (k_off_4*Hb4_O)-(k_on*O*Hb4);
reaction_4 = (k_off_3*Hb3_O)-(k_on*O*Hb3);
reaction_5 = (k_off_2*Hb2_O)-(k_on*O*Hb2);
reaction_6 = (k_off_1*Hb1_O)-(k_on*O*Hb1);
reaction_7 = (k_in*AIR)-(k_out*O);

%% 'ODEs:
% d(E)/dt = (reaction_2 - reaction_1)
% d(O)/dt = (reaction_3 + reaction_4 + reaction_5 + reaction_6 + reaction_7 - reaction_1)
% d(EO)/dt = (-reaction_2 + reaction_1)
% d(W)/dt = (reaction_2) % Not necessary
% d(Hb4_O)/dt = (-reaction_3)
% d(Hb3_O)/dt = (-reaction_4)
% d(Hb2_O)/dt = (-reaction_5)
% d(Hb1_O)/dt = (-reaction_6)
% d(Hb4)/dt = (reaction_3)
% d(Hb3)/dt = (reaction_4)
% d(Hb2)/dt = (reaction_5)
% d(Hb1)/dt = (reaction_6)
% d(AIR)/dt = (-reaction_7)
     
% Reactant_vector = [E  O  EO  Hb4_O  Hb3_O  Hb2_O  Hb1_O  Hb4  Hb3  Hb2  Hb1  AIR]';
          
du = [
         (reaction_2 - reaction_1)
         (reaction_3 + reaction_4 + reaction_5 + reaction_6 + reaction_7 - reaction_1)
         (-reaction_2 + reaction_1)
         (-reaction_3)
         (-reaction_4)
         (-reaction_5)
         (-reaction_6)
         (reaction_3)
         (reaction_4)
         (reaction_5)
         (reaction_6)
         0                
            ];

end
  
  

 
      

     
