function du = reaction(t,s) 

global reox_start_time k1 km1 k2 ko1 kom1 ko2 kom2 ...
       ko3 kom3 ko4 kom4 dk_in dk_ref ks1 ks2 ks3 ks4 sorc AIR
    
     if	t < reox_start_time                  	
         k1 = sorc;                      
         dk_in = 0;
         dk_out = 0;
         s1 = 0;
         s2 = 0;
         s3 = 0;
         s4 = 0;
     else                   	
         k1 = 0;                        
         dk_in = dk_ref;
         dk_out = dk_ref;
         s1 = ks1;
         s2 = ks2;
         s3 = ks3;
         s4 = ks4;         
     end

E = s(1);
O = s(2);
EO = s(3);
HO4 = s(4);
HO3 = s(5);
HO2 = s(6);
HO = s(7);
H = s(8);
% AIR = s(9);

%% Rate matrix representation.
% s = [E       O                 EO           HO4     HO3               HO2               HO                H          AIR]';

% K = [-k1*O   0                 (km1 + k2)   0       0                 0                 0                 0          0
%       0       (-k1*E -dk_out)   km1          ko1     (-kom1*O + ko2)   (-kom2*O + ko3)   (-kom3*O + ko4)   -kom4*O    dk_in 
%       0       k1*E              (-km1 -k2)   0       0                 0                 0                 0          0
%       0       0                 k2           0       0                 0                 0                 0          0
%       0       0                 0            -ko1    kom1*O            0                 0                 0          0
%       0       0                 0            ko1     (-kom1*O -ko2)    kom2*O            0                 0          0
%       0       0                 0            0       ko2               (-kom2*O -ko3)    kom3*O            0          0
%       0       0                 0            0       0                 ko3               (-kom3*O -ko4)    kom4*O     0
%       0       0                 0            0       0                 0                 ko4               -kom4*O    0 
%       0       0                 0            0       0                 0                 0                 0          0];

% du = K*s;

%% DIRECT 
% Fluxes:

reaction_1 = k1*E*O-km1*EO;
reaction_2 = k2*EO;
reaction_3 = (ko1 + s1)*HO4-kom1*O*HO3;
reaction_4 = (ko2 + s2)*HO3-kom2*O*HO2;
reaction_5 = (ko3 + s3)*HO2-kom3*O*HO;
reaction_6 = (ko4 + s4)*HO-kom4*O*H;
reaction_7 = dk_in*AIR-dk_out*O;

%% 'ODEs:
%      dEdt = (-reaction_1 + reaction_2)
%      dOdt = (-reaction_1 + reaction_3 + reaction_4 + reaction_5 + reaction_6 + reaction_7)
%      dEOdt = (reaction_1 - reaction_2)
%      dWdt = (reaction_2) % W is constant and not use anywhere else; this ODE can be removed
%      dHO4dt = (-reaction_3)
%      dHO3dt = (reaction_3 - reaction_4)
%      dHO2dt = (reaction_4 - reaction_5)
%      dHOdt = (reaction_5 - reaction_6)
%      dHdt = (reaction_6)
%      dAIR = 0;

du = [
        (-reaction_1 + reaction_2)
        (-reaction_1 + reaction_3 + reaction_4 + reaction_5 + reaction_6 + reaction_7)
        (reaction_1 - reaction_2)
        (-reaction_3)
        (reaction_3 - reaction_4)
        (reaction_4 - reaction_5)
        (reaction_5 - reaction_6)
        (reaction_6)
        0                
            ];


end
  
  

 
      

     
