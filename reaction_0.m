function du = reaction_0(t,s) 

global dk_in dk_out AIR    
   
O = s(1);

%% DIRECT 
% Fluxes:

reaction_1 = dk_in*AIR - dk_out*O;

%% 'ODEs:
%      dOdt = reaction_1
%      dAIR = 0;

du = [
        reaction_1 
        0                
            ];

end
  
  

 
      

     
