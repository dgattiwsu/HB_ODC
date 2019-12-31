 
 function loss = simul_4(EU,U,Time_vector,Product_vector)
 
global  k2 ko1 ko2 ko3 ko4 dk_ref AIR
         
k2 = U(1);
ko1 = U(3);    
ko2 = U(4);
ko3 = U(5);
ko4 = U(6);
dk_ref = U(7);

E = EU(1);
O = Product_vector(1);
EO = 0;
HO4 = EU(2);
HO3 = 0;
HO2 = 0;
HO = 0;
H = 0;
% AIR = Product_vector(end);
    
Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H AIR]';

options = odeset('AbsTol',1e-8);
[~,Product_matrix] = ode15s(@reaction,Time_vector,Reactant_vector,options);


loss = sum( ( (Product_matrix(:,2) - ...
    Product_vector ).^2) );
        
 end
 

