 
 function loss = simul_3(EU,U,Time_vector,Product_vector)
 
global km1 k2 kt_off kt_off_2 ...
       dk_ref AIR
         
k2 = U(1);
kt_off = U(3); % k_off T state
dk_ref = U(4);
km1 = U(5);
kt_off_2 = U(6); % k_off_2 T state

E = EU;
O = Product_vector(1);
EO = 0;
HtO4 = U(2);
HtO3 = 0;
HtO2 = 0;
HtO1 = 0;
Ht = 0;
% AIR = Product_vector(end);

Reactant_vector = [E  O  EO  HtO4  HtO3  HtO2  HtO1  Ht AIR]';

options = odeset('AbsTol',1e-8);
[~,Product_matrix] = ode15s(@reaction,Time_vector,Reactant_vector,options);


loss = sum( ( (Product_matrix(:,2) - ...
    Product_vector ).^2) );
        
 end
 

