 
 function loss = simul_0(U,Time_vector,Product_vector)
 
    global  dk_in dk_out AIR 

    dk_in = U(1);
    dk_out = U(2);

    O = Product_vector(1);

    Reactant_vector = [O  AIR]';

    options = odeset('AbsTol',1e-8);
    [~,Product_matrix] = ode15s(@reaction_0,Time_vector,Reactant_vector,options);

    loss = sum( ( (Product_matrix(:,1) - Product_vector ).^2) );
        
 end
 

