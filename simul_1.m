 
 function g = simul_1(K,Phi_vector,x)
     
% x = O2_conc;    
B_num = K(1)*x + 2*K(1)*K(2)*x.^2 + 3*K(1)*K(2)*K(3)*x.^3 + 4*K(1)*K(2)*K(3)*K(4)*x.^4;
B_den = 1 + K(1)*x + K(1)*K(2)*x.^2 + K(1)*K(2)*K(3)*x.^3 + K(1)*K(2)*K(3)*K(4)*x.^4;
Trial_vector = B_num./(4*B_den);    

g = (Trial_vector - Phi_vector);