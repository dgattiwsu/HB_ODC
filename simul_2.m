 
 function loss = simul_2(U,Time_vector,Product_vector,E_init,...
     phi_deox_ind,phi_reox_ind,weights,loss_type_1,loss_type_2,central)
 
global  km1 k2 ko1 ko2 ko3 ko4 dk_ref ks1 ks2 ks3 ks4 sorc AIR
         
k2 = U(1);
ko1 = U(3);    
ko2 = U(4);
ko3 = U(5);
ko4 = U(6);
dk_ref = U(7);
ks1 = U(8);
ks2 = U(9);
ks3 = U(10);
ks4 = U(11);
km1 = U(12);

E = E_init;
O = Product_vector(1);
EO = 0;
HO4 = U(2);
HO3 = 0;
HO2 = 0;
HO = 0;
H = 0;
% AIR = Product_vector(end);
    
Reactant_vector = [E  O  EO  HO4  HO3  HO2  HO  H AIR]';

options = odeset('AbsTol',1e-8);
[~,Product_matrix] = ode15s(@reaction,Time_vector,Reactant_vector,options);

%------------ORIGINAL CODE-------------------------------------------------
% Use this section below only if there is significant difference between
% the deoxygenation and reoxygenation ODC.
% % phi_deox = (Product_matrix(phi_deox_ind,4)*4+Product_matrix(phi_deox_ind,5)*3+...
% %     Product_matrix(phi_deox_ind,6)*2+Product_matrix(phi_deox_ind,7))/(Product_matrix(1,4)*4);
% % 
% % phi_reox = (Product_matrix(phi_reox_ind,4)*4+Product_matrix(phi_reox_ind,5)*3+...
% %     Product_matrix(phi_reox_ind,6)*2+Product_matrix(phi_reox_ind,7))/(Product_matrix(1,4)*4);
% % 
% % % find K
% % K = sorc./U(6:-1:3);
% % 
% % % solver options
% % options1 = optimoptions('lsqnonlin','FiniteDifferenceType','forward',...
% %     'MaxFunctionEvaluations',40000,'MaxIterations',40000,'Display','off');
% % 
% % % find K_deox
% % O2_conc_deox = Product_matrix(phi_deox_ind,2);
% % fun_d = @(K) simul_1(K,phi_deox,O2_conc_deox);
% % [K_deox,~,~,~,~,~,~] = lsqnonlin(fun_d,K,[0 0 0 0],[],options1);
% % 
% % % find K_reox
% % O2_conc_reox = Product_matrix(phi_reox_ind,2);
% % fun_r = @(K) simul_1(K,phi_reox,O2_conc_reox);
% % [K_reox,~,~,~,~,~,~] = lsqnonlin(fun_r,K,[0 0 0 0],[],options1);
% % 
% % Here we remove from the fit the part of the progress curve between deox
% % and reox
% phi_deox_reox_ind = [phi_deox_ind phi_reox_ind];
% weights_1 = [ones(size(Product_vector(phi_deox_reox_ind),1),1) ; weights];
% 
% switch loss_type
%     case 'sse'
%         loss = sum((([Product_matrix(phi_deox_reox_ind,2);K_deox'] - ...
%             [Product_vector(phi_deox_reox_ind);K_reox']).*weights_1).^2);
%     case 'res'
%         loss = ([Product_matrix(phi_deox_reox_ind,2);K_deox'] - ...
%             [Product_vector(phi_deox_reox_ind);K_reox']).*weights_1;
%         
% end
% 
% % Here we keep in the fit the central part of the progress curve between
% % deox and reox
% % weights_1 = [ones(size(Product_vector,1),1) ; weights];
% % 
% % switch loss_type
% %     case 'sse'
% %         loss = sum((([Product_matrix(:,2);K_deox'] - ...
% %             [Product_vector;K_reox']).*weights_1).^2);
% %     case 'res'
% %         loss = ([Product_matrix(:,2);K_deox'] - ...
% %             [Product_vector;K_reox']).*weights_1;
% %         
% % end
% 
% % Here we do not consider the contribution to the sum from the equilibrium
% % constants.
% %
% % switch loss_type
% %     case 'sse'
% %         loss = sum(((Product_matrix(:,2) - Product_vector)).^2);
% %     case 'res'
% %         loss = (Product_matrix(:,2) - Product_vector);
% % end
%--------------------------------------------------------------------------
switch loss_type_2
    
    case 'use_K'
    % Use this section if there is significant difference between
    % the deoxygenation and reoxygenation ODC.

    phi_deox = (Product_matrix(phi_deox_ind,4)*4+Product_matrix(phi_deox_ind,5)*3+...
        Product_matrix(phi_deox_ind,6)*2+Product_matrix(phi_deox_ind,7))/(Product_matrix(1,4)*4);

    phi_reox = (Product_matrix(phi_reox_ind,4)*4+Product_matrix(phi_reox_ind,5)*3+...
        Product_matrix(phi_reox_ind,6)*2+Product_matrix(phi_reox_ind,7))/(Product_matrix(1,4)*4);

    % find K
    K = sorc./U(6:-1:3);

    % solver options
    options1 = optimoptions('lsqnonlin','FiniteDifferenceType','forward',...
        'MaxFunctionEvaluations',40000,'MaxIterations',40000,'Display','off');

    % find K_deox
    O2_conc_deox = Product_matrix(phi_deox_ind,2);
    fun_d = @(K) simul_1(K,phi_deox,O2_conc_deox);
    [K_deox,~,~,~,~,~,~] = lsqnonlin(fun_d,K,[0 0 0 0],[],options1);

    % find K_reox
    O2_conc_reox = Product_matrix(phi_reox_ind,2);
    fun_r = @(K) simul_1(K,phi_reox,O2_conc_reox);
    [K_reox,~,~,~,~,~,~] = lsqnonlin(fun_r,K,[0 0 0 0],[],options1);


        switch central

        case 'discard_central'
        % Here we remove from the fit the part of the progress curve between deox
        % and reox
        phi_deox_reox_ind = [phi_deox_ind phi_reox_ind];
        weights_1 = [ones(size(Product_vector(phi_deox_reox_ind),1),1) ; weights];

            switch loss_type_1
                case 'sse'
                    loss = sum((([Product_matrix(phi_deox_reox_ind,2);K_deox'] - ...
                        [Product_vector(phi_deox_reox_ind);K_reox']).*weights_1).^2);
                case 'res'
                    loss = ([Product_matrix(phi_deox_reox_ind,2);K_deox'] - ...
                        [Product_vector(phi_deox_reox_ind);K_reox']).*weights_1;

            end

        case 'keep_central'
        % Here we keep in the fit the central part of the progress curve between
        % deox and reox
        weights_1 = [ones(size(Product_vector,1),1) ; weights];

            switch loss_type_1
                case 'sse'
                    loss = sum((([Product_matrix(:,2);K_deox'] - ...
                        [Product_vector;K_reox']).*weights_1).^2);
                case 'res'
                    loss = ([Product_matrix(:,2);K_deox'] - ...
                        [Product_vector;K_reox']).*weights_1;

            end
            
        end
        
        
    case 'do_not_use_K'
    % Use this section if there is not significant difference between
    % the deoxygenation and reoxygenation ODC.

        switch central

        case 'discard_central'
        % Here we remove from the fit the part of the progress curve between deox
        % and reox
        phi_deox_reox_ind = [phi_deox_ind phi_reox_ind];

            switch loss_type_1
                case 'sse'
                    loss = sum((Product_matrix(phi_deox_reox_ind,2) - ...
                        Product_vector(phi_deox_reox_ind)).^2);
                case 'res'
                    loss = Product_matrix(phi_deox_reox_ind,2) - ...
                        Product_vector(phi_deox_reox_ind);
            end

        case 'keep_central'
        % Here we keep in the fit the central part of the progress curve between
        % deox and reox
            switch loss_type_1
                case 'sse'
                    loss = sum((Product_matrix(:,2) - ...
                        Product_vector).^2);
                case 'res'
                    loss = Product_matrix(:,2) - ...
                        Product_vector;
            end            
        end

end

end

 
 function g = simul_1(K,Phi_vector,x)
 
 % x = O2_conc;    
 B_num = K(1)*x + 2*K(1)*K(2)*x.^2 + 3*K(1)*K(2)*K(3)*x.^3 + 4*K(1)*K(2)*K(3)*K(4)*x.^4;
 B_den = 1 + K(1)*x + K(1)*K(2)*x.^2 + K(1)*K(2)*K(3)*x.^3 + K(1)*K(2)*K(3)*K(4)*x.^4;
 Trial_vector = B_num./(4*B_den);    

 g = (Trial_vector - Phi_vector);

 end

