function [Rmtle, Rmlts, Rmle, Rols, Amtle, Amlts, Amle, Aols, Smtle, Smlts, Smle, Sols, R_sd, A_sd, S_sd, CI_R_Star_hat_mse, CI_A_hat_mse, CI_Sigma_hat_mse] = simulation(R_ast, A, Sigma, R0, T_array, nrep, out_pcent, ncp, bdp)

    p = size(A, 1);

    %iA = inv(A);
    iSigma = inv(Sigma);


    % Initialize error matrices for MLE, OLS, MLTS, and MTLE
    R_ast_errors_0 = zeros(length(T_array), nrep); % OLS
    A_errors_0     = zeros(length(T_array), nrep);
    Sigma_errors_0 = zeros(length(T_array), nrep);

    R_ast_errors_1 = zeros(length(T_array), nrep); % MLE
    A_errors_1     = zeros(length(T_array), nrep);
    Sigma_errors_1 = zeros(length(T_array), nrep);
    
    R_ast_errors_2 = zeros(length(T_array), nrep); % MLTS
    A_errors_2     = zeros(length(T_array), nrep);
    Sigma_errors_2 = zeros(length(T_array), nrep);
    
    R_ast_errors_3 = zeros(length(T_array), nrep); % MTLE
    A_errors_3     = zeros(length(T_array), nrep);
    Sigma_errors_3 = zeros(length(T_array), nrep);

    T_range = 1:length(T_array);

    parfor rep = 1:nrep    
        for i_T = T_range
            n  = T_array(i_T);
            dt = 1.0;

            rng(2*rep); 
            R = squeeze(forward_map_out(R_ast, A, Sigma, R0, dt, n, 1, out_pcent, ncp))';


            % MLE estimation
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map(R, dt, "ml", 1E-6, 5000);
            iSigma_hat = pinv(Sigma_hat);
            R_ast_errors_1(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_1(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_1(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');

            % OLS estimation
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map(R, dt, "ols", 1E-6, 5000);
            iSigma_hat = pinv(Sigma_hat);        
            R_ast_errors_0(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_0(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_0(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');

            % MLTS estimation
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map_rob(R, dt, bdp, "mlts", 1E-6, 5000);
            iSigma_hat = pinv(Sigma_hat);        
            R_ast_errors_2(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_2(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_2(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');

            % MTLE estimation
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map_rob(R, dt, bdp, "mtle", 1E-6, 5000); 
            iSigma_hat = pinv(Sigma_hat);        
            R_ast_errors_3(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_3(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_3(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');
        end
    end
    Rmtle = sqrt(mean(R_ast_errors_3, 2));
    Rmlts = sqrt(mean(R_ast_errors_2, 2));
    Rmle  = sqrt(mean(R_ast_errors_1, 2));
    Rols  = sqrt(mean(R_ast_errors_0, 2));

    Amtle = sqrt(mean(A_errors_3, 2));
    Amlts = sqrt(mean(A_errors_2, 2));
    Amle  = sqrt(mean(A_errors_1, 2));
    Aols  = sqrt(mean(A_errors_0, 2));

    Smtle = sqrt(mean(Sigma_errors_3, 2));
    Smlts = sqrt(mean(Sigma_errors_2, 2));
    Smle  = sqrt(mean(Sigma_errors_1, 2));
    Sols  = sqrt(mean(Sigma_errors_0, 2));

    Rmtle_sd = std(R_ast_errors_3, 0, 2);
    Rmlts_sd = std(R_ast_errors_2, 0, 2);
    Rmle_sd  = std(R_ast_errors_1, 0, 2);
    Rols_sd  = std(R_ast_errors_0, 0, 2);

    Amlte_sd = std(A_errors_3, 0, 2);
    Amlts_sd = std(A_errors_2, 0, 2);
    Amle_sd  = std(A_errors_1, 0, 2);
    Aols_sd  = std(A_errors_0, 0, 2);

    Smlte_sd = std(Sigma_errors_3, 0, 2);
    Smlts_sd = std(Sigma_errors_2, 0, 2);
    Smle_sd  = std(Sigma_errors_1, 0, 2);
    Sols_sd  = std(Sigma_errors_0, 0, 2);

    R_sd = [Rmtle_sd Rmlts_sd Rmle_sd Rols_sd];
    A_sd = [Amlte_sd Amlts_sd Amle_sd Aols_sd];
    S_sd = [Smlte_sd Smlts_sd Smle_sd Sols_sd];

    R_mean = [Rmtle Rmlts Rmle Rols];
    A_mean = [Amtle Amlts Amle Aols];
    S_mean = [Smtle Smlts Smle Sols];

    z0975 = norminv(0.975);

    llim_Rmtle = (Rmtle.^2 - z0975*R_sd(:,1)/sqrt(nrep));
    llim_Rmtle(llim_Rmtle < 0 ) = 0;
    ulim_Rmlte = (Rmtle.^2 + z0975*R_sd(:,1)/sqrt(nrep));

    llim_Rmlts = (Rmlts.^2 - z0975*R_sd(:,2)/sqrt(nrep));
    llim_Rmlts(llim_Rmlts < 0 ) = 0;
    ulim_Rmlts = (Rmlts.^2 + z0975*R_sd(:,2)/sqrt(nrep));

    llim_Rmle = (Rmle.^2  - z0975*R_sd(:,3)/sqrt(nrep));
    llim_Rmle(llim_Rmle < 0 ) = 0;
    ulim_Rmle = (Rmle.^2  + z0975*R_sd(:,3)/sqrt(nrep));

    llim_Rols = (Rols.^2  - z0975*R_sd(:,4)/sqrt(nrep));
    llim_Rols(llim_Rols < 0 ) = 0;
    ulim_Rols = (Rols.^2  + z0975*R_sd(:,4)/sqrt(nrep));

    CI_Rmtle = [sqrt(llim_Rmtle) sqrt(ulim_Rmlte)];
    CI_Rmlts = [sqrt(llim_Rmlts) sqrt(ulim_Rmlts)];
    CI_Rmle  = [sqrt(llim_Rmle)  sqrt(ulim_Rmle)];
    CI_Rols  = [sqrt(llim_Rols)  sqrt(ulim_Rols)];

    llim_Amtle = (Amtle.^2 - z0975*A_sd(:,1)/sqrt(nrep));
    llim_Amtle(llim_Amtle < 0 ) = 0;
    ulim_Amtle = (Amtle.^2 + z0975*A_sd(:,1)/sqrt(nrep));

    llim_Amlts = (Amlts.^2 - z0975*A_sd(:,2)/sqrt(nrep));
    llim_Amlts(llim_Amlts < 0 ) = 0;
    ulim_Amlts = (Amlts.^2 + z0975*A_sd(:,2)/sqrt(nrep));

    llim_Amle = (Amle.^2  - z0975*A_sd(:,3)/sqrt(nrep));
    llim_Amle(llim_Amle < 0) = 0;
    ulim_Amle = (Amle.^2  + z0975*A_sd(:,3)/sqrt(nrep));

    llim_Aols = (Aols.^2  - z0975*A_sd(:,4)/sqrt(nrep));
    llim_Aols(llim_Aols < 0) = 0;
    ulim_Aols = (Aols.^2  + z0975*A_sd(:,4)/sqrt(nrep));

    CI_Amtle = [sqrt(llim_Amtle) sqrt(ulim_Amtle)];
    CI_Amlts = [sqrt(llim_Amlts) sqrt(ulim_Amlts)];
    CI_Amle  = [sqrt(llim_Amle)  sqrt(ulim_Amle)];
    CI_Aols  = [sqrt(llim_Aols)  sqrt(ulim_Aols)];

    llim_Smtle = (Smtle.^2 - z0975*S_sd(:,1)/sqrt(nrep));
    llim_Smtle(llim_Smtle < 0) = 0;
    ulim_Smtle = (Smtle.^2 + z0975*S_sd(:,1)/sqrt(nrep));

    llim_Smlts = (Smlts.^2 - z0975*S_sd(:,2)/sqrt(nrep));
    llim_Smlts(llim_Smlts < 0) = 0;
    ulim_Smlts = (Smlts.^2 + z0975*S_sd(:,2)/sqrt(nrep));

    llim_Smle = (Smle.^2  - z0975*S_sd(:,3)/sqrt(nrep));
    llim_Smle(llim_Smle < 0) = 0;
    ulim_Smle = (Smle.^2  + z0975*S_sd(:,3)/sqrt(nrep));

    llim_Sols = (Sols.^2  - z0975*S_sd(:,4)/sqrt(nrep));
    llim_Sols(llim_Sols < 0) = 0;
    ulim_Sols = (Sols.^2  + z0975*S_sd(:,4)/sqrt(nrep));

    CI_Smtle = [sqrt(llim_Smtle) sqrt(ulim_Smtle)];
    CI_Smlts = [sqrt(llim_Smlts) sqrt(ulim_Smlts)];
    CI_Smle  = [sqrt(llim_Smle)  sqrt(ulim_Smle)];
    CI_Sols  = [sqrt(llim_Sols)  sqrt(ulim_Sols)];

    CI_R_Star_hat_mse = [CI_Rmtle CI_Rmlts CI_Rmle CI_Rols];
    CI_A_hat_mse      = [CI_Amtle CI_Amlts CI_Amle CI_Aols];
    CI_Sigma_hat_mse  = [CI_Smtle CI_Smlts CI_Smle CI_Sols]; 


    if not(isfolder("SimResults"))
        mkdir("SimResults")
    end

    cd ./SimResults/

    filename_R = ['R','_',num2str(bdp),'_',num2str(out_pcent),'_',num2str(ncp),'.mat'];
    filename_A = ['A','_',num2str(bdp),'_',num2str(out_pcent),'_',num2str(ncp),'.mat'];
    filename_S = ['S','_',num2str(bdp),'_',num2str(out_pcent),'_',num2str(ncp),'.mat'];
    save(filename_R,'R_mean','R_sd', 'CI_R_Star_hat_mse')
    save(filename_A,'A_mean','A_sd', 'CI_A_hat_mse')
    save(filename_S,'S_mean','S_sd', 'CI_Sigma_hat_mse')
    cd ../

end