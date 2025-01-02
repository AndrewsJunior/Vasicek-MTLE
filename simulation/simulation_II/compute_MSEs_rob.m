function compute_MSEs_rob(R_ast, A, Sigma, R0, T_array, nrep, Outlierfrac, out_pcent, ncp)
    
    p = size(A, 1);

    iSigma = inv(Sigma);

    
    R_ast_errors_0 = zeros(length(T_array), nrep);%OLS
    A_errors_0     = zeros(length(T_array), nrep);
    Sigma_errors_0 = zeros(length(T_array), nrep);

    R_ast_errors_1 = zeros(length(T_array), nrep);%MLE
    A_errors_1     = zeros(length(T_array), nrep);
    Sigma_errors_1 = zeros(length(T_array), nrep);

    R_ast_errors_2 = zeros(length(T_array), nrep);%MLTS
    A_errors_2     = zeros(length(T_array), nrep);
    Sigma_errors_2 = zeros(length(T_array), nrep);

    R_ast_errors_3 = zeros(length(T_array), nrep);%MTLE
    A_errors_3     = zeros(length(T_array), nrep);
    Sigma_errors_3 = zeros(length(T_array), nrep);

    T_range = 1:length(T_array);
    
    %start timing
    %tic;
    parfor rep = 1:nrep
        for i_T = T_range
            n  = T_array(i_T);
            dt = 1.0;

            rng(2*rep); 
            
            R = squeeze(forward_map_out(R_ast, A, Sigma, R0, dt, n, 1, out_pcent, ncp))';
    
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map(R, dt, "ml", 1E-6, 5000);
            iSigma_hat = pinv(Sigma_hat);
            R_ast_errors_1(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_1(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_1(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');%MLE
    
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map(R, dt, "ols", 1E-6, 5000);
            iSigma_hat = pinv(Sigma_hat);   
    
            R_ast_errors_0(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_0(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_0(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');%OLS
    
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map_rob(R, dt, Outlierfrac, "mlts", 1E-6, 5000);
            iSigma_hat = pinv(Sigma_hat);
    
            R_ast_errors_2(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_2(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_2(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');%MLTS
            
            [R_ast_hat, A_hat, Sigma_hat] = inverse_map_rob(R, dt, Outlierfrac, "mtle", 1E-6, 5000);
            iSigma_hat = pinv(Sigma_hat);  
    
            R_ast_errors_3(i_T, rep) = mean((R_ast_hat - R_ast).^2, 'all');
            A_errors_3(i_T, rep)     = mean((A - A_hat).^2, 'all');
            Sigma_errors_3(i_T, rep) = mean((eye(p) - iSigma*Sigma_hat).^2, 'all') + ...
                                       mean((eye(p) - Sigma*iSigma_hat).^2, 'all');%MTLE
    
        end
    end

    if not(isfolder("SimResultsMSEs"))
        mkdir("SimResultsMSEs")
    end

    cd ./SimResultsMSEs/

    filename_R = ['R_err','_',num2str(Outlierfrac),'_',num2str(out_pcent),'_',num2str(ncp),'.mat'];
    filename_A = ['A_err','_',num2str(Outlierfrac),'_',num2str(out_pcent),'_',num2str(ncp),'.mat'];
    filename_S = ['S_err','_',num2str(Outlierfrac),'_',num2str(out_pcent),'_',num2str(ncp),'.mat'];
    save(filename_R,'R_ast_errors_0','R_ast_errors_1','R_ast_errors_2','R_ast_errors_3')
    save(filename_A,'A_errors_0','A_errors_1', 'A_errors_2', 'A_errors_3')
    save(filename_S,'Sigma_errors_0','Sigma_errors_1', 'Sigma_errors_2', 'Sigma_errors_3')
    cd ..\

    
end