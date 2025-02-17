auxiliary = "../../auxiliary";

addpath(auxiliary);

R0 = [4; 7];

R_ast_true = [5; 6];
A_true     = (1/365)*[9 4; 2 3];
Sigma_true = (1/365^2)*[7 3; 3 7];

T_array = 60: 30: 360;
T_range = 1:length(T_array);

ncp = [25, 50]; 
bdp = [0.25, 0.35];

nrep = 5000;

if not(isfolder("SimResults"))
    mkdir("SimResults")
end

cd ./SimResults/

file_R_star_hat = 'R_Star_hat.txt';
file_A_hat = 'A_hat.txt';
file_Sigma_hat = 'Sigma_hat.txt';

if exist(file_R_star_hat, 'file') == 2
    delete(file_R_star_hat);
end

if exist(file_A_hat, 'file') == 2
    delete(file_A_hat);
end

if exist(file_Sigma_hat, 'file') == 2
    delete(file_Sigma_hat);
end

fileID = fopen(file_R_star_hat,'a');
fileIDA = fopen(file_A_hat,'a');
fileIDS = fopen(file_Sigma_hat,'a');

cd ../
%start timing
tic;
for i_bdp = 1:length(bdp)
    if i_bdp == 1
        alphavals = [0.00, 0.05, 0.10, 0.20];
    else
        alphavals = [0.0, 0.1, 0.2, 0.3, 0.4];
    end
    bdpval = bdp(i_bdp);
for i_outp = 1:length(alphavals)
    out_p = alphavals(i_outp);

    for i_ncp = 1:length(ncp)
        ncpval = ncp(i_ncp);
        [Rmtle, Rmlts, Rmle, Rols, Amtle, Amlts, Amle, Aols, Smtle, Smlts, Smle, Sols, R_sd, A_sd, S_sd, CI_R_Star_hat_mse, CI_A_hat_mse, CI_Sigma_hat_mse] = simulation(R_ast_true, A_true, Sigma_true, R0, T_array, nrep, out_p, ncpval, bdpval);

        %for R_star
 
        fprintf(fileID, '\n For bdp = %2.2f,  alpha = %2.2f and ncp =  %2d \n\n', bdpval, out_p, ncpval);
        fprintf(fileID, 'T \t sqrt hat err (MTLE) \t hat sd hat err (MTLE) \t 2.5 percent (MLTE) \t 97.5 percent (MTLE)  \t sqrt hat err (MLTS) \t hat sd hat err (MLTS) \t 2.5 percent (MLTS) \t 97.5 percent (MLTS)');
        fprintf(fileID, '\t sqrt hat err (MLE) \t hat sd hat err (MLE) \t 2.5 percent (MLE) \t 97.5 percent (MLE) \t sqrt hat err (OLS) \t hat sd hat err (OLS) \t 2.5 percent (OLS) \t 97.5 percent (OLS) \n');

        for j_T = 1:length(T_array)
            fprintf(fileID,'%d \t %5.7f \t\t %14.7f \t %5.5f \t\t %5.5f \t\t %5.7f \t\t %14.7f \t %5.5f \t\t %5.5f ',T_array(j_T), Rmtle(j_T), R_sd(j_T,1), CI_R_Star_hat_mse(j_T,1), CI_R_Star_hat_mse(j_T,2), Rmlts(j_T), R_sd(j_T,2), CI_R_Star_hat_mse(j_T,3), CI_R_Star_hat_mse(j_T,4));
            fprintf(fileID,   '\t\t %5.7f \t\t %14.7f \t %5.5f \t\t %5.5f \t\t %5.7f \t\t %14.7f \t %5.5f \t\t %5.5f  \n', Rmle(j_T), R_sd(j_T,3), CI_R_Star_hat_mse(j_T,5), CI_R_Star_hat_mse(j_T,6), Rols(j_T), R_sd(j_T,4), CI_R_Star_hat_mse(j_T,7), CI_R_Star_hat_mse(j_T,8));
        end

        %for A_hat
        fprintf(fileIDA, '\n For bdp = %2.2f, alpha = %2.2f and ncp =  %2d \n\n', bdpval, out_p, ncpval);
        fprintf(fileIDA, 'T \t sqrt hat err (MTLE) \t hat sd hat err (MTLE) \t 2.5 percent (MTLE) \t 97.5 percent (MTLE)  \t sqrt hat err (MLTS) \t hat sd hat err (MLTS) \t 2.5 percent (MLTS) \t 97.5 percent (MLTS)');
        fprintf(fileIDA, '\t sqrt hat err (MLE) \t hat sd hat err (MLE) \t 2.5 percent (MLE) \t 97.5 percent (MLE) \t sqrt hat err (OLS) \t hat sd hat err (OLS) \t 2.5 percent (OLS) \t 97.5 percent (OLS) \n');
  
        for j_T = 1:length(T_array)
            fprintf(fileIDA,'%d \t %5.7f \t\t %14.7f \t\t %5.7f \t  %14.7f \t\t %5.7f \t %14.7f \t\t %5.7f \t  %14.7f', T_array(j_T), Amtle(j_T), A_sd(j_T,1), CI_A_hat_mse(j_T,1), CI_A_hat_mse(j_T,2), Amlts(j_T),  A_sd(j_T,2), CI_A_hat_mse(j_T,3), CI_A_hat_mse(j_T,4));
            fprintf(fileIDA,  ' \t %5.7f \t %14.7f \t\t %5.7f \t  %14.7f \t\t %5.7f \t %14.7f \t\t %5.7f \t  %14.7f  \n', Amle(j_T), A_sd(j_T,3), CI_A_hat_mse(j_T,5), CI_A_hat_mse(j_T,6), Aols(j_T), A_sd(j_T,4), CI_A_hat_mse(j_T,7), CI_A_hat_mse(j_T,8));
        end

        %for S_hat
        
        fprintf(fileIDS, '\n For bdp = %2.2f, alpha = %2.2f and ncp =  %2d \n\n', bdpval, out_p, ncpval);
        fprintf(fileIDS, 'T \t sqrt hat err (MTLE) \t hat sd hat err (MTLE) \t 2.5 percent (MTLE) \t 97.5 percent (MTLE)  \t sqrt hat err (MLTS) \t hat sd hat err (MLTS) \t 2.5 percent (MLTS) \t 97.5 percent (MLTS)');
        fprintf(fileIDS, '\t sqrt hat err (MLE) \t hat sd hat err (MLE) \t 2.5 percent (MLE) \t 97.5 percent (MLE) \t sqrt hat err (OLS) \t hat sd hat err (OLS) \t 2.5 percent (OLS) \t 97.5 percent (OLS) \n');

        for j_T = 1:length(T_array)
            fprintf(fileIDS,'%d \t %5.7f \t\t %14.7f \t\t %5.7f \t  %14.7f \t\t %5.7f \t %14.7f \t\t %5.7f \t  %14.7f', T_array(j_T), Smtle(j_T), S_sd(j_T,1), CI_Sigma_hat_mse(j_T,1), CI_Sigma_hat_mse(j_T,2), Smlts(j_T), S_sd(j_T,2), CI_Sigma_hat_mse(j_T,3), CI_Sigma_hat_mse(j_T,4));
            fprintf(fileIDS,  ' \t %5.7f \t %14.7f \t\t %5.7f \t  %14.7f \t\t %5.7f \t %14.7f \t\t %5.7f \t  %14.7f  \n', Smle(j_T), S_sd(j_T,3), CI_Sigma_hat_mse(j_T,5), CI_Sigma_hat_mse(j_T,6), Sols(j_T), S_sd(j_T,4), CI_Sigma_hat_mse(j_T,7), CI_Sigma_hat_mse(j_T,8));
        end
        
    end
    
end
end

elapsedtime = toc;
% Display the elapsed time
fprintf(fileID,'Elapsed time for simulation: %.2f seconds\n', elapsedtime);
fprintf(fileIDA,'Elapsed time for simulation: %.2f seconds\n', elapsedtime);
fprintf(fileIDS,'Elapsed time for simulation: %.2f seconds\n', elapsedtime);

fclose(fileID);
fclose(fileIDA);
fclose(fileIDS);

fprintf('Elapsed time for simulation: %.2f seconds\n', elapsedtime);

rmpath(auxiliary);
