auxiliary = "../../auxiliary";

addpath(auxiliary);


bdp = [0.25 0.35];
ncp = [25, 50];

for i_bdp = 1:length(bdp)
    if i_bdp == 1
        alphavals = [0.0, 0.05, 0.10, 0.20];
    else
        alphavals = [0.0, 0.1, 0.2, 0.3, 0.4];
    end
    %alphavals = 0.2; 
    bdpval = bdp(i_bdp);
    for i_outp = 1:length(alphavals)
        out_p = alphavals(i_outp);
        for i_ncp = 1:length(ncp)
            ncpval = ncp(i_ncp);
            figname = ['RMSE6D_','vareps_', num2str(out_p),'_ncp_', num2str(ncpval),'_bdp_', num2str(bdpval)];
            filename_R = ['R_err','_',num2str(bdpval),'_',num2str(out_p),'_',num2str(ncpval),'.mat'];
            filename_A = ['A_err','_',num2str(bdpval),'_',num2str(out_p),'_',num2str(ncpval),'.mat'];
            filename_S = ['S_err','_',num2str(bdpval),'_',num2str(out_p),'_',num2str(ncpval),'.mat'];
            

            cd .\SimResultsMSEs\
            load(filename_R)
            load(filename_A)
            load(filename_S)
            cd ..\

            %makeplots_rmse2(R_mean, A_mean, S_mean, figname)
            plot_convergence_rates_rob_new(R_ast_errors_0, R_ast_errors_1, R_ast_errors_2, R_ast_errors_3, A_errors_0, ...
                A_errors_1, A_errors_2, A_errors_3, Sigma_errors_0, Sigma_errors_1, Sigma_errors_2, Sigma_errors_3, figname)
        end
    end
end

rmpath(auxiliary);