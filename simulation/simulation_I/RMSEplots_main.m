%For 2D
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
    bdpval = bdp(i_bdp);
    for i_outp = 1:length(alphavals)
        out_p = alphavals(i_outp);
        for i_ncp = 1:length(ncp)
            ncpval = ncp(i_ncp);
            figname = ['RMSE2D_','vareps_', num2str(out_p),'_ncp_', num2str(ncpval), '_bdp_', num2str(bdpval)];
            filename_R = ['R_',num2str(bdpval),'_',num2str(out_p),'_',num2str(ncpval),'.mat'];
            filename_A = ['A_',num2str(bdpval),'_',num2str(out_p),'_',num2str(ncpval),'.mat'];
            filename_S = ['S_',num2str(bdpval),'_',num2str(out_p),'_',num2str(ncpval),'.mat'];
            

            cd ./SimResults/
            load(filename_R)
            load(filename_A)
            load(filename_S)
            cd ..\

            makeplots_rmse2D(R_mean, A_mean, S_mean, figname)
        end
    end
end

rmpath(auxiliary);
