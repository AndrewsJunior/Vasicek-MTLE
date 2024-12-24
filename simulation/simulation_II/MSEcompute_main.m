addpath(auxiliary);

%% Simulation - MSE
%rng(1);

R0 = [-0.5; -0.05; 0.5; 0.6; -1.0; -1.3];

R_ast_true = [-0.5566; -0.0391; 0.3934; 0.4322; -1.0302; -1.0531];

A_true     = [0.6383   -0.0029    0.0503   -0.1402   -0.0285    0.0599;
              -0.0757    0.7442   -0.0803    0.0281   -0.0602    0.0973;
               0.0079   -0.0273    0.7957   -0.2526   -0.0333    0.0086;
              -0.0254   -0.0172   -0.0354    0.5785    0.0047   -0.0029;
               0.0584    0.0133   -0.1031    0.0685    0.6711   -0.1734;
               0.0902   -0.0054    0.0110   -0.0608   -0.1643    0.5691];

Sigma_true = [0.1246    0.0884    0.0297    0.0317    0.0458    0.0628;
              0.0884    0.1124    0.0252    0.0267    0.0344    0.0444;
              0.0297    0.0252    0.1197    0.1131    0.0118    0.0246;
              0.0317    0.0267    0.1131    0.1310    0.0106    0.0262;
              0.0458    0.0344    0.0118    0.0106    0.0935    0.0622;
              0.0628    0.0444    0.0246    0.0262    0.0622    0.1241];

T_array = 60:30:360;%60:30:360;

nrep = 5000; % Reduce to 500 or 50 for faster runs

bdp = [0.25, 0.35];
ncp = [25, 50]; 


tic;
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
            compute_MSEs_rob(R_ast_true, A_true, Sigma_true, R0, T_array, nrep, bdpval, out_p, ncpval);
        end
    end
end
 %compute_MSEs_rob(R_ast, A, Sigma, R0, T_array, nrep, Outlierfrac, out_pcent, ncp)

elapsedtime = toc;
% Display the elapsed time
fprintf('Elapsed time for simulation: %.2f seconds\n', elapsedtime);

rmpath(auxiliary);
