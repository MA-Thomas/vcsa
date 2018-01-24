function [] = inLoop2_BayesianOpt_GP(homeFolder)
% Loads from and Saves to /home/, not /scratch/

% NOTE: Contents of resampled_points_info.mat:
%       resampled_points ~ [ round, point#_in_round ]
%       close_points_to_remove ~ [ round, point#_in_round ]

display('In inLoop2_BayesianOpt_GP.m')
load([homeFolder,'paramFile.mat']);
load([homeFolder,'iterPrefix.mat']);

load([homeFolder,'xFile.mat']); 
load([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count'); 
display(['inLoop2.m: paramFile.mat,iterPrefix.mat,xFile.mat,',...
    'loop_eval_criteria.mat exist and are loaded.'])

load([homeFolder,'y_iter',num2str(count),'.mat'],'y');
display(['inLoop2_BayesianOpt_GP.m: y_iter',num2str(count),'.mat exists and is loaded.'])
%  --------------- DONE Loading -------------------------------------------    
addpath(genpath('gpml-matlab-v3.6-2015-07-07'))
display('inLoop2_BayesianOpt_GP.m: GP_folder added to path')
tic;
nPara = size(optInfo.paraList,2);    


if exist([homeFolder,'resampled_points_info.mat'])==2
    load([homeFolder,'resampled_points_info.mat'])

    % ALTERNATIVE BEHAVIOR.
    % If there were extra simulations run at previously evaluated points...
    if run_compObj_resampledPoints == 1
        
        % Reset run_compObj_resampledPoints and save the reset value.
        run_compObj_resampledPoints = 0    
        save_resampled_points_info = 1;
        
        % We need to know how many points were simulated in each round.
        % round 1 starts at index 1. Round 2 starts after init training data,
        % etc.
        start_index_each_round = 1;
        for yc = 1:count-1
            load([homeFolder,'y_iter',num2str(yc),'.mat'],'y');
            start_index_each_round(yc+1) = start_index_each_round(yc) + size(y,1);
            clear y
        end
        % Reload y from current round.
        load([homeFolder,'y_iter',num2str(count),'.mat'],'y');

        % resampled_points will tell which points in optInfo.offset and
        % optInfo.score need to be updated (rather than appending x,y).
        % NOTE: resampled_points ~ [ round, point#_in_round ]
        for ir = 1:size(resampled_points,1)
            round = resampled_points(ir,1);
            pInR = resampled_points(ir,2);
            index_to_alter = start_index_each_round(round) + (pInR-1);
            optInfo.offset(index_to_alter,:) = x(ir,:);
            optInfo.score(index_to_alter,:) = y(ir,:);
        end    

        % Remove points "close to" those we just updated, which weren't 
        % simulated sufficiently many times in their previous rounds.
        indices_to_alter = zeros(1,size(close_points_to_remove,1));
        for iclose = 1:size(close_points_to_remove,1)
            round = close_points_to_remove(iclose,1);
            pInR = close_points_to_remove(iclose,2);
            indices_to_alter(iclose) = start_index_each_round(round) + (pInR-1);
        end
        optInfo.offset(indices_to_alter,:) = [];
        optInfo.score(indices_to_alter,:) = [];
    else
        % DEFAULT BEHAVIOR.
        optInfo.offset = [optInfo.offset; x];
        optInfo.score = [optInfo.score; y];        
    end
    
else
    % DEFAULT BEHAVIOR.
    optInfo.offset = [optInfo.offset; x];
    optInfo.score = [optInfo.score; y];
end




s = 9; %0.5; % search box enclosing sphere the points were orig sampled from. %optInfo.initGridSize; 

display('InLoop2.m: Making prediction...');

% Skip the first index if this corresponds to ground-truth data.
start = 1; 
objVals = optInfo.score(start:end,1);        
objVals = objVals./(1*10^9); % smaller magn of obj vals works better (change of units).

objVals_std = optInfo.score(start:end,2);
objVals_std = objVals_std./(1*10^9); % smaller magn of obj vals works better (change of units).        
inputs = optInfo.offset(start:end,:);
Dim = size(inputs,2);

% Aquisition function, F_aq, has explore/exloit tradeoff parameter, kappa.
kappa = [0 1 3]; %[0 3 1000]; %[0 1 3] %[1/36, 1/12, 1/3, 1, 2]

% Gaussian Likelihood so exact inference can be done.    
likfunc = @likGauss;

lenkappa = length(kappa);
numKerns = 7;
X_new_kernels = []; % First <numKerns> rows correspond to first kappa val, etc.


for cov_iter = 1:numKerns
    
    % SET WHICH KERNEL FUNCTION TO USE.
    [covfunc,hyp2] = set_Kernel(cov_iter,Dim);    

    % TRAIN KERNEL
    % Initialize lenscale_matrix for the current <cov_iter>.
    lenscale_matrix = zeros(numIters,Dim);
    % Optimize the current kernel hyperparams with multiple random
    % restarts. Store them in <lenscale_matrix>.
    for iters = 1:numIters
    [ hyp2,lenscale_matrix ] = train_Kernel( cov_iter,...
        covfunc,hyp2,lenscale_matrix,iters,inputs,objVals );
    end
        
    nums = 500;
    X_store = zeros(nums,Dim,length(kappa));%[];
    parfor i = 1:nums
        % Using the new GP model (mean, cov, and updated hyperparameters, 
        % evaluate a large grid of test points (x_test, F_test), and feed
        % everything to the aquisition function optimization.
        % Test points consist of: 1) random points in search region.
        % 2) noisy approximations to the given parameters.    
        m = 100000; % number of points (xSize)
        r = s;  % radius of offset/parameter hypersphere.
        x_test = randsphere(m,nPara,r);
        noise_matrix = (.5/10)*randn(size(inputs));
        noise_matrix2 = (.5/10)*randn(size(inputs));
        x_test = [x_test; inputs+noise_matrix; inputs+noise_matrix2];
        [ymu ys2 fmu fs2] = gp(hyp2, @infExact, [], covfunc, likfunc, inputs, objVals, x_test);

        % Concatenate simulated data, and GP predictive data.
        means_f = [objVals; ymu];
        stds_f = [objVals_std; ys2.^(1/2)];
        x_total = [inputs; x_test];    

        % Evaluate Aquisition Functions, F.
        % Store minimizers.
        for j = 1:lenkappa
            F = means_f - kappa(j)*stds_f; % function values
            [M,I] = min(F);
            X_store(i,:,j) = x_total(I,:);
        end            

    end

    % Calc avg minimizer for current kernel, at each kappa val.
    X_kappa = zeros(lenkappa,Dim);
    for j = 1:lenkappa
        X_kappa(j,:) = mean( X_store(:,:,j) );
    end

    % Store the avg minimizers across kappa values for current kernel.
    X_new_kernels = [X_new_kernels; X_kappa];
end


% % Save everything for next round.
fprintf('\n Argmin(s) of Aquisition Function:\n');
x = X_new_kernels
xSize = size(x,1); 
toc
%%{
fprintf('Iteration %d finished.\n\n',count);

save([homeFolder,'xFile.mat'],'x','xSize');
save([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count')

if exist('save_resampled_points_info') && save_resampled_points_info == 1
    save([homeFolder,'resampled_points_info.mat'],'resampled_points',...
        'close_points_to_remove','run_compObj_resampledPoints')
    display('In inLoop2. Just saved updated variable,run_compObj_resampledPoints=0, to resampled_points_info.mat')
end

display('first save done (xFile and loop...mat).')
% Save another copy of loop_eval to be used in PI_script.sh (source)
system(['[ -f ',homeFolder,'loop_eval_criteria.txt ] && rm -r ',...
    homeFolder,'loop_eval_criteria.txt']);
formatspec = '%s \n';
fid_loop = fopen([homeFolder,'loop_eval_criteria.txt'],'w');
fprintf(fid_loop,formatspec,['minScore=',num2str(minScore)]);
fprintf(fid_loop,formatspec,['noImprove=',num2str(noImprove)]);
fprintf(fid_loop,formatspec,['count=',num2str(count)]);
fclose(fid_loop);

pF = ['paramFile',num2str(count),'_post_inLoop2.mat'];
save([homeFolder,'storage/',pF],'fileInfo','optInfo','qInfo','simInfo');

save([homeFolder,'paramFile.mat'],'fileInfo','qInfo','simInfo','optInfo');
%}   
    
    
% % OLD CODE    
%{    
    tempFile = [fileInfo.homeWorkFolder,fileInfo.prefix,'.mat']
% %     tempFile = 'temp.mat';    

    if ~isequal( size(y,1),size(x,1) )
        display('problem with sizes of x and y. Quiting now.')
        assert( isequal( size(y,1),size(x,1)) ) 
    end
    optInfo.score = [optInfo.score;y]
    optInfo.offset = [optInfo.offset;x]
% % %     save([homeFolder,'paramFile.mat'],'fileInfo','qInfo','simInfo','optInfo');


    display('InLoop2.m: Making prediction...');

    lb = ones(nPara,1) * optInfo.lb * s;
    ub = ones(nPara,1) * optInfo.ub * s;
    dx = (ub - lb) * optInfo.eps;
    
    %MT SNOBFIT is evaluating the aquisition function. We want it to return argmin(s).
    %MT Sometimes, it returns slightly more points. Just be aware of this.
    numPointsToReturn = 2 
    params = struct('bounds',{lb,ub},'nreq',numPointsToReturn,'p',0.5);
    

    % Evaluate the Aquisition Function for input into SNOBFIT.
    % Set F_aq equal to a matrix of values.
    % The first col should be F_aq values, the second col is uncertainties.
    % Set the second col (as in p7 of SNOBFIT paper) to nonpositive vals.
    % This will set the uncertainties to machine precision epsilon which 
    % is fine bc the aquisition function already captures the true 
    % uncertainty.
    F_aq = {};
    newPoints = [];    
    
    % Skip the first index if this corresponds to ground-truth data.
    start = 1        
    objVals = optInfo.score(start:end,1)        
    objVals = objVals./(1*10^9) % smaller magn of obj vals works better (change of units).

    objVals_std = optInfo.score(start:end,2)
    objVals_std = objVals_std./(1*10^9) % smaller magn of obj vals works better (change of units).        
    inputs = optInfo.offset(start:end,:)
    Dim = size(inputs,2);

    % Gaussian Likelihood so exact inference can be used.
    % Kernel is Matern 5/2.
    % Optimize hyperparameters w.r.t. simulated data.
    likfunc = @likGauss;
    covfunc = {@covMaternard,5};
    hyp2.cov = zeros(Dim+1,1);
    hyp2.lik = log(0.1);
    hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, inputs, objVals);
    hyp2.cov
    
    % Using the new GP model (mean, cov, and updated hyperparameters, 
    % evaluate a large grid of test points (x_test, F_test), and feed
    % everything to the aquisition function optimization.
    % Test points consist of: 1) random points in search region.
    % 2) noisy approximations to the given parameters.    
    m = 200 % number of points (xSize)
    r = s  % radius of offset/parameter hypersphere.
    x_test = randsphere(m,nPara,r);
    noise_matrix = (.5/10)*randn(size(inputs));
    noise_matrix2 = (.5/10)*randn(size(inputs));
    x_test = [x_test; inputs+noise_matrix; inputs+noise_matrix2];
    [ymu ys2 fmu fs2] = gp(hyp2, @infExact, [], covfunc, likfunc, inputs, objVals, x_test);

    % Concatenate simulated data, and GP predictive data.
    means_f = [objVals; ymu];
    stds_f = [objVals_std; ys2.^(1/2)];
    x_total = [inputs; x_test];
    
    % Aquisition function, F_aq, has explore/exloit tradeoff parameter, kappa.
    kappa = [1/36, 1/12, 1/3, 1, 2]
    
    for j = 1:length(kappa)

        display(['j is ',num2str(j)])
        
        % Evaluate Aquisition Function.
        F = means_f - kappa(j)*stds_f; % function values
        F_aq{j} = [F,-1*ones(size(F))]; % snobfit wants uncertainties too

        try            
            % Column nPara+2 of n will show the point class (1,2, or 3 mean 
            % local. 4,or 5 is more global. See p.7 of snobfit paper.)
            % Column nPara+1 will show an estimated obj val at the point. 
            [n,~,~] = snobfit(tempFile, x_total, F_aq{j}, params, dx)
            newPoints = [newPoints; n(:,1:nPara)];
            delete(tempFile)
        catch ME
        end
        if ~exist('newPoints','var')
           ME.identifier
           ME.message
           display('QUITING inLoop2.m NOW')
           quit
        end

    end


    fprintf('\n Argmin(s) of Aquisition Function:\n');
    x = newPoints
    xSize = size(x,1);    
%%{
    fprintf('Iteration %d finished.\n\n',count);

    save([homeFolder,'xFile.mat'],'x','xSize');
    save([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count')
    
    display('first save done (xFile and loop...mat).')
    % Save another copy of loop_eval to be used in PI_script.sh (source)
    system(['[ -f ',homeFolder,'loop_eval_criteria.txt ] && rm -r ',...
        homeFolder,'loop_eval_criteria.txt']);
    formatspec = '%s \n';
    fid_loop = fopen([homeFolder,'loop_eval_criteria.txt'],'w');
    fprintf(fid_loop,formatspec,['minScore=',num2str(minScore)]);
    fprintf(fid_loop,formatspec,['noImprove=',num2str(noImprove)]);
    fprintf(fid_loop,formatspec,['count=',num2str(count)]);
    fclose(fid_loop);
    
    save([homeFolder,'paramFile.mat'],'fileInfo','qInfo','simInfo','optInfo');
%}
display('Post inLoop2_BayesianOpt_GP.m')
end