function lenscale_min_cell = get_hyperparams(x_cell_to_Matrix, y_cell_to_Matrix)

    % Skip the first index if this corresponds to ground-truth data.
    start = 1; 
    objVals = y_cell_to_Matrix(start:end,1);        
    objVals = objVals./(1*10^9); % smaller magn of obj vals works better (change of units).

    objVals_std = y_cell_to_Matrix(start:end,2);
    objVals_std = objVals_std./(1*10^9); % smaller magn of obj vals works better (change of units).        
    inputs = x_cell_to_Matrix(start:end,:);
    Dim = size(inputs,2);

    % Gaussian Likelihood so exact inference can be done.    
    likfunc = @likGauss;
    
    numKerns = 7;
    numIters = 5; % number of random restarts for hyperparam opt
    lenscale_min_cell = cell(1,numKerns);
    lenscale_mean_cell = cell(1,numKerns);
    lenscale_std_cell = cell(1,numKerns);
    lenscale_mode_cell = cell(1,numKerns);
% ------------------------------------------------------------------------
    for cov_iter = 1:numKerns
        
        % SET WHICH KERNEL FUNCTION TO USE.
        [covfunc,hyp2] = set_Kernel(cov_iter,Dim);
        
        % Initialize lenscale_matrix for the current <cov_iter>.
        lenscale_matrix = zeros(numIters,Dim);
        
        % Optimize the current kernel hyperparams with multiple random
        % restarts. Store them in <lenscale_matrix>.
        for iters = 1:numIters
        [ hyp2,lenscale_matrix ] = train_Kernel( cov_iter,...
            covfunc,hyp2,lenscale_matrix,iters,inputs,objVals );
        end
        
        % Update lenscale_cell{cov_iter} with the min lengthscale 
        % (of each dimension) from the repeated hyperparameter 
        % optimizations - each repeated set of hyperparameter length 
        % scales is a row of <lenscale_matrix_updated>.
        lenscale_min_cell{cov_iter} = min(lenscale_matrix,[],1);
        lenscale_mean_cell{cov_iter} = mean(lenscale_matrix,1);
        lenscale_std_cell{cov_iter} = std(lenscale_matrix,1,1);
        lenscale_mode_cell{cov_iter} = mode(lenscale_matrix,1);

    end

end
