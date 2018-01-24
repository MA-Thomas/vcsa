function [ hyp2,lenscale_matrix ] = train_Kernel( cov_iter,covfunc,hyp2,lenscale_matrix,iter,inputs,objVals )

        % --------      Training.   ---------------------------------------
        % --------(i.e. optimizing kernel hyperparameters) ----------------
        
        % Gaussian Likelihood so exact inference can be done.    
        likfunc = @likGauss;
        Dim = size(inputs,2);
        try_again = 0;
        optIters = 100;
        while try_again == 0 && optIters > 0
            try
                pstruct.method = 'CG';
                pstruct.length = -260;
                pstruct.verbosity = 2;
                pstruct.mem = 4;
                [hyp2,fx] = minimize(hyp2, @gp, pstruct, @infExact, [], covfunc, likfunc, inputs, objVals);
%                 hyp2 = minimize(hyp2, @gp, -optIters, @infExact, [], covfunc, likfunc, inputs, objVals);
                try_again = 1;
                if length(hyp2.cov) >= Dim; lenscales = exp(hyp2.cov(1:Dim));else lenscales = exp(hyp2.cov(1));end
                display(['iter is: ', num2str(iter)])
                lenscale_matrix(iter,:) = lenscales'; % save as row vector
                
                % If hyper-parameter optimization fails by returning NaNs,
                % reduce the number of function evaluations and RESET the
                % hyperparameters.
                if any(isnan(hyp2.cov))
                    display('Try again')
                    optIters = optIters - 5;
                    pstruct.length = -optIters;
                    
                    [covfunc,hyp2] = set_Kernel(cov_iter,Dim);
                    try_again = 0;
                end
                
                % If function evaluations start becoming NaN, reduce the
                % number of function evaluations and RESET the hyp-params.
                if any(isnan(fx))
                    optIters = optIters - 5;
                    pstruct.length = -optIters;
                    
                    [covfunc,hyp2] = set_Kernel(cov_iter,Dim);
                    try_again = 0;
                end
                
            catch ME
                warning('MT: problem likely that PD no longer applicable, reducing optIters')
                optIters = optIters - 5;
                pstruct.length = -optIters;
                
                display('rethrowing error message from TRY block')
                rethrow(ME)
                
                [covfunc,hyp2] = set_Kernel(cov_iter,Dim);
                try_again = 0;                
            end
        end
    


end

