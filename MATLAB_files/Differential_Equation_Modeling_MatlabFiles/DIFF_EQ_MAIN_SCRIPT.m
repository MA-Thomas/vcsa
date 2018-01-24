%{
Notes:

Using parfor only over the diff params at the current iteration works.

%}

% RS1 version of this code
addpath(genpath('/home/marcust/Capsid_Differential_Eq_FILES'))


% define the constant matrices a,b,O.
% a: forward reaction degeneracies
% b: backward reaction degeneracies
% etc.
dodecahedron_monomer_addition
save('dodecahedron_a_b_n_mers_cell_face_centers.mat','a','b','n_mers_cell','face_centers')
%
% ---------------------------------------------------
% ---------------------------------------------------
% First add the monomer/zero-mer row to 'a' and to 'b'. Degeneracy is of
% course 0 for each.
load('dodecahedron_a_b_n_mers_cell_face_centers.mat')

b = [[1 1 0 0 0];b];
a = [[0 0 1 1 0];a];
 
% % Compute s(j,k,m,n)
delG = 3.5 % Units: kcal/mole
RT = 0.593 % Units: kcal/mole. See https://en.wikipedia.org/wiki/KT_(energy) 
s = zeros(size(a));
for v = 2:size(a,1)
    s(v,[1:4]) = a(v,[1:4]);
    
    m_ind = s(v,3);
    n_ind = s(v,4);
    j_ind = s(v,1);
    k_ind = s(v,2);
    
    c_diff = c(m_ind,n_ind) - c(j_ind,k_ind);
    s(v,5) = exp( -delG * c_diff / (RT) );
    
end
%%
tic;
% ------------- Ground Truth settings, V.1 ------------
% % parameters = [10^-1, 10^2];
% init_conc(1) = 1*10^-1;
% tspan = [0 100]
% k_on = 10^2;

% Perhaps the parameter inference can assume a model 
% in which there are multiple oligomer dependent reaction rates
% (which, when translated into the GT model, are all equal to k_on.)
% ------------- Ground Truth settings, V.1 ------------

% How many species are there? Count the unique (m,n) in 'b' (or in 'a').
[unique_mn_rows,ia,ic] = unique( b(:,[1,2]),'rows' );
init_conc = zeros(length(ia),1);


init_conc(1) = 1 * 10^-1; % 10uM in units of M is 1*10^-5.
tspan = [0 1]
k_on = 10^2;

% NOTE: In my model, k_on is completely independent of the concentration.
%       This, coupled with the fact that monomer concentration never
%       actually reaches zero, can lead to unrealistic behavior. (such as
%       ending up with more complete capsids at the final time point than
%       the initial monomer concentration allows for.
options = odeset('Stats','on')
time_points = [0:1e-4:1e-3,2e-3:1e-3:1];
sol = ode45(@(t,conc) dodec_diff_eq(t,conc,a,b,s,k_on),time_points,...
    init_conc,options);
t = time_points;
conc = deval(sol,t)'
% % % [t,conc] = ode45(@(t,conc) dodec_diff_eq(t,conc,a,b,s,k_on),tspan,...
% % %     init_conc);

% Because the diff eqs deal in concentrations and not directly in number,
% we need a way to translate the conc matrix to a number matrix.

% STRATEGY 1: enforce that all simulations produce a single complete
% capsid at the final time point. This provides a way to convert conc to
% numMatrix.
% One issue with this approach is that the numbers of many of the
% intermediates may vary only in non-significant digits.
% % last = conc(end,end);
% % f = 1/last;
% % NUMmatrix = f*conc;

% STRATEGY 2: Determine the median value of each col in conc. Choose the 
% smallest such median value. Determine f by enforing that this value
% be at least 1.

% STRATEGY 3: Since we are dealing with monomer addition, enforce that
% there be at least 1 monomer left at the last time point.
% % last = conc(end,1);
% % f = 1/last;
% % NUMmatrix = f*conc;

% STRATEGY 4: Enforce initial monomer count to be exactly enough for 100
% complete capsids to form. Use this to compute f. Convert conc to
% NUMmatrix and perform floor(NUMmatrix).
% With this strategy, we can truncate all rows of NUMmatrix which have more
% than 100 complete capsids.
numCapsids = 100;
last = conc(1,1);
f = 12*numCapsids/last;
NUMmatrix = f*conc;

NUMmatrix = floor(NUMmatrix);
NUMmatrix(NUMmatrix<0) = 0;

% % % truncate_at = find(NUMmatrix(:,end) > numCapsids)
% % % if ~isempty(truncate_at)
% % %     truncate_at = truncate_at(1)
% % %     NUMmatrix = NUMmatrix(1:truncate_at, :);
% % %     t = t(1:truncate_at);
% % % end

% Save Ground Truth Data --------
% t_gt = t(1:2:end);
t_gt = t;
load('crysol_Iq.mat');crysolIq = crysol_Iq_cc;
SAXScurves = VirusGenSAXS_fromDiffEqs(t,NUMmatrix,t_gt,b,...
    n_mers_cell,face_centers,crysolIq);
SAXScurves_GT = SAXScurves;
save('GroundTruth_SAXS_diffeq.mat','SAXScurves_GT','t','NUMmatrix','t_gt','b',...
    'n_mers_cell','face_centers','crysolIq')
% Save Ground Truth Data --------
toc
%%
% -----------------------------------------------------
% ---------- Parameter Inference ----------------------
% -----------------------------------------------------

%{
Assumptions:
1. The concentration is not a parameter to be searched. It will not change.
2. The parameters to be searched will be the vector of microscopic rates,
k_oligomer.
3. The ground truth assumed a single microscopic rate k_on. Thus, our
parameter inference is successful when it returns a vector with each
element taking the value k_on. 


%}

load('crysol_Iq.mat')
crysolIq = crysol_Iq_cc;

load('GroundTruth_SAXS_diffeq.mat');
expCurve = SAXScurves_GT;

display('yo')
% Initialize search with set of parameter vectors
r = 100;%50;
m = 200;
nPara = 12; 
initial_params = abs(randsphere(m,nPara,r) + k_on*ones(m,nPara));

% Create grid of training data at these regularly spaced points.
% % % range = 99:.1:101;
% % % [grid_a,grid_b] = meshgrid(range, range);
% % % initial_params = [ grid_a(:), grid_b(:) ];

numRounds = 20; % number of search rounds
allRMSDs = [];
inputs = [];
clear new_params
best_input = [];
best_objVal = [];
best_objVal_storage = [];
best_input_storage = [];
numFailedKerns = 0;
next_for = 0;

numKerns = 9;% When using the 2 additive kernels in add. to the 7 normal kernels %7;
%numKerns = 1; % temporary. setKernel() also temp altered.
tic;

%%
% Loop
tic;
lastwarn('')
for nR = 1:numRounds
    
    
    if exist('new_params','var')
        param_matrix = new_params
        
        % There are 12 params, one for each oligomer size. But since there
        % may be multiple oligomer types per size, some params need to be
        % repeated in the rows of k_matrix.
        k_matrix = formatParams(param_matrix,b);
    else
        
        param_matrix = initial_params
        k_matrix = formatParams(param_matrix,b);
    end
    
    %  ----- ITERATE OVER ALL PARAMETER VECTORS ------
    RMSDs = zeros(size(k_matrix,1),1);
    clear conc
    parfor p = 1:size(k_matrix,1)
%     for p = 1:size(k_matrix,1)
        
        % ---- Begin Run Simulations and Compute Objective Values ---------        
        
        k = k_matrix(p,:);     
    
        [unique_mn_rows,ia,ic] = unique( b(:,[1,2]),'rows' );
        init_conc = zeros(length(ia),1);
        init_conc(1) = 1 * 1*10^-1; % Ground Truth value.
        
%         clear conc 
        
        % Simulate assembly.
        options = odeset('Stats','on');
        try
            sol = ode45(@(t,conc) dodec_diff_eq(t,conc,a,b,s,k),time_points,...
                init_conc,options);
            
            % If integration does not proceed normally (warning issued), 
            % manually create error.
            [msgstr, msgid] = lastwarn
            if isequal(msgid,'MATLAB:ode45:IntegrationTolNotMet')                
                %reset last warning
                lastwarn('')                
                % manually create error.
                eye(4) + eye(5)
            end
            
        catch
            display('catch block entered')
            sol = ode23s(@(t,conc) dodec_diff_eq(t,conc,a,b,s,k),time_points,...
                init_conc,options);
            % If integration does not proceed normally (warning issued), 
            % manually create error.
            [msgstr, msgid] = lastwarn
            if isequal(msgid,'MATLAB:ode23s:IntegrationTolNotMet')                
                % Integration again failed due to too small step size.
                % deval will therefore fail if it's called below. 
                % As a workaround, I will manually set this p-th RMSD to
                % inf, continue to next parameter, and then remove this 
                % parameter from consideration below.
                RMSD(p) = inf
                display(['RMSD parameter ', num2str(p),' is: ',num2str(RMSDs(p))])
                continue
            end            
        end        
        
        % Now that simulation is done, interpolate solution at designated
        % time points.
        t = time_points;
        conc = deval(sol,t)';
        
        % Reset init_conc.
        init_conc(:) = 0;        
        
        numCapsids = 100;
        last = conc(1,1);
        f = 12*numCapsids/last;
        NUMmatrix = f*conc;

        NUMmatrix = floor(NUMmatrix);
        NUMmatrix(NUMmatrix<0) = 0;

% % % % %         truncate_at = find(NUMmatrix(:,end) > numCapsids);
% % % % %         if ~isempty(truncate_at)
% % % % %             truncate_at = truncate_at(1);
% % % % %             NUMmatrix = NUMmatrix(1:truncate_at, :);
% % % % %             t = t(1:truncate_at);
% % % % %         end
        
        SAXScurves = VirusGenSAXS_fromDiffEqs(t,NUMmatrix,t_gt,b,...
            n_mers_cell,face_centers,crysolIq);

        RMSDs(p) = calcRMSD(SAXScurves,expCurve);
        display(['RMSD parameter ', num2str(p),' is: ',num2str(RMSDs(p))])
    end
    
    % If any parameter points failed to integrate, remove them from
    % consideration.
    inf_indices = find(isinf(RMSDs));
    RMSDs(inf_indices) = [];
    param_matrix(inf_indices,:) = [];    
    
    %toc
    % --------- End Run Simulations and Compute Objective Values ----------
    % ---------------------------------------------------------------------
    %tic;
    
    
    
    
    % -------------- Begin Gaussian Process Regression --------------------
    allRMSDs = [allRMSDs;RMSDs];
    objVals = allRMSDs;        
    objVals = objVals./(1*10^9); % smaller magn of obj vals works better (change of units).
    
    objVals_std = zeros(size(objVals)); % smaller magn of obj vals works better (change of units).        
    inputs = [inputs;param_matrix];
    Dim = size(inputs,2);

    [best_objVal,indy] = min(objVals)
    best_input = inputs(indy,:)
    best_objVal_storage(end+1) = best_objVal;
    best_input_storage(end+1,:) = best_input;
    
    
    % Aquisition function, F_aq, has explore/exloit tradeoff parameter, kappa.
    kappa = 0;%[0 1 3] Note: There is no noise, so tradeoff param is 0.

    % Gaussian Likelihood so exact inference can be done.    
    likfunc = @likGauss;

    lenkappa = length(kappa);    
    new_params = []; % First <numKerns> rows correspond to first kappa val, etc.

    for cov_iter = 1:numKerns
        
        % SET WHICH KERNEL FUNCTION TO USE.
        [covfunc,hyp2] = set_Kernel(cov_iter,Dim)

        % --------      Training.   ---------------------------------------
        % --------(i.e. optimizing kernel hyperparameters) ----------------
        try_again = 0;
        optIters = 100;
        while try_again == 0 && optIters > 0
            try
                pstruct.method = 'CG';
                pstruct.length = -optIters;
                pstruct.verbosity = 2;
                pstruct.mem = 4;
                [hyp2,fx] = minimize(hyp2, @gp, pstruct, @infExact, [], covfunc, likfunc, inputs, objVals);
%                 hyp2 = minimize(hyp2, @gp, -optIters, @infExact, [], covfunc, likfunc, inputs, objVals);
                try_again = 1;
                
                % If hyper-parameter optimization fails by returning NaNs,
                % reduce the number of function evaluations and RESET the
                % hyperparameters.
                if any(isnan(hyp2.cov))
                    optIters = optIters - 10;
                    pstruct.length = -optIters;
                    
                    [covfunc,hyp2] = set_Kernel(cov_iter,Dim);
                    try_again = 0;
                end
                
                % If function evaluations start becoming NaN, reduce the
                % number of function evaluations and RESET the hyp-params.
                if any(isnan(fx))
                    optIters = optIters - 10;
                    pstruct.length = -optIters;
                    
                    [covfunc,hyp2] = set_Kernel(cov_iter,Dim);
                    try_again = 0;
                end
                
            catch ME
                warning('MT: problem likely that PD no longer applicable, reducing optIters')
                optIters = optIters - 10;
                pstruct.length = -optIters;
                
                display('rethrowing error message from TRY block')
                rethrow(ME)
                
                [covfunc,hyp2] = set_Kernel(cov_iter,Dim);
                try_again = 0;                
            end
        end
        
        if next_for == 1
            next_for = 0;            
            continue
        end
        % --------   End Training.  ---------------------------------------
        % --------(i.e. optimizing kernel hyperparameters) ----------------
        
        
        nums = 20;%100;
        X_store = zeros(nums,Dim,length(kappa));%[];
        F_store = zeros(nums,1,length(kappa));
%         parfor i = 1:nums
        for i = 1:nums %Using parfor here leads to errors in 'minimize.m'
            
            % Using the new GP model (mean, cov, and updated hyperparameters, 
            % evaluate a large grid of test points (x_test, F_test), and feed
            % everything to the aquisition function optimization.
            
            % Include as part of test set randomly chosen points obtained 
            % similarly to the initial parameters.                       
            % 'r' (defined above) radius of offset/parameter hypersphere. 
            m_gp = 10000 % number of points (xSize) 
            x_test = abs(randsphere(m_gp,nPara,r) + k_on*ones(m_gp,nPara));
            
            % Include as part of test set points near the current best
            % parameter. Near means within hypersphere surrounding x_best.
            m_near = 10000
            r_near = sqrt(min(best_input));
            near_best_input = abs(repmat(best_input,m_near,1) ...
                + randsphere(m_near,nPara,r_near));
            
            % Include as part of test set points near the current best
            % parameter. Near means identical in all but 1 dimension, with
            % noise introduced to the remaining coordinate.                      
            near_best_coord_wise = [];
            for le = 1:size(best_input,2)
                temp = repmat(best_input,m_near,1);
                temp(:,le) = near_best_input(:,le);
                
                near_best_coord_wise = [near_best_coord_wise; abs(temp)];
            end
            
            % Include as part of the test set points near the 'new_points'
            % from the last round (if new_params obtained from last round).
            near_new_params = []
            if ~isequal(param_matrix, initial_params)
            for newp = 1:size(new_params,1)
                near_new_params = [near_new_params;...
                    abs( repmat(new_params(newp,:),m_near,1) ...
                    + randsphere(m_near,nPara,r_near) )];
            end
            end
            
            x_test = [x_test; near_best_input; near_new_params; near_best_coord_wise];
            [ymu ys2 fmu fs2] = gp(hyp2, @infExact, [], covfunc, likfunc, inputs, objVals, x_test);

            % Concatenate simulated data, and GP predictive data.
            means_f = [objVals; ymu];
            stds_f = [objVals_std; ys2.^(1/2)];
            x_total = [inputs; x_test];    

            % Evaluate Aquisition Functions, F.
            % Store minimizers. Store minima.
            for j = 1:lenkappa
                F = means_f - kappa(j)*stds_f; % function values
                [M,I] = min(F);
                X_store(i,:,j) = x_total(I,:);
                F_store(i,1,j) = M;
            end            

        end

        % Calc avg minimizer for current kernel, at each kappa val.
        X_kappa = zeros(lenkappa,Dim);
        for j = 1:lenkappa
%             X_kappa(j,:) = mean( X_store(:,:,j) );
            
            % Alternate: instead of taking the avg minimizer, take the 
            % smallest minimizer yet seen.
            [M,I] = min( F_store(:,1,j) );
            X_kappa(j,:) = X_store(I,:,j);
        end

        % Store the avg minimizers across kappa values for current kernel.
        new_params = [new_params; X_kappa];
    end
    % End Gaussian Process Regression --------------------
    
    
end
%toc
% Display results
allRMSDs
inputs

 
m
RMSDs_kernel_1 = allRMSDs(m+1:numKerns:end)
RMSDs_kernel_2 = allRMSDs(m+2:numKerns:end)
RMSDs_kernel_3 = allRMSDs(m+3:numKerns:end)
RMSDs_kernel_4 = allRMSDs(m+4:numKerns:end)
RMSDs_kernel_5 = allRMSDs(m+5:numKerns:end)
RMSDs_kernel_6 = allRMSDs(m+6:numKerns:end)
RMSDs_kernel_7 = allRMSDs(m+7:numKerns:end)
RMSDs_kernel_8 = allRMSDs(m+8:numKerns:end)
RMSDs_kernel_9 = allRMSDs(m+9:numKerns:end)
inputs_kernel_1 = inputs(m+1:numKerns:end,:)
inputs_kernel_2 = inputs(m+2:numKerns:end,:)
inputs_kernel_3 = inputs(m+3:numKerns:end,:)
inputs_kernel_4 = inputs(m+4:numKerns:end,:)
inputs_kernel_5 = inputs(m+5:numKerns:end,:)
inputs_kernel_6 = inputs(m+6:numKerns:end,:)
inputs_kernel_7 = inputs(m+7:numKerns:end,:)
inputs_kernel_8 = inputs(m+8:numKerns:end,:)
inputs_kernel_9 = inputs(m+9:numKerns:end,:)


% ------------- Plot using subplot -------------
subplot(2,5,1)
plot(RMSDs_kernel_1,'b')
hold on
lgd1 = legend('Matern3-ARD')
lgd1.FontSize = 14;

subplot(2,5,2)
plot(RMSDs_kernel_2,'b')
hold on
lgd2 = legend('Matern5-ARD')
lgd2.FontSize = 14;

subplot(2,5,3)
plot(RMSDs_kernel_3,'Color',[ 0.9100 0.4100 0.1700],'LineWidth',1.1)
hold on
lgd3 = legend('RationalQuad-ARD')
lgd3.FontSize = 14;

subplot(2,5,4)
plot(RMSDs_kernel_4,'Color',[ 0.9100 0.4100 0.1700],'LineWidth',.6)
hold on
lgd4 = legend('RationalQuad-ISO')
lgd4.FontSize = 14;

subplot(2,5,5)
plot(RMSDs_kernel_5,'c','LineWidth',1)
hold on
lgd5 = legend('Gabor-ARD')
lgd5.FontSize = 14;

subplot(2,5,6)
plot(RMSDs_kernel_6,'Color',[1 .85 .1],'LineWidth',1.3)
hold on
lgd6 = legend('NN')
lgd6.FontSize = 14;

subplot(2,5,7)
plot(RMSDs_kernel_7,'k','LineWidth',1.3)
hold on
lgd7 = legend('SqExp-ARD')
lgd7.FontSize = 14;

subplot(2,5,8)
plot(RMSDs_kernel_8,'k','LineWidth',1.3)
hold on
lgd7 = legend('Add.-NN-ISO')
lgd7.FontSize = 14;

subplot(2,5,9)
plot(RMSDs_kernel_9,'k','LineWidth',1.3)
hold on
lgd7 = legend('Add.-SqExp-ISO')
lgd7.FontSize = 14;


ax1 = subplot(2,5,1)
scatter(ax1,1:length(RMSDs_kernel_1),RMSDs_kernel_1,11,'b','filled')
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_1); minInputK1 = inputs_kernel_1(I,:);
hold on
scatter(ax1,I,RMSDs_kernel_1(I),40,'g','filled'); % plot best RMSD in green

ax2 = subplot(2,5,2)
scatter(ax2,1:length(RMSDs_kernel_2),RMSDs_kernel_2,11,'b')
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_2); minInputK2 = inputs_kernel_2(I,:);
hold on
scatter(ax2,I,RMSDs_kernel_2(I),40,'g','filled'); % plot best RMSD in green

ax3 = subplot(2,5,3)
scatter(ax3,1:length(RMSDs_kernel_3),RMSDs_kernel_3,11,'MarkerFaceColor',[ 0.9100 0.4100 0.1700],'MarkerEdgeColor',[ 0.9100 0.4100 0.1700])
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_3); minInputK3 = inputs_kernel_3(I,:);
hold on
scatter(ax3,I,RMSDs_kernel_3(I),40,'g','filled'); % plot best RMSD in green

ax4 = subplot(2,5,4)
scatter(ax4,1:length(RMSDs_kernel_4),RMSDs_kernel_4,11,'MarkerEdgeColor',[ 0.9100 0.4100 0.1700])
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_4); minInputK4 = inputs_kernel_4(I,:);
hold on
scatter(ax4,I,RMSDs_kernel_4(I),40,'g','filled'); % plot best RMSD in green

ax5 = subplot(2,5,5)
scatter(ax5,1:length(RMSDs_kernel_5),RMSDs_kernel_5,11,'c','filled','s')
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_5); minInputK5 = inputs_kernel_5(I,:);
hold on
scatter(ax5,I,RMSDs_kernel_5(I),40,'g','filled'); % plot best RMSD in green

ax6 = subplot(2,5,6)
scatter(ax6,1:length(RMSDs_kernel_6),RMSDs_kernel_6,11,'d','MarkerEdgeColor',[1 .8 .1])
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_6); minInputK6 = inputs_kernel_6(I,:);
hold on
scatter(ax6,I,RMSDs_kernel_6(I),40,'g','filled'); % plot best RMSD in green

ax7 = subplot(2,5,7)
scatter(ax7,1:length(RMSDs_kernel_7),RMSDs_kernel_7,11,'k','*')
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_7); minInputK7 = inputs_kernel_7(I,:);
hold on
scatter(ax7,I,RMSDs_kernel_7(I),40,'g','filled'); % plot best RMSD in green

ax8 = subplot(2,5,8)
scatter(ax8,1:length(RMSDs_kernel_8),RMSDs_kernel_8,11,'m','*')
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_8); minInputK8 = inputs_kernel_8(I,:);
hold on
scatter(ax8,I,RMSDs_kernel_8(I),40,'g','filled'); % plot best RMSD in green

ax9 = subplot(2,5,9)
scatter(ax9,1:length(RMSDs_kernel_9),RMSDs_kernel_9,11,'m','*')
xlabel('Iterations','FontSize',15)
ylabel('RMSD','FontSize',15)
[~,I] = min(RMSDs_kernel_9); minInputK9 = inputs_kernel_9(I,:);
hold on
scatter(ax9,I,RMSDs_kernel_9(I),40,'g','filled'); % plot best RMSD in green

axis([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9],[0 inf 0 7e9])

% Plot best input returned by acq func in last round
ax10 = subplot(2,5,10)
X_new = inputs(m+1:end,:);
grouplabels = {'Matern3-ARD','Matern5-ARD','RationalQuad-ARD','RationalQuad-ISO','Gabor-ARD','NN','SqExp-ARD','Add.-NN','Add.-SqExp-ISO'};
varNames2 = {'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12'};
% varNames2 = {'p1','p2','p3','p4','p5','p6'};
% varNames2 = {'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10'};
Kappa = 0;
lenKappa = 1; % i.e. Kappa = [0 1 3]

start = 1;
X_best_diff_rounds = [minInputK1;minInputK2;minInputK3;minInputK4;minInputK5;minInputK6;minInputK7;minInputK8;minInputK9];

X = X_best_diff_rounds
h = parallelcoords(X,'group',grouplabels,'labels',varNames2)
title('best input for each kernel across all rounds')
h(1).LineWidth = .1;
h(2).LineWidth = .1; h(2).LineStyle = '--'; 
h(3).LineWidth = .1;
h(4).LineWidth = .1; h(4).LineStyle = '--'; 
h(5).LineWidth = .1;
h(6).LineWidth = .1; h(6).LineStyle = '--'; 
h(7).LineWidth = .1;
h(8).LineWidth = .3; h(8).LineStyle = '--'; 
h(9).LineWidth = .3;
% Make figure window larger using format [lowerleftcornerx lowerleftcornery width height]
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2)/4 pos(3:4)*1.5])

savename = ['workspace_DIFF_EQ_',num2str(nPara),...
    'param_Init_r_',num2str(r),'_',num2str(length(allRMSDs)),...
    'points','_',num2str(m),'init_',num2str(numRounds),...
    '_More_rounds_Oct14.mat']
save(savename)
toc
%{
Junk Code

% % % %                 if any(isnan(hyp2.cov))
% % % %                     infoString = ['kernel ',num2str(cov_iter),...
% % % %                         ' hyperparam optimization failed in round ',...
% % % %                         num2str(nR),'.mat']
% % % %                     save(infoString)
% % % %                     numFailedKerns = numFailedKerns + 1;
% % % %                     next_for = 1;
% % % %                 end
% % % %                 if any(isnan(hyp2.lik))
% % % %                     infoString = ['kernel ',num2str(cov_iter),...
% % % %                         ' hyperparam optimization failed in round ',...
% % % %                         num2str(nR),'.mat']
% % % %                     save(infoString)
% % % %                     numFailedKerns = numFailedKerns + 1;
% % % %                     next_for = 1;
% % % %                 end

%}