%%%%function [ store_reEval_info,reEvalPoints,toDeletePoints ] = decide_points_to_resample(homeFolder) % DONT USE THIS VERSION
% function [ store_reEval_info,store_toDelete_info ] = decide_points_to_resample(homeFolder) % USE THIS VERSION

%NOTE: this function requires the following auxialliary .m files:
% set_Kernel, train_Kernel,
% returnBinaryNotCloseVecsISO,returnBinaryNotCloseVecsARD,get_hyperparams

% load('all_x_20rounds.mat')
% load('all_y_20rounds.mat')  % % %  inputs and function evaluations (x and y) for all prev rounds should be stored in .mat files
% load([homeFolder,'paramFile.mat']);

% DESSA search. Assumes all_x_20rounds.mat has been loaded from the
% lanec1 folder '/CCMV_syn'. Also all_y_20rounds.mat.
% The former contains cell arrays x_cell (each cell contains a matrix of 
% size [#params,dim] at corresponding round), and x_cell_kernels_kappas (each cell
% corresponds to the cells of x_cell, i.e.  a matrix of size [#params,2] where 
% first col specifies the kernel used, second col specifies the kappa
% value. Cells in y_cell contain corresponding RMSDs and stds.

% -------------------------------------------------------------------------
% Major User-Defined Variables:
% -------------------------------------------------------------------------

% Number of suitably distinct parameter vectors to plot initially.
% Also, the points to resample will come from among the top scoring
%       <numParamsToShow> of all parameters. 
numParamsToShow = 50; 

% Maximum Number of points we will run extra simulations at.
numReEval = 16;

% -------------------------------------------------------------------------
% Non-User-Defined Parameters:
% numSurviving - The number of cluster centers in the <numParamsToShow>
%                best scoring parameters. There is no predefined number of
%                clusters, rather clustering is done using a distance
%                measure (i.e., there are multiple distance measures and
%                "closeness" is done by majority vote).
%                As defined below, numReEval = min(numReEval,numSurviving).
% -------------------------------------------------------------------------
%%
% PLOT THE <numParamsToShow> BEST POINTS

x_cell_kernels_kappas
x_cell
numRounds = length(x_cell);
currplot = 0;
current_best_RMSD = [inf inf];
best_Kernel = [];
best_kappa = [];
best_RMSDs = [inf];
for i = 1:numRounds 
    if i == 1; continue;end;
    if isempty(x_cell_kernels_kappas{i});continue;end;

    % Maintain vector of best RMSDs found as of the current round, and which
    % round/kernel/kappa values were used to get it.
    best_RMSDs(end+1,1) = min([current_best_RMSD(1,1), min(y_cell{i}(:,1))])

    if ~isequal(current_best_RMSD(1,1),best_RMSDs(end))

       [M,I] = min(y_cell{i}(:,1))

       best_Round = i;
       best_Kernel = x_cell_kernels_kappas{i}(I,1);
       best_kappa = x_cell_kernels_kappas{i}(I,2);
       current_best_RMSD = [best_RMSDs(end,1) y_cell{i}(I,2)];
    end
end


% The plot will show the parameter vectors whose RMSDs are closest to
% the minimum RMSD together with the minimum RMSD.
ncols = 6;
nrows = ceil(numParamsToShow/ncols);

numInitTrain = size(x_cell{1},1);%71;

x_cell_to_Matrix = cell2mat(x_cell(:));
x_cell_to_Matrix = x_cell_to_Matrix(numInitTrain+1:end,:); 
y_cell_to_Matrix = cell2mat(y_cell(:));

globalMin = best_RMSDs(end);
[~,Inds] = sort(y_cell_to_Matrix(:,1) - globalMin)
sortedRMSDs = [y_cell_to_Matrix(Inds,1), y_cell_to_Matrix(Inds,2)];

sortedParams = x_cell_to_Matrix(Inds,:);

[sortedParamsUnique,ia,ic] = unique(sortedParams,'rows','stable');
sortedRMSDsUnique = [y_cell_to_Matrix(Inds(ia),1),...
    y_cell_to_Matrix(Inds(ia),2)];

% The first row is the global min so far.
closestParams_to_globalMin = sortedParamsUnique(1:numParamsToShow,:);

%%
color_mat = rand(numParamsToShow,3);
grouplabels = {};
for gl = 1:numParamsToShow
    % OPTION 1:
grouplabels{end+1} = ['RMSD = ',num2str(sortedRMSDsUnique(gl,1),'%1.3e'),...
    ' (+- ',num2str(sortedRMSDsUnique(gl,2),'%1.3e'),')']

end
figure
for i = 1:numParamsToShow
    ax = subplot(nrows,ncols,i) %numParamsToShow==9
    
    varNames2 = {'p1','p2','p3','p4','p5','p6','p7','p8','p9',' p10'};
%     varNames2 = {'p1',' ',' ',' ',' ',' ',' ',' ',' ','p10'};
    X = closestParams_to_globalMin(i,:)

    h = parallelcoords(ax,X,'group',grouplabels{i},'labels',varNames2)
    ylabel('Value')

    axis(ax,[0 inf -8 6])
    set(gca,'fontsize',13)

end

% Make sure Documents/MATLAB/mtit is on path
maintitle = [num2str(numParamsToShow),' Best Scoring Parameters (excluding exact duplicate parameters)'];
p=mtit(maintitle,...
     'fontsize',20,'color',[1 0 0],...
     'xoff',-.1,'yoff',.025);
% refine title using its handle <p.th>
set(p.th,'edgecolor',.5*[1 1 1]);
%%
% DETERMINE POINTS TO KEEP.

% Here is where we remove points "close" to the current globalMin.
% And then iteratively only keep points far from the next best point.
params_to_keep(1,:) = closestParams_to_globalMin(1,:);
remaining_params = closestParams_to_globalMin(2:end,:);

lengthscale_cell = get_hyperparams(x_cell_to_Matrix, y_cell_to_Matrix);

numKerns = length(lengthscale_cell);
dim = size(params_to_keep,2);


% While there are still parameters to consider (i.e. while length of
% <remaining_params> is greater than 0), determine which of the remaining
% params will be added to <params_to_keep> based on proximity to the latest
% entry of <params_to_keep>.
% When <remaining_params> is empty, we're done. We can simply search for
% each entry of the <params_to_keep> from within x_cell.
while ~isempty(remaining_params)
    
    % Determine all points NOT "close" to params_to_keep(end) from among 
    % <remaining_params>. Considering these points, add the best scoring
    % one to <params_to_keep>, and remove it from <remaining_params>.
    
    % binary_not_close_vecs will be a matrix. Cols correspond to kernels.
    % Rows correspond to points.
    binary_not_close_vecs = [];
    for k = 1:numKerns
        display(['Kernel number: ',num2str(k)])
        lengthscale = lengthscale_cell{k};
        
        % Convert real space lengthscale to log-scaled (offset) space
        lengthscale = log(lengthscale)/log(10);
        
        switch length(lengthscale)
            case 1 % Isotropic kernel, single length scale
                binary_not_close_vecs = returnBinaryNotCloseVecsISO(...
                    binary_not_close_vecs,lengthscale,...
                    remaining_params,params_to_keep(end,:));
                
            case dim % ARD kernel, <dim> length scales
                binary_not_close_vecs = returnBinaryNotCloseVecsARD(...
                    binary_not_close_vecs,lengthscale,...
                    remaining_params,params_to_keep(end,:));
        end
    end
    
    % Examine all remaining points. If a point is considered "close" by a
    % majority of kernels (more 0s than 1s), we will consider it close.
    close_points = mean(binary_not_close_vecs,2) < 0.5;
    
    % Alternate def of "close": If at least 1 kernel considers a point
    % close, its close. 
    % THIS WOULD ULTIMATELY LEAD TO FEWER SURVIVING PARAMS.
% % %     close_points = sum(~binary_not_close_vecs,2) >= 1;
    
    indices_to_remove = find(close_points);
    indices_to_keep = find(~close_points);
    
    % Remove the "close" params from <remaining_params>
    remaining_params(indices_to_remove,:) = [];
    
    % Append the best of the remaining entries to <params_to_keep>
    if ~isempty(remaining_params)
        params_to_keep = [params_to_keep; remaining_params(1,:)];
    end
    
    % Compute necessary info and write it to the filename:
    % [homeFolder,'resampled_points_info.mat'] 
    % ('resampled_points_info.mat' as well, for non cluster use)
end
%%
% PLOT surviving points
numSurviving = 16;%size(params_to_keep,1)
ncols = 4;
nrows = ceil(numSurviving/ncols);

color_mat = rand(numSurviving,3);
grouplabels = {};
for gl = 1:numSurviving
    grouplabels{end+1} = num2str(gl);
end

figure
for i = 1:numSurviving
    ax = subplot(nrows,ncols,i) 

    varNames2 = {'p1','p2','p3','p4','p5','p6','p7','p8','p9',' p10'};
    X = params_to_keep(i,:);

    h = parallelcoords(ax,X,'group',grouplabels{i},'labels',varNames2)
    ylabel('Value')

    axis(ax,[0 inf -8 10]);
%     axis(ax,[0 inf -7 9]);
    set(gca,'fontsize',13);

end

% Make sure Documents/MATLAB/mtit is on path
maintitle = [num2str(numSurviving),' Distinct Points to be Re-Evaluated at Lower Noise Level'];
p=mtit(maintitle,...
     'fontsize',20,'color',[1 0 0],...
     'xoff',-.1,'yoff',.025);
% refine title using its handle <p.th>
set(p.th,'edgecolor',.5*[1 1 1]);

%%
% RETURN the follwoing:
% 1. re-eval points (subset of <params_to_keep> to be re-evaluated 
%    with more simulations)
% 2. previously evaluated points "close" to the re-eval points. 
%    Remember to use the same measure of "close" as used above.
%    These will be deleted from the GP training set used in future rounds.


cell_start_indices = 1;
for cs = 1:length(x_cell)
    cell_start_indices = [cell_start_indices,...
        cell_start_indices(end) + size(x_cell{cs},1)] 
end

% % % Re-Eval points
numReEval = min(numSurviving,numReEval); 
reEvalPoints = params_to_keep(1:numReEval,:);
all_x_cell_to_Matrix = cell2mat(x_cell(:));

% Store which cells/position-in-cell the entries of <reEvalPoints>
% correspond to. I.e. store the round/position-in-round
store_reEval_info = cell(1,numReEval);
for ir = 1:numReEval
    currPointAsMatrix = repmat(reEvalPoints(ir,:),size(all_x_cell_to_Matrix,1),1);
    temp = find(sum(currPointAsMatrix == all_x_cell_to_Matrix,2) == dim);

    % If temp contains more than 1 index, then reEvalPoints(ir,:)
    % was evaluted multiple times in previous rounds.
    store_reEval_info{ir} = zeros(length(temp),2);
    for t = 1:length(temp)            
        prev_to_curr_rounds = find((temp(t) - cell_start_indices) >= 0)
        store_reEval_info{ir}(t,1) = prev_to_curr_rounds(end); % store round
        store_reEval_info{ir}(t,2) = temp(t) - cell_start_indices(prev_to_curr_rounds(end)) + 1; % store position in round
    end
end


% % % To-Delete points
remaining_params = sortedParamsUnique; % Initial set to be considered for deletion.
[C,ia,ib] = intersect(reEvalPoints,sortedParamsUnique,'rows')
remaining_params(ib,:) = []; % Get rid of <reEvalPoints> from initial consideration.   
toDeletePoints = []; % To be appended to.

% Append the "close to <reEvalPoints>" to <toDeletePoints>.
for r = 1:size(remaining_params,1)

    % Iterate over <reEvalPoints>
    for irEP = 1:size(reEvalPoints,1)

        binary_not_close_vecs = [];
        for k = 1:numKerns
            display(['Kernel number: ',num2str(k)])
            lengthscale = lengthscale_cell{k};

            % Convert real space lengthscale to log-scaled (offset) space
            lengthscale = log(lengthscale)/log(10);

            switch length(lengthscale)
                case 1 % Isotropic kernel, single length scale
                    binary_not_close_vecs = returnBinaryNotCloseVecsISO(...
                        binary_not_close_vecs,lengthscale,...
                        remaining_params,reEvalPoints(irEP,:));

                case dim % ARD kernel, <dim> length scales
                    binary_not_close_vecs = returnBinaryNotCloseVecsARD(...
                        binary_not_close_vecs,lengthscale,...
                        remaining_params,reEvalPoints(irEP,:));
            end
        end    

        % Examine all remaining points. If a point is considered "close" by a
        % majority of kernels (more 0s than 1s), we will consider it close.
        close_points = mean(binary_not_close_vecs,2) < 0.5;

        % Alternate def of "close": If at least 1 kernel considers a point
        % close, its close. 
        % THIS WOULD ULTIMATELY LEAD TO FEWER SURVIVING PARAMS.
    % % %     close_points = sum(~binary_not_close_vecs,2) >= 1;
        indices_to_delete = find(close_points);
        toDeletePoints = [toDeletePoints;...
            remaining_params(indices_to_delete,:)];            
    end
end
toDeletePoints = unique(toDeletePoints,'rows');    

% Store which cells/position-in-cell the entries of <toDeletePoints>
% correspond to. I.e. store the round/position-in-round
numToDelete = size(toDeletePoints,1);
store_toDelete_info = cell(1,numToDelete);
for ir = 1:numToDelete
    currPointAsMatrix = repmat(toDeletePoints(ir,:),size(all_x_cell_to_Matrix,1),1);
    temp = find(sum(currPointAsMatrix == all_x_cell_to_Matrix,2) == dim);

    % If temp contains more than 1 index, then reEvalPoints(ir,:)
    % was evaluted multiple times in previous rounds.
    store_toDelete_info{ir} = zeros(length(temp),2);
    for t = 1:length(temp)            
        prev_to_curr_rounds = find((temp(t) - cell_start_indices) >= 0)
        store_toDelete_info{ir}(t,1) = prev_to_curr_rounds(end); % store round
        store_toDelete_info{ir}(t,2) = temp(t) - cell_start_indices(prev_to_curr_rounds(end)) + 1; % store position in round
    end
end

resampled_points = cell2mat(store_reEval_info(:));
close_points_to_remove = cell2mat(store_toDelete_info(:));
run_compObj_resampledPoints = 1; % tells compObj.m to use the alternate code in the resampling case.
save('resampled_points_info.mat','resampled_points',...
    'close_points_to_remove','run_compObj_resampledPoints') % For use in compObj.m and inLoop2_BayesianOpt_GP.m
x = reEvalPoints;
xSize = size(x,1);
save('xFile_resamplePoints.mat','x','xSize') % For use in inLoop1.m
% save([homeFolder,'resampled_points_info.mat'],'resampled_points','close_points_to_remove')
% Return
%end