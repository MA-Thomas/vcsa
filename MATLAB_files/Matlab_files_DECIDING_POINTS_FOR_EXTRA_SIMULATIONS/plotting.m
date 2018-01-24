% Number of suitably distinct parameter vectors to plot initially.
% Also, the points to resample will come from among the top scoring
%       <numParamsToShow> of all parameters. 
numParamsToShow = 50; 


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


% PLOTTING
ax1 = subplot(1,2,1)
ax1 = errorbar(y_cell_to_Matrix(:,1),y_cell_to_Matrix(:,2),'-s','MarkerSize',7,...
'MarkerEdgeColor','red','MarkerFaceColor','red')
xlabel(''); ylabel('RMSD');
addpath('C:\Users\Marcus\Documents\MATLAB\mtit')
maintitle = 'RMSDs of Successive AF Minimizers';
p=mtit(maintitle,...
     'fontsize',15,'color',[1 0 0],...
     'xoff',-.1,'yoff',.025);
% refine title using its handle <p.th>
set(p.th,'edgecolor',.5*[1 1 1]);
% errorbar(sortedRMSDsUnique(:,1),sortedRMSDsUnique(:,2),'-s','MarkerSize',10,...
% 'MarkerEdgeColor','red','MarkerFaceColor','red')

% ------------------------------------------------------
grouplabel = ['RMSD = ',num2str(y_cell_to_Matrix(end,1),'%1.3e'),...
    ' (+- ',num2str(y_cell_to_Matrix(end,2),'%1.3e'),')']
ax2 = subplot(1,2,2) 

varNames2 = {'p1','p2','p3','p4','p5','p6','p7','p8','p9',' p10'};
X = x_cell_to_Matrix(end,:)
h = parallelcoords(ax2,X,'group',grouplabel,'labels',varNames2)
ylabel('Value')

axis(ax2,[0 inf -8 6])
set(gca,'fontsize',13)
% ------------------------------------------------------

