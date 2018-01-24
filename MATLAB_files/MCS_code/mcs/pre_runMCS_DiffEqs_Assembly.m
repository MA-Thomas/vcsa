%MT  pre_runMCS_DiffEqs_Assembly
% This script will run the code from DIFF_EQ_MAIN_SCRIPT.m that's
% necessary for the generation of RMSDs by the MCS based function
% saxs_Diff_Eq_MT.m used by the MCS caller. 


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
global global_a
global global_b
global global_s 
global global_time_points
global global_t_gt
global global_n_mers_cell
global global_face_centers
global global_crysolIq
global global_expCurve

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
% r = 100;%50;
% m = 200;
nPara = 12; 
% initial_params = abs(randsphere(m,nPara,r) + k_on*ones(m,nPara));

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
next_for = 0;

% numKerns = 9;% When using the 2 additive kernels in add. to the 7 normal kernels %7;
%numKerns = 1; % temporary. setKernel() also temp altered.

u = 1*ones(nPara,1)
v = 200*ones(nPara,1)

global_a = a;
global_b = b;
global_s = s;
global_time_points = time_points;
global_t_gt = t_gt;
global_n_mers_cell = n_mers_cell;
global_face_centers = face_centers;
global_crysolIq = crysolIq;
global_expCurve = expCurve;
