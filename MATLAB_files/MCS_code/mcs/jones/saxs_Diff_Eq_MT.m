function [ RMSDs ] = saxs_Diff_Eq_MT( new_params )

global global_a
global global_b
global global_s 
global global_time_points
global global_t_gt
global global_n_mers_cell
global global_face_centers
global global_crysolIq
global global_expCurve
a = global_a;
b = global_b;
s = global_s;
time_points = global_time_points;
t_gt = global_t_gt;
n_mers_cell = global_n_mers_cell;
face_centers = global_face_centers;
crysolIq = global_crysolIq;
expCurve = global_expCurve;
    
    new_params = new_params'; % since MCS likes to use cols but I use rows.
    param_matrix = new_params;
    % There are 12 params, one for each oligomer size. But since there
    % may be multiple oligomer types per size, some params need to be
    % repeated in the rows of k_matrix.
    k_matrix = formatParams(param_matrix,b);

    RMSDs = zeros(size(k_matrix,1),1);
    clear conc
%     parfor p = 1:size(k_matrix,1)
% % %     for p = 1:size(k_matrix,1)
        p = 1; % Needed with for loop commented out, MCS version.
        
        % ---- Begin Run Simulations and Compute Objective Values ---------        
        
        k = k_matrix % for MCS, k_matrix contains a single param vector.     
    
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
                RMSDs(p) = inf
                display(['RMSD parameter ', num2str(p),' is: ',num2str(RMSDs(p))])
                return; %MT MCS needs no for loop, so 'continue' not appropriate. 
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
% % %     end


end

