function [isisomorphic,indices] = test_isomorphic(C, C_cell,...
    subunits_cart)
% Outputs and Inputs:
% isisomorphic is a logical. 1 means the input, connection matrix C, is
% isomorphic to some connection matrix within input cellarray C_cell.
% indices identifies the location(s) within C_cell of the isomorphic structure.
% If no matrix is isomorphic, index==-1.

% DODECAHEDRON

% TERMINOLOGY: 'tumbling orientation' - a rotation resulting in a new
% face in the fixed location.
% 'rotation/orientation' - a rotation about axis perpendicular to the fixed face.
% Because I oriend the dodecahedron so that the fixed face is parallel to
% the x-y plane, at the highest z location, 'rotations' are w.r.t. the
% coordinate axes, {x,y,z}.

indices = -1;
isisomorphic = 0;

% --------------------------
% How this function works
% --------------------------
% Compute all orientations of C.
% Logic: with the location of a face fixed, there are 5 orientations it can
% take (achieved by successive rotations of 360/5 degrees around an axis
% perpendicular to the face), and therefore 5 orientations for the
% structrure.

% Each of the 12 faces has 5 orientations (rotations in plane).

% So, start by generating the initial 5 orientations (the original, plus 4
% reached by rotation). If none of these are found in C_cell, move on.
% Next, tumble the structure so that a new face becomes the fixed face.
% (The fixed face is always located at the same place)
% Explore its 5 orientations. 
% After tumbling through all 12 faces, if none of the structures are found
% in C_cell, exit.

% Precision. Don't set this too large or numerical errors may prevent 
% sortrows from sorting consistently. Even 5 sign figs is too many.
decimals = 4; 
subunits_cart = round( subunits_cart, decimals);

subunits_in_C = find(sum(C,1));
C_coords = round( subunits_cart(subunits_in_C,:), decimals);

% Replace the pesky -0.0000 entries with 0.
C_coords( abs(C_coords - zeros(size(C_coords)))<0.0005 ) = 0;

location_of_fixed_face = round( subunits_cart(1,:), decimals);


% When creating row vector of coordinates [x...,y...,z...] for a given
% structure, we need a consistent ordering scheme.
C_coords_vector = reshape(sortrows(C_coords),1,numel(C_coords));


% Construct a matrix in which each row contains all of the coordinates
% [x...,y...,z...] of a particular structure in C_cell. 
% WE CAN THEN IDENTIFY EQUIVALENT STRUCTURES AS HAVING IDENTICAL ROWS.
C_cell_coords_matrix = zeros( length(C_cell),numel(C_coords) );
for j = 1:length(C_cell)
    subunits = find(sum(C_cell{j},1));
    coords = round( subunits_cart(subunits,:), decimals);
    
    % Replace any pesky -0.0000 entries with 0.
    coords( abs(coords - zeros(size(coords)))<0.005 ) = 0;

    C_cell_coords_matrix(j,:) = reshape(sortrows(coords),1,numel(coords));    
    
end

C_cell_coords_matrix = round(C_cell_coords_matrix, decimals);

% ----------------------------------------------------------------
% Explore initial 4 orientations of the 1st 'tumbling orientation'
% (rotations in the plane of the fixed face)
% ----------------------------------------------------------------
initial_rotd_coords = zeros(5,numel(C_coords));
initial_rotd_coords(1,:) = C_coords_vector;
alpha = 360/5;  

% Angle btw successive orientations in plane.
twoDim_rots = [0,alpha,-alpha,2*alpha,-2*alpha];
for i = 2:length(twoDim_rots)    
    theta = twoDim_rots(i);
    Crot = round( (rotz(theta)*C_coords')', decimals);
    
    % Replace any pesky -0.0000 entries with 0.
    Crot( abs(Crot - zeros(size(Crot)))<0.0005 ) = 0;   
    initial_rotd_coords(i,:) = reshape(sortrows(Crot),1,numel(Crot));

end



% If any row in initial_rotd_coords is found in C_cell_coords_matrix, 
% set isisomorphic = 1.
% *** Because there still may be differences related to numerical
%     precision, round initial_rotd_coords and C_cell_coords_matrix to 2
%     sigfigs before comparison. 
% ***
initial_rotd_coords = round( initial_rotd_coords,2 );
C_cell_coords_matrix = round( C_cell_coords_matrix,2 );
[CC,ia,ib] = intersect( initial_rotd_coords, C_cell_coords_matrix, 'rows' );

% % % for tt = 1:size(initial_rotd_coords,1)
% % %    r = initial_rotd_coords(tt,:);
% % %    r_repmat = repmat(r,size(C_cell_coords_matrix,1),1);
% % %    diffm = abs(r_repmat - C_cell_coords_matrix);
% % %    diffs = sum(diffm,2);
% % %    ind_below = find(diffs<1);
% % %    if ~isempty(ind_below)
% % %       r' 
% % %       C_cell_coords_matrix(ind_below,:)'
% % %    end
% % % end

if ~isempty(CC)
    indices = ib;
    isisomorphic = 1; 
%     display('is isomorphic, pre tumble')
    return
end








% -----------------------------------------------------------------
% Tumble the structure so that a new face is in the fixed location
% Explore this structure's 5 orientations. Etc.
% -----------------------------------------------------------------

% tumble_list is the list of rotation-tuples used to transform the
% 'tumbling orientation' of the implicit dodecahedron into all of its
% (12-1) remaining tumbling orientations.
% Each orientation is defined by which face is on top in the 
% 'fixed location'.
% Each row of tumble_list contains [angle_x, angle_y, angle_z], the angles
% about the x,y,z axes (in radians) needed to transform the original 
% tumbling orientation into each of the (12-1) others. 


phi = rad2deg(acos(-1/sqrt(5))); % Dihedral angle in degrees.
alpha = 360/5;          % 5 rotatns by alpha returns pentagon to original

% In the fifth col of tubmle_list:
% '1' means use order x-x-y-z.
% '2' means use order x-z-y-x.
tumble_list = ...
    [ % Bottom pentagon. Rotation Order is x-x-y-z
      180, 0, 0, 0, 1;...
      
      % Place each element of upper ring of 5 pentagons in fixed location.       
      0, -(180-phi), 0,   alpha/2,    1;...% Rotation Order is x-x-y-z.    
      0, (180-phi),  0,   alpha/2,    2;...% Rotation Order is x-z-y-x.
      0, (180-phi),  0,   -alpha/2,   2;...% Rotation Order is x-z-y-x.
      0, (180-phi),  0,   3*alpha/2,  2;...% Rotation Order is x-z-y-x.
      0, (180-phi),  0,   -3*alpha/2, 2;...% Rotation Order is x-z-y-x.
      
      % Lower ring of 5 pentagons. Flip the dodecah upside-down first.
      % Then apply same operations as upper ring.
      180, -(180-phi), 0,   alpha/2,    1;...% Rotation Order is x-x-y-z.    
      180, (180-phi),  0,   alpha/2,    2;...% Rotation Order is x-z-y-x.
      180, (180-phi),  0,   -alpha/2,   2;...% Rotation Order is x-z-y-x.
      180, (180-phi),  0,   3*alpha/2,  2;...% Rotation Order is x-z-y-x.
      180, (180-phi),  0,   -3*alpha/2, 2];% Rotation Order is x-z-y-x.
      
for t = 1:size(tumble_list,1)
    
    angle_x_pre = tumble_list(t,1);
    angle_x = tumble_list(t,2);
    angle_y = tumble_list(t,3);
    angle_z = tumble_list(t,4);
    
    if tumble_list(t,5) == 1
        tumbled_coords = (rotx(angle_x_pre)*C_coords')';
        tumbled_coords = (rotx(angle_x)*tumbled_coords')';
        tumbled_coords = (roty(angle_y)*tumbled_coords')';
        tumbled_coords = round( (rotz(angle_z)*tumbled_coords')', decimals);
    else
        tumbled_coords = (rotx(angle_x_pre)*C_coords')';
        tumbled_coords = (rotz(angle_z)*tumbled_coords')';
        tumbled_coords = (roty(angle_y)*tumbled_coords')';
        tumbled_coords = round( (rotx(angle_x)*tumbled_coords')', decimals);
    end
    
    
%     % If new tumbling orientation does not result in a face being present
%     % in the 'fixed location' (with less than a complete dodecahedron
%     % structure, this will happen), CONTINUE.
%     diffMatrix = abs(tumbled_coords - ...
%         repmat(location_of_fixed_face, size(tumbled_coords,1), 1));
%     if ~any( abs(sum(diffMatrix,2)) < 0.001 )
%         continue
%     end
    
    
    
    % Explore the 5 rotations of this 'tumbling orientation'.
    rotd_coords = zeros(5,numel(tumbled_coords));
    rotd_coords(1,:) = reshape(sortrows(tumbled_coords),1,numel(tumbled_coords));
    twoDim_rots = [0,alpha,-alpha,2*alpha,-2*alpha]; 
    for i = 2:length(twoDim_rots)
        theta = twoDim_rots(i);
        Crot = round( (rotz(theta)*tumbled_coords')', decimals);
        rotd_coords(i,:) = reshape(sortrows(Crot),1,numel(Crot));

    end
    
    % Replace any pesky -0.0000 entries with 0.
    rotd_coords(...
        abs(rotd_coords - zeros(size(rotd_coords)))<0.005 )...
        = 0;
   
    % If any row in rotd_coords is found in C_cell_coords_matrix, 
    % set isisomorphic = 1. 
    % *** Because there still may be differences related to numerical
    %     precision, round rotd_coords and C_cell_coords_matrix to 2
    %     sigfigs before comparison. 
    % ***
    rotd_coords = round( rotd_coords,2 );
    C_cell_coords_matrix = round( C_cell_coords_matrix,2 );    
    [CC,ia,ib] = intersect( rotd_coords, C_cell_coords_matrix, 'rows' );

% % %     for tt = 1:size(rotd_coords,1)
% % %        r = rotd_coords(tt,:);
% % %        r_repmat = repmat(r,size(C_cell_coords_matrix,1),1);
% % %        diffm = abs(r_repmat - C_cell_coords_matrix);
% % %        diffs = sum(diffm,2);
% % %        ind_below = find(diffs<1);
% % %        if ~isempty(ind_below)
% % %           r' 
% % %           C_cell_coords_matrix(ind_below,:)'
% % %        end
% % %     end
    
    
    if ~isempty(CC)
        indices = ib;
        isisomorphic = 1;
%         display('is isomorphic, tumble')
        return
    end


end



end


% Old code
%{
function [mean_azimuth1,mean_elevation1,mean_azimuth2,mean_elevation2] = ...
    calc_mean_angles(subunits_sph,subunits_in_C, subunits_in_C_other)

   numer1 = sum( sin( subunits_sph(subunits_in_C,2) ).*sin( subunits_sph(subunits_in_C,1) ) );
   denom1 = sum( sin( subunits_sph(subunits_in_C,2) ).*cos( subunits_sph(subunits_in_C,1) ) );
   numer2 = sum( sin( subunits_sph(subunits_in_C_other,2) ).*sin( subunits_sph(subunits_in_C_other,1) ) );
   denom2 = sum( sin( subunits_sph(subunits_in_C_other,2) ).*cos( subunits_sph(subunits_in_C_other,1) ) );   
   mean_azimuth1 = atan( numer1 / denom1 );
   mean_azimuth2 = atan( numer2 / denom2 );
   
   N = length(subunits_in_C);
   mean_elevation1 = acos( sum( cos( subunits_sph(subunits_in_C,2) ) ) / N );
   mean_elevation2 = acos( sum( cos( subunits_sph(subunits_in_C_other,2) ) ) / N );


end
%}

%{
   % Check that mean(r) is the same for C and C_other. step2 in Misra Algo.
   mean_r1 = mean( subunits_sph(subunits_in_C,3) );
   mean_r2 = mean( subunits_sph(subunits_in_C_other,3) );
   if ~isequal(mean_r1,mean_r2)
       isisomorphic = 0;
       index = -1;
       continue
   end
   
   % Calculate mean phi and theta for each graph. Step1 in Misra Algorithm.  
    [mean_azimuth1,mean_elevation1,mean_azimuth2,mean_elevation2] = ...
        calc_mean_angles(subunits_sph,subunits_in_C, subunits_in_C_other)   
   
   if isequal(mean_azimuth1,mean_azimuth2) &...
           isequal(mean_elevation1,mean_elevation1)
       isisomorphic = 0;
       index = -1;
       continue       
   end
   
   % Use rodrigues formula to rotate the vectors defining the subunits 
   % of C_other so that centroid(C_other') is the same as centroid(C).
   % See wiki on Rodrigues Rotation Formula
   [x1,y1,z1] = sph2cart(mean_azimuth1,mean_elevation1,mean_r1);
   [x2,y2,z2] = sph2cart(mean_azimuth2,mean_elevation2,mean_r2);
   a = [x1,y1,z1]
   b = [x2,y2,z2]
   angle_btw_centroid_vectors = atan2(norm(cross(a,b)), dot(a,b))
   k = cross(a,b) / ( norm(cross(a,b)) )
   
   v = subunits_cart(subunits_in_C_other,:)
   v_rot = zeros(size(v));
   for j = 1:size(v,1)
   v_rot(j,:) = v(j,:)*cos(angle_btw_centroid_vectors) + ...
           cross(k,v(j,:))*sin(angle_btw_centroid_vectors) + ...
           k*dot(k,v(j,:))*(1-cos(angle_btw_centroid_vectors))
   end

   
   
   % Check that mean angles are same for both structures
   subunits_cart(subunits_in_C_other,:) = v_rot
   [x,y,z] = cart2sph(subunits_cart(:,1),subunits_cart(:,2),subunits_cart(:,3));
   subunits_sph = [x,y,z]
   [mean_azimuth1,mean_elevation1,mean_azimuth2,mean_elevation2] = ...
    calc_mean_angles(subunits_sph,subunits_in_C, subunits_in_C_other)   

   
   % Step3 in Misra Algo.
   R_z_azimuth1 = [cos(mean_elevation1) -sin(mean_elevation1) 0;...
                   sin(mean_elevation1) cos(mean_elevation1) 0;...
                   0 0 1]   
   rot1_toNorthPole = R_z_azimuth1;
   
%-------------------------------------------------   
   % Perform Euler rotations on each such that each centroid coincides with
   % azimuth==0 and elevation==0.
   subunits_cart_new = (rot1_toNorthPole*subunits_cart')'
   [az_new,elev_new,r_new] = cart2sph(subunits_cart_new(:,1),...
       subunits_cart_new(:,2),subunits_cart_new(:,3));
   subunits_sph_new = [az_new,elev_new,r_new]
   % Calculate mean phi and theta for each graph. Step1 in Misra Algorithm.  
   N = length(subunits_in_C);
   mean_elevation1 = acos( sum( cos( subunits_sph_new(subunits_in_C,2) ) ) / N )
   
   numer1 = sum( sin( subunits_sph_new(subunits_in_C,2) ).*sin( subunits_sph_new(subunits_in_C,1) ) );
   denom1 = sum( sin( subunits_sph_new(subunits_in_C,2) ).*cos( subunits_sph_new(subunits_in_C,1) ) );
   mean_azimuth1 = atan( numer1 / (denom1+10^-15) ) 
   
   % -------------
   
   isisomorphic = 1
   
   %}