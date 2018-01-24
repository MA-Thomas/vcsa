%{
This script generates an origin centered dodecahedron, computes the 
centroid of each face, and implements an algorithm equivalent to 
Misra & Schwartz to compute all distinct intermediate structures 
of the dodecahedron, as well as degeneracies.

%}

% This is the 2D forward reaction degeneracy matrix.
% 'a' row format: [j,k,m,n,deg]. 
% deg is the degeneracy of the forward reaction btw [j,k] & [m,n].
% SEE Misra and Schwartz 2008 for details.
% Because this file only considers monomer addition, m==j+1. 
a = zeros(1,5); 


% From the geo3d package on MathWorks FileExchange.
[V, E, F] = createDodecahedron;
numFaces = size(F,1);
numNeighbors = size(F,2);

% s(j,k,m,n) is c(j,k) - c(m,n). SEE misra and Schwartz 2008 for details.
% 50 is an overestimate for the number of sister species at a given size.
maxNumSisSpecies = 50;
c = zeros(numFaces,maxNumSisSpecies);

% Rotate around x-axis so that face1 is on top when viewed 
% in the X-Z plane. This will be the standard orientation from which
% all others are derived.
phi_inDegrees = 116.56505; % Dihedral angle
theta_tumble = (180 - phi_inDegrees)/2;
V = (rotx(theta_tumble)*V')';
% % figure
% % drawMesh(V,F)
% % hold on

% Get the representative (centroid) of each face.
face_centers = zeros( numFaces,3 );
for f = 1:numFaces
   verts_curr_face = V( F(f,:),: );
   face_centers(f,:) = mean(verts_curr_face);    
    
end
      

% % % Plot on same figure the center of each face.
% % scatter3(face_centers(1,1),face_centers(1,2),face_centers(1,3),'filled','c')
% % scatter3(face_centers(2,1),face_centers(2,2),face_centers(2,3),'filled','b')
% % scatter3(face_centers(5,1),face_centers(5,2),face_centers(5,3),'filled','y')
% % scatter3(face_centers(6,1),face_centers(6,2),face_centers(6,3),'filled','k')
% % scatter3(face_centers(9,1),face_centers(9,2),face_centers(9,3),'filled','g')
% % scatter3(face_centers(2:end,1),face_centers(2:end,2),face_centers(2:end,3))
% % xlabel('X')
% % ylabel('Y')
% % zlabel('Z')
% % hold off

% -----  TESTING: DISPLAY ROTATED STRUCTURES --------------------
%{
% Testing other rotatations. Tumble next:
phi = rad2deg(acos(-1/sqrt(5)));
alpha = 360/5;
angle_x1 = 180
angle_x = (180-phi);%-(180);
angle_y = 0%90-(180-phi);%(180-phi)
angle_z = -3*alpha/2;%5*72/2;

V = (rotx(angle_x1)*V')'
V = (rotz(angle_z)*V')'
V = (roty(angle_y)*V')'
V = (rotx(angle_x)*V')'

% V = (rotz(angle_z)*V')'
% Get the representative (centroid) of each face.
numFaces = size(F,1);
numNeighbors = size(F,2);
face_centers_rot = zeros( numFaces,3 );
for f = 1:numFaces
   verts_curr_face_rot = V( F(f,:),: );
   face_centers_rot(f,:) = mean(verts_curr_face_rot);    
    
end  

% Plot on same figure the center of each face.
figure
drawMesh(V,F)
hold on
scatter3(face_centers_rot(1,1),face_centers_rot(1,2),face_centers_rot(1,3),'filled','c')
scatter3(face_centers_rot(2,1),face_centers_rot(2,2),face_centers_rot(2,3),'filled','b')
scatter3(face_centers_rot(5,1),face_centers_rot(5,2),face_centers_rot(5,3),'filled','y')
scatter3(face_centers_rot(6,1),face_centers_rot(6,2),face_centers_rot(6,3),'filled','k')
scatter3(face_centers_rot(9,1),face_centers_rot(9,2),face_centers_rot(9,3),'filled','g')
scatter3(face_centers_rot(3,1),face_centers_rot(3,2),face_centers_rot(3,3),'filled','k')
scatter3(face_centers_rot(2:end,1),face_centers_rot(2:end,2),face_centers_rot(2:end,3))
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off
%}
% ------------------ End TESTING: DISPLAY ROTATED STRUCTURES ---------
%%

% Calc 5 nearest neighbors for each face (pentagons have 5 neighbors).
% nearestNeighbor_indices(j,:) will give the 5 faces adjacent to the 
% j-th face.
nearestNeighbor_indices = zeros(numFaces,numNeighbors);
distMat_cartesian = squareform( pdist(face_centers) );
for f = 1:numFaces
    curr_neighbor_dists = sort(distMat_cartesian(f,:));
    curr_bond_dist = curr_neighbor_dists(2); %First elem gives self dist 
    nearestNeighbor_indices(f,:) = ...
        sort(find( abs(distMat_cartesian(f,:)-curr_bond_dist) < 0.001 ));
    
end

% ------ ALGORITHM BELOW ----------
% Iterate through all structures of sizes 1:n.
% Start with a monomer. There are 5 possible dimers that could be formed.
% Form them all.
% Next, check pairs for isomorphism.
% Keep the unique ones dimers (There should be a single unique dimer).
% Next, check all possible trimers, keep the unique ones. Etc.
% ---------------------------------
n_mers_cell = cell(1,numFaces);
unique_face_list = cell(1,numFaces);

% Iterate through all structures of sizes 2:numFaces.
for n = 2:numFaces
    tic;
    display(['n is ',num2str(n)])    
    
    if n == 2
        partners = nearestNeighbor_indices(1,:);
        for i = 1:length(partners)
            connection_matrix = zeros(numFaces,numFaces); 
            connection_matrix(1,partners(i)) = 1;
            connection_matrix(partners(i),1) = 1;
            if i == 1; n_mers_cell{2}{i} = connection_matrix; else
            [isisomorphic,ind] = test_isomorphic(connection_matrix,n_mers_cell{2},...
                face_centers);
            
            if ~isisomorphic % --------------------------------------------
                % Not isomorphic. Add this novel structure to n_mers_cell.
                n_mers_cell{2}{i} = connection_matrix;                
                
                row = [n-1 1 n i];
                diffMat_a = abs(a(:,1:4) - repmat(row, size(a,1), 1));
                if all(sum(diffMat_a,2))
                    % If none of the rows of 'a' match 'row', append to 'a', and set deg=1.
                    a(end+1,:) = [row,1];
                else
                    % Otherwise, increment the approrpiate element of 'a'.
                    row_ind = find(~sum(diffMat_a,2));
                    a(row_ind,5) = a(row_ind,5) + 1;
                end 

                
            else % --------------------------------------------------------
                % Isomorphic. Increment forward reaction degeneracy for 
                % corresponding reaction.
                row = [n-1 1 n ind];
                diffMat_a = abs(a(:,1:4) - repmat(row, size(a,1), 1));               
                if all(sum(diffMat_a,2))
                    % If none of the rows of a match 'row', add it to 'a'.
                    a(end+1,:) = [row,1];
                else
                    % Otherwise, increment the approrpiate element of 'a'.
                    row_ind = find(~sum(diffMat_a,2));
                    a(row_ind,5) = a(row_ind,5) + 1;
                end
                
            end      
            end
        end
        
    else
        
       count = 1;
       
       
       % Look at unique connection matrices for structures of size (n-1).
       for i = 1:length(n_mers_cell{n-1})
           connect_mat = n_mers_cell{n-1}{i};
           faces_in_structure = [];          
           
           % connect_mat is a structure of size (n-1). Identify all the
           % ways its size can be incremented by 1. Create the new
           % connection matrices and store in n_mers_cell{n}.           
           
           % faces (i.e. rows) with at most numNeighbors-1 free edges
           % and at least 1 free edge.
           rows_to_alter = find((0<sum(connect_mat,2))&...
               (5>sum(connect_mat,2)));
           
           if ~isempty(rows_to_alter)
              
              for j = 1:length(rows_to_alter)                  
                 face = rows_to_alter(j);
                 faces_in_structure = [faces_in_structure, ...
                     find(connect_mat(face,:))];
                 
                 
                 % Other faces close enough to be added to 'face'.
                 inds_add = find( ~ismember(...
                 nearestNeighbor_indices(face,:),...
                 find(connect_mat(face,:)) ) );             
                 faces_to_add = nearestNeighbor_indices(face,inds_add);
                 
                 % So as not to create duplicates, only keep entries of
                 % 'faces_to_add' which are not in already in_structure.
                 faces_to_add = faces_to_add( ~ismember(faces_to_add,...
                     faces_in_structure) );              
                 
                 % Add each face (in turn) and, for each unique structure,
                 % save its connection matrix in n_mers_cell{n}.
                 if ~isempty(faces_to_add)
                 for jj = 1:length(faces_to_add)
                     temp = connect_mat;
                     temp(face,faces_to_add(jj)) = 1;
                     temp(faces_to_add(jj),face) = 1;                     
                       
                     list = find(sum(temp,1));
                     b = find(sum(temp,2))';
                     assert( isequal(list,b) )
                     
                     % First, only include structures with n faces. 
                     if length(list)==n                         
                         
                         if isempty(unique_face_list{n})
                             unique_face_list{n} = list;            
                             n_mers_cell{n}{count} = temp;
                             count = count+1;                             
                            
                             row = [n-1,i,n,i];
                             a(end+1,:) = [row,1];

                         else  
                             % If this list appears in the set of lists for
                             % structures of size n, don't add to n_mers_cell{n}
                             pre_z = abs(unique_face_list{n} -...
                                 repmat(list,size(unique_face_list{n},1),1));
                             z = sum(pre_z,2);
                             if any(sum(z,2)==0) && n < numFaces % list exists already and capsid not complete.
                             else

                                 % If 'temp' connection matrix makes it to
                                 % this point, its particular set of faces
                                 % does not appear in another structure of
                                 % the same size. It may still, however, be
                                 % isomorphic to another saved structure.
                                 [isisomorphic,ind] = test_isomorphic(temp,...
                                     n_mers_cell{n},face_centers);
% %                                 if length(ind)>1;T = temp; Cell = n_mers_cell{n};return;end
                                 if length(ind)>1;return;end                                 
                                 if ~isisomorphic %------------------------
                                    % Not isomorphic. Add structure to n_mers_cell{n}.
                                    unique_face_list{n} = [unique_face_list{n};list];            
                                    n_mers_cell{n}{count} = temp;                                                                     
                                    
                                    row = [n-1,i,n,count]; 
                                    diffMat_a = abs(a(:,1:4) - repmat(row, size(a,1), 1));
                                    if all(sum(diffMat_a,2))
                                        % If none of the rows of 'a' match 'row', append to 'a', and set deg=1.
                                        a(end+1,:) = [row,1];
                                    else
                                        % Otherwise, increment the approrpiate element of 'a'.
                                        row_ind = find(~sum(diffMat_a,2));
                                        a(row_ind,5) = a(row_ind,5) + 1;
                                    end
                                    count = count+1;
                                    
                                 else % -----------------------------------
                                    % Isomorphic. Increment forward reaction degeneracy for 
                                    % corresponding reaction.
                                    row = [n-1,i,n,ind];
                                    diffMat_a = abs(a(:,1:4) - repmat(row, size(a,1), 1));
                                    if all(sum(diffMat_a,2))
                                        % If none of the rows of 'a' match 'row', append to 'a', and set deg=1.
                                        a(end+1,:) = [row,1];
                                    else
                                        % Otherwise, increment the approrpiate element of 'a'.
                                        row_ind = find(~sum(diffMat_a,2));
                                        a(row_ind,5) = a(row_ind,5) + 1;
                                    end                                
                                                                        
                                 end

                                
                                 
                             end
                             
                         end
                     end
                         
                 end
                 end
                 
              end
              
           else
               % DONE!
               display('no free edges in structure.')
               
           end           
           
       end       
       
        % Feed pairs of connection graphs to isomorphism algorithm.
        % Remove duplicate (i.e. equivalent) structures identified by 
        % the algorithm.
        
        
    end
    


    
    % Finally, compute c(m,n) using the connection matrices.
    % To do this, we need to first alter the connection matrix for [m,n] to
    % include all bonds, not merely the subset of bonds originally included to 
    % add each subsequent face into the structure.
    for mm = 2:numFaces
        len = length(n_mers_cell{mm});
        
        for nn = 1:maxNumSisSpecies           
            if nn <= len
               connect_mat = n_mers_cell{mm}{nn};
               
               % Alter connect_mat to reflect all bonds present.
               faces_present = find(sum(connect_mat,1));
               for fp = 1:length(faces_present)
                  partners = nearestNeighbor_indices( faces_present(fp),: );
                  partners_present = partners( ismember(partners, faces_present) );
                  connect_mat(faces_present(fp), partners_present) = 1;
               end
               
               % The connection matrix will tell us the number of bonds.
               % Since the matrix is a symmetric matrix, it double counts,
               % and thus, we can simply count the number of '1's on the
               % main diagonal and above.
               c(mm,nn) = sum(sum(triu(connect_mat)));
               
            end
        end
    end
    toc
end
    


% Get rid of first row of 'a' used to initialize it.
a = a(2:end,:);   

% ---------------------------------------------------
% ---------------------------------------------------
% % Determine O(m,n), the symmetry of oligomer (m,n).
a_ = a;
a_(:,2) = 0;
a_(:,end) = [];
[C,ia,ic] = unique(a_,'rows');
 
% Order of elements in rows of C is [j,-,m,n].
% accumarray(ic,1) is a col vector. Each entry tells how many times
% [j,-,m,n] appears in a_.

% Order of elements in rows of O is [m,n,O(m,n)].
% O = [C(:,[3,4]) accumarray(ic,1)];

% % Determine the backward degeneracies.
%   It turns out that 'a' ~ [j,k,m,n,deg] contains the 
%   backward degeneracy info also. 
b = [a(:,[3,4]) a(:,[1,2]) a(:,5)];



