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
figure
drawMesh(V,F)
hold on

% Get the representative (centroid) of each face.
face_centers = zeros( numFaces,3 );
for f = 1:numFaces
   verts_curr_face = V( F(f,:),: );
   face_centers(f,:) = mean(verts_curr_face);    
    
end
      

% Plot on same figure the center of each face.
scatter3(face_centers(1,1),face_centers(1,2),face_centers(1,3),'filled','c')
scatter3(face_centers(2,1),face_centers(2,2),face_centers(2,3),'filled','b')
% scatter3(face_centers(5,1),face_centers(5,2),face_centers(5,3),'filled','y')
% scatter3(face_centers(6,1),face_centers(6,2),face_centers(6,3),'filled','k')
scatter3(face_centers(7,1),face_centers(7,2),face_centers(7,3),'filled','g')
% scatter3(face_centers(8,1),face_centers(8,2),face_centers(8,3),'filled','g')
% scatter3(face_centers(9,1),face_centers(9,2),face_centers(9,3),'filled','g')
scatter3(face_centers(2:end,1),face_centers(2:end,2),face_centers(2:end,3))
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off

% -----  TESTING: DISPLAY ROTATED STRUCTURES --------------------
%%{
% Testing other rotatations. Tumble next:
phi = rad2deg(acos(-1/sqrt(5)));
alpha = 360/5;
angle_x1 = 0;%180
angle_x = -(180-phi);%-(180);
angle_y = 0%90-(180-phi);%(180-phi)
angle_z = alpha/2;%5*72/2;

V = (rotx(angle_x1)*V')'
V = (rotx(angle_x)*V')'
V = (roty(angle_y)*V')'
V = (rotz(angle_z)*V')'

% 2D rots
V = (rotz(2*alpha)*V')'


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
% scatter3(face_centers_rot(5,1),face_centers_rot(5,2),face_centers_rot(5,3),'filled','y')
% scatter3(face_centers_rot(6,1),face_centers_rot(6,2),face_centers_rot(6,3),'filled','k')
scatter3(face_centers_rot(7,1),face_centers_rot(7,2),face_centers_rot(7,3),'filled','g')
% scatter3(face_centers_rot(9,1),face_centers_rot(9,2),face_centers_rot(9,3),'filled','g')
% scatter3(face_centers_rot(3,1),face_centers_rot(3,2),face_centers_rot(3,3),'filled','k')
scatter3(face_centers_rot(2:end,1),face_centers_rot(2:end,2),face_centers_rot(2:end,3))
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off
%}