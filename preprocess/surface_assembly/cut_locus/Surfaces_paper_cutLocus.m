% This script produces the surfaces we need for the paper on OT and Cut
% Locus
clear;
close all;
clc;
addpath('../distmesh');
% % NOT IN THE REFERENCE - NOT TO BE USED - GENERAL CODE TO BE COPIED AND
% ADAPTED TO THE DESIRED SURFACE
% % Uniform mesh on an asymmetrical ellipsoid x^2/a^2+y^2/b^2+z^2/c^2=1
% % Umbilic points are (\pm a \sqrt(\frac{a^2-b^2}{a^2-c^2}),0,\pm c
% % \sqrt{\frac{b^2-c^2}{a^2-c^2}}
% % a^2=1
% % b^2=2
% % c^2=4
% fd=@(p) p(:,1).^2+p(:,2).^2/2+p(:,3).^2/4-1;
% [p,t]=distmeshsurface(fd,@huniform,0.2,[-1,-1.5,-2; 1,1.5,2]);
% % Fixed points we want: we substitute them to the closest point on the
% % ellipsoid mesh, using the Euclidean distance function in R^3
% % Here we input the umbilic point and its antipodal
% fixed_points = [1/sqrt(3) 0 2*sqrt(2/3); -1/sqrt(3) 0 -2*sqrt(2/3);
%     1 0 0; -1 0 0];
% p=substitute_point_in_mesh(p,fixed_points);
% fid=fopen('mesh_asymmetric_ellipsoid.dat','w');
% fprintf(fid,'%g\t',size(p,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(t,1));
% fprintf(fid,'\n');
% for ii = 1:size(p,1)
%     fprintf(fid,'%g\t',p(ii,:));
%     fprintf(fid,'\n');
% end
% t = [t ones(size(t,1),1)];
% for ii = 1:size(t,1)
%     fprintf(fid,'%g\t',t(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid)

% Uniform mesh on an asymmetrical ellipsoid x^2/a^2+y^2/b^2+z^2/c^2=1
% Umbilic points are (\pm a \sqrt(\frac{a^2-b^2}{a^2-c^2}),0,\pm c
% \sqrt{\frac{b^2-c^2}{a^2-c^2}}
% a^2=0.2^2
% b^2=0.6^2
% c^2=1
% REFERENCE: "Thaw: A Tool for 
% Approximating Cut Loci on a Triangulation of a Surface
% (4-1) WITH UMBILIC POINT AND RANDOM POINT
fd=@(p) p(:,1).^2/0.2^2+p(:,2).^2/0.6^2+p(:,3).^2-1;
[p,t]=distmeshsurface(fd,@huniform,0.1,[-0.2,-0.6,-1; 0.2,0.6,1]);
% Fixed points we want: we substitute them to the closest point on the
% ellipsoid mesh, using the Euclidean distance function in R^3
% Here we input the umbilic point and its antipodal
fixed_points = [-0.115470 0 0.816497; 0.115470  0 -0.816497; % Umbilic point below (4-1)
    -0.151128  -0.350718 0.295520; 0.151128  0.350718 -0.295520]; %Random point below (5-1)
p=substitute_point_in_mesh(p,fixed_points);
fid=fopen('ellipsoid_grid.dat','w');
fprintf(fid,'%g\t',size(p,1));
fprintf(fid,'\n');
fprintf(fid,'%g\t',size(t,1));
fprintf(fid,'\n');
for ii = 1:size(p,1)
    fprintf(fid,'%g\t',p(ii,:));
    fprintf(fid,'\n');
end
t = [t ones(size(t,1),1)];
for ii = 1:size(t,1)
    fprintf(fid,'%g\t',t(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);

% Uniform mesh on a torus
% fd=@(p) (sqrt(p(:,1).^2+p(:,2).^2)-2).^2+p(:,3).^2-1;
% [p,t]=distmeshsurface(fd,@huniform,0.2,[-3.5,-3.5,-1; 3.5,3.5,1]);
% fixed_points=[-3 0 0; 3 0 0; 0 3 0; 0 -3 0; 1 0 0; -1 0 0; 0 1 0; 0 -1 0];
% p=substitute_point_in_mesh(p,fixed_points);
% fid=fopen('torus_grid.dat','w');
% fprintf(fid,'%g\t',size(p,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(t,1));
% fprintf(fid,'\n');
% for ii = 1:size(p,1)
%     fprintf(fid,'%g\t',p(ii,:));
%     fprintf(fid,'\n');
% end
% t = [t ones(size(t,1),1)];
% for ii = 1:size(t,1)
%     fprintf(fid,'%g\t',t(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid);

% % Uniform mesh on a unit sphere
% fd=@(p) dsphere(p,0,0,0,1);
% [p,t]=distmeshsurface(fd,@huniform,0.2,1.1*[-1,-1,-1;1,1,1]);
% fixed_points=[0 0 1; 0 0 -1];
% p=substitute_point_in_mesh(p,fixed_points);
% fid=fopen('sphere_grid.dat','w');
% fprintf(fid,'%g\t',size(p,1));
% fprintf(fid,'\n');
% fprintf(fid,'%g\t',size(t,1));
% fprintf(fid,'\n');
% for ii = 1:size(p,1)
%     fprintf(fid,'%g\t',p(ii,:));
%     fprintf(fid,'\n');
% end
% t = [t ones(size(t,1),1)];
% for ii = 1:size(t,1)
%     fprintf(fid,'%g\t',t(ii,:));
%     fprintf(fid,'\n');
% end
% fclose(fid);

%Uniform mesh on a quartic surface (figure 18 paper "Thaw: A Tool for 
%Approximating Cut Loci on a Triangulation of a Surface")
fd=@(p) p(:,1).^4+p(:,2).^4+p(:,3).^4-1;
[p,t]=distmeshsurface(fd,@huniform,0.2,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
fixed_points=[0.533843 0.800764 0.844080];
p=substitute_point_in_mesh(p,fixed_points);
fid=fopen('quartic_grid.dat','w');
fprintf(fid,'%g\t',size(p,1));
fprintf(fid,'\n');
fprintf(fid,'%g\t',size(t,1));
fprintf(fid,'\n');
for ii = 1:size(p,1)
    fprintf(fid,'%g\t',p(ii,:));
    fprintf(fid,'\n');
end
t = [t ones(size(t,1),1)];
for ii = 1:size(t,1)
    fprintf(fid,'%g\t',t(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);