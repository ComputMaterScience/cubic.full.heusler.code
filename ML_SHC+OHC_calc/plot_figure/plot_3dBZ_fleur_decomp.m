clear all;clc;figure; hold on;

% code for BZ boundary
%view-source:http://lampx.tugraz.at/~hadley/ss1/bzones/drawing_BZ.php

% load and input data

prefix = 'Al\pbe';

kdiv = 25;

% n_type = 1 <px|Lz|py>
% n_type = 2 <py|Lz|px>
% n_type = 3 <dxz|Lz|dyz>
% n_type = 4 <dyz|Lz|dxz>
% n_type = 5 <dx2-y2|Lz|dxy>
% n_type = 6 <dxy|Lz|dx2-y2>

n_type = 1:6;

zoom_factor = 1.7;

echg = 1.602176634e-19;

hbar = 1.054571817e-34;

scale_factor = (echg^2/hbar)*1.0e8;

winfile = fopen(strcat('thi\',prefix,'\data\wannier90.win'),'rt');
lat_vec = zeros(3,3);
while ~feof(winfile)
    line = fgetl(winfile);
    if contains(line,'begin unit_cell_cart')
        line = fgetl(winfile);
        for i = 1:3
            line = fgetl(winfile);
            lat_vec(i,:) = sscanf(line, '%f %f %f');
        end
    end
    if ~isempty(regexp(line, '\<fermi_energy\>', 'once'))
        efermi = split(line); efermi = str2double(efermi{3});
    end
end

fclose(winfile);

inv_lat = 2*pi*inv(lat_vec)';

% lat_vec = [pi,pi,0;pi,0,pi;0,pi,pi];

% define variables

G = zeros(1,4); % the reciprocal lattice vectors, the fourth index checks
% if it is in the 1 Bz

G1 = zeros(1,3); % the reciprocal lattice vectors used to draw the 1st Bz

% hkl plane
hkl = [[0,0,1];[0,0,-1];[0,1,0];[0,-1,0];[0,1,1];[0,-1,-1];[1,0,0];...
    [-1,0,0];[1,0,1];[-1,0,-1];[1,1,0];[-1,-1,0];[1,1,1];[-1,-1,-1];...
    [1,1,-1];[-1,-1,1];[1,-1,1];[-1,1,-1];[-1,1,1];[1,-1,-1];[-1,1,0];...
    [1,-1,0];[1,0,-1];[-1,0,1];[0,1,-1];[0,-1,1]];

corners = zeros(1,4); %the first elements are the coordinates,
%the next three are the planes that intersect
%at that point, the last is the Bz index.

cornersbz = zeros(1,3); % the corners of the Bz zone

edges = zeros(1,6); % the edges of the Bz zone

% lattice vectors
a1x = lat_vec(1,1); a1y = lat_vec(1,2); a1z = lat_vec(1,3);
a2x = lat_vec(2,1); a2y = lat_vec(2,2); a2z = lat_vec(2,3);
a3x = lat_vec(3,1); a3y = lat_vec(3,2); a3z = lat_vec(3,3);

Np = 0; % number of planes

Nc = 0; % number of corners

Ne = 0; % number of edges

% volume
v = a1x*(a2y*a3z-a2z*a3y)+a1y*(a2z*a3x-a2x*a3z)+a1z*(a2x*a3y-a2y*a3x);

% reciprocal lattice vectors
b1x = 2*pi*(a2y*a3z-a2z*a3y)/v;
b1y = 2*pi*(a2z*a3x-a2x*a3z)/v;
b1z = 2*pi*(a2x*a3y-a2y*a3x)/v;
b2x = 2*pi*(a3y*a1z-a3z*a1y)/v;
b2y = 2*pi*(a3z*a1x-a3x*a1z)/v;
b2z = 2*pi*(a3x*a1y-a3y*a1x)/v;
b3x = 2*pi*(a1y*a2z-a1z*a2y)/v;
b3y = 2*pi*(a1z*a2x-a1x*a2z)/v;
b3z = 2*pi*(a1x*a2y-a1y*a2x)/v;

% define the reciprocal lattice vectors
for i = 1:size(hkl,1)
    G(i,:) = [hkl(i,1)*b1x+hkl(i,2)*b2x+hkl(i,3)*b3x,...
        hkl(i,1)*b1y+hkl(i,2)*b2y+hkl(i,3)*b3y,...
        hkl(i,1)*b1z+hkl(i,2)*b2z+hkl(i,3)*b3z,1];
end

% find the planes that form the boundaries of the first Bz
for i = 1:size(hkl,1)
    Dgamma = sqrt((0.5*G(i,1))^2 + (0.5*G(i,2))^2 + (0.5*G(i,3))^2);
    for j = 1:size(hkl,1)/2
        if (j~=i)
            GG = sqrt((G(i,1)/2-G(j,1))^2+(G(i,2)/2-G(j,2))^2+...
                (G(i,3)/2-G(j,3))^2);
            if (GG <= Dgamma)
                G(i,4) = 0; % this G is not part of the Bz boundary
                break
            end
        end
    end
end

fprintf('The first Brillouin zone boundary consists:\n');

fprintf('Planes:\n');
G_hkl = zeros(1,3);
for i = 1:size(hkl,1)
    if (G(i,4) == 1)
        Np = Np + 1;
        G1(Np,:) = [G(i,1),G(i,2),G(i,3)];
        fprintf('(%d,%d,%d)\n',hkl(i,1),hkl(i,2),hkl(i,3));
        G_hkl(Np,:) = [hkl(i,1),hkl(i,2),hkl(i,3)];
    end
end
fprintf('Number of planes: %d\n',Np);

% This section considers all distinct combinations of three planes
% and determines the corners where they intersect
n = 0;
for i = 3:Np
    for j = 2:i-1
        for k = 1:j-1
            A = [[G1(i,1), G1(i,2), G1(i,3)];...
                [G1(j,1), G1(j,2), G1(j,3)];...
                [G1(k,1), G1(k,2), G1(k,3)]];
            b(1) = (G1(i,1)*G1(i,1)+G1(i,2)*G1(i,2)+G1(i,3)*G1(i,3))/2;
            b(2) = (G1(j,1)*G1(j,1)+G1(j,2)*G1(j,2)+G1(j,3)*G1(j,3))/2;
            b(3) = (G1(k,1)*G1(k,1)+G1(k,2)*G1(k,2)+G1(k,3)*G1(k,3))/2;
            if (det(A) ~= 0)
                n = n+1;
                corners(n,:) = [linsolve(A,b')',1];
            end
        end
    end
end

% find the corners of the Bz by choosing only the corners closer to Gamma
% than another G
for i = 1:size(corners,1)
    for j = 1:size(hkl,1)
        Dgamma = sqrt(corners(i,1)^2+corners(i,2)^2+corners(i,3)^2);
        if (j~=i)
            DGG = sqrt((corners(i,1)-G(j,1))^2+...
                (corners(i,2)-G(j,2))^2+(corners(i,3)-G(j,3))^2);
            if DGG < Dgamma
                corners(i,4) = 0;
                break
            end
        end
    end
end

% check to see if the corners are unique
for i = 2:size(corners,1)
    for j = 1:i-1
        if (corners(i,1) == corners(j,1)) &&...
                (corners(i,2) == corners(j,2)) &&...
                (corners(i,3) == corners(j,3))
            corners(i,4) = 0;
        end
    end
end

fprintf('Corners:\n');
dmax = 0;
for i = 1:size(corners,1)
    if corners(i,4) == 1
        d = sqrt(corners(i,1)^2+corners(i,2)^2+corners(i,3)^2);
        dmax = max(dmax,d);
        Nc = Nc + 1;
        cornersbz(Nc,:) = [corners(i,1), corners(i,2), corners(i,3)];
        fprintf('(%.2f,%.2f,%.2f) d=%f\n',cornersbz(Nc,1),...
            cornersbz(Nc,2),cornersbz(Nc,3), dmax);
    end
end
fprintf('Number of corners: %d\n',Nc);

% for every pair of corners check every pair of planes to see
% if the corners both lie in those planes
delta = 1E-20;
for i = 1:Nc
    for j = 1:i-1
        for k = 2:Np
            for l = 1:k-1
                d_ki = G1(k,1)*cornersbz(i,1)+G1(k,2)*cornersbz(i,2)+...
                    G1(k,3)*cornersbz(i,3)-(G1(k,1)*G1(k,1)+G1(k,2)*G1(k,2)+...
                    G1(k,3)*G1(k,3))/2;
                d_li = G1(l,1)*cornersbz(i,1)+G1(l,2)*cornersbz(i,2)+...
                    G1(l,3)*cornersbz(i,3)-(G1(l,1)*G1(l,1)+G1(l,2)*G1(l,2)+...
                    G1(l,3)*G1(l,3))/2;
                d_kj = G1(k,1)*cornersbz(j,1)+G1(k,2)*cornersbz(j,2)+...
                    G1(k,3)*cornersbz(j,3)-(G1(k,1)*G1(k,1)+G1(k,2)*G1(k,2)+...
                    G1(k,3)*G1(k,3))/2;
                d_lj = G1(l,1)*cornersbz(j,1)+G1(l,2)*cornersbz(j,2)+...
                    G1(l,3)*cornersbz(j,3)-(G1(l,1)*G1(l,1)+G1(l,2)*G1(l,2)+...
                    G1(l,3)*G1(l,3))/2;
                if (abs(d_ki)<delta) && (abs(d_li)<delta)...
                        && (abs(d_kj)<delta)&& (abs(d_lj)<delta)
                    Ne = Ne+1;
                    edges(Ne,:) = [cornersbz(i,1),cornersbz(i,2),...
                        cornersbz(i,3),cornersbz(j,1),cornersbz(j,2),...
                        cornersbz(j,3)];
                    break
                end
            end
        end
    end
end

fprintf('Edges:\n');
for i = 1:size(edges,1)
    fprintf('(%.2f,%.2f,%.2f) to (%.2f,%.2f,%.2f)\n',edges(i,1),...
        edges(i,2),edges(i,3),edges(i,4),edges(i,5),edges(i,6))
end
fprintf('Number of edges: %d\n',size(edges,1));

for i = 1:size(edges,1)
    plot3([edges(i,1),edges(i,4)],[edges(i,2),edges(i,5)],...
        [edges(i,3),edges(i,6)],'k-','LineWidth',3.0);
end

% for patch of BZ faces
% k = boundary(cornersbz(:,1),cornersbz(:,2),cornersbz(:,3),0.001);
%
% bz_patch = trisurf(k,cornersbz(:,1),cornersbz(:,2),cornersbz(:,3),...
%     'Facecolor','w','FaceAlpha',0.05);
%
% bz_patch.LineStyle = 'none';


% import data
if isfile(strcat('thi\',prefix,'\berry\fermi_data.mat'))
    load(strcat('thi\', prefix,'\berry\fermi_data.mat'));
else
    fid = fopen(strcat('thi\',prefix,'\berry\fermi_surface_ohc_decomp.dat'),'r');
    % get number of kgrid
    line = fgetl(fid); line = split(line,':');
    line = split(line{2});
    ngrid = str2double(line(2:end));
    nk = ngrid(1)*ngrid(2)*ngrid(3);
    % get number of bands
    line = fgetl(fid); line = split(line,':');
    nband = str2double(line{2});
    % get e-fermi
    line = fgetl(fid); line = split(line,':');
    efermi = str2double(line{2});
    % get orbital projection
    line = fgetl(fid); line = split(line,':');
    proj = split(line{2}); proj = proj(2:end-1);
    % get eigenvalues
    eig = zeros(nk,nband);
    fgetl(fid);
    for i = 1:nband
        fgetl(fid);
        for k = 1:nk
            eig(k,i) = str2double(fgetl(fid));
        end
    end
    
    % get conductivity
    if size(proj,2) == 1
        num_proj = 2;
    else
        num_proj = 6;
    end
    conduct = zeros(nk,nband,num_proj);
    fgetl(fid);
    for i = 1:nband
        fgetl(fid);
        for j = 1:num_proj
            fgetl(fid);
            for k = 1:nk
                conduct(k,i,j) = str2double(fgetl(fid));
            end
        end
    end
    
    % get berry
    berry = zeros(nk,nband,num_proj);
    fgetl(fid);
    for i = 1:nband
        fgetl(fid);
        for j = 1:num_proj
            fgetl(fid);
            for k = 1:nk
                berry(k,i,j) = str2double(fgetl(fid));
            end
        end
    end
    
    % get kpoints
    kpoints = zeros(nk,3);
    fgetl(fid);
    for i = 1:nk
        kpoints(i,:) = sscanf(fgetl(fid),'%f %f %f');
    end
    fclose(fid);
    
    save(strcat('thi\',prefix,'\berry\fermi_data.mat'),'ngrid','nk',...
        'nband','efermi','proj','num_proj','eig','conduct','berry','kpoints');
    
end

% processing data

if isfile(strcat('thi\',prefix,'\berry\processed_data.mat'))
    load(strcat('thi\',prefix,'\berry\processed_data.mat'))
else

% duplicate data
kx = kpoints(:,1);
ky = kpoints(:,2);
kz = kpoints(:,3);

kpoints = [kx,ky,kz,ones(size(kx,1),1)];
conduct = sum(sum(conduct(:,:,n_type),3),2)/scale_factor;

[id,~] = get_inside_BZ(kpoints,hkl,G,inv_lat,true);
k0 = kpoints(id,:);
conduct_new = conduct(id,:);

for i = 1:size(G_hkl,1)
    kpoints = [kx+G_hkl(i,1),ky+G_hkl(i,2),...
        kz+G_hkl(i,3),ones(size(kx,1),1)];
    [id,~] = get_inside_BZ(kpoints,hkl,G,inv_lat,true);
    k0 = [k0;kpoints(id,:)];
    conduct_new = [conduct_new;conduct(id,:)];
end

[id,k0] = get_inside_BZ(k0,hkl,G,inv_lat,true);
conduct_new = conduct_new(id,:);

% convert to log scale
conduct_new = sign(conduct_new).*log10(abs(conduct_new)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpolate data
d1 = linspace(min(k0(:,1)),max(k0(:,1)),kdiv);
d2 = linspace(min(k0(:,2)),max(k0(:,2)),kdiv);
d3 = linspace(min(k0(:,3)),max(k0(:,3)),kdiv);
[x0,y0,z0] = ndgrid(d1,d2,d3);
XI = [x0(:) y0(:) z0(:)];
YI = griddatan(k0(:,1:3),conduct_new,XI);
YI = reshape(YI, size(x0));

save(strcat('thi\',prefix,'\berry\processed_data.mat'),'YI','x0','y0',...
    'z0','k0','conduct_new','n_type','kdiv');

end

% plot positive contribution
cmax = max(abs(conduct_new));
step = 0.05;
for i = 0:step:cmax
    p = patch(isosurface(x0,y0,z0,YI,i),'LineWidth',0.5);
    isonormals(x0,y0,z0,YI,p)
    set(p,'FaceColor','r','EdgeColor','none','FaceAlpha',i/cmax,...
        'AmbientStrength',0.8,'DiffuseStrength',0.8,'SpecularStrength',0.9,...
        'SpecularExponent',50,'SpecularColorReflectance',1);
end

for i = 0:-step:-cmax
    p = patch(isosurface(x0,y0,z0,YI,i),'LineWidth',0.5);
    isonormals(x0,y0,z0,YI,p)
    set(p,'FaceColor','b','EdgeColor','none','FaceAlpha',abs(i/cmax),...
        'AmbientStrength',0.8,'DiffuseStrength',0.8,'SpecularStrength',0.9,...
        'SpecularExponent',50,'SpecularColorReflectance',1);
end

% read high symmetry points
% read kpoints info
kpoints = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/kpts/coordinates');
kid = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/kpts/specialPointIndices');
klb = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/kpts/specialPointLabels');

% for i = 1:size(klb,1)
%    klb{i} = strtrim(klb{i});
%    if strcmp(klb{i},'gamma')
%       klb{i} = '\Gamma'; 
%    else
%        klb{i} = upper(klb{i}); 
%    end
%    tmp = kpoints(:,kid(i))'*inv_lat;
%    scatter3(tmp(1),tmp(2),tmp(3),30,'k','filled');
%    text(tmp(1),tmp(2),tmp(3)+0.1,klb{i},'FontSize',18)
% end

% plot custom points
custom_points = [0 0.5 0.5; 0.5 0.5 0;0 -0.5 -0.5;0.5,0,0.5];
custom_label = {'X','X1','X2','X3'};

for i = 1:size(custom_points,1)
    tmp = custom_points(i,:)*inv_lat;
    scatter3(tmp(1),tmp(2),tmp(3),30,'c','filled');
    text(tmp(1),tmp(2),tmp(3)+0.1,custom_label{i},'FontSize',18)
end

% % plot kx, ky, and kz axes
% 
% px1 = [-1,0,0]; px1 = px1*inv_lat;
% px2 = [1,0,0]; px2 = px2*inv_lat;
% quiver3(0,0,0,px2(1),px2(2),px2(3),0,'LineWidth',2);
% text(px2(1),px2(2),px2(3)+0.1,'kx','FontSize',18)
% 
% py1 = [0,-1,0]; py1 = py1*inv_lat;
% py2 = [0,1,0]; py2 = py2*inv_lat;
% quiver3(0,0,0,py2(1),py2(2),py2(3),0,'LineWidth',2)
% text(py2(1),py2(2),py2(3)+0.1,'ky','FontSize',18)
% 
% pz1 = [0,0,-1]; pz1 = pz1*inv_lat;
% pz2 = [0,0,1]; pz2 = pz2*inv_lat;
% quiver3(0,0,0,pz2(1),pz2(2),pz2(3),0,'LineWidth',2)
% text(pz2(1),pz2(2),pz2(3)+0.1,'kz','FontSize',18)

% plot x, y, and z axes

quiver3(0,0,0,1.5,0,0,0,'k','LineWidth',2);
text(1.5,0,0+0.1,'x','FontSize',18)

quiver3(0,0,0,0,1.5,0,0,'k','LineWidth',2)
text(0,1.5,0+0.1,'y','FontSize',18)

quiver3(0,0,0,0,0,1.5,0,'k','LineWidth',2)
text(0,0,1.5+0.1,'z','FontSize',18)


xlim([-zoom_factor zoom_factor])
ylim([-zoom_factor zoom_factor])
zlim([-zoom_factor zoom_factor])

axis vis3d

view(3)
camlight;
lighting gouraud;
% daspect([1,1,1])
axis off;
set(gcf,'color','w');
load('colormap_ud.mat');
colormap(SpinUpDown);
colorbar;
caxis([-cmax cmax]);

%slice along x axis
figure; hold on;

% find corners points
id = zeros(1,1);
cnt = 0;
for i = 1:size(cornersbz,1)
    if (abs(cornersbz(i,3)) < 0.01)
       cnt = cnt + 1;
       id(cnt) = i;
    end
end

slice_corners = cornersbz(id,:);

id = [1,7,8,4,2,5,6,3,1];

scatter3(slice_corners(:,1),slice_corners(:,2),slice_corners(:,3),10,'k','filled');

plot3(slice_corners(id,1),slice_corners(id,2),slice_corners(id,3),'k-','LineWidth',2.5);

%plot axis

quiver3(0,0,0,1.5,0,0,0,'k','LineWidth',2)
text(1.6,0,0.2,'x','FontSize',18)

quiver3(0,0,0,0.0,0,1.5,0,'k','LineWidth',2)
text(0,0,1.7,'z','FontSize',18)

quiver3(0,0,0,0.0,1.5,0,0,'k','LineWidth',2)
text(0,1.6,0.0,'y','FontSize',18)

xlim([-zoom_factor zoom_factor])
ylim([-zoom_factor zoom_factor])
zlim([-zoom_factor zoom_factor])

s = slice(x0,y0,z0,YI,[],[],0);
s.LineStyle = 'none';
view(0,0);

axis off;
axis vis3d;
set(gcf,'color','w');
load('colormap_ud.mat');
colormap(SpinUpDown);
colorbar;
caxis([-cmax cmax]);






