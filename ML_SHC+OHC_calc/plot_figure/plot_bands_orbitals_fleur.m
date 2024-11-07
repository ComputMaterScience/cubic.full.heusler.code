clear all;clc; figure;

prefix = 'Al\pbe';

factor = 70;

% Hatree to eV
har2eV = 27.211396641308;

ang2bohr = 1.8897259886;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PBE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read kpoints info
kpoints = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/kpts/coordinates');
kid = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/kpts/specialPointIndices');
klb = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/kpts/specialPointLabels');

for i = 1:size(klb,1)
   klb{i} = strtrim(klb{i});
   if strcmp(klb{i},'gamma')
      klb{i} = '\Gamma'; 
   else
       klb{i} = upper(klb{i}); 
   end
end

% read lattice info
lat_rec = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/cell/reciprocalCell')'*ang2bohr;

% read eigenvalues
eig = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Local/BS/eigenvalues')*har2eV;

% read E-fermi
efermi = h5readatt(strcat('thi\',prefix,'\band\banddos.hdf'),'/general','lastFermiEnergy')*har2eV;

% read number of atom
natom = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/atoms/atomicNumbers');
natom = size(natom,1);

% 1: s
orb_s = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:1');
if natom == 2
    orb_s = orb_s + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:1');
end

% 2: px
orb_px = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:2');
if natom == 2
    orb_px = orb_px + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:2');
end

% 3: py
orb_py = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:3');
if natom == 2
    orb_py = orb_py + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:3');
end

% 4: pz
orb_pz = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:4');
if natom == 2
    orb_pz = orb_pz + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:4');
end

% 5: dxy
orb_dxy = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:5');
if natom == 2
    orb_dxy = orb_dxy + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:5');
end

% 6: dyz
orb_dyz = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:6');
if natom == 2
    orb_dyz = orb_dyz + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:6');
end

% 7: dxz
orb_dxz = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:7');
if natom == 2
    orb_dxz = orb_dxz + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:7');
end

% 8: dx2-y2
orb_dx2y2 = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:8');
if natom == 2
    orb_dx2y2 = orb_dx2y2 + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:8');
end

% 9: dz2
orb_dz2 = h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:1,ind:9');
if natom == 2
    orb_dz2 = orb_dz2 + h5read(strcat('thi\',prefix,'\band\banddos.hdf'),'/Orbcomp/BS/ORB:2,ind:9');
end

% scale fermi to zero
eig = eig - efermi;

% calculate k-path
kpath = zeros(size(kpoints,2),1);
cnt = 0;

for i = 2:size(kpoints,2)
    kpath(i) = kpath(i-1) + 2*pi*norm(lat_rec*(kpoints(:,i-1)-kpoints(:,i)));
end

% plot bands
hold on;

for i = 1:size(eig,1)
   plot(kpath,eig(i,:,1)','k-','LineWidth',1.5); 
end

% plot_wline2(kpath,eig(3,:,1)',orb_px(i,:)*factor,0.2,0.1,'r');

for i = 1:size(eig,1)
%     scatter(kpath,eig(i,:,1),orb_px(i,:)*scale_f+0.001,'r','filled');
     plot_wline(kpath,eig(i,:,1)',orb_px(i,:)*factor,0.2,0.1,'r');
     plot_wline(kpath,eig(i,:,1)',orb_py(i,:)*factor,0.2,0.1,'b');
     plot_wline(kpath,eig(i,:,1)',orb_pz(i,:)*factor,0.2,0.1,'g');
end

% plot k-zones
for i = 1:size(kid,1)
   plot([kpath(kid(i)) kpath(kid(i))],[-100 100],'k--'); 
end

% plot E-fermi line
plot([0 kpath(end)],[0 0],'k--','LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% format figure
xlim([0 kpath(end)]);
ylim([-15 11]);
xticks(kpath(kid))
xticklabels(klb);
ylabel('Energy (eV)')
box on
set(gca,'FontSize',18)
set(gca,'LineWidth',2.5);




