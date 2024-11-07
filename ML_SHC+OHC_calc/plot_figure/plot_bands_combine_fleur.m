clear all;clc; figure;

prefix = 'Al\pbe';

scale_f = 1;

% Hatree to eV
har2eV = 27.211396641308;

ang2bohr = 1.8897259886;

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

% % shift E-fermi to zero
if scale_f == 1
    eig = eig - efermi;
end

% calculate k-path
kpath = zeros(1,1);
cnt = 0;

for i = 2:size(kpoints,2)
    kpath(i,:) = kpath(i-1,:) + norm(lat_rec*kpoints(:,i-1)-lat_rec*kpoints(:,i));
end

% plot bands
hold on;

for i = 1:size(eig,1)
   plot(kpath,eig(i,:,1),'k-','LineWidth',2.5); 
end

% plot k-zones
for i = 1:size(kid,1)
   plot([kpath(kid(i)) kpath(kid(i))],[-100 100],'k--'); 
end

% plot E-fermi line
if scale_f == 1
    plot([0 kpath(end)],[0 0],'k--','LineWidth',1.5);
else
    plot([0 kpath(end)],[efermi efermi],'k--','LineWidth',1.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load k-points
kfile = fopen(strcat('thi\',prefix,'\band\WF1-path.kpt'),'rt');
nk = str2num(fgetl(kfile));
kpoints = zeros(1,4);
for i = 1:nk
    line = fgetl(kfile);
    kpoints(i,:) = sscanf(line, '%f %f %f %f');
end
fclose(kfile);

% load band strucutre data
data_band = load(strcat('thi\',prefix,'\band\WF1-bands.dat'));

% processing data
data_band(:,2) = data_band(:,2) - efermi;

% plot
scatter(data_band(1:20:end,1),data_band(1:20:end,2),20,'r');

% format figure
xlim([0 kpath(end)]);
ylim([-12 10]);
xticks(kpath(kid))
xticklabels(klb);
ylabel('Energy (eV)')
box on
set(gca,'FontSize',18)
set(gca,'LineWidth',2.5);




