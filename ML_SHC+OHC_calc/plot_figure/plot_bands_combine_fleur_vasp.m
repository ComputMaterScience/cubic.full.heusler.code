clear all;clc; figure;

trans_factor = 0.3;

% fleur data
prefix = 'Fe\lda';

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

for i = 2:size(kpoints,2)
    kpath(i,:) = kpath(i-1,:) + 2*pi*norm(lat_rec*kpoints(:,i-1)-lat_rec*kpoints(:,i));
end

% plot bands
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FLEUR
for i = 1:size(eig,1)
   plot(kpath,eig(i,:,1),'r-','LineWidth',3.3); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scatter(2*pi*data_band(1:60:end,1),data_band(1:60:end,2),28,'r','filled');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xticks(kpath(kid))
xticklabels(klb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vasp data

% load k-points
kfile = fopen(strcat('qanh\',prefix,'\ohc\band\wannier90-path.kpt'),'rt');
nk = str2num(fgetl(kfile));
kpoints = zeros(1,4);
for i = 1:nk
    line = fgetl(kfile);
    kpoints(i,:) = sscanf(line, '%f %f %f %f');
end
fclose(kfile);
% load high-symmetry points
winfile = fopen(strcat('qanh\',prefix,'\ohc\band\wannier90.win'),'rt');
lat_vec = zeros(3,3);
k_path = cell(1,1);
cnt = 0;
while ~feof(winfile)
    line = fgetl(winfile);
    if contains(line,'begin unit_cell_cart')
        line = fgetl(winfile);
        for i = 1:3
            line = fgetl(winfile);
            lat_vec(i,:) = sscanf(line, '%f %f %f');
        end
    end
    if contains(line,'num_wann')
        num_wann = split(line,'='); num_wann = str2double(num_wann{2});
    end
    if contains(line,'begin kpoint_path')
        while ~contains(line,'end kpoint_path')
            line = fgetl(winfile);
            cnt = cnt + 1;
            k_path{cnt} = split(line);
        end
    end
    if ~isempty(regexp(line, '\<fermi_energy\>', 'once'))
        efermi = split(line); efermi = str2double(efermi{3});
    end
    if contains(line,'kpath_num_points')
        kdiv = split(line); kdiv = str2double(kdiv{3});
    end
end

rec_vec = inv(lat_vec)';

fclose(winfile);

% load band strucutre data
data_band = load(strcat('qanh\',prefix,'\ohc\band\wannier90-bands.dat'));

% processing data
data_band(:,2) = data_band(:,2) - efermi;

w90_band = reshape(data_band(:,2),[nk,num_wann]);

% calculate k-path
kpath = zeros(1,1);
for i = 2:size(kpoints,1)
    kpath(i) = kpath(i-1) + 2*pi*norm(kpoints(i-1,1:3)*lat_rec'-kpoints(i,1:3)*lat_rec');
end

%%% plot band structure from vasp

% A = import_outcar(strcat('qanh\',prefix,'\shc\data\OUTCAR'),'lat-basis'); 
A = lat_rec';

[eigen, kpoints_scf, occvalues, nelect ] = import_eigenval(strcat('qanh\',prefix,'\ohc\band\EIGENVAL'));

nk_scf = size(eigen,1);
nband = size(eigen,2);

s = zeros(nk_scf,1);
bplot = zeros(nk_scf,nband);
for i = 2:nk_scf
    s(i) = s(i-1) + 2*pi*norm((kpoints_scf(i,1:3)-kpoints_scf(i-1,1:3))*A);
end

eigen = eigen - efermi;

%plot###########################################################VASP
for i = 1:nband
    plot(s,eigen(:,i),'-','Color',[0 0 1 trans_factor],'LineWidth',3.5);
end

for i = 1:num_wann
   scatter(kpath(1:28:end),w90_band(1:28:end,i),35,'bs','filled'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% format figure
xlim([0 kpath(end)]);
ylim([-13 5]);
ylabel('Energy (eV)')
box on
set(gca,'FontSize',18,'TickLength',[0.025,0.002])
set(gca,'LineWidth',2);
axis square
yticks(-12:2:4)




