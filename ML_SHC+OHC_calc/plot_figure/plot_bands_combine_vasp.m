clear all;clc; figure; hold on;

prefix = 'Fe\lda';

% load k-points
kfile = fopen(strcat('qanh\',prefix,'\shc\band\wannier90-path.kpt'),'rt');
nk = str2num(fgetl(kfile));
kpoints = zeros(1,4);
for i = 1:nk
    line = fgetl(kfile);
    kpoints(i,:) = sscanf(line, '%f %f %f %f');
end
fclose(kfile);
% load high-symmetry points
winfile = fopen(strcat('qanh\',prefix,'\shc\band\wannier90.win'),'rt');
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
high_k = zeros(1,3);
label_k = cell(1,1);
for i = 1:size(k_path,2)-1
    high_k(i,:) = [str2double(k_path{i}{2}) str2double(k_path{i}{3}) str2double(k_path{i}{4})];
    label_k{i} = k_path{i}{1};
end
high_k(i+1,:) = [str2double(k_path{i}{6}) str2double(k_path{i}{7}) str2double(k_path{i}{8})];
label_k{i+1} = k_path{i}{5};

for i = 1:size(label_k,2)
    if strcmp(label_k{i},'G')
        label_k{i} = '\Gamma';
    end
    if strcmp(label_k{i},'K2')
        label_k{i} = 'K''';
    end
end

rec_vec = inv(lat_vec)';

fclose(winfile);

% load band strucutre data
data_band = load(strcat('qanh\',prefix,'\shc\band\wannier90-bands.dat'));

% processing data
data_band(:,2) = data_band(:,2) - efermi;

% find position of high symmetry points

id_high_k = zeros(1,1);
for i = 1:size(high_k,1)
    for j = 1:size(kpoints,1)
        if kpoints(j,1) == high_k(i,1) && kpoints(j,2) == high_k(i,2) && kpoints(j,3) == high_k(i,3)
            id_high_k(i) = j;
            kpoints(j,:) = [-10 -10 -10 -10];
            break
        end
    end
end

% plot data
plot([0 max(data_band(:,1))],[0 0],'-k','LineWidth',1.5);

for i = 2:size(high_k,1)-1
    plot([data_band(id_high_k(i),1) data_band(id_high_k(i),1)],[-50 50],'-k','LineWidth',1.5);
end

% scatter(data_band(:,1),data_band(:,2),5,'k','filled');

%%% plot band structure from vasp

A = import_outcar(strcat('qanh\',prefix,'\shc\data\OUTCAR'),'lat-basis'); A = A.rec;

[eigen, kpoints, occvalues, nelect ] = import_eigenval(strcat('qanh\',prefix,'\shc\band\EIGENVAL'));

nk = size(eigen,1);
nband = size(eigen,2);

s = zeros(nk,1);
bplot = zeros(nk,nband);
for i = 2:nk
    s(i) = s(i-1) + 2*pi*norm((kpoints(i,1:3)-kpoints(i-1,1:3))*A);
end

eigen = eigen - efermi;

%plot
for i = 1:nband
    plot(s,eigen(:,i),'k-','LineWidth',1.5);
end

scatter(data_band(1:20:end,1),data_band(1:20:end,2),20,'r');

% figure options
box on
set(gca,'LineWidth',3.5)
set(gca,'FontSize',20);
ylabel('Energy (eV)')
xlim([min(data_band(:,1)) max(data_band(:,1))]);
% ylim([-2 1]);
ylim([-15 11]);
xticks(data_band(id_high_k,1));
xticklabels(label_k);

