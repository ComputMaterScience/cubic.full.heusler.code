clear all;clc; figure; hold on;

prefix = 'Al';

% PBE: red
% LDA: blue


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PBE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot band structure from vasp

% load high-symmetry points
winfile = fopen(strcat('qanh\',prefix,'\pbe\shc\band\wannier90.win'),'rt');
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

A = import_outcar(strcat('qanh\',prefix,'\pbe\shc\data\OUTCAR'),'lat-basis'); A = A.rec;

[eigen, kpoints, occvalues, nelect ] = import_eigenval(strcat('qanh\',prefix,'\pbe\shc\band\EIGENVAL'));

nk = size(eigen,1);
nband = size(eigen,2);

s = zeros(nk,1);
bplot = zeros(nk,nband);
for i = 2:nk
    s(i) = s(i-1) + 2*pi*norm((kpoints(i,1:3)-kpoints(i-1,1:3))*A);
end

% find position of high symmetry points

kdiv = size(kpoints,1)/(size(label_k,2)-1);

x_pos = zeros(1,1);
cnt = 0;
for i = 1:kdiv:size(s,1)
    cnt = cnt + 1;
    x_pos(cnt) = s(i);
    plot([s(i) s(i)],[-40 40],'k--');
end

x_pos(cnt+1) = s(end);

plot([0 s(end)],[0 0],'k--','LineWidth',1.5);

eigen = eigen - efermi;

%plot
for i = 1:nband
    plot(s,eigen(:,i),'r-','LineWidth',1.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot band structure from vasp

% load high-symmetry points
winfile = fopen(strcat('qanh\',prefix,'\lda\shc\band\wannier90.win'),'rt');
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

A = import_outcar(strcat('qanh\',prefix,'\pbe\shc\data\OUTCAR'),'lat-basis'); A = A.rec;

[eigen, kpoints, occvalues, nelect ] = import_eigenval(strcat('qanh\',prefix,'\lda\shc\band\EIGENVAL'));

nk = size(eigen,1);
nband = size(eigen,2);

s = zeros(nk,1);
bplot = zeros(nk,nband);
for i = 2:nk
    s(i) = s(i-1) + 2*pi*norm((kpoints(i,1:3)-kpoints(i-1,1:3))*A);
end

% find position of high symmetry points

% kdiv = size(kpoints,1)/(size(label_k,2)-1);
% 
% x_pos = zeros(1,1);
% cnt = 0;
% for i = 1:kdiv:size(s,1)
%     cnt = cnt + 1;
%     x_pos(cnt) = s(i);
%     plot([s(i) s(i)],[-40 40],'k--');
% end
% 
% x_pos(cnt+1) = s(end);
% 
% plot([0 s(end)],[0 0],'k--','LineWidth',1.5);

eigen = eigen - efermi;

%plot
for i = 1:nband
    plot(s,eigen(:,i),'b-','LineWidth',1.5);
end

% figure options
box on
set(gca,'LineWidth',3.5)
set(gca,'FontSize',20);
ylabel('Energy (eV)')
xlim([s(1) s(end)]);
% ylim([-2 1]);
ylim([-15 11]);
xticks(x_pos);
xticklabels(label_k);

