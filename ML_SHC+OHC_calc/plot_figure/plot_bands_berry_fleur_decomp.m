clear all;clc;figure; hold on;

prefix = 'Al\pbe';

calc_type = 'ohc';

factor = 50;

% n = 1 <px|Lz|py>
% n = 2 <py|Lz|px>
% n = 3 <dxz|Lz|dyz>
% n = 4 <dyz|Lz|dxz>
% n = 5 <dx2-y2|Lz|dxy>
% n = 6 <dxy|Lz|dx2-y2>

n = 1;

plot_mode = 2; % 1:berry, 2:conductivity

% load data
fid = fopen(strcat('thi\',prefix,'\berry\bandstructure_',calc_type,'_decomp.dat'),'r');

% get number of kpoints
line = fgetl(fid); line = split(line,':');
nk = str2double(line{2});

% get number of bands
line = fgetl(fid); line = split(line,':');
nband = str2double(line{2});

% get e-fermi
line = fgetl(fid); line = split(line,':');
efermi = str2double(line{2});

% get orbital projection
line = fgetl(fid); line = split(line,':');
proj = split(line{2}); proj = proj(2:end-1);

% get k-labels
line = fgetl(fid); line = split(line,':');
klabels = split(line{2}); klabels = klabels(2:end-1);

% get k-path
kpath = zeros(nk,1);
fgetl(fid);
for i = 1:nk
    kpath(i) = str2double(fgetl(fid));
end

for i = 1:size(klabels,1)
    if klabels{i} == 'G'
       klabels{i} = '\Gamma';
    end
end

% get k-label position
kpos = zeros(size(klabels,1),1);
fgetl(fid);
for i = 1:size(klabels,1)
    kpos(i) = str2double(fgetl(fid));
end

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

fclose(fid);

% processing data
eig = eig - efermi;

% scale in log10(|x|+1)
berry = sign(berry).*log10(abs(berry)+1);
conduct = sign(conduct).*log10(abs(conduct)+1);

%scatter berry values
if plot_mode == 1
    berry_n = zeros(size(berry,1),size(berry,2));
    berry_p = zeros(size(berry,1),size(berry,2));
    for i = 1:nband
        id_n = find(berry(:,i,n) < 0);
        id_p = find(berry(:,i,n) > 0);
        berry_n(id_n,i) = berry(id_n,i,n);
        berry_p(id_p,i) = berry(id_p,i,n);
    end
    for i = 1:nband
%         scatter(kpath,eig(:,i),factor,berry(:,i),'filled');
        plot_wline(kpath,eig(:,i),berry_n(:,i),0.2,1.0,'b');
        plot_wline(kpath,eig(:,i),berry_p(:,i),0.2,1.0,'r');
    end
else
    conduct_n = zeros(size(conduct,1),size(conduct,2));
    conduct_p = zeros(size(conduct,1),size(conduct,2));
    for i = 1:nband
        id_n = find(conduct(:,i,n) < 0);
        id_p = find(conduct(:,i,n) > 0);
        conduct_n(id_n,i) = conduct(id_n,i,n);
        conduct_p(id_p,i) = conduct(id_p,i,n);
    end
    for i = 1:nband
%         scatter(kpath,eig(:,i),factor,berry(:,i),'filled');
        plot_wline(kpath,eig(:,i),conduct_n(:,i),0.2,1.0,'b');
        plot_wline(kpath,eig(:,i),conduct_p(:,i),0.2,1.0,'r');
    end
end

% if plot_mode == 1
%     for i = 1:nband
%         scatter(kpath,eig(:,i),factor,berry(:,i,n),'filled');
%     end
% else
%     for i = 1:nband
%         scatter(kpath,eig(:,i),factor,conduct(:,i,n),'filled');
%     end
% end


% plot band structure
for i = 1:nband
    plot(kpath,eig(:,i),'k-','LineWidth',1.5);
end

plot([0 kpath(end)],[0 0],'k--','LineWidth',1.2);
for i = 2:size(kpos,1)-1
    plot([kpos(i) kpos(i)],[-20 20],'k--','LineWidth',1);
end

% format figure
load('colormap_ud.mat');
colormap(SpinUpDown);
xlim([0 kpath(end)]);
ylim([-15 6]);
caxis([-3 3]);
xticks(kpos)
xticklabels(klabels);
ylabel('Energy (eV)')
box on
set(gca,'FontSize',18)
set(gca,'LineWidth',2.5);
% set(gca,'color',[251, 248, 241]/255);






