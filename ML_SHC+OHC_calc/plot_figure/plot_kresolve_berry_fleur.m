clear all;clc;figure; hold on;

prefix = 'Al\pbe';

calc_type = 'ohc';

factor = 50;

echg = 1.602176634e-19;

hbar = 1.054571817e-34;

scale_factor = (echg^2/hbar)*1.0e8;

% load data
fid = fopen(strcat('thi\',prefix,'\berry\bandstructure_',calc_type,'.dat'),'r');

% get number of kpoints
line = fgetl(fid); line = split(line,':');
nk = str2double(line{2});

% get number of bands
line = fgetl(fid); line = split(line,':');
nband = str2double(line{2});

% get e-fermi
line = fgetl(fid); line = split(line,':');
efermi = str2double(line{2});

% get k-labels
line = fgetl(fid); line = split(line,':');
klabels = split(line{2}); klabels = klabels(2:end-1);

for i = 1:size(klabels,1)
    if klabels{i} == 'G'
       klabels{i} = '\Gamma';
    end
end

% get k-path
kpath = zeros(nk,1);
fgetl(fid);
for i = 1:nk
    kpath(i) = str2double(fgetl(fid));
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
conduct = zeros(nk,nband);
fgetl(fid);
for i = 1:nband
    fgetl(fid);
    for k = 1:nk
        conduct(k,i) = str2double(fgetl(fid));
    end
end

% get berry
berry = zeros(nk,nband);
fgetl(fid);
for i = 1:nband
    fgetl(fid);
    for k = 1:nk
        berry(k,i) = str2double(fgetl(fid));
    end
end

fclose(fid);

% processing data
eig = eig - efermi;

% scale in log10(|x|+1)
% berry = sign(berry).*log10(abs(berry)+1);
% conduct = sign(conduct).*log10(abs(conduct)+1);

%scatter berry values
plot(kpath,sum(conduct,2)/scale_factor,'k-','LineWidth',1.5);

plot([0 kpath(end)],[0 0],'k--','LineWidth',1.2);
% for i = 2:size(kpos,1)-1
%     plot([kpos(i) kpos(i)],[-20 20],'k--','LineWidth',1);
% end

% format figure
xlim([0 kpath(end)]);
% ylim([-15 6]);
xticks(kpos)
xticklabels(klabels);
ylabel('Berry Curvature')

box on
set(gca,'FontSize',18)
set(gca,'LineWidth',2.5);

y_l = ylim;
for i = 2:size(kpos,1)-1
    plot([kpos(i) kpos(i)],[y_l(1) y_l(2)],'k--','LineWidth',1);
end

set(gcf,'Position',[441.8000  249.8000  641.6000  420.0000])




