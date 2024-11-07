clear all;clc;figure; hold on;

prefix = 'RhCoAl';

kdiv = 150;

pos = import_poscar(strcat(prefix,'/POSCAR'));
[param,data] = import_procar(strcat(prefix,'/PROCAR'),2);
efermi = import_outcar(strcat(prefix,'/OUTCAR'),'efermi');

% processing data
% kpath
d = zeros(1,1);
for i = 2:param.num_kpoint
    d(i) = d(i-1) + norm((data.kpoint1.pos{i}-data.kpoint1.pos{i-1})*pos.lattice);
end

% eigen
up = zeros(param.num_kpoint,param.num_band);
down = zeros(param.num_kpoint,param.num_band);
for i = 1:param.num_kpoint
    for j = 1:param.num_band
        up(i,j) = data.kpoint1.data{i}.energy{j} - efermi;
        down(i,j) = data.kpoint2.data{i}.energy{j} - efermi;
    end
end

% plot bands
% up
for i = 1:param.num_band
    plot(d,up(:,i),'LineWidth',2.0,'Color','b','LineStyle','-');
end

% down
for i = 1:param.num_band
    plot(d,down(:,i),'LineWidth',2.0,'Color','r','LineStyle','-');
end

% Figure option
xlim([d(1) d(end)]);
ylim([-10 5]);
ylabel('Energy (eV)');
set(gca,'fontsize',16);
box('on');
set(gca,'linewidth',3.0);