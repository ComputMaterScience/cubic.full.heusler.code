clear all;clc;figure;hold on;

elements = {'Ta','W','Pt'};

pot_type = 'pbe';

% load fleur data

data_fleur = zeros(1,1);

for i = 1:size(elements,2)
    fid = fopen(strcat('thi\',elements{i},'\',pot_type,'\data\wannier90.w4hout'),'r');
    while ~feof(fid)
        line = fgetl(fid);
        if strfind(line,'# Orbital hall conductivity')
            fgetl(fid);line = fgetl(fid);
            line = split(line);
            data_fleur(i) = str2double(line(2))/1000;
            break
        end
    end
    fclose(fid);
end

% load vasp data

data_vasp = zeros(1,1);

for i = 1:size(elements,2)
    fid = fopen(strcat('qanh\',elements{i},'\',pot_type,'\ohc\data\wannier90.w4hout'),'r');
    while ~feof(fid)
        line = fgetl(fid);
        if strfind(line,'# Orbital hall conductivity')
            fgetl(fid);line = fgetl(fid);
            line = split(line);
            data_vasp(i) = str2double(line(2))/1000;
            break
        end
    end
    fclose(fid);
end

% plot OHC

plot([0 size(elements,2)+1],[0 0],'k--','LineWidth',1.0);

% for i = 1:size(elements,2)
%     plot([i i],[-20 20],'k:','LineWidth',1.0);
% end

plot(1:size(elements,2),data_fleur,'ro-','LineWidth',2.5,'MarkerFaceColor','r','MarkerSize',18);
plot(1:size(elements,2),data_vasp,'bs--','LineWidth',2.5,'MarkerFaceColor','b','MarkerSize',15);

% format figure
ylabel('$\sigma_{OH} [10^3(\hbar/e)(\Omega\cdot{cm})^{-1}]$', 'Interpreter', 'latex')
box on
set(gca,'FontSize',18,'TickLength',[0.025,0.002])
set(gca,'LineWidth',3);
xticks(1:size(elements,2))
xticklabels(elements)
axis square

xlim([0.5 size(elements,2)+0.5]);
ylim([-2,11]);
