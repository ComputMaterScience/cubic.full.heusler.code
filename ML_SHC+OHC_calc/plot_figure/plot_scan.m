clear all;clc;figure;hold on;

prefix = 'Ta\lda';

% fleur data

num_smooth = 3;

har2eV = 27.211396641308;

% read E-fermi
efermi = h5readatt(strcat('thi\',prefix,'\band\banddos.hdf'),'/general','lastFermiEnergy')*har2eV;

% load scan files

fid = fopen(strcat('thi\',prefix,'\scan\ohc_values.dat'),'r');

fgetl(fid); fgetl(fid); % ignore lines

data = zeros(1,2);
cnt = 0;

while ~feof(fid)
    try
        tmp = sscanf(fgetl(fid),'%f %f %f %f');
        cnt = cnt + 1;
        data(cnt,:) = [tmp(1)-efermi,tmp(2)/1000];
    catch
    end
end

fclose(fid);

if num_smooth ~= 0
    for i = 1:num_smooth
        data(:,2) = smooth(data(:,2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
efermi2 = import_outcar(strcat('qanh\',prefix,'\shc\data\OUTCAR'),'efermi'); 

fid = fopen(strcat('qanh\',prefix,'\ohc\scan\ohc_values.dat'),'r');

fgetl(fid); fgetl(fid); % ignore lines

data2 = zeros(1,2);
cnt = 0;

while ~feof(fid)
    try
        tmp = sscanf(fgetl(fid),'%f %f %f %f');
        cnt = cnt + 1;
        data2(cnt,:) = [tmp(1)-efermi2,tmp(2)/1000];
    catch
    end
end

fclose(fid);

if num_smooth ~= 0
    for i = 1:num_smooth
        data2(:,2) = smooth(data2(:,2));
    end
end

% plot
plot(data(:,2),data(:,1),'r-','LineWidth',2.5);
plot(data2(:,2),data2(:,1),'b-','LineWidth',2.5);
x_lim = xlim;

% format figure
ylabel('Energy (eV)')
xlabel('$\sigma_{OH} [10^3(\hbar/e)(\Omega\cdot{cm})^{-1}]$', 'Interpreter', 'latex')
box on
set(gca,'FontSize',18,'TickLength',[0.025,0.002])
set(gca,'LineWidth',3);
% ylim(x_lim);
plot([-12 12],[0,0],'k--','LineWidth',0.5);
plot([0 0],[-1,1],'k--','LineWidth',0.5);
xticks(-2:2:12)
yticks(-1:0.5:1)
ylim([-1,1]);
xlim([-4 12]);
