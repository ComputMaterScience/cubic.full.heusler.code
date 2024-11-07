clear all;clc;figure; hold on;

prefix = 'new2/MnMnSn/xa';

[energy,total_dos,efermi,pdos] = import_doscar(strcat(prefix,'/DOSCAR'));

% processing data
energy = energy - efermi;

% plot(energy,total_dos(:,1),'k-','LineWidth',2.0);
% plot(energy,-total_dos(:,2),'k-','LineWidth',2.0);

% atom resolved
dos_X_up = sum(sum(pdos(:,1:2:17,[1,2]),2),3);
dos_X_down = sum(sum(pdos(:,2:2:18,[1,2]),2),3);

dos_Y_up = sum(sum(pdos(:,1:2:17,3),2),3);
dos_Y_down = sum(sum(pdos(:,2:2:18,3),2),3);

dos_Z_up = sum(sum(pdos(:,1:2:17,4),2),3);
dos_Z_down = sum(sum(pdos(:,2:2:18,4),2),3);

% type1
plot(energy,dos_X_up,'-','LineWidth',2.0,'Color','#FF0000');
plot(energy,-dos_X_down,'-','LineWidth',2.0,'Color','#FF0000');

plot(energy,dos_Y_up,'-','LineWidth',2.0,'Color','#4c89d2');
plot(energy,-dos_Y_down,'-','LineWidth',2.0,'Color','#4c89d2');

plot(energy,dos_Z_up,'-','LineWidth',2.0,'Color','#71b578');
plot(energy,-dos_Z_down,'-','LineWidth',2.0,'Color','#71b578');

% % type2
% plot(energy,-dos_X_up,'-','LineWidth',2.0,'Color','#FF0000');
% plot(energy,dos_X_down,'-','LineWidth',2.0,'Color','#FF0000');
% 
% plot(energy,-dos_Y_up,'-','LineWidth',2.0,'Color','#4c89d2');
% plot(energy,dos_Y_down,'-','LineWidth',2.0,'Color','#4c89d2');
% 
% plot(energy,-dos_Z_up,'-','LineWidth',2.0,'Color','#71b578');
% plot(energy,dos_Z_down,'-','LineWidth',2.0,'Color','#71b578');

plot([-8,6],[0 0],'k-','LineWidth',1.0);
plot([0,0],[-10 10],'k--','LineWidth',1.0);

% Figure option
xlim([-6 2]);
ylim([-8 8]);
xlabel('Energy (eV)');
ylabel('DOS (eV/states/atom)');
set(gca,'fontsize',16);
box('on');
xticks(-6:2:2)
yticks(-8:4:8)
set(gca,'linewidth',3.0);
set(gcf,'Position',[481.0000  289.8000  560.0000  286.4000]);