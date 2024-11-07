clear all;clc;figure; hold on;

prefix = 'RhCoAl';

[energy,total_dos,efermi,pdos] = import_doscar(strcat(prefix,'/DOSCAR'));

% processing data
energy = energy - efermi;

% plot(energy,total_dos(:,1),'k-','LineWidth',2.0);
% plot(energy,-total_dos(:,2),'k-','LineWidth',2.0);

% atom resolved
dos_X_up_xy = sum(sum(pdos(:,9,3),2),3);
dos_X_up_yz = sum(sum(pdos(:,11,3),2),3);
dos_X_up_z2 = sum(sum(pdos(:,13,3),2),3);
dos_X_up_xz = sum(sum(pdos(:,15,3),2),3);
dos_X_up_x2y2 = sum(sum(pdos(:,17,3),2),3);

dos_X_down_xy = sum(sum(pdos(:,10,3),2),3);
dos_X_down_yz = sum(sum(pdos(:,12,3),2),3);
dos_X_down_z2 = sum(sum(pdos(:,14,3),2),3);
dos_X_down_xz = sum(sum(pdos(:,16,3),2),3);
dos_X_down_x2y2 = sum(sum(pdos(:,18,3),2),3);

% Type1
% plot(energy,dos_X_up_xy+dos_X_up_yz+dos_X_up_xz,'-','LineWidth',2.0,'Color','#0A6847');
% plot(energy,-(dos_X_down_xy+dos_X_down_yz+dos_X_down_xz),'-','LineWidth',2.0,'Color','#0A6847');
% 
% plot(energy,dos_X_up_z2+dos_X_up_x2y2,'-','LineWidth',2.0,'Color','#FB6D48');
% plot(energy,-(dos_X_down_z2+dos_X_down_x2y2),'-','LineWidth',2.0,'Color','#FB6D48');
ylim([-4 4]);

% Type2
plot(energy,-(dos_X_up_xy+dos_X_up_yz+dos_X_up_xz),'-','LineWidth',2.0,'Color','#0A6847');
plot(energy,(dos_X_down_xy+dos_X_down_yz+dos_X_down_xz),'-','LineWidth',2.0,'Color','#0A6847');

plot(energy,-(dos_X_up_z2+dos_X_up_x2y2),'-','LineWidth',2.0,'Color','#FB6D48');
plot(energy,(dos_X_down_z2+dos_X_down_x2y2),'-','LineWidth',2.0,'Color','#FB6D48');
ylim([-10 10]);

plot([-8,6],[0 0],'k-','LineWidth',1.0);
plot([0,0],[-10 10],'k--','LineWidth',1.0);

% Figure option
xlim([-6 2]);
xlabel('Energy (eV)');
ylabel('DOS (eV/states/atom)');
set(gca,'fontsize',16);
box('on');
xticks(-6:2:2)
set(gca,'linewidth',3.0);
set(gcf,'Position',[481.0000  289.8000  560.0000  286.4000]);