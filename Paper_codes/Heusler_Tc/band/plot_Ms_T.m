clear all;clc;figure; hold on;

prefix = 'new2/MnMnSn/xa';

data = import_vampire_Tc(strcat(prefix,'/output'));

mag = import_outcar(strcat(prefix,'/OUTCAR'),'mag_x');

% Tc = max([0 max(data(ischange(data(:,2),'linear'),1)),max(data(ischange(data(:,3),'linear'),1)),...
%     max(data(ischange(data(:,4),'linear'),1)),max(data(ischange(data(:,5),'linear'),1))]);

mag1 = data(:,2)*mag(1);
mag2 = data(:,3)*mag(2);
mag3 = data(:,4)*mag(3);
mag4 = data(:,5)*mag(4);

tot = mag1 + mag2 + mag3 + mag4;

if sum(tot) > 0
    type = 1;
else
    type = -1;
end



plot([0 max(data(:,1))],[0 0],'k:')

plot(data(:,1),tot*type,'-',Color='k',LineWidth = 2.0);

plot(data(:,1),mag2*type,'-','LineWidth',4.0,Color='#1679AB');
plot(data(:,1),mag3*type,'-','LineWidth',1.0,Color='#219C90');
plot(data(:,1),mag4*type,'-','LineWidth',2.0,Color='#850F8D');
plot(data(:,1),mag1*type,'-','LineWidth',2.0,Color='#EE4E4E'); % #C1

% plot([Tc Tc],[-5 5],'k--','LineWidth',1.5)

% Figure option
ylim([-4 4]);
xlabel('Temperature (K)');
ylabel('Ms (\mu_B)');
set(gca,'fontsize',16);
box('on');
set(gca,'linewidth',3.0);
