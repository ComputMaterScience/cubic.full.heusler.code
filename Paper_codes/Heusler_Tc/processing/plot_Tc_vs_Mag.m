clear all;clc; figure; hold on;

load('redmap.mat');

data = readtable("data_full_stable2.csv");

elements = {'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Ru','Rh','Pd'};

% get stable struct
id = table2array(data(:,'form_en'))<0 & table2array(data(:,'Tc'))>300;%& table2array(data(:,'magsum'))>0.1 & table2array(data(:,'Tc'))>0;
data = data(id,:);

% for i = 1:size(elements,2)
% 
%     
% 
%     % test
% %     data_test = data(strcmp(data.X,elements{i}),:);
% % %     h = heatmap(data_test,"Z","Y",'ColorVariable','Tc','ColorMethod','mean','MissingDataColor',[0.8 0.8 0.8]);
% % 
% %     scatter(table2array(data_test(:,'magsum')),table2array(data_test(:,'Tc')))
% 
% %     data_test = data(strcmp(data.X,'Sc'),:);
% %     scatter(table2array(data_test(:,'magsum')),table2array(data_test(:,'Tc')),10,'b')
% % 
% %     data_test = data(strcmp(data.X,'Fe'),:);
% %     scatter(table2array(data_test(:,'magsum')),table2array(data_test(:,'Tc')),10,'r')
% 
%     data_test = data(strcmp(data.X,elements{i}),:);
%     x = table2array(data_test(:,'ne'));
%     y = table2array(data_test(:,'form_en'));
%     scatter(x,y,30,'b')
% 
%     % Figure option
% %     xlim([20 36]);
% %     ylim([0 1200]);
% 
% %     data_test = data(strcmp(data.X,'Ni'),:);
% %     scatter(table2array(data_test(:,'magsum')),table2array(data_test(:,'Tc')),10,'m')
% 
% %     % format figure
% %     h.Title = elements{i};
% %     clim(h,[0,1000]);
% %     colormap(gca,redmap);
% %     set(gca,'FontSize', 12);
% %     set(gcf,'Position',[507.4000   93.0000  450.4000  548.0000]);
% %     saveas(gcf,strcat('Tc_',elements{i},'.eps'),'epsc')
% %     %close(gcf);
% end
