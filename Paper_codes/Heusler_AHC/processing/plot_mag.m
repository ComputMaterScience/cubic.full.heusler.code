clear all;clc;

load('redmap.mat');

data = readtable("data_full.csv");

elements = {'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Ru','Rh','Pd'};

% get stable struct
id = table2array(data(:,'form_en'))<0;
data = data(id,:);

for i = 1:size(elements,2)

    figure;

    % test
    data_test = data(strcmp(data.X,elements{i}),:);
    h = heatmap(data_test,"Z","Y",'ColorVariable','magsum','ColorMethod','mean','MissingDataColor',[0.8 0.8 0.8]);

    % format figure
    h.Title = elements{i};
    clim(h,[0,4]);
    colormap(gca,redmap);
    set(gca,'FontSize', 12);
    set(gcf,'Position',[529.0000  151.4000  280.0000  546.4000]);
    saveas(gcf,strcat('MAG_',elements{i},'.eps'),'epsc')
    %close(gcf);
end
