clear all;clc;

load('redbluemap.mat');

data = readtable("data_full.csv");

elements = {'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Ru','Rh','Pd'};

for i = 1:size(elements,2)

    figure;

    % test
    data_test = data(strcmp(data.X,elements{i}),:);
    h = heatmap(data_test,"Z","Y",'ColorVariable','form_en','ColorMethod','mean');

    % format figure
    h.Title = elements{i};
    clim(h,[-0.5,0.5]);
    colormap(gca,redbluemap);
    set(gca,'FontSize', 12);
    set(gcf,'Position',[529.0000  151.4000  280.0000  546.4000]);
    saveas(gcf,strcat(elements{i},'.eps'),'epsc')
    close(gcf);
end
