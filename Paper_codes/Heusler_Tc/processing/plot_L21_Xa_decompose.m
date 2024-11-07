clear all;clc;figure;

load('redmap.mat');

data = readtable("data_full.csv");

elements = {'Co','Cr','Fe','Mn','Ni','Pd','Rh','Ru','Sc','Ti','V'};

% get stable struct
id = table2array(data(:,'form_en'))<0;
data = data(id,:);

num_struct = zeros(1,1);
num_struct_p = zeros(1,1);
num_struct_n = zeros(1,1);
for i = 1:size(elements,2)

    % test
    data_test = data(strcmp(data.X,elements{i}),:);
    
    
    num_struct_p(i) = sum(table2array(data_test(:,'dE'))>0);
    num_struct_n(i) = sum(table2array(data_test(:,'dE'))<0);
    num_struct(i) = size(data_test,1);

end

h = bar(1:11,[num_struct;num_struct_p;num_struct_n]);
xticklabels(elements)
