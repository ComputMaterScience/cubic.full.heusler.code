function result = import_vampire_Tc(filename)

if nargin == 0
    filename='output';
end

fid = fopen(filename);
if fid==-1
    result = 0;
else
    % ignore comment lines
    for i = 1:8
        fgetl(fid);
    end

    data = zeros(1,18);
    cnt = 0;
    while ~feof(fid)
        line = fgetl(fid);
        line = split(line);
        line = str2double(line(~cellfun('isempty',line)));
        cnt = cnt + 1;
        data(cnt,:) = line;
    end

    % Tc = max([0 max(data(ischange(data(:,6),'linear'),1)),max(data(ischange(data(:,10),'linear'),1)),...
        % max(data(ischange(data(:,14),'linear'),1)),max(data(ischange(data(:,18),'linear'),1))]);

    result = [data(:,1), data(:,6), data(:,10), data(:,14), data(:,18)];
    
    fclose(fid);
end

end