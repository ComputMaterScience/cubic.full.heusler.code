function [param,data] = import_procar_ohc(filename)

if nargin == 0
    filename='PROCAR';
end

fid = fopen(filename);
if fid==-1
    error(['File ' filename ' not found']);
end

fgetl(fid); %jump first line

%get data
s = fgetl(fid);
param.num_kpoint = sscanf(s(15:25),'%d');
param.num_band = sscanf(s(41:50),'%d');
param.num_ion = sscanf(s(65:end),'%d');
for i = 1:param.num_kpoint
    fgetl(fid); %jump line
    s = fgetl(fid);
    data.kpoint1.pos{i} = sscanf(s(17:54),'%lf %lf %lf')';
    data.kpoint1.weight{i} = sscanf(s(66:end),'%lf');
    for j = 1:param.num_band
        band_atom_tmp = zeros(1,10);
        fgetl(fid); %jump line
        s = fgetl(fid);
        band_atom.energy{j} = sscanf(s(20:33),'%lf');
        band_atom.occupation{j} = sscanf(s(41:end),'%lf');
        s = fgetl(fid); %jump line
        s = fgetl(fid); %jump line
        for k = 1:param.num_ion
            s = fgetl(fid);
            band_atom_tmp(k,:) = sscanf(s(7:end),'%f %lf %lf %lf %lf %lf %lf %lf %lf %lf')';
        end
        if param.num_ion == 1
            for n = 1:3
                s = fgetl(fid);
            end
        else
            for n = 1:param.num_ion*3+3
                s = fgetl(fid);
            end
        end
        band_atom.band{j} = band_atom_tmp;
        if param.num_ion > 1
            s = fgetl(fid); %jump line
        end
    end
    data.kpoint1.data{i} = band_atom;
    s = fgetl(fid); %jump line
end
end
