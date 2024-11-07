function [id,knew] = get_inside_BZ(kpoints,hkl,G,inv_lat,conv)
% convert to reciprocal space
if conv == true
    knew = zeros(size(kpoints,1),3);
    for i = 1:size(kpoints,1)
        knew(i,:) = kpoints(i,1:3)*inv_lat;
    end
else
    knew = kpoints(:,1:3);
end

% find the kpoints inside the Bz by choosing only the kpoints closer
% to Gamma than another G
for i = 1:size(knew,1)
    Dgamma = sqrt(knew(i,1)^2+knew(i,2)^2+knew(i,3)^2);
    for j = 1:size(hkl,1)
        if (j~=i)
            DGG = sqrt((knew(i,1)-G(j,1))^2+...
                (knew(i,2)-G(j,2))^2+(knew(i,3)-G(j,3))^2);
            if DGG < Dgamma
                kpoints(i,4) = 0;
                break
            end
        end
    end
end
id = find(kpoints(:,4)==1);
knew = knew(id,:);
end