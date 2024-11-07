clear all;clc; figure;

%red: PBE
%blue: LDA

prefix = 'V';

scale_f = 1;

% Hatree to eV
har2eV = 27.211396641308;

ang2bohr = 1.8897259886;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PBE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read kpoints info
kpoints = h5read(strcat('thi\',prefix,'\pbe\band\banddos.hdf'),'/kpts/coordinates');
kid = h5read(strcat('thi\',prefix,'\pbe\band\banddos.hdf'),'/kpts/specialPointIndices');
klb = h5read(strcat('thi\',prefix,'\pbe\band\banddos.hdf'),'/kpts/specialPointLabels');

for i = 1:size(klb,1)
   klb{i} = strtrim(klb{i});
   if strcmp(klb{i},'gamma')
      klb{i} = '\Gamma'; 
   else
       klb{i} = upper(klb{i}); 
   end
end

% read lattice info
lat_rec = h5read(strcat('thi\',prefix,'\pbe\band\banddos.hdf'),'/cell/reciprocalCell')'*ang2bohr;

% read eigenvalues
eig = h5read(strcat('thi\',prefix,'\pbe\band\banddos.hdf'),'/Local/BS/eigenvalues')*har2eV;

% read E-fermi
efermi = h5readatt(strcat('thi\',prefix,'\pbe\band\banddos.hdf'),'/general','lastFermiEnergy')*har2eV;

% % shift E-fermi to zero
if scale_f == 1
    eig = eig - efermi;
end

% calculate k-path
kpath = zeros(1,1);
cnt = 0;

for i = 2:size(kpoints,2)
    kpath(i,:) = kpath(i-1,:) + 2*pi*norm(lat_rec*(kpoints(:,i-1)-kpoints(:,i)));
end

% plot bands
hold on;

for i = 1:size(eig,1)
   plot(kpath,eig(i,:,1),'r-','LineWidth',2.5); 
end

% plot k-zones
for i = 1:size(kid,1)
   plot([kpath(kid(i)) kpath(kid(i))],[-100 100],'k--'); 
end

% plot E-fermi line
if scale_f == 1
    plot([0 kpath(end)],[0 0],'k--','LineWidth',1.5);
else
    plot([0 kpath(end)],[efermi efermi],'k--','LineWidth',1.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read kpoints info
kpoints = h5read(strcat('thi\',prefix,'\lda\band\banddos.hdf'),'/kpts/coordinates');
kid = h5read(strcat('thi\',prefix,'\lda\band\banddos.hdf'),'/kpts/specialPointIndices');
klb = h5read(strcat('thi\',prefix,'\lda\band\banddos.hdf'),'/kpts/specialPointLabels');

for i = 1:size(klb,1)
   klb{i} = strtrim(klb{i});
   if strcmp(klb{i},'gamma')
      klb{i} = '\Gamma'; 
   else
       klb{i} = upper(klb{i}); 
   end
end

% read lattice info
lat_rec = h5read(strcat('thi\',prefix,'\pbe\band\banddos.hdf'),'/cell/reciprocalCell')'*ang2bohr;

% read eigenvalues
eig = h5read(strcat('thi\',prefix,'\lda\band\banddos.hdf'),'/Local/BS/eigenvalues')*har2eV;

% read E-fermi
efermi = h5readatt(strcat('thi\',prefix,'\lda\band\banddos.hdf'),'/general','lastFermiEnergy')*har2eV;

% % shift E-fermi to zero
if scale_f == 1
    eig = eig - efermi;
end

% calculate k-path
kpath = zeros(1,1);
cnt = 0;

for i = 2:size(kpoints,2)
    kpath(i,:) = kpath(i-1,:) + 2*pi*norm(lat_rec*(kpoints(:,i-1)-kpoints(:,i)));
end

% plot bands
hold on;

for i = 1:size(eig,1)
   plot(kpath,eig(i,:,1),'b-','LineWidth',2.5); 
end
% 
% % plot k-zones
% for i = 1:size(kid,1)
%    plot([kpath(kid(i)) kpath(kid(i))],[-100 100],'k--'); 
% end

% % plot E-fermi line
% if scale_f == 1
%     plot([0 kpath(end)],[0 0],'k--','LineWidth',1.5);
% else
%     plot([0 kpath(end)],[efermi efermi],'k--','LineWidth',1.5);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% format figure
xlim([0 kpath(end)]);
ylim([-12 10]);
xticks(kpath(kid))
xticklabels(klb);
ylabel('Energy (eV)')
box on
set(gca,'FontSize',18)
set(gca,'LineWidth',2.5);




