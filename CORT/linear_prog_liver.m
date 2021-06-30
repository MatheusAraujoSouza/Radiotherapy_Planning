%Gantry and couch angles to use:
ga = [58 106 212 328 216 226 296];
ca = [0 0 0 0 32 -13 17];
Dij = [];

%Form the Dij matrix
% num2string - convert a number to string
for i=1:length(ga)
    fname = ['Gantry' num2str(ga(i)) '_Couch' num2str(ca(i)) '_D.mat'];
    load(fname)
    
    %number of beamlets at the i-th angle
    % size(D,2) takes the size of the second dimension; the number os
    % columns
    nba(i) = size(D,2);
    
    if i==1
        Dij = D;
    else
        Dij = [Dij D];
    end
end

%Load structures. These are the prescriptions for each voxel of the
%respective structures
load('PTV_VOILIST.mat');
V{1} = v;
load('Liver_VOILIST.mat');
V{2} = v;
load('Heart_VOILIST.mat');
V{3} = v;
load('entrance_VOILIST.mat');
V{4} = v;
load('Skin_VOILIST.mat');
V{5} = v;

% mean doses contributions of all beamlets
Dlivermean = mean(Dij(V{2},:));
Dheartmean = mean(Dij(V{3},:));
Dentrancemean = mean(Dij(V{4},:));

%construct the linear inequality constraints
%to enforce a minimum dose of 1 to the PTV
A = -Dij(V{1},:);
b = -1*ones(size(A,1),1);

%cost vector
c = Dlivermean + Dheartmean + 0.6*Dentrancemean;

%bounds in beamlet intensity
nb = size(Dij,2); % total number of beamlets
lb = zeros(nb,1); % lower bound of zero
ub = 25*ones(nb,1); % upper bound of 25 MU

%optimization options
opt = optimset('Display','iter');
opt.LargeScale = 'on';

%solve problem using matlab's LP solver
[x, fval, eflag] = linprog(c,A,b,[],[],lb,ub,[],opt);

%save solution in our recommended format
ctr = 1;
for i=1:length(ga)
    fname = ['Gantry' num2str(ga(i)) '_Couch' ...
    num2str(ca(i)) '_beamletSol.mat'];
    %num beamlets at angle i = nba(i)
    beamx = x(ctr:ctr+nba(i)-1);
    save(fname,'beamx');
    ctr = ctr+nba(i);
end

%report mean doses to structs:

Dlivermean*x
Dheartmean*x
Dentrancemean*x

%calculate dose distribution
d = Dij*x;
%reshape dose vector to 3-dimensional array

%dose grid dimensions
dim = [217 217 168];
%total number of voxels
nVoxels = 217*217*168;
%reshape dose vector
dose = reshape(d,dim);
%create 3-dimensional masks for structures
for s=1:length(V)
    mask{s} = zeros(nVoxels,1);
    for i=1:length(V{s})
        mask{s}(V{s}(i))=1;
    end
    mask{s} = reshape(mask{s},dim);
end

%select axial slice to plot
slice = 50;

%plot dose and structures
figure;
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'YDir','rev');
axis([30 190 40 160]);
hold on;

%plot colorwash dose
imagesc(dose(:,:,slice));

%plot contours for PTV, Liver, Skin
for s=[1 2 5]
    contour(mask{s}(:,:,slice),0.5,'k','LineWidth',2);
end
hold off;