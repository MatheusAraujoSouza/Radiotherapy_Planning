%About "dir" command, read the "output arguments" part
%https://www.mathworks.com/help/matlab/ref/dir.html#bup_1_c-3

beamdir = './';
allFiles = dir([beamdir '*D.mat']);


%Regarding the allFiles structure, ".name" gets the name columns. If we had
%used ".folder" instead, we would have taken the "folder" column.
allNames = { allFiles.name };
d = [];

for i=1:length(allNames)
    f = allNames{i};

    %load the matrix D for the gantry/couch pair
    %stored in the file f (but I don't know it takes the name "D").
    load(f)

    %nv = number of voxels
    %nb = nuber of beams
    [nv, nb]=size(D);
    
    %Computation of d
    if i==1
        d = D*ones(nb,1);
    else
        d = d + D*ones(nb,1);
    end
end

fname = 'OuterTarget_VOILIST.mat';
%load a vector v of the target voxels indices:
load(fname)

%get dose distribution for just those voxels:
dstruct = d(v);

%compute and display dose stats

%gets the lowest value of the vector dstruct
dmin = min(dstruct);
%gets the weighted average of the vector dstruct
dmean = mean(dstruct);
%gets the highest value of the vector dstruct
dmax = max(dstruct);

%The solution mean that the entire range of the dose is between the minimum
%value and the maximum value. And the average dose is given by "dmean".
disp(['min, mean, max = ' num2str(dmin) ', ' ...
num2str(dmean) ', ' num2str(dmax)]);
