function M = estMotionIter3(vol1,vol2,numIters,Minitial,rotFlag,robustFlag,crop,CB,SC)
%
% function M = estMotionIter3(vol1,vol2,numIters,Minitial,[rotFlag],[robustFlag],[crop])
%
% vol1 and vol2 are volumes, 3d arrays
% numIters is number of iterations to run
% Minitial is initial guess for M.  Default is 3x3 identity matrix.
% crop specifies border size to crop/ignore  around all sides of the volume. 
%     Should be of the form [ymin xmin zmin; ymax xmax zmax]
%     Default crops 2 pixel border: [2 2 2; (size(vol1) - [2 2 2])].
%
% M is 4x4 translation+rotation or affine transform matrix: X' = M X
% where X=(x,y,1) is starting position in homogeneous coords
% and X'=(x',y'',1) is ending position
%
% where X=(x,y,z,1) is starting position in homogeneous coords
% and X'=(x',y',z',1) is ending position
%
% Each iteration warps the volumes according to the previous
% estimate, and estimates the residual motion.
%
% robustFlag is passed to estMotion3 (if activated, uses robust
% M-estimator).

% default values
if ieNotDefined('robustFlag')
    robustFlag = 0;
end
if ieNotDefined('rotFlag')
    rotFlag = 1;
end
if ieNotDefined('crop')
    crop = [2 2 2; (size(vol1) - [2 2 2])];
end
if ieNotDefined('CB')
    CB = [];
end
if ieNotDefined('SC')
    SC = [];
end
if ieNotDefined('numIters')
    numIters=3;
end

if ieNotDefined('Minitial')
    M=eye(4);
else
    M=Minitial;
end

for iter=1:numIters
  Mhalf2=real(sqrtm(M));
  Mhalf1=real(sqrtm(inv(M)));
  volWarp1=warpAffine3(vol1,Mhalf1);
  volWarp2=warpAffine3(vol2,Mhalf2);
  deltaM=estMotion3(volWarp1,volWarp2,rotFlag,robustFlag,crop,CB,SC);
  M=deltaM*M;
end

return;


%%%%%%%%%
% Debug %
%%%%%%%%%
Niter=3;

% input
filter = [0.03504 0.24878 0.43234 0.24878 0.03504];
in = convXYsep(convZ(rand(30,40,20),filter),filter,filter);

% translation
A= [1 0 0 .2;
    0 1 0 .3;
    0 0 1 .4;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
% default - rot and LS
ArotLS = estMotionIter3(vol1,vol2,Niter);
% rot and robust
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1);
% affine and LS
AaffLS = estMotionIter3(vol1,vol2,Niter,[],0,0);
% affine and robust
AaffRob = estMotionIter3(vol1,vol2,Niter,[],0,1);
A*A
ArotLS
ArotRob
AaffLS
AaffRob

% test crop
crop = [2 2 2; (size(in) - [2 2 2])];
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1,crop)
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1,[])

% rotation in x,y
theta=.03;
A=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
% default - rot and LS
ArotLS = estMotionIter3(vol1,vol2,Niter);
% rot and robust
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1);
% affine and LS
AaffLS = estMotionIter3(vol1,vol2,Niter,[],0,0);
% affine and robust
AaffRob = estMotionIter3(vol1,vol2,Niter,[],0,1);
A*A
ArotLS
ArotRob
AaffLS
AaffRob

% rotation in x,z
theta=.03;
A=[cos(theta) 0 sin(theta) 0;
    0 1 0 0;
    -sin(theta) 0 cos(theta) 0;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
% default - rot and LS
ArotLS = estMotionIter3(vol1,vol2,Niter);
% rot and robust
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1);
% affine and LS
AaffLS = estMotionIter3(vol1,vol2,Niter,[],0,0);
% affine and robust
AaffRob = estMotionIter3(vol1,vol2,Niter,[],0,1);
A*A
ArotLS
ArotRob
AaffLS
AaffRob

%%%%%%%%%%%%%%%%%%%%
% test with outliers
%%%%%%%%%%%%%%%%%%%%

% translation
A= [1 0 0 .2;
    0 1 0 .3;
    0 0 1 .4;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
% putting inconsistent information in upper left corner of vol2
Nc=3;
dims=size(vol2);
vol2(1:round(dims(1)/Nc), 1:round(dims(2)/Nc), :) = ...
			rand(round(dims(1)/Nc), round(dims(2)/Nc), dims(3));
% default - rot and LS
ArotLS = estMotionIter3(vol1,vol2,Niter);
% rot and robust
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1);
% affine and LS
AaffLS = estMotionIter3(vol1,vol2,Niter,[],0,0);
% affine and robust
AaffRob = estMotionIter3(vol1,vol2,Niter,[],0,1);
A*A
ArotLS
ArotRob
AaffLS
AaffRob

% rotation in x,y
theta=.03;
A=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
% putting inconsistent information in upper left corner of vol2
Nc=3;
dims=size(vol2);
vol2(1:round(dims(1)/Nc), 1:round(dims(2)/Nc), :) = ...
			rand(round(dims(1)/Nc), round(dims(2)/Nc), dims(3));
% default - rot and LS
ArotLS = estMotionIter3(vol1,vol2,Niter);
% rot and robust
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1);
% affine and LS
AaffLS = estMotionIter3(vol1,vol2,Niter,[],0,0);
% affine and robust
AaffRob = estMotionIter3(vol1,vol2,Niter,[],0,1);
A*A
ArotLS
ArotRob
AaffLS
AaffRob

% rotation in x,z
theta=.03;
A=[cos(theta) 0 sin(theta) 0;
    0 1 0 0;
    -sin(theta) 0 cos(theta) 0;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
% putting inconsistent information in upper left corner of vol2
Nc=3;
dims=size(vol2);
vol2(1:round(dims(1)/Nc), 1:round(dims(2)/Nc), :) = ...
			rand(round(dims(1)/Nc), round(dims(2)/Nc), dims(3));
% default - rot and LS
ArotLS = estMotionIter3(vol1,vol2,Niter);
% rot and robust
ArotRob = estMotionIter3(vol1,vol2,Niter,[],1,1);
% affine and LS
AaffLS = estMotionIter3(vol1,vol2,Niter,[],0,0);
% affine and robust
AaffRob = estMotionIter3(vol1,vol2,Niter,[],0,1);
A*A
ArotLS
ArotRob
AaffLS
AaffRob
