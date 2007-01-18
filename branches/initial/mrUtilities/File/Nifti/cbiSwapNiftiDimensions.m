function varargout=cbiSwapNiftiDimensions(varargin);
% [newdata,newhdr,swapMatrix]=cbiSwapNiftiDimensions(data,hdr,swapvect);
% - OR -
% swapMatrix=cbiSwapNiftiDimensions(inputfile,outputfile,swapvect,[outputdatatype]);
% 
% Swaps the dimensions of a Nifti file or data set.
% Modifies the qform44 and quaternions accordingly.
% If sform44 is non-identity, also modifies the sform.
% Equivalent to avwswapdim except it works - including modifying qform44/sform44 matrices
% NOTE: Only modifies the first 3 dimensions of a file - but works with data with 4 or more dimensions.
%
% swapvect: a 3-vector with entries +/-(1,2,3) indicating axis changes/flips, e.g.:
% - flip x (dimension 1):      swapvect=[-1 2 3]
% - exchange x and z:          swapvect=[3 2 1]
% - exchange x with flipped y: swapvect=[-2 1 3]
%
% Optionally returns swap matrix
  
% File or data input?
if (nargin<3) 
  error('not enough inputs');
end

isfile=0;
isarray=0;
if (isstr(varargin{1}) & isstr(varargin{2}))
  % Input is a file
  isfile=1;
  [data,hdr]=cbiReadNifti(varargin{1});
  datasize=size(data);
else
  % Input is a Matlab data array
  if (nargout<2) 
    error('Not enough outputs')
  end
  if (~isstruct(varargin{2}))
    error('Second input must be a header struct')
  end
  hdr=varargin{2};
  isarray=1;
  datasize=size(varargin{1});
end

swapvect=varargin{3};

if (length(swapvect)~=3)
  error('Permutation (swap) vector must have 3 dimensions!');
end

if (~find(abs(swapvect)==1) | ~find(abs(swapvect)==2) | ~find(abs(swapvect)==3 ))
  error('swapvect MUST contain 1, 2, and 3')
end

% Generate swap matrix P
% Old coordinates: X
% New coordinates after flipping: Y = P*X
% Scanner coordinates: S = Q*X where Q is the qform44
% Hence the new qform R is given by
% S = Q*X = Q*inv(P)*Y 
% R = Q*inv(P)
% and similarly, the new sform Z is given by
% Z = W*inv(P) where W is the old sform

swapMatrix=eye(4);
id44=eye(4);
for n=1:3
  % axis changes correspond to exchanging columns
  srccol=swapvect(n);
  sgn=sign(srccol);
  srccol=abs(srccol);
  swapMatrix(:,n)=sgn*id44(:,srccol);
  if (sgn<0)
    % axis inversions correspond to inverting sign and adding axis dimension-1  
    swapMatrix(srccol,4)=datasize(srccol)-1;
  end
end

% Calculate new qform & sform. Qform is needed for pixdim etc.
if (~isfield(hdr,'qform44') | isempty(hdr.qform44))
  hdr.qform44=eye(4);
end

newqform44=hdr.qform44*inv(swapMatrix);
% Set new qform and update quaternions, qoffset, and pixdim
fliphdr=cbiSetNiftiQform( hdr, newqform44 );

if (isfield(hdr,'sform44') & ~isempty(hdr.sform44))
  newsform44=hdr.sform44*inv(swapMatrix);
  fliphdr=cbiSetNiftiSform( fliphdr, newsform44 );
end

% Flip data
permutevect=1:length(datasize);
permutevect(1:3)=abs(swapvect);
if (isarray)
  flipdata=permute(varargin{1},permutevect);
else
  flipdata=permute(data,permutevect);  
end
for n=1:3
  if (swapvect(n)<0)
    flipdata=flipdim(flipdata,n);
  end
end


if (isarray)
  varargout{1}=flipdata;
  varargout{2}=fliphdr;
  if (nargout==3) 
    varargout{3}=swapMatrix;
  end
else
  % Save if desired
  if (nargin==4)
    [b,h]=cbiWriteNifti(varargin{2},flipdata,fliphdr,varargin{4});    
  else
    [b,h]=cbiWriteNifti(varargin{2},flipdata,fliphdr);
  end
  if (nargout==1) 
    varargout{1}=swapMatrix;
  end
end


