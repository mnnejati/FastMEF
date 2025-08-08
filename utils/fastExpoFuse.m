function FI = fastExpoFuse(Io, r)
%===============================================================
% Very Fast Exposure Fusion Method Based On Local Mean Intensity 
% and Exposedness Function
%
% Reference: 
% M. Nejati, M. Karimi, S. M. Reza Soroushmehr, N. Karimi, S. Samavi, and Kayvan Najarian, 
% "Fast Exposure Fusion Using Exposedness Function", IEEE International Conference on Image 
% Processing (ICIP), Sept. 2017, pp. 2234-2238. 
%
%
% Copyright (C) 2017, Mansour Nejati
%===============================================================
[Hei,Wid,~,N] = size(Io);

%---< Parameters >---%
alpha = 1.1; % gain of structure in the final construction
gSig = 0.2;
lSig = 0.5;
SigD = 0.12;
%---< Parameters >---%

% r = (wSize-1)/2;
epsil = 0.25;%(0.5*1)^2;
np = Hei*Wid;
H = ones(7); % small average filter is enough
H = H/sum(H(:));
L = zeros(Hei,Wid,N); 
gMu = zeros(Hei, Wid, N); % global mean intensity
lMu   = zeros(Hei, Wid, N); % local mean intensity
Iones = ones(Hei, Wid);

for i = 1:N
    
    %---- luminance component
    Ig = rgb2grey(Io(:,:,:,i)); 
    
    %--- filtered luminance component
    IgPad = padimage(Ig,[3,3]);
    L(:,:,i) = conv2(IgPad, H, 'valid');

    %--- global mean intensity
    gMu(:,:,i) = Iones * sum(Ig(:))/np;
    
    %--- local mean intensity (Base layer)
    lMu(:,:,i) = fastGF(Ig, r, epsil, 2.5);
    
end


%============< Computing Weight Maps  >================%
%--- Detail layer's blending weights
Sig2 = 2*SigD.^2;
sMap = exp(-1*(L - .5).^2 /Sig2)+1e-6; 
normalizer = sum(sMap, 3);
sMap = sMap ./ repmat(normalizer,[1, 1, N]); 
%--- Base layer's blending weights
muMap =  exp( -.5 * ( (gMu - .5).^2 /gSig.^2 +  (lMu - .5).^2 /lSig.^2 ) ); % mean intensity weighting map
normalizer = sum(muMap, 3);
muMap = muMap ./ repmat(normalizer,[1, 1, N]);


%=====================< Fusion >======================%
FI  = zeros(Hei, Wid, 3);
% sMap = alpha*sMap;
for j=1:3
    Ist = (squeeze(Io(:,:,j,:))-lMu);
    FI(:,:,j) = sum((alpha*sMap.*Ist + muMap.*lMu),3);
end
FI(FI > 1) = 1;
FI(FI < 0) = 0;
return


function I = rgb2grey(X)

origSize = size(X);

% Determine if input includes a 3-D array 
threeD = (ndims(X)==3);

% Calculate transformation matrix
coef = [0.2,0.6,0.2];

if threeD
  %RGB
  
  % Do transformation
  if isa(X, 'double') || isa(X, 'single')

    % Shape input matrix so that it is a n x 3 array and initialize output matrix  
    X = reshape(X(:),origSize(1)*origSize(2),3);
    sizeOutput = [origSize(1), origSize(2)];
    I = X * coef';
    I = min(max(I,0),1);

    %Make sure that the output matrix has the right size
    I = reshape(I,sizeOutput);
    
  else
    %uint8 or uint16
    I = imapplymatrix(coef, X, class(X));
  end

else
  I = X * coef';
  I = min(max(I,0),1);
  I = [I,I,I];
end
