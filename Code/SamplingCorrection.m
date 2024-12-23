%%   Sampling correction for 3p datasets.
% 
% See github wiki for documentation. 
%
% Dependencies:
% -
% 
% Input
% frames: Frames of the 3p dataset (x: fast axis).
%
% Output
% out: Corrected frames.
%
%
% --SW, last modified: 12/14/2018.


function [out] = SamplingCorrection(frames)

% Data properties:
[sizex,sizey,sizet] = size(frames);

% Extract line modulus information:
pulseInfo(2,sizey,sizet) = uint8(0);
for kk = 1:sizet
pulseInfo(:,:,kk) = reshape(typecast(frames(1,:,kk),'uint8'),2,[]);
end

% Generate mask matrix:
MaskMatrix(sizex,sizey,sizet) = uint8(0);
OddLine = uint8(0:sizex-1)';
EvenLine = uint8(sizex-1:-1:0)';
for kk = 1:sizet
  MaskMatrix(:,1:2:end,kk) = mod(repmat(OddLine,[1,sizey/2]) + repmat(pulseInfo(2,1:2:end,kk),[sizex,1]),5);
  MaskMatrix(:,2:2:end,kk) = mod(repmat(EvenLine,[1,sizey/2]) + repmat(pulseInfo(2,2:2:end,kk),[sizex,1]),5);
end

% Determine correct modulus:
PixelSums = [sum(double(frames(MaskMatrix==0))), sum(double(frames(MaskMatrix==1))), ...
  sum(double(frames(MaskMatrix==2))),sum(double(frames(MaskMatrix==3))),sum(double(frames(MaskMatrix==4)))];
[~,Modulus] = sort(PixelSums);

% Mask frames:
out = frames;
out(MaskMatrix==(Modulus(1)-1)) = 0;
out(MaskMatrix==(Modulus(2)-1)) = 0;
out(MaskMatrix==(Modulus(3)-1)) = 0;
out(MaskMatrix==(Modulus(4)-1)) = 0;

% Remove line modulus information from the dataset:
out = out(2:end,:,:);

end