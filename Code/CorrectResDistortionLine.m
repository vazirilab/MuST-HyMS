%%   Correct resonant scanner distortion.
% 
% See github wiki for documentation. 
%
% Dependencies:
% -
% 
% Input
% frame: Input frame (fast axis: x).
% FillFracSpatial: Spatial fill fraction of the recording.
% SamplingRate: Sampling rate of the recording.
% ResFreq: Resonant scanner scan frequency.
% 
%
% Output
% frame_corr: Corrected frame.
%
%
% --SW, last modified: 12/14/2018.


function [frame_corr] = CorrectResDistortionLine(frame,FillFracSpatial,SamplingRate,ResFreq)

% Extract dimensions
x = size(frame,1);
y = size(frame,2);

% Calculate correction
pixelBoundaries = linspace(-FillFracSpatial,FillFracSpatial,x+1)';
pixelBoundariesTime = sin(pixelBoundaries)/(2*pi*ResFreq);
pixelBoundariesSamples = pixelBoundariesTime*SamplingRate;
mask = round(diff(pixelBoundariesSamples)*6.25);
x2 = [];
for kk = 1:length(mask)
x2 = [x2 kk*ones(1,mask(kk))];
end

% Resampling grid
X2 = 1:length(x2);
X = 1:x;
Xp = linspace(min(x2),length(x2),x);

% Resample data
frame_corr = zeros(x,y);
for kk = 1:y
tmp = interp1(X,frame(:,kk),x2,'linear')';
frame_corr(:,kk) = interp1(X2,tmp,Xp,'linear');
end

end