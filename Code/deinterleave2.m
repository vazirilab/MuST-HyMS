%%   Deinterleave volumetric data channels.
% 
% See github wiki for documentation. 
%
% Dependencies:
% -
% 
% Input
% in: Volumetric data.
% ch: Number of multiplexed channels.
% z: Number of axial (z) planes.
%
% Output
% out: Deinterleaved volumetric data.
%
%
% --SW, last modified: 12/14/2018.

function out = deinterleave2(in,ch,z)

x = size(in,1);
y = size(in,2);
t = size(in,3);
out = zeros(x,y,t,z*ch,'single');

for kk = 1:ch
    tmp = in(:,:,:,kk:ch:ch*(z-1)+ch);
    out(:,:,:,(1:z)+(z*(kk-1))) = squeeze(tmp);
end

end