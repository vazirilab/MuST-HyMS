%%   Deinterleave volumetric stack.
% 
% See github wiki for documentation. 
%
% Dependencies:
% -
% 
% Input
% in: Volumetric stack.
% z: Number of axial (z) planes.
%
% Output
% out: Deinterleaved volumetric stack.
%
%
% --SW, last modified: 12/14/2018.

function out = deinterleave(in,z)

x = size(in,1);
y = size(in,2);
t = size(in,3)/z;
out = zeros(x,y,t,z,'single');

for kk = 1:z
    out(:,:,:,kk) = in(:,:,kk:z:z*(t-1)+z);
end

end