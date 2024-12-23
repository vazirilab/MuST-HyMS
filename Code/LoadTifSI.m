%%   Read data and metadata from ScanImage TIF file.
% 
% See github wiki for documentation. 
%
% Dependencies:
% ScanImage Tiff Reader package
% 
% Input
% filename: ScanImage TIF file.
%
% Output
% frames: Dataset.
% metadata: Meta data.
% descriptions: Descriptions.
%
%
% --SW, last modified: 12/14/2018.


function [frames,metadata,descriptions] = LoadTifSI(filename)

frames = ScanImageTiffReader(filename).data();
metadata = ScanImageTiffReader(filename).metadata();
descriptions = ScanImageTiffReader(filename).descriptions();

end