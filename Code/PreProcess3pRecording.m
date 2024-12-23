%%   Pre-process 3p single plane and volumetric movies.
% 
% See github wiki for documentation. 
%
% Dependencies:
% LoadTifSI.m
% SamplingCorrection.m
% GetMetaDataNumber.m
% deinterleave.m
% 
% Input
% files: List of *.tif files of 3p recordings to be pre-processed.
%
% Output
% Writes *.mat files containing pre-processed 3p datasets.
%
%
% --SW, last modified: 12/14/2018.


function [] = PreProcess3pRecording(files)

if ~exist('files','var') || isempty(files)
  % Load all files in the current directory and process
    files = dir('*.tif');
end

for file = files'
    [~,name,~] = fileparts(file.name);
    disp(['Processing ',name,'..']);

    % Load ScanImage produced TIF file:
    [frames,meta,~] = LoadTifSI(file.name);
    frames = SamplingCorrection(frames);    
    
    % Assess plane vs volumetric from meta data:
    z = GetMetaDataNumber(meta,'SI.hFastZ.numFramesPerVolume',2);
    if z>1
        type = '3p-volumetric';        
    else
        type = '3p-plane';
    end
    
    switch type
        
%% =========================================================================
        
        case '3p-plane'
            % Convert to single
            Y = single(frames);
            
            % Save data in .m format
            save([name,'.mat'],'Y','z','meta','-v7.3');
            
%% =========================================================================
            
        case '3p-volumetric'
            % Convert to single
            Y = single(frames);
            
            % Deinterleave planes and time points
            Y = deinterleave(Y,z);            

            % Remove flyback frames
            FlybackFrames = GetMetaDataNumber(meta,'SI.hFastZ.numDiscardFlybackFrames',1);
            Y = Y(:,:,:,1:z-FlybackFrames);
            
            % Save data in .m format
            z = z-FlybackFrames;
            save([name,'.mat'],'Y','z','meta','-v7.3');
    end
end

end