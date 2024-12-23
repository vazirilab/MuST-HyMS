%%   Pre-process 2p single plane and volumetric movies.
% 
% See github wiki for documentation. 
%
% Dependencies:
% LoadTifSI.m
% GetMetaDataNumber.m
% deinterleave.m
% deinterleave2.m
% 
% Input
% files: List of *.tif files of 2p recordings to be pre-processed.
%
% Output
% Writes *.mat files containing pre-processed 2p datasets.
%
%
% --SW, last modified: 12/14/2018.


function [] = PreProcess2pRecording(files)

if ~exist('files','var') || isempty(files)
    % Load all files in the current directory and process
    files = dir('*.tif');
end

for file = files'
    [~,name,~] = fileparts(file.name);
    disp(['Processing ',name,'..']);

    % Load ScanImage produced TIF file:
    [Y,meta,~] = LoadTifSI(file.name);
    
    % Assess plane vs volumetric from meta data:
    z = GetMetaDataNumber(meta,'SI.hFastZ.numFramesPerVolume',2);
    if z>1
        type = '2p-volumetric';        
    else
        type = '2p-plane';
    end
    
    % Assess number of channels (for multiplexing) from meta data:
    channels = numel(GetMetaDataNumber(meta,'SI.hChannels.channelSave',10));
    
    switch type
        
%% =========================================================================
        
        case '2p-plane'            
            % Convert to single
            Y = single(Y);
            
            % Deinterleave multi-plane recording
            planes = channels;
            Y = deinterleave(Y,planes);
                        
            % Save data in .m format
            save([name,'.mat'],'Y','z','meta','planes','-v7.3');

%% =========================================================================
            
        case '2p-volumetric'
            % Convert to single
            Y = single(Y);
            
            % Deinterleave planes and time points, stitch volumes
            Y = deinterleave(Y,z*channels);

            % Remove flyback frames
            FlybackFrames = GetMetaDataNumber(meta,'SI.hFastZ.numDiscardFlybackFrames',1);
            keep = 1:size(Y,4) - (channels*FlybackFrames);
            Y = Y(:,:,:,keep);

            % Deinterleave channels (for multiplexing) and planes
            Y = deinterleave2(Y,channels,z-FlybackFrames);

            % Save data in .m format            
            z = (z-FlybackFrames)*channels;            
            save([name,'.mat'],'Y','z','meta','-v7.3');
            clear Y;            
    end
end

end