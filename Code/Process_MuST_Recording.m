%%   Process 2p/3p single plane or volumetric movie.
% 
% See github wiki for documentation. 
%
% System requirements:
% Matlab Base
% Matlab Parallel Computing Toolbox
% Matlab Image Processing Toolbox
% Matlab Statistics and Machine Learning Toolbox
%
% Tested with Matlab R2017b on Red Hat Enterprise Linux 7.2
%
%
% Dependencies:
% GetMetaDataNumber.m
% NoRMCorre package
% CaImAn Matlab package
% 
% Input
% flag3p: Set "0" for 2p datasets and "1" for 3p datasets.
% MultFiles: Set to number of files the dataset is spread over.
% files: List of pre-processed *.mat files of recordings to be processed.
%
% Output
% Writes *.mat files containing processed datasets.
%
%
% --SW, last modified: 12/14/2018.


function [] = Process_MuST_Recording(flag3p,MultFiles,files)

gcp;

if flag3p==1
    disp('3-photon dataset');
else
    disp('2-photon MuST dataset');
end

if MultFiles > 1
    disp('Multiple files');
end

if ~exist('files','var') || isempty(files)
    % List of all files in the current directory:
    files = dir('*.mat');
end

% Process file / list of files:
for file = files'
    timer1 = tic;
        
    if MultFiles > 1
        % Concatenate if dataset is spread over multiple files:
        [folder_name,name,~] = fileparts(file.name);
        disp(['Processing ',name(1:end-6),'..']);
        
        % Load pre-processed *.mat file:
        load([name,'.mat']);
        Y_c = Y; clear Y;        
        for kk = 2:MultFiles
            name_next = [name(1:end-1),num2str(kk)];
            load([name_next,'.mat']);
            Y_c = cat(3,Y_c,Y);
            clear Y;
        end
        Y = Y_c;
        clear Y_c;
        name = name(1:end-6);
    else  
        [folder_name,name,~] = fileparts(file.name);
        disp(['Processing ',name,'..']);        

        % Load pre-processed *.mat file:
        load([name,'.mat']);        
    end
    
    % Remove sample timing information for 2p datasets:
    if flag3p~=1
        Y(1,:,:,:) = [];
    end
    
    % Assess plane vs volumetric from meta data:        
    if z>1
        type = 'volumetric';        
    else
        type = 'plane';
    end
    
    % Process data accordingly:
    switch type

 
% =========================================================================


%% SINGLE PLANE DATA:
        case 'plane'            
            
            disp('Identified single plane recording.');
            
            % Get frame rate of recording:
            frate = GetMetaDataNumber(meta,'SI.hRoiManager.scanFrameRate',8); % frame rate          
                        
            % Check if motion corrected file already exists:
            mc_filename = fullfile(folder_name,[name,'_mc.h5']);
            if exist(mc_filename, 'file') == 2
                disp('Motion corrected file exists. Skipping.');
            else
                
            % -------------------------------------------------------------
                            
                % Rigid motion correction using NoRMCorre algorithm:    
                options_rigid = NoRMCorreSetParms(...
                    'd1',size(Y,1),...
                    'd2',size(Y,2),...
                    'bin_width',24,...
                    'max_shift',8,...
                    'us_fac',20,...
                    'init_batch',120,...
                    'correct_bidir',false...
                    );
                [M1,shifts1,~,~] = normcorre_batch(Y,options_rigid);
                
                % Compute template from the best frames:
                shifts_r = squeeze(cat(3,shifts1(:).shifts));
                shifts_v = movvar(shifts_r,24,1);
                [~,minv_idx] = mink(shifts_v,120,1);
                best_idx = unique(reshape(minv_idx,1,[]));
                template_good = mean(M1(:,:,best_idx),3);
                
                % Piecewise non-rigid motion correction using NoRMCorre algorithm:                
                options_nonrigid = NoRMCorreSetParms(...                    
                    'd1',size(Y,1),...
                    'd2',size(Y,2),...
                    'grid_size',[32,32],...
                    'mot_uf',4,...
                    'bin_width',24,...
                    'max_shift',8,...
                    'max_dev',3,...
                    'us_fac',20,...
                    'init_batch',120,...
                    'correct_bidir',false,...
                    'output_type','h5'...
                    );
                options_nonrigid.h5_filename = fullfile(folder_name,[name,'_mc.h5']);                    
                [~,~,~,~] = normcorre_batch(Y,options_nonrigid,template_good);
                
                clear M1;                
             end
                
            % -------------------------------------------------------------
            
            % Check if CaImAn results already exist:
            cnmf_results_filename = fullfile(folder_name,[name,'_cnmf_results.mat']);
            if exist(cnmf_results_filename, 'file') == 2
                disp('CaImAn results exist. Loading results...');
                load([name,'_cnmf_results.mat']);
            else 
                % Handle for motion-corrected file:
                h5_file = subdir(fullfile(folder_name,[name,'_mc.h5'])); 
                
                % Source extraction using CaImAn on motion-corrected dataset:            
                fr = frate;
                tsub = round(fr/5);
                ds_filename = [folder_name,'ds_data.mat'];
                data_type = class(read_file(h5_file(1).name,1,1));
                data = matfile(ds_filename,'Writable',true);
                FOV = size(read_file(h5_file(1).name,1,1));
                data.Y  = zeros([FOV,0],data_type);
                data.Yr = zeros([prod(FOV),0],data_type);
                data.sizY = [FOV,0];
                F_dark = Inf;
                batch_size = 500;
                batch_size = round(batch_size/tsub)*tsub;       
                cnt = 0;
                tt1 = tic;            
                h5name = h5_file.name;
                info = h5info(h5name);
                dims = info.Datasets.Dataspace.Size;
                ndimsY = length(dims);
                Ts = dims(end);
                Ysub = zeros(FOV(1),FOV(2),floor(Ts/tsub),data_type);
                data.Y(FOV(1),FOV(2),sum(floor(Ts/tsub))) = zeros(1,data_type);
                data.Yr(prod(FOV),sum(floor(Ts/tsub))) = zeros(1,data_type);
                cnt_sub = 0;
                for t = 1:batch_size:Ts
                    Y = bigread2(h5name,t,min(batch_size,Ts-t+1));    
                    F_dark = min(nanmin(Y(:)),F_dark);
                    ln = size(Y,ndimsY);
                    Y = reshape(Y,[FOV,ln]);
                    Y = cast(downsample_data(Y,'time',tsub),data_type);
                    ln = size(Y,3);
                    Ysub(:,:,cnt_sub+1:cnt_sub+ln) = Y;
                    cnt_sub = cnt_sub + ln;
                end
                data.Y(:,:,cnt+1:cnt+cnt_sub) = Ysub;
                data.Yr(:,cnt+1:cnt+cnt_sub) = reshape(Ysub,[],cnt_sub);
                toc(tt1);
                cnt = cnt + cnt_sub;
                data.sizY(1,3) = cnt;
                data.F_dark = F_dark;

                sizY = data.sizY;
                patch_size = [60,60];
                overlap = [5,5];

                patches = construct_patches(sizY(1:end-1),patch_size,overlap);
                decay_time = 0.14;                    % GCaMP6f
                
                if flag3p==1
                    ndens = 9.2e4;                    % Number of neurons per mm^3 in mouse cortex
                    vol_fov = 0.33*0.33*0.025;        % Volume FOV in mm^3
                    ntotal = ndens*vol_fov;           % Number of expected neurons in the volume
                else
                    ndens = 9.2e4;                    % Number of neurons per mm^3 in mouse cortex
                    vol_fov = 0.5*0.5*0.025;          % Volume FOV in mm^3
                    ntotal = ndens*vol_fov;           % Number of expected neurons in the volume
                end
                npatches = length(patches);           % Number of patches
                K = round(ntotal/npatches);           % Number of expected components to be found
                disp(['Initializing CaImAn on ',num2str(npatches),' patches with K = ',num2str(K),' components...']);
                
                if flag3p==1
                    tau = [3,3];                      % STD of gaussian kernel (half size of neuron)
                else
                    tau = [2,2];                      % STD of gaussian kernel (half size of neuron)
                end                
                p = 2;                                % Order of autoregressive system
                merge_thr = 0.8;                      % Merging threshold
                sizY = data.sizY;

                options = CNMFSetParms(...
                    'd1',sizY(1),'d2',sizY(2),...
                    'deconv_method','constrained_foopsi',...    % Neural activity deconvolution method
                    'search_method','ellipse',...               % Method for determining footprint of spatial components
                    'min_size',ceil(tau(1)),...                 % Minimum size of ellipse axis
                    'max_size',ceil(3*tau(1)),...               % Maximum size of ellipse axis
                    'dist',2,...                                % Expansion factor of ellipse
                    'p',p,...                                   % Order of AR dynamics
                    'temporal_iter',3,...                       % Number of block-coordinate descent steps 
                    'maxIter',15,...                            % Number of NMF iterations during initialization
                    'ssub',1,...                                % Spatial downsampling when processing
                    'tsub',1,...                                % Further temporal downsampling when processing
                    'merge_thr',merge_thr,...                   % Merging threshold
                    'gSig',tau,...                              % STD of gaussian kernel (half size of neuron) 
                    'max_size_thr',ceil(pi*(2*tau(1)-1)^2),...  % Maximum acceptable size for each component
                    'min_size_thr',ceil(pi*(tau(1)-1)^2),...    % Minimum acceptable size for each component
                    'spatial_method','regularized',...          % Method for updating spatial components
                    'df_prctile',20,...                         % Take the median of background fluorescence to compute baseline fluorescence 
                    'fr',fr/tsub,...                            % Downsamples
                    'space_thresh',0.1,...                      % Spatial correlation acceptance threshold
                    'min_SNR',1.5,...                           % Trace SNR acceptance threshold                    
                    'nb',1,...                                  % Number of background components per patch
                    'gnb',3,...                                 % Number of global background components
                    'decay_time',decay_time...                         
                    );

                % Run CNMF on patches:
                [A,b,C,f,~,P,~,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

                % Classify components:
                rval_space = classify_comp_corr(data,A,C,b,f,options);
                ind_corr = rval_space > options.space_thresh;

                % Further classification with cnn_classifier
                try  % Matlab 2017b or later is needed
                    [ind_cnn,~] = cnn_classifier(A,FOV,'cnn_model',options.cnn_thr);
                catch
                    ind_cnn = true(size(A,2),1); % Components that pass the CNN classifier
                end     

                % Event exceptionality:
                fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
                ind_exc = (fitness < options.min_fitness);

                % Select components:
                keep = (ind_corr | ind_cnn) & ind_exc;

                % Keep only the active components:  
                A_keep = A(:,keep);
                C_keep = C(keep,:);                

                % Extract fluorescence on native temporal resolution:
                options.fr = options.fr*tsub;                   % Revert to original frame rate
                N = size(C_keep,1);                             % Total number of components
                T = sum(Ts);                                    % Total number of timesteps
                C_full = imresize(C_keep,[N,T]);                % Upsample to original frame rate
                f_full = imresize(f,[size(f,1),T]);             % Upsample temporal background
                S_full = zeros(N,T);
                P.p = 0;
                options.nb = options.gnb;            
                [C_full,f_full,~,~,R_full] = update_temporal_components_fast(h5_file.name,A_keep,b,C_full,f_full,P,options);
                disp('Extracting raw fluorescence at native frame rate.');

                % Order components:
                [A_or,C_or,~,~] = order_ROIs(A_keep,C_full,S_full,P);

                % Calculate dF/F traces:
                [F_dff,F0] = detrend_df_f(A_or,[b,ones(prod(FOV),1)],C_or,[f_full;-double(F_dark)*ones(1,T)],R_full,options);
                if flag3p==1
                    w_mean = 3;                    
                    C_df = movmean(F_dff,w_mean,2,'omitnan');
                else
                    C_df = full(F_dff);
                end 
                
                % Trace deconvolution: 
                C_dec = zeros(N,T);         % Deconvolved dF/F traces
                S_dec = zeros(N,T);         % Deconvolved neural activity
                bl = zeros(N,1);            % Baseline for each trace (should be close to zero since traces are dF/F)
                neuron_sn = zeros(N,1);     % Noise level at each trace
                g = cell(N,1);              % Discrete time constants for each trace
                if p == 1; model_ar = 'ar1'; elseif p == 2; model_ar = 'ar2'; else; error('This order of dynamics is not supported'); end

                disp('Performing deconvolution.');
                for i = 1:N
                    spkmin = options.spk_SNR*GetSn(F_dff(i,:));
                    lam = choose_lambda(exp(-1/(options.fr*options.decay_time)),GetSn(F_dff(i,:)),options.lam_pr);
                    [cc,spk,opts_oasis] = deconvolveCa(F_dff(i,:),model_ar,'method','thresholded','optimize_pars',true,'maxIter',20,...
                        'window',150,'lambda',lam,'smin',spkmin);
                    bl(i) = opts_oasis.b;
                    C_dec(i,:) = cc(:)' + bl(i);
                    S_dec(i,:) = spk(:);
                    neuron_sn(i) = opts_oasis.sn;
                    g{i} = opts_oasis.pars(:)';                    
                end
                
                % Center of mass positions of spatial components:
                sizY = FOV;
                center = com(A_or,sizY(1),sizY(2));
                
                % Save CaImAn results to Matlab:
                C_df_dec = full(C_dec);
                time = linspace(0,size(C_df,1)/frate,size(C_df,1));
                try
                    save([name,'_cnmf_results.mat'],'C_df','F0','C_df_dec','S_dec','A_or',...
                        'center','neuron_sn','bl','g','time','frate','options','-v7.3');
                catch
                    disp('Writing to Matlab failed.');
                end
                
            end 
            
            
% =========================================================================


%% VOLUMETRIC DATA:
        case 'volumetric'
            
            disp('Identified volumetric recording.');
            
            % Get frame rate and number of channels of recording:
            frate = GetMetaDataNumber(meta,'SI.hRoiManager.scanVolumeRate',8);
            
            % -------------------------------------------------------------            
            
            % Check if motion-correction file already exists:
            mc_filename = fullfile(folder_name,[name,'_mc.h5']);
            if exist(mc_filename, 'file') == 2
                disp('Motion corrected file exists. Skipping...');
            else
                Y = permute(Y,[1,2,4,3]);
                disp(['Size of stack: ',num2str(size(Y))]);
            
            % -------------------------------------------------------------
                
                % Rigid motion correction using NoRMCorre algorithm:
                options_rigid = NoRMCorreSetParms(...
                    'd1',size(Y,1),...
                    'd2',size(Y,2),...
                    'd3',size(Y,3),...
                    'bin_width',24,...
                    'max_shift',[8,8,2],...
                    'us_fac',20,...
                    'init_batch',120,...
                    'correct_bidir',false...
                    );
                [M1,shifts1,~,~] = normcorre_batch(Y,options_rigid);

                % Compute template from the lowest motion frames:
                shifts_r = squeeze(cat(3,shifts1(:).shifts));
                shifts_v = movvar(shifts_r,24,1);
                [~,minv_idx] = mink(shifts_v,120,1);
                best_idx = unique(reshape(minv_idx,1,[]));
                template_good = mean(M1(:,:,:,best_idx),4);

                % Non-rigid motion correction using NoRMCorre algorithm:    
                options_nr = NoRMCorreSetParms(...
                    'd1',size(Y,1),...
                    'd2',size(Y,2),...
                    'd3',size(Y,3),...
                    'bin_width',24,...
                    'max_shift',[8,8,2],...
                    'us_fac',20,...
                    'init_batch',120,...
                    'correct_bidir',false,...
                    'mem_batch_size',100,...
                    'output_type','h5'...
                    );
                options_nr.h5_filename = fullfile(folder_name,[name,'_mc.h5']);
                [~,~,~,~] = normcorre_batch(Y,options_nr,template_good);
                
                clear M1;
            end
                        
            % -------------------------------------------------------------
            
            % Check if CaImAn results already exist:
            cnmf_results_filename = fullfile(folder_name,[name,'_cnmf_results.mat']);
            if exist(cnmf_results_filename, 'file') == 2
                disp('CNMF results exist. Loading results...');
                load([name,'_cnmf_results.mat']);
            else            
                % Handle for motion-corrected file:
                h5_file = subdir(fullfile(folder_name,[name,'_mc.h5']));            

                % Source extraction using CaImAn on motion-corrected dataset:  
                disp('Memory mapping motion corrected dataset...');
                fr = frate; % volume rate              
                if fr > 5
                    tsub = round(fr/5);
                else
                    tsub = 1;
                end                
                ds_filename = [folder_name,'ds_data.mat'];
                data_type = class(read_file(h5_file(1).name,1,1));
                data = matfile(ds_filename,'Writable',true);
                FOV = size(read_file(h5_file(1).name,1,1));
                data.Y  = zeros([FOV,0],data_type);
                data.Yr = zeros([prod(FOV),0],data_type);
                data.sizY = [FOV,0];
                F_dark = Inf;
                batch_size = 500;
                batch_size = round(batch_size/tsub)*tsub;         
                cnt = 0;
                tt1 = tic;            
                h5name = h5_file.name;
                info = h5info(h5name);
                dims = info.Datasets.Dataspace.Size;
                ndimsY = length(dims);
                Ts = dims(end);
                Ysub = zeros(FOV(1),FOV(2),FOV(3),floor(Ts/tsub),data_type);
                data.Y(FOV(1),FOV(2),FOV(3),sum(floor(Ts/tsub))) = zeros(1,data_type);
                data.Yr(prod(FOV),sum(floor(Ts/tsub))) = zeros(1,data_type);
                cnt_sub = 0;
                for t = 1:batch_size:Ts
                    Y = bigread2(h5name,t,min(batch_size,Ts-t+1));    
                    F_dark = min(nanmin(Y(:)),F_dark);
                    ln = size(Y,ndimsY);
                    Y = reshape(Y,[FOV,ln]);
                    Y = cast(downsample_data(Y,'time',tsub),data_type);
                    ln = size(Y,4);
                    Ysub(:,:,:,cnt_sub+1:cnt_sub+ln) = Y;
                    cnt_sub = cnt_sub + ln;
                end
                data.Y(:,:,:,cnt+1:cnt+cnt_sub) = Ysub;
                data.Yr(:,cnt+1:cnt+cnt_sub) = reshape(Ysub,[],cnt_sub);
                disp('Done.');
                toc(tt1);                
                cnt = cnt + cnt_sub;
                data.sizY(1,4) = cnt;
                data.F_dark = F_dark;

                sizY = data.sizY;
                if flagLat==1
                    patch_size = [ceil(sizY(1)/6),...
                        ceil(sizY(2)/6),50];
                else                    
                    patch_size = [ceil(sizY(1)/4),...
                        ceil(sizY(2)/3),50];
                end
                overlap = [4,4,0];
                
                patches = construct_patches(sizY(1:end-1),patch_size,overlap);
                decay_time = 0.14;                       % GCaMP6f                
                
                if flag3p==1
                    ndens = 9.2e4;                      % Number of neurons per mm^3 in mouse cortex
                    vol_fov = 0.350*0.350*0.250;        % Volume FOV in mm^3
                    ntotal = ndens*vol_fov;             % Number of neurons in the volume
                else
                    ndens = 9.2e4;                      % Number of neurons per mm^3 in mouse cortex
                    vol_fov = 0.500*0.500*0.600;        % Volume FOV in mm^3
                    ntotal = ndens*vol_fov;             % Number of neurons in the volume
                end
                npatches = length(patches);             % Number of patches
                K = round(ntotal/npatches);             % Number of components to be found
                disp(['Initializing CNMF on ',num2str(npatches),' patches with K = ',num2str(K),' components...']);
                
                if flag3p==1
                    tau = [3,3,1];                      % STD of gaussian kernel (half size of neuron) 
                    p = 1;                              % Order of autoregressive system
                else
                    tau = [2,2,2];                      % STD of gaussian kernel (half size of neuron)
                    p = 2;                              % Order of autoregressive system
                end                
                merge_thr = 0.8;                        % Merging threshold
                sizY = data.sizY;

                options = CNMFSetParms(...
                    'd1',sizY(1),'d2',sizY(2),'d3',sizY(3),...                    
                    'deconv_method','constrained_foopsi',...        % Neural activity deconvolution method
                    'p',p,...                                       % Order of AR dynamics
                    'search_method','ellipse',...                   % Method for determining footprint of spatial components
                    'min_size',ceil(0.5*tau(1)),...                 % Minimum size of ellipse axis
                    'max_size',ceil(3.0*tau(1)),...                 % Maximum size of ellipse axis
                    'dist',2,...                                    % Expansion factor of ellipse
                    'temporal_iter',3,...                           % Number of block-coordinate descent steps 
                    'maxIter',15,...                                % Number of NMF iterations during initialization
                    'ssub',1,...                                    % Spatial downsampling when processing
                    'tsub',1,...                                    % Further temporal downsampling when processing
                    'merge_thr',merge_thr,...                       % Merging threshold                    
                    'gSig',tau,...                                  % STD of gaussian kernel (half size of neuron)
                    'min_size_thr',0.5*ceil(pi*(tau(1)-1)^2),...    % Minimum acceptable size for each component
                    'max_size_thr',3.0*ceil(pi*(2*tau(1)-1)^2),...  % Maximum acceptable size for each component
                    'spatial_method','regularized',...              % Method for updating spatial components
                    'df_prctile',10,...                             % Take the median of background fluorescence to compute baseline fluorescence 
                    'cl_thr',0.8,...                                % Overlap threshold for energy for a component to be classified as true
                    'fr',fr/tsub,...                                % Downsamples
                    'space_thresh',0.10,...                         % Spatial correlation acceptance threshold
                    'min_SNR',1.5,...                               % Trace SNR acceptance threshold
                    'nb',1,...                                      % Number of background components per patch
                    'gnb',3,...                                     % Number of global background components
                    'decay_time',decay_time...                      % Length of typical transient for the indicator used
                    );

                % Run CNMF on patches:
                [A,b,C,f,~,P,~,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

                % Classify components:
                rval_space = classify_comp_corr(data,A,C,b,f,options);
                ind_corr = rval_space > options.space_thresh;                     

                % Event exceptionality:
                fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
                ind_exc = (fitness < options.min_fitness);

                % Select components:
                keep = ind_corr & ind_exc;

                % Keep only the active components:  
                A_keep = A(:,keep);
                C_keep = C(keep,:);

                % Extract fluorescence on native temporal resolution:
                options.fr = options.fr*tsub;                   % Revert to original frame rate
                N = size(C_keep,1);                             % Total number of components
                T = sum(Ts);                                    % Total number of timesteps
                C_full = imresize(C_keep,[N,T]);                % Upsample to original frame rate                
                f_full = imresize(f,[size(f,1),T]);             % Upsample temporal background
                S_full = zeros(N,T);
                P.p = 0;
                options.nb = options.gnb;
                [C_full,f_full,~,~,R_full] = update_temporal_components_fast(h5_file.name,A_keep,b,C_full,f_full,P,options);
                disp('Extracting raw fluorescence at native frame rate.');               

                % Order components:
                [A_or,C_or,~,~] = order_ROIs(A_keep,C_full,S_full,P);

                % Calculate dF/F traces:
                [F_dff,F0] = detrend_df_f(A_or,[b,ones(prod(FOV),1)],C_or,[f_full;-double(F_dark)*ones(1,T)],R_full,options);
                if flag3p==1
                    w_mean = 3;                    
                    C_df = movmean(F_dff,w_mean,2,'omitnan');
                else
                    C_df = full(F_dff);
                end 

                % Trace deconvolution: 
                disp('Performing deconvolution.');
                C_dec = zeros(N,T);         % Deconvolved dF/F traces
                S_dec = zeros(N,T);         % Deconvolved neural activity
                bl = zeros(N,1);            % Baseline for each trace (should be close to zero since traces are DF/F)
                neuron_sn = zeros(N,1);     % Noise level at each trace
                g = cell(N,1);              % Discrete time constants for each trace
                if p == 1; model_ar = 'ar1'; elseif p == 2; model_ar = 'ar2'; else; error('This order of dynamics is not supported'); end

                for i = 1:N
                    spkmin = options.spk_SNR*GetSn(F_dff(i,:));
                    lam = choose_lambda(exp(-1/(options.fr*options.decay_time)),GetSn(F_dff(i,:)),options.lam_pr);
                    [cc,spk,opts_oasis] = deconvolveCa(F_dff(i,:),model_ar,'method','thresholded','optimize_pars',true,'maxIter',20,...
                                                'window',150,'lambda',lam,'smin',spkmin);
                    bl(i) = opts_oasis.b;
                    C_dec(i,:) = cc(:)' + bl(i);
                    S_dec(i,:) = spk(:);
                    neuron_sn(i) = opts_oasis.sn;
                    g{i} = opts_oasis.pars(:)';                    
                end
                
                % Center of mass positions of spatial components:                
                sizY = FOV;
                center = com(A_or,sizY(1),sizY(2),sizY(3));                

                % Save CNMF results to Matlab:
                C_df_dec = full(C_dec);
                time = linspace(0,size(C_df,1)/frate,size(C_df,1));
                try
                    save([name,'_cnmf_results.mat'],'C_df','F0','C_df_dec','S_dec','A_or',...
                        'center','neuron_sn','bl','g','time','frate','options','-v7.3');
                catch
                    disp('Writing to Matlab failed.');
                end
            
            end           
                     
    end
    
    tEnd = toc(timer1);
    disp(['Total runtime: ',num2str(floor(tEnd/60)),' minutes and ',num2str(rem(tEnd,60)),' seconds.']);
    
end

end