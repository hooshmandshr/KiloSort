function [ spt, cluster_id ] = kilosortDeconvolution( recording, templates, ops )
%KILOSORTDECONVOLUTION Summary of this function goes here
%   INPUT:
%       recording: C x T matrix
%       templates: input templates, should have the form NT x C x K
%       spikeTimes: 1 x N
%       cluster_id: 1 x N
%       ops: options, after loading kilosort env load deconvolution_options
%   OUTPUT:
%       spt:
%       cluster_id

C = size(recording, 1);
ops.NchanTOT = C;
ops.Nchan = C;
ops.Nfilt = size(templates, 3);

%% convert recording in mat
output_file = strcat(ops.temp_path, '\temp_bin.dat');
% fid = fopen(output_file, 'w');
% %%recording = int16([recording; zeros(2, size(recording, 2))]*(10^4));
% recording = int16(recording*(ops.scaleproc));
% fwrite(fid, recording,'*int16');
% fclose(fid);
% clear fid rec_file

%%  create a channel map file

Nchannels = C;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = ops.fs; % sampling frequency
save(strcat(ops.temp_path, '\chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
ops.chanMap = strcat(ops.temp_path, '\chanMap.mat');

clear Nchannels connected chanMap chanMap0ind xcoords ycoords kcoords fs

%%
ops.fbinary             = output_file; % will be created for 'openEphys'		
ops.fproc               = strcat(ops.temp_path, 'temp_wh.dat'); % residual from RAM of preprocessed data		
ops.root                = ops.temp_path;
clear output_file
%%
tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

Nbatch = ceil(size(recording,2) /(ops.NT-ops.ntbuff));
rez.ops = ops;
rez.temp.Nbatch = Nbatch;
rez.temp.Nbatch_buff = Nbatch;
rez.init_templates = templates;
rez.ops.kcoords = ones(ops.Nchan);
rez = setTemplates(rez);

DATA = zeros(ops.NT,ops.Nchan,Nbatch,'single');
for j = 1:Nbatch
    startpoint = max(0, 2*ops.Nchan*((ops.NT - ops.ntbuff) * (j-1) - 2*ops.ntbuff))/(2*ops.Nchan);
    if j == 1
        ioffset = 0;
    else
        ioffset = ops.ntbuff;
    end
    
    if j < Nbatch
        DATA(:,:,j) = recording(:,startpoint+ioffset+(1:ops.NT))';
    else
        DATA(1:(size(recording,2)-startpoint-ioffset),:,j) = recording(:,(startpoint+ioffset+1):end)';
    end
end

rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% remove temporary file
%delete(ops.fproc);
%%

spt = rez.st3(:, 1);
cluster_id = rez.st3(:, 2);

end


% function [ spt, cluster_id ] = kilosortDeconvolution( recording, templates, ops )
% %KILOSORTDECONVOLUTION Summary of this function goes here
% %   INPUT:
% %       recording: C x T matrix
% %       templates: input templates, should have the form NT x C x K
% %       spikeTimes: 1 x N
% %       cluster_id: 1 x N
% %       ops: options, after loading kilosort env load deconvolution_options
% %   OUTPUT:
% %       spt:
% %       cluster_id
% 
% C = size(recording, 1);
% ops.NchanTOT = C;
% ops.Nchan = C;
% ops.Nfilt = size(templates, 3);
% 
% %% convert recording in mat
% output_file = strcat(ops.temp_path, '\temp_bin.dat');
% fid = fopen(output_file, 'w');
% %%recording = int16([recording; zeros(2, size(recording, 2))]*(10^4));
% recording = int16(recording*(ops.scaleproc));
% fwrite(fid, recording,'*int16');
% fclose(fid);
% clear fid rec_file
% 
% %%  create a channel map file
% 
% Nchannels = C;
% connected = true(Nchannels, 1);
% chanMap   = 1:Nchannels;
% chanMap0ind = chanMap - 1;
% xcoords   = ones(Nchannels,1);
% ycoords   = [1:Nchannels]';
% kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
% 
% fs = ops.fs; % sampling frequency
% save(strcat(ops.temp_path, '\chanMap.mat'), ...
%     'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
% ops.chanMap = strcat(ops.temp_path, '\chanMap.mat');
% 
% clear Nchannels connected chanMap chanMap0ind xcoords ycoords kcoords fs
% 
% %%
% ops.fbinary             = output_file; % will be created for 'openEphys'		
% ops.fproc               = strcat(ops.temp_path, 'temp_wh.dat'); % residual from RAM of preprocessed data		
% ops.root                = ops.temp_path;
% clear output_file
% %%
% 
% tic; % start timer
% %
% if ops.GPU     
%     gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
% end
% 
% if strcmp(ops.datatype , 'openEphys')
%    ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
% end
% %
% [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
% rez.init_templates = templates;
% rez                = setTemplates(rez);  % fit templates iteratively
% rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
% %rez.init_clusterid = cluster_id;
% %rez.init_spiketime = spiketimes;
% 
% % remove temporary file
% %delete(ops.fproc);
% %%
% 
% spt = rez.st3(:, 1);
% cluster_id = rez.st3(:, 2);
% 
% end
% 
