
function acs_filemerge_1minbin(dir_acs)

% This function 1) merges all the de-spiked Wetlabs AC-S data
% (raw_acs_rmspikes*.mat) into one .mat file and 2) bins the merged data
% into 1-minute interval.

% Detailed in:
% 1) Liu et al. (2018) Underway spectrophotometry in the Fram Strait
%(European Arctic Ocean): a highly resolved chlorophyll a data source for
%complementing satellite ocean color. Optics Express, 26(14), A678-A696. 
% 2) Liu et al. (2019) Retrieval of Phytoplankton Pigments from Underway
%Spectrophotometry in the Fram Strait. Remote Sensing, 11(3), 318.

% Input:
% dir_file - the directory in string format
% (e.g.dir_file=['/Users/yliu/Data/test/99/']) that contains all
% "raw_acs_rmspikes*.mat" files. Each "raw_acs_rmspikes*.mat" file contains
% a structure variable "raw_acs_rmspikes*" with 5 fields - a, c,t, t2 and
% wl which stand for absorption and attenuation coefficients, time in
% numeric and vector formats and wavelength, respectively.

% Output:
% The merged de-spiked data "merged_raw_acs_rmspikes.mat" and the 1-minute
% binned data "merged_raw_acs_rmspikes_1min.mat". They contain structure
% variables "merged_raw_acs_rmspikes" and "merged_raw_acs_rmspikes_1min"
% respectively, with 5 fields - a, c,t, t2 and wl in the above mentioned
% directory.

% Author:Yangyang Liu (yangyang.liu@awi.de), March 2018.


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% merge all raw_acs_rmspikes*.mat files
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

filenames = dir(cat(2,dir_acs, 'raw_acs_rmspikes*.mat'));
aux=numel(filenames);
filenames = struct2cell(filenames);
filenames = filenames(1,1:aux);

if aux>0
    a=[];c=[];t=[];t2=[];wl=[];
    for fn = 1:aux
        clear matFile* VarName
        matFilePath = strcat(dir_acs,char(filenames(fn)));
        matFileName = char(filenames(fn));
        
        % load raw_acs_rmspikes*.mat files.
        load(matFilePath);
        matFileName(end-3:end)=[];
        VarName=eval(matFileName);
        a0=VarName.a; c0=VarName.c;
        t0=VarName.t; t20=VarName.t2; wl=VarName.wl;
        
        tsize=size(t0);
        if tsize(1)<tsize(2)
            t0=t0';
        end
        
        wlsize=size(wl);
        if wlsize(1)<wlsize(2)
            wl=wl';
        end
        
        asize=size(a0);
        if asize(2)~=length(wl)
            a0=a';
        end
        
        csize=size(c0);
        if csize(2)~=length(wl)
            c0=c';
        end
        
        a=[a;a0]; c=[c;c0]; t=[t;t0]; t2=[t2;t20];
        clear *0*
        
    end
    % sort the data in time-ascending order.
    [t,pos]=sort(t);
    a=a(pos,:); c=c(pos,:); t2=t2(pos,:);
    
    clear merged_raw_acs_rmspikes
    merged_raw_acs_rmspikes.a=a;
    merged_raw_acs_rmspikes.c=c;
    merged_raw_acs_rmspikes.t=t;
    merged_raw_acs_rmspikes.t2=t2;
    merged_raw_acs_rmspikes.wl=wl;
    
    str1=['merged_raw_acs_rmspikes.mat'];
    save (strcat(dir_acs,(''),str1),str1(1:end-4),'-v7.3');
    
    fprintf('File merging finished!\n');
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    % bin the merged file into 1-minute interval
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    t2(:,6)=0;
    t=datenum(t2);
    [t_1min,row]=unique(t,'rows');
    t_1min(:)=t(row);
    t2_1min=t2(row,:);
    
    for i=1:length(t_1min)
        pos=find(t==t_1min(i)); 
        a_1min(i,:)=nanmedian(a(pos,:));% Take the median over one minute of measurements.
        c_1min(i,:)=nanmedian(c(pos,:));
    end
    
    merged_raw_acs_rmspikes_1min.a=a_1min;
    merged_raw_acs_rmspikes_1min.c=c_1min;
    merged_raw_acs_rmspikes_1min.t=t_1min;
    merged_raw_acs_rmspikes_1min.t2=t2_1min;
    merged_raw_acs_rmspikes_1min.wl=wl;
    
    str2=['merged_raw_acs_rmspikes_1min.mat'];
    save (strcat(dir_acs,(''),str2),str2(1:end-4));
    
    fprintf('Data 1-minute binning finished!\n');
    
else
    fprintf('Path or .mat file does not exist.\n');
end
