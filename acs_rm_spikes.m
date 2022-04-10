
function acs_rm_spikes(dir_acs)

% This function removes the spikes of Wetlabs AC-S data. The spikes are
% defined either as the relative median standard deviation of the data in
% visible band measured within 2s are greater than a threshold (in this
% function the threshold is 2), or the magnitude at ~440nm of each spectrum
% is greater than a threshold (in this function the threshhold is 5). Data
% that are before and after 10s with respect to the spikes are also removed.

% Detailed in:
% 1) Liu et al. (2018) Underway spectrophotometry in the Fram Strait
%(European Arctic Ocean): a highly resolved chlorophyll a data source for
%complementing satellite ocean color. Optics Express, 26(14), A678-A696. 
% 2) Liu et al. (2019) Retrieval of Phytoplankton Pigments from Underway
%Spectrophotometry in the Fram Strait. Remote Sensing, 11(3), 318.

% Input: 
% dir_acs - the directory in string format
% (e.g.dir_file=['/Users/yliu/Data/test/99/']) that contains all
% "raw_acs_data*.mat" files. Each raw_acs_data*.mat file contains a
% structure variable "raw_acs_data*" with 5 fields - a, c,t, t2 and wl
% which stand for absorption and attenuation coefficients, time in numeric
% and vector formats and wavelength, respectively.

% Output: 
% The de-spiked data which are saved as raw_acs_rmspikes*.mat files that
% contain structure variables "raw_acs_rmspikes*" with 5 fields - a, c,t, t2
% and wl in the above mentioned directory.

% Author:Yangyang Liu (yangyang.liu@awi.de), March 2018.


k=4; % 4 measurements per second.
k_s=10*4; % 10 seconds measurements.

% get the names of all raw_acs_data*.mat files in the designated directory.
filenames = dir(cat(2,dir_acs, 'raw_acs_data*.mat'));
aux=numel(filenames);
filenames = struct2cell(filenames);
filenames = filenames(1,1:aux);

if aux>0
    a=[];c=[];t=[];t2=[];wl=[];
    for fn = 1:aux
        
        clear matFile* VarName
        matFilePath = strcat(dir_acs,char(filenames(fn)));
        matFileName = char(filenames(fn));
        
        % load raw_acs_data*.mat files.
        load(matFilePath);
        matFileName(end-3:end)=[];
        VarName=eval(matFileName);
        a=VarName.a; c=VarName.c;
        t=VarName.t; t2=VarName.t2; wl=VarName.wl;
        
        tsize=size(t);
        if tsize(1)<tsize(2)
            t=t';
        end
        
        wlsize=size(wl);
        if wlsize(1)<wlsize(2)
            wl=wl';
        end
        
        asize=size(a);
        if asize(2)~=length(wl)
            a=a';
        end
        
        csize=size(c);
        if csize(2)~=length(wl)
            c=c';
        end
        
        pos_NIR=find(wl>700);
        n=length(t); 
        for i = 1:n
            if i-k>=1 && i+k<=n
                nn=i-k:i+k;
            elseif i-k<1
                nn=1:i+k;
            else
                nn=i-k:n;
            end
            a_rsd=rsd_median(a(nn,1:pos_NIR(1)-1));
            c_rsd=rsd_median(c(nn,1:pos_NIR(1)-1));
            
            % a,c(nn,12) thresholds for: PS107[0.7 1.5];PS932[0.5
            % 2];PS991[2 6];PS992_acs213[2 7];PS992_acs219[0.8 1.8];
            if any(a_rsd>0.5)|| any(c_rsd>0.5)|| any(a(nn,12)>0.7) || any(c(nn,12)>1.5) %thresholds
                a(nn,:)=nan;
                c(nn,:)=nan;
                if i-k_s+1>=1 && i+k_s<=n
                    a(i-k_s+1:i+k_s,:)=nan;
                    c(i-k_s+1:i+k_s,:)=nan;
                end
            end            
        end
        
        clear raw_acs_rmspikes
        raw_acs_rmspikes.a=a;
        raw_acs_rmspikes.c=c;
        raw_acs_rmspikes.t=t;
        raw_acs_rmspikes.t2=t2;
        raw_acs_rmspikes.wl=wl;
        
        eval(['raw_acs_rmspikes',num2str(matFileName(13:end)),'=','raw_acs_rmspikes',';']);
        str=['raw_acs_rmspikes' num2str(matFileName(13:end)) '.mat'];
        save (strcat(dir_acs,(''),str),str(1:end-4));
        
        fprintf('%s Spikes removal finished!\n', matFileName)
    end
else
    fprintf('Path or .mat file does not exist.\n');
end
