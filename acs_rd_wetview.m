function acs_rd_wetview(dir_acs)

% This function extracts Wetlabs AC-S data from the Wetview software output
% (.dat files).

% Detailed in:
% 1) Liu et al. (2018) Underway spectrophotometry in the Fram Strait
%(European Arctic Ocean): a highly resolved chlorophyll a data source for
%complementing satellite ocean color. Optics Express, 26(14), A678-A696. 
% 2) Liu et al. (2019) Retrieval of Phytoplankton Pigments from Underway
%Spectrophotometry in the Fram Strait. Remote Sensing, 11(3), 318.

% Before using this function, modify manually the AC-S .dat files names as:
% acs_01_ ... acs_20_... (i.e. just write a 0 before the numbers from 1 in
% order to read the files in a sorted way).

% Input: dir_acs - the directory in string format
% (e.g.dir_acs=['/Users/yliu/Data/test/99/']) that contains all the AC-S
% "acs*.dat" files.

% Output: The extracted AC-S data which are saved as raw_acs_data*.mat
% files that contain structure variables "raw_acs_data*" with 5 fields - a,
% c,t, t2 and wl that stand for absorption and attenuation coefficients,
% time in numeric and vector formats and wavelength, respectively, in the
% above mentioned directory.

% Author:Yangyang Liu (yangyang.liu@awi.de), March 2018.


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% get wavelenghts from the first .dat file
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% get the names of all acs*.dat files in the designated directory.
filenames = dir(cat(2,dir_acs, 'acs*.dat'));
aux=numel(filenames);
filenames = struct2cell(filenames);
filenames = filenames(1,1:aux);

fid=fopen(filenames{1},'r');

%skip the first 8 lines and get into the 9th line "83; output wavelengths."
%This number could vary. Check it!
for i=1:9
    line=fgetl(fid);
end

n_wl=str2num(line(1:2)); % get number of wavelengths

%skip the next 2 lines
for i=1:2
    line=fgetl(fid);
end

wla=[];wlc=[];
for i=1:n_wl
    line=fgetl(fid);
    wla0=str2num(line(9:13));
    wlc0=str2num(line(2:6));
    wla=[wla;wla0];
    wlc=[wlc;wlc0];
end
fclose(fid);
clear fid i


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% get data from the .dat files
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

for n=1:aux
    
    fid=fopen(filenames{n},'r');
    line1=fgetl(fid);
    datetime=datenum(line1(17:35),'mm-dd-yyyy HH:MM:SS');
    
    data=dlmread(filenames{n},'\t',97,0); % skip the first 97 lines. This number could vary. Check!
    
    t=data(:,1)/1000/86400+datetime;
    c=data(:,2:n_wl+1);
    a=data(:,n_wl+2:2*n_wl+1);
    c=interp1(wlc,c',wla,'linear')';
    t2=datevec(t);
    
    clear raw_acs_data
    raw_acs_data.a=a;
    raw_acs_data.c=c;
    raw_acs_data.t=t;
    raw_acs_data.t2=t2;
    raw_acs_data.wl=wla;
    
    eval(['raw_acs_data',num2str(n),'=','raw_acs_data',';']);
    str=['raw_acs_data' num2str(n) '.mat'];
    save (strcat(dir_acs,(''),str),str(1:end-4));
    fclose(fid);
    clear a c t* data
    
end
