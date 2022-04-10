function acs_TempSalCorr(dir_acs,acs,dev,TSG,TScoeff)

% This function does the temperature and salinity correction of Wetlabs
% AC-S data. The temperature and salinity correction coefficients are from
% Sullivan et al, 2006.

% Detailed in:
% 1) Liu et al. (2018) Underway spectrophotometry in the Fram Strait
%(European Arctic Ocean): a highly resolved chlorophyll a data source for
%complementing satellite ocean color. Optics Express, 26(14), A678-A696. 
% 2) Liu et al. (2019) Retrieval of Phytoplankton Pigments from Underway
%Spectrophotometry in the Fram Strait. Remote Sensing, 11(3), 318.

% Input: 
% dir_acs - the directory of AC-S data file (filetype:.mat) in string
% format.
% acs - the absolute path of AC-S data file (filetype:.mat).This file
% contains a structure variable with 5 fields - a, c,t, t2 and wl which
% stand for absorption and attenuation coefficients, time in numeric and
% vector formats and wavelength, respectively.
% dev - the absolute path of AC-S device file (filetype:.dev).
% TSG - the absolute path of temperature and salinity data file
% (filetype:.mat). This file contains a structure variable with 5 fields -
% a, c,t, t2 and wl which stand for absorption and attenuation
% coefficients, time in numeric and vector formats and wavelength,
% respectively.
% TScoeff - the absolute path of temperature and salinity correction
% coefficients (Sullivan et al, 2006) file (filetype:.xls).

% e.g. 
% dir_acs=['/Users/yliu/Data/test/99/'];
% acs=['/Users/yliu/Data/test/99/merged_raw_acs_rmspikes_1min.mat'];
% TSG=['/Users/yliu/Data/test/TSG_keel_PS99.mat'];
% dev=['/Users/yliu/Data/test/acs219.dev'];
% TScoeff=['/Users/yliu/Data/test/Sullivan_etal_2006_instrumentspecific.xls'];

% Output: 
% The temperature and salinity corrected AC-S data which is saved as
% "acs_tscorr.mat".

% Author:Yangyang Liu (yangyang.liu@awi.de), March 2018.

load (acs)
VarName1=strsplit(acs,'/');
VarName1=char(VarName1(end));
VarName1=VarName1(1:end-4);
VarName1=eval(VarName1);

a=VarName1.a; c=VarName1.c;
t_acs=VarName1.t; t2_acs=VarName1.t2; wl=VarName1.wl;


load(TSG)
VarName2=strsplit(TSG,'/');
VarName2=char(VarName2(end));
VarName2=VarName2(1:end-4);
VarName2=eval(VarName2);

temp=VarName2.temp; sal=VarName2.sal;
t_TSG=VarName2.t; t2_TSG=VarName2.t2;

[t_overlap,pos_acs,pos_TSG]=intersect(t_acs,t_TSG);
temp_overlap=temp(pos_TSG);
sal_overlap=sal(pos_TSG);

% read temparature value during factory calibration (tcal) from device file.
% (salinity value during factory calibration is 0, i.e. scal=0)

dev=fopen(dev,'r');
%skip first 3 lines
for i=1:4
    tcal=fgetl(dev);
end
tcal=str2num(tcal(7:11)); %Line 4 .dev file.

% T-S corrections (Sullivan et al, 2006)
TScoeff = xlsread(TScoeff); % Correction coefficients

psi_T=interp1(TScoeff(:,1),TScoeff(:,2),wl','linear','extrap');
psi_Sc=interp1(TScoeff(:,1),TScoeff(:,4),wl','linear','extrap');
psi_Sa=interp1(TScoeff(:,1),TScoeff(:,6),wl','linear','extrap');

for j=1:size(a,1)
    a_tscorr(j,:)=a(j,:)-(temp_overlap(j)-tcal)*psi_T-sal_overlap(j)*psi_Sa;
    c_tscorr(j,:)=c(j,:)-(temp_overlap(j)-tcal)*psi_T-sal_overlap(j)*psi_Sc;
end

acs_tscorr.a=a_tscorr;
acs_tscorr.c=c_tscorr;
acs_tscorr.wl=wl;
acs_tscorr.t=t_acs;
acs_tscorr.t2=t2_acs;

namestr=[strcat(dir_acs,'acs_tscorr'),'.mat'];
save (namestr,'acs_tscorr')


