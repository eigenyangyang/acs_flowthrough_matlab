function acs_calc_apcp(dir_acs, acs, minute1,minute2)

% This function calculates ap (particulate absorption coefficient) and cp
% (particulate beam attenuation coefficient) from AC-S inline system.

% Input:
% dir_acs - the directory of AC-S data file (filetype:.mat) in string
% format.
% acs - the absolute path (in string format) of temperature and salinity corrected AC-S data
% file (filetype:.mat). This file contains a structure variable with 5
% fields - a, c,t, t2 and wl which stand for absorption and attenuation
% coefficients, time in numeric and vector formats and wavelength,
% respectively.
% minute1 - the starting minute of filtering period.
% minute2 - the ending minute of filtering period.
% minute1<minute2.

% Output:
% The a and c data of particulate and dissovled materials saved as
% "acs_p.mat" and "acs_filter.mat", respectively.

% Author:
% Emmanuel Boss (emmanuel.boss@maine.edu), March 2018.
% Yangyang Liu (yangyang.liu@awi.de), March 2018.

load (acs)
VarName=strsplit(acs,'/');
VarName=char(VarName(end));
VarName=VarName(1:end-4);
VarName=eval(VarName);

a=VarName.a; c=VarName.c;
t=VarName.t; t2=VarName.t2; wl=VarName.wl;

tmp=t2;
tmp(:,[5 6]) = 0;
tmp = datenum(tmp);
unique_hr = unique(tmp);
datestr(unique_hr);

acs_filter.a_median=[];
acs_filter.t_a=[];
% acs_filter.c_median=[];
% acs_filter.t_c=[];
acs_filter.rsd_a=[]; % relative median standard deviation
% acs_filter.rsd_c=[];

pos_NIR=find(wl>700);
acs_filter.wl_for_rsd=wl(1:pos_NIR-1);
acs_filter.wl=wl;

for kh = 1:length(unique_hr)
    %datestr(unique_hr(kh),'yyyy-mm-dd HH:MM:SS');
    ix = find( t>=unique_hr(kh)+datenum(0,0,0,0,minute1,0) ...
        & t<=unique_hr(kh)+datenum(0,0,0,0,minute2,0) );
    
    if ~isempty(ix)
        rsd_a=rsd_median(a(ix,1:pos_NIR(1)-1));
        %rsd_c=rsd_median(c(ix,1:pos_NIR(1)-1));
        
        % a,c(ix,12) thresholds for:
        % PS107: row(1:end-800)->[0.6 0.6],row(end-801:end)->[0.13 0.4];
        % PS932:row(1:1300)->[0.6 0.6],row(1301:end)->[0.3 0.3];
        % PS991:[1 4]; PS992_acs213[2 7],PS992_acs219[0.25 0.4].
        if any(rsd_a>0.2)
            a(ix,:)=nan;
        end
%         if any(rsd_c>0.2)
%             c(ix,:)=nan;
%         end
        
        tmp_a = nanmedian(a(ix,:),1);
        %tmp_c = nanmedian(c(ix,:),1);
        
        if all(isfinite(tmp_a))
            
            acs_filter.t_a=[acs_filter.t_a; unique_hr(kh)+datenum(0,0,0,0,minute1,0)];
            acs_filter.t2_a=datevec(acs_filter.t_a);
            acs_filter.a_median=[acs_filter.a_median; tmp_a];
            acs_filter.rsd_a=[acs_filter.rsd_a;rsd_a];
        end
        
        %if all(isfinite(tmp_c))
            
            acs_filter.t_c=acs_filter.t_a;
           % acs_filter.t2_c=acs_filter.t2_a;
            acs_filter.c_median=acs_filter.a_median;
            acs_filter.rsd_c=acs_filter.rsd_a;
        %end
        namestr1=[strcat(dir_acs,'acs_filter'),'.mat'];
        save (namestr1,'acs_filter')
    end
end

acs_filter.a_interp=interp1(acs_filter.t_a,acs_filter.a_median,t,'linear','extrap');
acs_filter.c_interp=acs_filter.a_interp;
acs_filter.t=t;
acs_filter.t2=t2;
save (namestr1,'acs_filter')

% Calculate particle properties using linear interpolation between
% filtered measurements

acs_p.t=t;
acs_p.t2=t2;
acs_p.wl=wl;
acs_p.ap=a-interp1(acs_filter.t_a,acs_filter.a_median,t,'linear','extrap');
acs_p.cp=c-interp1(acs_filter.t_c,acs_filter.c_median,t,'linear','extrap');
acs_p.bp=acs_p.cp -acs_p.ap;

tmp_min=t2(:,5);
tmp_filter=tmp_min>=minute1 & tmp_min<=minute2;

acs_p.ap(tmp_filter,:)=nan;
acs_p.cp(tmp_filter,:)=nan;
acs_p.bp(tmp_filter,:)=nan;

namestr2=[strcat(dir_acs,'acs_p'),'.mat'];
save (namestr2,'acs_p')
