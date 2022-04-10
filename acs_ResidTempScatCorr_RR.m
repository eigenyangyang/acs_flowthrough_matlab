function acs_ResidTempScatCorr_RR(dir_acs,acs,TScoeff)

% This function corrects the residual temperature and scattering errors of
% ap and residual temperature errors of cp from AC-S inline system.
% Residual temperature correction of ap and cp follows Slade et al (2010),
% and scattering errors correction of ap follows Röttgers et al (2013).

% Detailed in:
% 1) Liu et al. (2018) Underway spectrophotometry in the Fram Strait
%(European Arctic Ocean): a highly resolved chlorophyll a data source for
%complementing satellite ocean color. Optics Express, 26(14), A678-A696. 
% 2) Liu et al. (2019) Retrieval of Phytoplankton Pigments from Underway
%Spectrophotometry in the Fram Strait. Remote Sensing, 11(3), 318.

% Input:
% dir_acs - the directory of AC-S data file (filetype:.mat) in string
% format.
% acs - the absolute path of temperature and salinity corrected AC-S ap and
% cp data file (filetype:.mat). This file contains a structure variable
% with 5 fields - ap, cp,t, t2 and wl which stand for particulate
% absorption and attenuation coefficients, time in numeric and vector
% formats and wavelength, respectively.
% TScoeff - the absolute path of temperature and salinity correction
% coefficients (Sullivan et al, 2006) file (filetype:.xls).

% e.g.
% dir_acs=['/Users/yliu/Data/test/99/'];
% acs=['/Users/yliu/Data/test/99/acs_p.mat'];
% TScoeff=['/Users/yliu/Data/test/Sullivan_etal_2006_instrumentspecific.xls'];

% Output: The particulate absorption coefficient corrected by residual
% temperature and scattering effect and attenuation coefficient saved as
% "acs_TSalScatCorr.mat"

% Author:
% Yangyang Liu (yangyang.liu@awi.de), August 2018.

load (acs)
VarName=strsplit(acs,'/');
VarName=char(VarName(end));
VarName=VarName(1:end-4);
VarName=eval(VarName);

ap=VarName.ap; cp=VarName.cp;
t=VarName.t; t2=VarName.t2; wl=VarName.wl;

pos_NIR=find(wl>700 & wl<750);
pos_ref=find(wl>=715, 1,'first');
pos_440=find(wl>440);
pos_675=find(wl>675);


% Temperature corrections coefficients (Sullivan et al, 2006)
TScoeff=xlsread(TScoeff);
psi_T=interp1(TScoeff(:,1),TScoeff(:,2),wl','linear','extrap');

opts=optimset('fminsearch');
opts=optimset(opts,'MaxIter',20000000);
opts=optimset(opts,'MaxFunEvals',20000);
opts=optimset(opts,'TolX',1e-8);
opts=optimset(opts,'TolFun',1e-8);

ap_TSalScatCorr=nan(size(ap));
cp_ResidTCorr=nan(size(cp));
deltaT=nan(size(ap,1),1);
fiterr=nan(size(ap,1),1);

for k=1:size(ap,1)
    
    if all(isfinite(ap(k,:)))
        % guess for paramters (beamc at lambda0, beamc slope)
        delT=0;
        
        % minimization routine (Slade et al., 2010)
        [deltaT(k), fiterr(k)]=fminsearch(@f_TS,0,opts,ap(k,:),cp(k,:),psi_T,pos_NIR,pos_ref);
        
        bp=cp(k,:)-ap(k,:);
        
        tmp=ap(k,:);
        
        ap_residT_ref=tmp(:,pos_ref)-psi_T(pos_ref).*deltaT(k);
        
        ap_ref=ap_residT_ref-0.212*sign(ap_residT_ref).*abs(ap_residT_ref).^1.135; %(Roettgers et al., 2013)
        
        ap_TSalScatCorr(k,:)=ap(k,:)-psi_T.*deltaT(k)- ...
            (ap_ref./bp(pos_ref)).*bp;
      
        
        cp_ResidTCorr(k,:)=cp(k,:)-psi_T.*deltaT(k);
        
        %disp(['  Residual T,S,Scat correction: k=' ...
        %    num2str(k) ' deltaT=' num2str(deltaT(k)) ' error=' num2str(fiterr(k))]);
        
    end
    
end

bp=cp_ResidTCorr-ap_TSalScatCorr;
for i=1:size(ap_TSalScatCorr,1)
    pos1(i,:)=any(ap_TSalScatCorr(i,1:pos_NIR(1)-1)<=0);
    pos2(i,:)=ap_TSalScatCorr(i,pos_440(1))/ap_TSalScatCorr(i,pos_675(1))<=0.95;
    pos3(i,:)=any(bp(i,1:pos_NIR(1)-1)<=0);
    pos4(i,:)=any(ap_TSalScatCorr(i,1:5)<=0.001);
end
pos=find(pos1==1 | pos2==1 | pos3==1 | pos4==1);

ap_TSalScatCorr(pos,:)=[]; cp_ResidTCorr(pos,:)=[];
t(pos,:)=[]; t2(pos,:)=[];
deltaT(pos,:)=[]; fiterr(pos,:)=[];

acs_TSalScatCorr_RR.ap=ap_TSalScatCorr;
acs_TSalScatCorr_RR.cp=cp_ResidTCorr;
acs_TSalScatCorr_RR.t=t;
acs_TSalScatCorr_RR.t2=t2;
acs_TSalScatCorr_RR.wl=wl;
acs_TSalScatCorr_RR.deltaT=deltaT;
acs_TSalScatCorr_RR.fit_error=fiterr;

namestr2=[strcat(dir_acs,'acs_TSalScatCorr_RR'),'.mat'];
save (namestr2,'acs_TSalScatCorr_RR')

return


function costf=f_TS(delT,ap,cp,psi_T,pos_NIR,pos_ref) %(Slade et al., 2010)

costf=sum( abs(ap(pos_NIR) - psi_T(pos_NIR).*delT - ((ap(pos_ref)-psi_T(pos_ref).*delT)) ));

return
