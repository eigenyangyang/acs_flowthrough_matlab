function rsd=rsd_median(a)
% it returns the relative median standard deviation (i.e. median coefficient of variation)

% Input:
% x: a matrix

% Output:
% rsd: the relative median standard deviation for several rows. dim=1.

% Author:Yangyang Liu (yangyang.liu@awi.de), March 2018.

s=size(a);
row=s(1,1);
a_median=nanmedian(a,1);
abs_median=abs(a_median);
a_median_match=repmat(a_median,row,1);
a_diff_sq=(a-a_median_match).^2;

if length(a_diff_sq(~isnan(a_diff_sq)))==0
    summa=nan;
else
    summa=nansum(a_diff_sq,1);
end

if row>1
    v=summa/(row-1);
else
    v=summa/row;
end

rsd=sqrt(v)./abs_median;

end
