function pv = calcmulti(f,VS,H,VP,den)
%   Summary of this function goes here.
%   pv= calcmulti(f,VS,H,VP,den)
%   Detailed explanation goes here.
%   The function is for calculating multi-mode of Rayleigh wave.
%
%   IN    f: row vector of frequency.
%        VS: row vector of shear wave velocity.
%         H: row vector of thickness of layer.
%        VP: row vector of primary wave velocity.
%       den: row vector of density.
%
%  OUT   pv: the phase velocity of multi-mode Rayleigh wave.
%
%  Example:
%  f=5:100;VS=[200 400];VP=[400 800];H=[10];den=[2000 2000];
%  pv=calcmulti(f,VS,H,VP,den);
%
%  Author(s): Yan Yingwei
%  Copyright: 2017-2019 
%  Revision: 1.0  Date: 2/27/2017
%
%  Department of Geophysics, Jilin University.

if(nargin<=4)
    [~,c]=size(VS);
    den=2*ones(1,c);
end
if(nargin<=3)
    VP=2*VS;
end
global n_mode mode_base 
cmax=max(VS);
cmin=0.88*min(VS);
dc1=(cmax-cmin)/400;
cc=cmin:dc1:cmax;

[~,nc]=size(cc);
[~,nf]=size(f);
mode_base=[1,VS,H,VP,den];
for i=1:nf 
    mode_base(1)=f(i);
    r(i,1)=fastcalc(cc(1));
    n=1;
    if(r(i,1)==0)
      ccc(i,n)=cc(1);
      n=n+1;
    end
    for j=2:nc
       r(i,j)=fastcalc(cc(j));
       if r(i,j)*r(i,j-1)<0
           ccd1=cc(j-1);ccd2=cc(j);rrd1=r(i,j-1);rrd2=r(i,j);
           for icd=1:4
                ccd=(ccd2+ccd1)/2;
                rrd=fastcalc(ccd);
                if rrd1*rrd<0
                    ccd2=ccd;
                    rrd2=rrd;
                end
                if rrd2*rrd<0
                    ccd1=ccd;
                    rrd1=rrd;
                end
            end
           ccc(i,n)=(ccd1+ccd2)/2;
           n=n+1;
       end
      if(r(i,j)==0)
           ccc(i,n)=cc(j);
           n=n+1;
       end
      if n==n_mode+1
            break
      end
    end
end
pv=ccc;
end
