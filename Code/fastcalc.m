function y=fastcalc(x)
%   Summary of this function goes here
%   function y=fastcalc(x)
%   The function is the key function for calculating the disperive function
%   of Rayleigh wave.
%   Detailed explanation goes here
%   The function would search the zero solution of fastcalc.
%    x: The solution of phase velocity.
%
%   References: 
%   Fan, Y. & Liu, J. 2001. Research on the dispersion of Rayleigh waves in 
%   multilayered media, Journal of HarBin Institute of Technology. (in Chinese), 
%   33(5), 577-581, https://doi.org/10.3321/j.issn:0367-6234.2001.05.001.
%
%  Author(s): Yan Yingwei
%  Revision: 1.0  Date: 2/25/2016.
%
%  Department of Geophysics, Jilin University.

%%
global mode_base
f=mode_base(1);
[~,n]=size(mode_base);
n=n/4;
VS=mode_base(2:n+1);
H=mode_base(n+2:2*n);
VP=mode_base(2*n+1:3*n);
den=mode_base(3*n+1:4*n);
% k is the circular wavenumber.
k=2*pi*f/x;  
xs=x^2;
%%
% nth Layer. 
btv=VS(n)^2;
bt=den(n)*btv;
af=VP(n)^2;
g=xs/(btv+btv);
t=1-g;
r=xs/af;
r=r-1;
s=g+g-1;
r1=sqrt(r);
s1=sqrt(s);
ps=r1*s1;
bt1=bt;
x1=1+ps;
x2=t+ps;
x4=1i*s1*g;
x5=-1i*r1*g;
x3=-t^2-ps;
%%
% (n-1)th ~ first Layer.
for ii=n-1:-1:1
   btv=VS(ii)^2;
   bt=den(ii)*btv;
   af=VP(ii)^2;
   g=xs/(btv+btv);
   t=1-g;
   r=xs/af;
   r=r-1;
   s=g+g-1;
   r1=sqrt(r);
   s1=sqrt(s);
   l=bt1/bt;
   bt1=bt;
   x1=x1/l;
   x3=x3*l;
   k1=k*H(ii);
   p=k1*r1;
   q=k1*s1;
   
   tx1=t*x1;
   ttx1=t*tx1;
   tx2=t*x2;
   p1=x1-x2-x2-x3;
   p2=-ttx1+tx2+tx2+x3;
   p3=g*x4;
   p4=g*x5;
   p5=-tx1+x2+tx2+x3;
   if x>=VP(ii)
     if r1==0
        c=k1;
     else
        c=sin(p)/r1;
     end
     d=sin(q)/s1;
     a=cos(p);
     b=cos(q);
     ab=a*b;
     ad=a*d;
     cd=c*d;
     bc=b*c;
     ads=ad*s;
     bcr=bc*r;
     cds=cd*s;
     cdr=cd*r;
     cdrs=cdr*s;
  	 q1=ab*p1+cd*p2-ad*p3+bc*p4;
     q2=cdrs*p1+ab*p2+bcr*p3-ads*p4;
  	 q3=ads*p1-bc*p2+ab*p3+cds*p4;
     q4=-bcr*p1+ad*p2+cdr*p3+ab*p4;
  	 tq1=t*q1;
     ttq1=t*tq1;
     tp5=t*p5;
  	 x1=q1-q2+p5+p5;
     x2=tq1-q2+p5+tp5;
     x3=-ttq1+q2-tp5-tp5;
     x4=g*q3;
     x5=g*q4;
   end
%%
% delete the exp-growing part.
   if x<VP(ii)&&x>=VS(ii) 
      if s1==0 
        d=k1;
      else
        d=sin(q)/s1;
      end
      ar1=abs(r1);
      ark1=ar1*k1;
      ee=exp(-2*ark1);
  	  b=cos(q);
      c=(1-ee)/(1+ee)/ar1; 
  	  ds=d*s;
      br=b*r;
      dr=d*r;
      drs=dr*s;
      q1=b*p1+c*d*p2-d*p3+b*c*p4;
      q2=c*drs*p1+b*p2+c*br*p3-ds*p4;
  	  q3=ds*p1-b*c*p2+b*p3+c*ds*p4;
      q4=-c*br*p1+d*p2+c*dr*p3+b*p4;
  	  tq1=t*q1;
      ttq1=t*tq1;
  	  if ark1<20
     	p5=p5/cos(p);
  	  else
     	p5=0;
  	  end
      tp5=t*p5;
      x1=q1-q2+p5+p5;
      x2=tq1-q2+p5+tp5;
      x3=-ttq1+q2-tp5-tp5;
      x4=g*q3;
      x5=g*q4;
   end
%%
% delete the exp-growing part
   if x<VS(ii)      
     ar1=abs(r1);
     ark1=ar1*k1;
     ee=exp(-2*ark1);
     
     as1=abs(s1);
     ask1=as1*k1;
     ees=exp(-2*ask1);
     c=(1-ee)/(1+ee)/ar1;
     d=(1-ees)/(1+ees)/as1;
     
     ds=d*s;
     cr=c*r;
     cd=c*d;
     cds=c*ds;
     cdr=cr*d;
     cdrs=ds*cr;
     
     q1=p1+cd*p2-d*p3+c*p4;
     q2=cdrs*p1+p2+cr*p3-ds*p4;
     q3=ds*p1-c*p2+p3+cds*p4;
     q4=-cr*p1+d*p2+cdr*p3+p4;
     tq1=t*q1;
     ttq1=t*tq1;
     if (ask1+ark1)<20
        p5=p5/cos(p)/cos(q);
     else
        p5=0;
     end
    	tp5=t*p5;
    	x1=q1-q2+p5+p5;
        x2=tq1-q2+p5+tp5;
        x3=-ttq1+q2-tp5-tp5;
        x4=g*q3;
        x5=g*q4;
   end
end
y=real(x3);
end

