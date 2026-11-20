function [f,Z]=co2_fugacity_pr(p,T)
% 纯 CO2 Peng–Robinson 逸度；Pa、°C 输入。稳定根由最低残余 Gibbs 能选择。
% 逸度用 PR，密度和焓用 HEOS。
% 先数值选择根，再隐式微分 cubic，保留对 MRST ADI 的压力/温度导数。
R=8.314462618; Tc=304.1282; Pc=7.3773e6; omega=.22394;
Tk=T+273.15;
kappa=.37464+1.54226*omega-.26992*omega^2;
aa=.45724*R^2*Tc^2/Pc*(1+kappa*(1-(Tk/Tc).^.5)).^2;
bb=.07780*R*Tc/Pc;
A=aa.*p./(R^2*Tk.^2); B=bb*p./(R*Tk);
av=value(A); bv=value(B);
c2=bv-1; c1=av-3*bv.^2-2*bv; c0=-(av.*bv-bv.^2-bv.^3);
pp=c1-c2.^2/3; qq=2*c2.^3/27-c2.*c1/3+c0;
disc=(qq/2).^2+(pp/3).^3;
z=zeros(numel(av),3); z(:)=NaN;
one=disc>=0;
cube=@(x) sign(x).*abs(x).^(1/3);
z(one,1)=cube(-qq(one)/2+sqrt(disc(one)))+cube(-qq(one)/2-sqrt(disc(one)))-c2(one)/3;
three=~one;
theta=acos(max(-1,min(1,(-qq(three)/2)./sqrt(-(pp(three)/3).^3))));
for j=1:3
    z(three,j)=2*sqrt(-pp(three)/3).*cos((theta+2*pi*(j-1))/3)-c2(three)/3;
end
cost=Inf(size(z));
for j=1:3
    valid=z(:,j)>bv;
    zz=z(valid,j); b=bv(valid); a=av(valid);
    cost(valid,j)=zz-1-log(zz-b)-a./(2*sqrt(2)*b).*log((zz+(1+sqrt(2))*b)./(zz+(1-sqrt(2))*b));
end
[best,index]=min(cost,[],2); assert(all(isfinite(best)),'PR 未找到可行稳定根。');
zn=z(sub2ind(size(z),(1:numel(av))',index));
if isa(A,'ADI')
    residual=zn.^3+(B-1).*zn.^2+(A-3*B.^2-2*B).*zn-(A.*B-B.^2-B.^3);
    deriv=3*zn.^2+2*c2.*zn+c1;
    Z=zn-(residual-value(residual))./deriv;
else
    Z=zn;
end
lnphi=Z-1-log(Z-B)-A./(2*sqrt(2)*B).*log((Z+(1+sqrt(2))*B)./(Z+(1-sqrt(2))*B));
f=p.*exp(lnphi);
end
