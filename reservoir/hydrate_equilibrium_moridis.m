function [peq,deltaT]=hydrate_equilibrium_moridis(T,s,variant)
% Moridis (2025), Tables 4/6/10；T 为 °C，s 为盐水盐质量分数，peq 为 Pa。
% 盐度温降采用 DL 或 DB 关联式。
if nargin<3, variant='DL'; end
switch variant
    case 'DL', a=[0,41.960134480071723,-49.772417348714633,763.49669277683938];
    case 'DB', a=[0,47.372104577621734,96.991483627559859,-438.35370775310838,3391.110070660055];
    otherwise, error('盐度关联只允许已核查的 DL 或 DB。');
end
deltaT=horner(a,s);
tc=T+deltaT;
av=[1.048,.109562462881416,.025880176219489209,-.0078300297683361261,...
    .001559436970817869,-.00013146338302432161,.0000042930042382365193];
% 同一温度在不同分支中的中间求值必须有限；最终逐单元选用相应分支。
v=horner(av,min(tc,10.05)+1.35)*1e6;
th=max(tc,10.05)-10.05;
c=[2.502372010423192,-.2792100157839153,.2594950206586341,-1.569081314601355,5.719452199762003];
% 两条平衡曲线在 Q2 处连续。
pq=horner(av,11.4)*1e6;
l=pq.*10.^(c(2)*(exp(-th/c(3))-1)+c(4)*(exp(-th/c(5))-1));
mask=double(value(tc)>10.05);
peq=(1-mask).*v+mask.*l;
end

function y=horner(a,x)
y=0*x+a(end);
for j=numel(a)-1:-1:1, y=y.*x+a(j); end
end
