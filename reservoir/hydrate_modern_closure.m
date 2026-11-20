function [reaction,formation,dissociation,peq]=hydrate_modern_closure(pg,T,s,sg,sh,phi,vol,c)
% Li (2024), Table 3：qH 单位为 kg 水合物/(m³ 地层·s)。
[peq,~]=hydrate_equilibrium_moridis(T-c.eq_shift,s,c.salt_equilibrium);
f=co2_fugacity_pr(pg,T); feq=co2_fugacity_pr(peq,T);
delta=f-feq;
sa=max(0,1-sg-sh); g=max(0,sg); h=max(0,sh);
small=c.area_regularization;
garea=(g+small).^(2/3)-small^(2/3);
harea=(h+small).^(2/3)-small^(2/3);
area=3*(1-phi)/c.grain_radius;
arr=exp(-c.activation_energy./(8.314*(T+273.15)));
qh_form=c.kinetic_multiplier*c.formation_prefactor*arr.*area.*...
    sa.^c.area_beta.*(1-h).^c.area_beta.*garea.*max(delta,0);
% 分解面积随 Sh^(2/3) 缩小，在 Sh=0 时为零。
qh_diss=c.kinetic_multiplier*c.dissociation_prefactor*arr.*area.*harea.*max(-delta,0);
formation=c.wc*vol.*qh_form; dissociation=c.wc*vol.*qh_diss;
reaction=double(c.reaction)*(formation-dissociation);
end
