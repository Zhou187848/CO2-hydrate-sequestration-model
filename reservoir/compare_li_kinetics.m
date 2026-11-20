function [report,data]=compare_li_kinetics(nrings,write_outputs)
% Run 1 轴线温度和恒温水浴约束下的反应釜动力学重建。
if nargin<1, nrings=8; end
if nargin<2, write_outputs=true; end
validateattributes(nrings,{'numeric'},{'scalar','integer','positive'});
root=fileparts(mfilename('fullpath')); addpath(root);
assert(~isempty(which('initVariablesADI')),'Start MRST 2025b first.');
source=load(fullfile(root,'li_run1_observations.mat'));
obs=source.observed; time=obs(:,2); pp=obs(:,3);
axis_temp=obs(:,4:6)-273.15;
c=struct('salt_equilibrium','DL','eq_shift',0,'grain_radius',1.645e-4,...
    'activation_energy',102880,'formation_prefactor',8.4e11,'area_beta',5.3,...
    'dissociation_prefactor',2.783e7,'kinetic_multiplier',1,'area_regularization',1e-10,...
    'reaction',true,'rho_h',1100);
c.wc=.044010/.152100;
bulk=1.155e-3; pore=3.60e-4; porosity=pore/bulk; water0=.09;
% Li 图 1 和第 4.1 节：三个测点均在轴线上，水浴边界为 1 °C。
% 假设每个等厚轴向层的径向剖面为圆柱均匀放热稳态解的抛物线形状；
% 环带在 (r/R)^2 上等距，因此等体积。各环带独立更新含水量。
wall=1; area_midpoints=((1:nrings)-.5)/nrings;
weights=repmat(1-area_midpoints,3,1);
ncells=numel(weights);
pressure=griddedInterpolant(time,pp,'linear');
temperatures=cell(1,3);
for z=1:3, temperatures{z}=griddedInterpolant(time,axis_temp(:,z),'linear'); end
opt=odeset('RelTol',1e-7,'AbsTol',1e-11,'MaxStep',30,'NonNegative',1:ncells,...
    'JPattern',speye(ncells));
[~,mass]=ode15s(@rhs,time,zeros(ncells,1),opt);
% 局部质量以相同饱和度下的整床等效质量计，取体积平均即总质量。
mh=mean(mass,2); ref=obs(:,7); residual=mh-ref;
assert(all(mass(:)>=-1e-10 & mass(:)<=water0/(1-c.wc)+1e-10));
report=struct('test','Run 1, axis thermometry and fixed bath; parabolic radial reconstruction',...
    'n',numel(time),'final_predicted_hydrate_kg',mh(end),'final_source_inferred_hydrate_kg',ref(end),...
    'endpoint_relative_error',residual(end)/ref(end),'RMSE_kg',sqrt(mean(residual.^2)),...
    'normalized_RMSE',sqrt(mean(residual.^2))/ref(end),'volume_m3',bulk,'porosity',porosity,...
    'water0_kg',water0,'rings_per_layer',nrings,'wall_temperature_C',wall,...
    'limitations','Prescribed radial temperature shape and measured pressure. Hydrate mass inferred from experiment.');
data=table(time,pp,axis_temp(:,2),ref,mh,residual,'VariableNames',...
    {'time_s','observed_P_Pa','observed_T_C','source_hydrate_kg','predicted_hydrate_kg','error_kg'});
if write_outputs
    if ~isfolder(fullfile(root,'results_modern')), mkdir(fullfile(root,'results_modern')); end
    writetable(data,fullfile(root,'results_modern','li_kinetic_comparison.csv'));
    fid=fopen(fullfile(root,'results_modern','li_kinetic_comparison.json'),'w');
    closer=onCleanup(@() fclose(fid));
    fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
end
    function rate=rhs(t,m)
        axis_now=cellfun(@(f) f(t),temperatures);
        temp=wall+(axis_now(:)-wall).*weights;
        sh=max(0,m)/(c.rho_h*pore);
        sw=max(0,water0-(1-c.wc)*max(m,0))/(1000*pore);
        sg=max(0,1-sh-sw);
        rate=hydrate_modern_closure(repmat(pressure(t),ncells,1),temp(:),0,sg,sh,porosity,bulk,c)/c.wc;
    end
end
