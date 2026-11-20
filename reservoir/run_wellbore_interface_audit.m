function audit_summary=run_wellbore_interface_audit(audit_bhp,audit_tin,audit_rate,cfg)
% 给定井底压力、井口温度和注入速率，计算出口温度。
if nargin<1, audit_bhp=20.4; end
if nargin<2, audit_tin=-15; end
if nargin<3, audit_rate=10; end
audit_root=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(audit_root),'wellbore'));
assert(nargin>=4 && isstruct(cfg),'Supply a wellbore configuration as the fourth argument.');
cfg.injection.mass_flow_rate_kg_s=audit_rate*1000/3600;
cfg.injection.free_water_mass_ratio=0; cfg.injection.free_water_ppm=0;
cfg.injection.temp_inlet=audit_tin;
cfg.injection.bh_target_pressure=audit_bhp;
cfg.injection.solve_inlet_pressure=true;
cfg.injection.fracture_gradient_g_cm3=NaN;
cfg.injection.fracture_factor=NaN;
cfg.injection.liquid_transport.enabled=false;
cfg.injection.liquid_transport.enforce=false;
audit_log=evalc('audit_result=CO2_2(cfg);');
audit_folder=fullfile(audit_root,'results_interface_audit');
if ~isfolder(audit_folder), mkdir(audit_folder); end
audit_tag=sprintf('wellbore_P%05d_T%+03d_Q%02d',round(audit_bhp*1000),audit_tin,audit_rate);
save(fullfile(audit_folder,[audit_tag '.mat']),'audit_result','cfg','-v7');
fid=fopen(fullfile(audit_folder,[audit_tag '.log']),'w'); fprintf(fid,'%s',audit_log); fclose(fid);
rr=audit_result.results;
audit_summary=struct('target_BHP_MPa',audit_bhp,'actual_BHP_MPa',rr.P_co2(end),...
    'outlet_T_C',rr.T_co2(end),'inlet_P_MPa',rr.P_co2(1),...
    'mass_rate_kg_s',cfg.injection.mass_flow_rate_kg_s,'inlet_T_C',audit_tin,...
    'reference_depth_mRT',cfg.well.total_depth,'max_vapor_fraction',max(rr.beta_vapor));
audit_summary.pressure_residual_MPa=rr.P_co2(end)-audit_bhp;
audit_summary.pressure_matched=abs(audit_summary.pressure_residual_MPa)<=.02;
fid=fopen(fullfile(audit_folder,[audit_tag '.json']),'w'); fprintf(fid,'%s',jsonencode(audit_summary,PrettyPrint=true)); fclose(fid);
disp(audit_summary);
end
