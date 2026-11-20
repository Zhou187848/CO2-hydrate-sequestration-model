function run_c1(workbook,parameter_file)
% Read authorized local forcing and layer inputs.
here=fileparts(mfilename('fullpath')); addpath(fullfile(fileparts(here),'reservoir'));
base=load(parameter_file,'cfg');
p=readtable(workbook,'Sheet','C1_重建排程','VariableNamingRule','preserve');
t=readtable(workbook,'Sheet','C1_入口温度','VariableNamingRule','preserve');
layers=readtable(workbook,'Sheet','C1_分层渗透率','VariableNamingRule','preserve');
targets=readtable(workbook,'Sheet','C1_温度与库存','VariableNamingRule','preserve');
temps=targets.T123_initial_C;
assert(height(targets)==3);

for i=1:3
 c=base.cfg;c.name=sprintf('c1_temperature_%d',i);c.save=true;c.results_folder='results_c1';
 c.layer_profile=[layers.top_mBML,layers.bottom_mBML,layers.model_mD];
 c.field_schedule=[p.start_s,p.end_s,p.rate_t_h*1000/3600];
 c.field_temp=[t.seconds_since_0700,t.inlet_proxy_C];
 c.geothermal_gradient=(temps(i)-3)/123;
 out=profile_case(c);
 assert(abs(out.history(end,5)-1000*targets.injection_t(i))<1e-4);
end
end
