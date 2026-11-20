function run_fig6()
root=fileparts(mfilename('fullpath'));
cases=readtable(fullfile(fileparts(root),'data','Supplementary_Data_1_public.xlsx'),'Sheet','F6 Operating cases','TextType','string','VariableNamingRule','preserve');
for j=1:height(cases)
  cfg=modern_configuration(struct('name',char(string(cases{j,1})),'well_control','common_bhp',...
    'sand_k',140*9.869233e-16,'rate_kg_s',cases.rate_t_h(j)*1000/3600,...
    'inject_temperature',cases.entry_C(j),'results_folder','results_fig6',...
    'snapshot_days',[.5,1,2,3,5,7,10]));
  out=axisymmetric_co2(cfg);
  assert(abs(out.summary.injected_kg/1000-cases.rate_t_h(j)*239.5)<1e-4);
  assert(abs(out.summary.hydrate_fraction_pct-cases.hydrate_percent(j))<1e-4);
  assert(abs(out.summary.overpressure_MPa-cases.pressure_MPa(j))<1e-3);
end
end
