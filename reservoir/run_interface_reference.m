function run_interface_reference()
% 十天注入基准算例，入口温度取井筒计算值。
c=modern_configuration(struct('name','baseline_interface','well_control','common_bhp',...
    'inject_temperature',5.50,'sand_k',140*9.869233e-16,'results_folder','results_interface_final',...
    'snapshot_days',[.5,1,2,3,5,7,10]));
o=axisymmetric_co2(c);
assert(abs(o.summary.CO2_defect_kg)/o.summary.injected_kg<1e-6);
assert(max(abs(o.diagnostics(:,3)-diff([0;o.history(:,5)])./o.history(:,2)))<1e-8);
end
