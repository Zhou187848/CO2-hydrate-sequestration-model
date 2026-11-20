function run_interface_grid(partition,parts)
% 以同一井压连接、入口热条件重算地层参数扫描。
c=struct('well_control','common_bhp','inject_temperature',5.50,...
    'results_folder','results_interface_final','snapshot_days',[.5,1,2,3,5,7,10]);
run_modern_suite('grid',partition,parts,c);
end
