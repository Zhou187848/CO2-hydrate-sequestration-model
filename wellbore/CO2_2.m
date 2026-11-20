function out = CO2_2(cfg)
% 井筒准稳态多相流与传热。沿井深输出压力、温度、相态及沉积厚度。

if nargin < 1 || ~isstruct(cfg)
    error('CO2_2: cfg 必须为 struct');
end
if ~isfield(cfg, 'options') || ~isstruct(cfg.options)
    error('CO2_2: cfg.options 缺失或不是 struct');
end
if ~isfield(cfg, 'model') || ~isstruct(cfg.model)
    error('CO2_2: cfg.model 缺失或不是 struct');
end
if ~isfield(cfg, 'well') || ~isstruct(cfg.well)
    error('CO2_2: cfg.well 缺失或不是 struct');
end
if ~isfield(cfg, 'injection') || ~isstruct(cfg.injection)
    error('CO2_2: cfg.injection 缺失或不是 struct');
end
if ~isfield(cfg, 'thermal') || ~isstruct(cfg.thermal)
    error('CO2_2: cfg.thermal 缺失或不是 struct');
end
if isfield(cfg.options,'storage') && isfield(cfg.options.storage,'enabled') && cfg.options.storage.enabled
    error('CO2_2:Storage','Run the reservoir model separately with axisymmetric_co2.');
end

options = cfg.options;
model = cfg.model;
well = cfg.well;
injection = cfg.injection;
thermal = cfg.thermal;

out = struct();
tic;

plot_level = "summary";
if isfield(options, 'plot_level')
    plot_level = lower(strtrim(string(options.plot_level)));
end
if plot_level ~= "summary" && plot_level ~= "full"
    error('options.plot_level 仅支持 summary / full');
end

gcfg = setup_graphics_runtime(options);
if isstruct(gcfg) && isfield(gcfg, 'auto_export') && gcfg.auto_export
    fprintf('[可视化] auto_export=ON -> %s\n', char(getfield_or_default(gcfg, 'export_dir_abs', gcfg.export_dir)));
end

environment = options.environment;
environment.source = lower(string(getfield_or_default(environment, 'source', 'woa18')));
environment.lat = getfield_or_default(environment, 'lat', NaN);
environment.lon = getfield_or_default(environment, 'lon', NaN);
environment.month = round(getfield_or_default(environment, 'month', 11));
environment.resolution = char(getfield_or_default(environment, 'resolution', '1.00'));
environment.cache_mat = char(getfield_or_default(environment, 'cache_mat', 'woa18_env_cache.mat'));
environment.enable_cache = logical(getfield_or_default(environment, 'enable_cache', true));
environment.timeout_s = max(5.0, getfield_or_default(environment, 'timeout_s', 30));
environment.max_extrap_m = max(0.0, getfield_or_default(environment, 'max_extrap_m', 300));
environment.verbose = logical(getfield_or_default(environment, 'verbose', true));
environment.air_temp_C = getfield_or_default(environment, 'air_temp_C', 25.0);

env_profile = struct('ok', false);
if environment.source == "woa18"
    assert(isfinite(environment.lat) && isfinite(environment.lon), "Supply a permitted WOA18 sampling point in cfg.options.environment.");
    env_profile = woa18_profile_ncss(environment.lat, environment.lon, environment.month, environment.resolution, ...
        environment.timeout_s, environment.cache_mat, environment.enable_cache, environment.verbose);
    if env_profile.ok
        depth_max = max(env_profile.depth_m);
        depth_ok = isfield(env_profile, 'depth_m') && ~isempty(env_profile.depth_m) && (depth_max + environment.max_extrap_m) >= well.seawater_depth;
        if ~depth_ok
            warning('[环境场] WOA18 点位的有效水深不足以覆盖 seawater_depth=%.0fm，回退到参数化温度分层。', well.seawater_depth);
            env_profile.ok = false;
            environment.source = "param";
        else
            if depth_max < well.seawater_depth
                warning('[环境场] WOA18 剖面最大深度=%.0fm，海底深度=%.0fm，将对深层做线性外推（max_extrap_m=%.0fm）。', ...
                    depth_max, well.seawater_depth, environment.max_extrap_m);
            end
            thermal.sea_surface_temp = interp1(env_profile.depth_m, env_profile.T_C, 0.0, 'linear', 'extrap');
            thermal.seabed_temp = interp1(env_profile.depth_m, env_profile.T_C, well.seawater_depth, 'linear', 'extrap');
            s_mud_psu = interp1(env_profile.depth_m, env_profile.S_psu, well.seawater_depth, 'linear', 'extrap');
            thermal.seawater_salinity = max(0.0, s_mud_psu) / 10.0; % 近似：PSU≈g/kg => wt%≈PSU/10
        end
    else
        warning('[环境场] WOA18 数据获取失败，回退到参数化温度分层（thermal.sea_surface_temp/seabed_temp）。');
        environment.source = "param";
    end
end

% 按给定泥线温度校正海水剖面。
mudline_temp_target_C = getfield_or_default(environment, 'mudline_temp_C', NaN);
if isfinite(mudline_temp_target_C)
    z_mud_sea = getfield_or_default(well, 'seawater_depth', NaN); % 相对海平面的水深（m）
    if ~isfinite(z_mud_sea) && isfield(well, 'mudline_depth') && isfield(well, 'platform_elevation')
        z_mud_sea = well.mudline_depth - well.platform_elevation;
    end
    if isfinite(z_mud_sea) && z_mud_sea >= 0
        if isfield(env_profile, 'ok') && env_profile.ok ...
                && isfield(env_profile, 'depth_m') && isfield(env_profile, 'T_C') ...
                && ~isempty(env_profile.depth_m) && ~isempty(env_profile.T_C)
            T_mud_curr_C = interp1(env_profile.depth_m, env_profile.T_C, z_mud_sea, 'linear', 'extrap');
            dT_C = mudline_temp_target_C - T_mud_curr_C;
            env_profile.T_C = env_profile.T_C + dT_C;
            thermal.sea_surface_temp = thermal.sea_surface_temp + dT_C;
            thermal.seabed_temp = mudline_temp_target_C;
        else
            thermal.seabed_temp = mudline_temp_target_C;
        end
    end
end

%% CO2物性常数
CO2 = struct();
CO2.Tc = 31.1;              % 临界温度（°C）
CO2.Pc = 7.38;              % 临界压力（MPa）
CO2.Mw = 44.01;             % 分子量（g/mol）
CO2.R = 8.314;              % 气体常数（J/(mol·K)）
CO2.H_vap = 235e3;          % 蒸发潜热（J/kg）

%% 数值网格
dz = options.dz;  % 网格步长（m）
depth_grid = (0:dz:well.total_depth)';
n_nodes = length(depth_grid);
fprintf('[网格生成] %d 计算节点, 步长 %.2f m\n', n_nodes, dz);

%% 海洋/地层温度剖面
fprintf('[环境场] 构建海洋/地层温度剖面...\n');
T_env = zeros(n_nodes, 1);

idx_sea_surface = find(depth_grid >= well.platform_elevation, 1);
if isempty(idx_sea_surface)
    idx_sea_surface = 1;
end

if env_profile.ok
    fprintf('  环境源: WOA18 decav | 月=%02d | 点位(%.3f°N, %.3f°E) | 网格=%s\n', ...
        environment.month, environment.lat, environment.lon, environment.resolution);
end

T_200_param_C = thermal.sea_surface_temp - 20.0; % 参数化温跃层底界温度（°C）

for i = 1:n_nodes
    z = depth_grid(i);

    if z <= well.mudline_depth
        % 海水段：基于实际海洋温度剖面
        z_sea = z - well.platform_elevation;  % 相对海平面的深度

        if z_sea < 0
            % 平台段（暴露在空气中）
            T_env(i) = environment.air_temp_C;
        elseif z_sea < 50
            if env_profile.ok
                T_env(i) = interp1(env_profile.depth_m, env_profile.T_C, z_sea, 'linear', 'extrap');
            else
                % 混合层（参数化）
                T_env(i) = thermal.sea_surface_temp;
            end
        elseif z_sea < 200
            if env_profile.ok
                T_env(i) = interp1(env_profile.depth_m, env_profile.T_C, z_sea, 'linear', 'extrap');
            else
                % 温跃层（参数化）
                T_env(i) = thermal.sea_surface_temp - 20*(z_sea-50)/150;
            end
        elseif z_sea < 1000
            if env_profile.ok
                T_env(i) = interp1(env_profile.depth_m, env_profile.T_C, z_sea, 'linear', 'extrap');
            else
                % 中层水（参数化）
                T_env(i) = T_200_param_C + (thermal.seabed_temp - T_200_param_C) * (z_sea - 200) / 800;
            end
        else
            if env_profile.ok
                T_env(i) = interp1(env_profile.depth_m, env_profile.T_C, z_sea, 'linear', 'extrap');
            else
                % 深层水（参数化）
                T_env(i) = thermal.seabed_temp;
            end
        end
    else
        % 地层段
        z_formation = z - well.mudline_depth;
        T_env(i) = thermal.seabed_temp + thermal.geothermal_gradient * z_formation;
    end
end

fprintf('  海表温度: %.1f°C\n', T_env(idx_sea_surface));
idx_mudline_sea = find(depth_grid <= well.mudline_depth, 1, 'last');
if isempty(idx_mudline_sea)
    idx_mudline_sea = 1;
end
fprintf('  泥线海水温度: %.1f°C (水深%.0fm)\n', T_env(idx_mudline_sea), well.seawater_depth);
if env_profile.ok
    s_mud_psu = interp1(env_profile.depth_m, env_profile.S_psu, well.seawater_depth, 'linear', 'extrap');
    fprintf('  泥线盐度: %.2f wt%% (WOA18 S=%.2f)\n', thermal.seawater_salinity, s_mud_psu);
end
fprintf('  井底温度: %.1f°C (深度%dm)\n', T_env(end), well.total_depth);

%% 主计算：多相水力 + 多层传热 + 水合物动力学（OLGA风格闭合）
fprintf('\n[主计算] 多相-传热-水合物耦合求解...\n');

const = struct('g', 9.80665, 'R', 8.314462618);

% 约束：破裂压力上限（可选）与井底目标压力
[P_frac_MPa, fracinfo] = estimate_fracture_pressure_MPa(well, injection, const);
P_bh_max_MPa = NaN;
if isfinite(P_frac_MPa)
    P_bh_max_MPa = P_frac_MPa - max(getfield_or_default(injection, 'fracture_margin_MPa', 0.0), 0.0);
end
if ~isfield(injection, 'bh_target_pressure') || ~isfinite(injection.bh_target_pressure)
    injection.bh_target_pressure = injection.formation_pressure;
end
if isfinite(P_bh_max_MPa) && injection.bh_target_pressure > P_bh_max_MPa
    warning('[约束] 井底目标压力 %.2f MPa 超过破裂上限 %.2f MPa，已截断。', injection.bh_target_pressure, P_bh_max_MPa);
    injection.bh_target_pressure = P_bh_max_MPa;
end
fprintf('[约束] 井底目标=%.2f MPa', injection.bh_target_pressure);
if isfinite(P_frac_MPa)
    frac_method_str = '';
    if isstruct(fracinfo) && isfield(fracinfo, 'method_str') && ~isempty(fracinfo.method_str)
        frac_method_str = [' | ' char(fracinfo.method_str)];
    end
    fprintf(' | 破裂≈%.2f MPa%s（裕度%.2f -> 上限%.2f）', ...
        P_frac_MPa, frac_method_str, getfield_or_default(injection, 'fracture_margin_MPa', 0.0), P_bh_max_MPa);
end
fprintf('\n');

% 井口T/P调参扫掠（可选）：固定 T_inlet，扫描 P_inlet（表压），检查全程液态 & 井底不超破裂
if isfield(options, 'inlet_tuning') && isstruct(options.inlet_tuning) ...
        && isfield(options.inlet_tuning, 'enabled') && options.inlet_tuning.enabled
    tune_cfg = options.inlet_tuning;
    if isfield(tune_cfg, 'temp_inlet_C') && isfinite(tune_cfg.temp_inlet_C)
        injection.temp_inlet = tune_cfg.temp_inlet_C;
        if isfield(injection, 'liquid_transport') && isstruct(injection.liquid_transport)
            injection.liquid_transport.min_temp_inlet_C = injection.temp_inlet;
        end
    end

	    sweep = sweep_inlet_pressure_constraints(depth_grid, dz, T_env, injection, model, well, thermal, CO2, const, P_bh_max_MPa, tune_cfg);
	    if isfield(tune_cfg, 'stop_after_scan') && tune_cfg.stop_after_scan
	        if isstruct(sweep) && isfield(sweep, 'ok') && sweep.ok
	            fprintf('[调参] 扫掠完成（stop_after_scan=true）。\n');
	        end
	        fprintf('[调参] 提示：stop_after_scan=true 将在扫掠后提前返回，不会生成可视化图窗；如需出图请设 stop_after_scan=false 或关闭 inlet_tuning。\n');
	        try
	            out.options = options;
	            out.model = model;
	            out.well = well;
	            out.injection = injection;
            out.thermal = thermal;
            out.environment = environment;
            out.environment_profile = env_profile;
            out.CO2 = CO2;
            out.const = const;
            out.sweep = sweep;
        catch ME
            warning('[调参] stop_after_scan 输出组装失败：%s', ME.message);
        end
        return;
    end
end

simopt = struct();
simopt.verbose = true;
simopt.progress_every = 500;
simopt.show_progress = isfield(options, 'progress') && isstruct(options.progress) ...
    && isfield(options.progress, 'enabled') && options.progress.enabled; % 是否显示总进度（-）
simopt.progress_update_interval_s = getfield_or_default(getfield_or_default(options, 'progress', struct()), 'update_interval_s', 0.5); % 刷新间隔（s）
% 等效注入持续时间（s）：用于把瞬时沉积通量外推为累计沉积厚度（准稳态筛查参数，可改）。
simopt.injection_duration_s = max(getfield_or_default(getfield_or_default(options, 'storage', struct()), 'injection_duration_s', ...
    getfield_or_default(getfield_or_default(model, 'hydrate', struct()), 'simulation_time_s', 0.0)), 0.0);
% 工程约束为全程液态下输 + 井底不超破裂上限。自动调整顺序：
% 先在破裂上限内提高井底目标压力（超压），若仍不足再降低井口温度（平台冷却）。
% 井口压力由边界匹配反求。
use_boundary_solver = isfield(injection, 'solve_inlet_pressure') && logical(injection.solve_inlet_pressure);
dz_solve = dz;
depth_grid_solve = depth_grid;
T_env_solve = T_env;
if use_boundary_solver
    dz_solve = dz;
    if isfield(injection, 'solve_dz_m')
        dz_solve = max(dz, injection.solve_dz_m);
    end
    if dz_solve > dz
        n_seg_solve = max(1, ceil(well.total_depth / dz_solve));
        dz_solve = well.total_depth / n_seg_solve;
        depth_grid_solve = (0:dz_solve:well.total_depth)';
        T_env_solve = interp1(depth_grid, T_env, depth_grid_solve, 'linear', 'extrap');
        fprintf('  [边界匹配] 求解网格步长=%.2f m（主网格=%.2f m）\n', dz_solve, dz);
    end
end

simopt_solver = simopt;
simopt_solver.verbose = false;
simopt_solver.progress_every = 1000000000;
simopt_solver.show_progress = simopt.show_progress;

% 先在求解网格上完成边界匹配/液态校核
injection_work = injection;
if use_boundary_solver
    fprintf('  [提示] solve_inlet_pressure=true：pressure_inlet=%.3f MPa(g) 仅作初值，后续将被反求结果覆盖；若要固定井口表压请设 solve_inlet_pressure=false。\n', ...
        max(injection_work.pressure_inlet, 0.0));
    fprintf('\n[边界匹配] 反求井口压力，使井底压力≈目标压力 %.2f MPa...\n', injection_work.bh_target_pressure);
    [injection_work, sim_solve] = solve_inlet_pressure_to_match_formation(...
        depth_grid_solve, dz_solve, T_env_solve, injection_work, model, well, thermal, CO2, const, simopt_solver);
else
    sim_solve = simulate_olga_coupled(depth_grid_solve, dz_solve, T_env_solve, injection_work, model, well, thermal, CO2, const, simopt_solver);
end

% 液态下输：工程约束全程液态
margin_req = NaN;
if isfield(injection_work, 'liquid_transport') && isstruct(injection_work.liquid_transport) ...
        && isfield(injection_work.liquid_transport, 'enabled') && injection_work.liquid_transport.enabled
    margin_req = max(getfield_or_default(injection_work.liquid_transport, 'margin_MPa', 0.0), 0.0);
    lt = check_liquid_transport(depth_grid_solve, sim_solve, CO2, margin_req);
    fprintf('[液态下输] 校核(求解网格): min(P-Psat)=%.2f MPa @ %.0fm | max(β_v)=%.3f | %s\n', ...
        lt.min_margin_MPa, lt.depth_min_m, lt.max_beta_vapor, ternary(lt.ok, 'OK', 'WARN'));

    enforce = logical(getfield_or_default(injection_work.liquid_transport, 'enforce', false));
    if enforce && ~lt.ok
        max_iter_cfg = max(1, round(getfield_or_default(injection_work.liquid_transport, 'max_iter', 30)));

        % 1) 优先提高井底目标压力（通过超压；受破裂上限约束）
        dP_step = max(getfield_or_default(injection_work.liquid_transport, 'overpressure_step_MPa', 0.2), 0.05);
        P_bh_limit = P_bh_max_MPa;
        if isfinite(P_bh_limit)
            fprintf('[液态下输] 自动调整：优先提高井底目标压力（步长=%.2f MPa，上限=%.2f MPa）...\n', dP_step, P_bh_limit);
        else
            fprintf('[液态下输] 自动调整：优先提高井底目标压力（步长=%.2f MPa，无破裂上限）...\n', dP_step);
        end

        for it = 1:max_iter_cfg
            if lt.ok
                break;
            end
            if isfinite(P_bh_limit) && injection_work.bh_target_pressure >= P_bh_limit - 1e-9
                break;
            end

            P_next = injection_work.bh_target_pressure + dP_step;
            if isfinite(P_bh_limit)
                P_next = min(P_next, P_bh_limit);
            end
            if P_next <= injection_work.bh_target_pressure + 1e-12
                break;
            end
            injection_work.bh_target_pressure = P_next;
            if isfield(injection_work, 'formation_pressure') && isfinite(injection_work.formation_pressure)
                injection_work.injection_overpressure_MPa = max(P_next - injection_work.formation_pressure, 0.0);
            end

            if use_boundary_solver
                [injection_work, sim_solve] = solve_inlet_pressure_to_match_formation(...
                    depth_grid_solve, dz_solve, T_env_solve, injection_work, model, well, thermal, CO2, const, simopt_solver);
            else
                sim_solve = simulate_olga_coupled(depth_grid_solve, dz_solve, T_env_solve, injection_work, model, well, thermal, CO2, const, simopt_solver);
            end

            lt = check_liquid_transport(depth_grid_solve, sim_solve, CO2, margin_req);
            fprintf('  [调压 %d/%d] P_bh_target=%.2f MPa -> min(P-Psat)=%.2f MPa @ %.0fm | %s\n', ...
                it, max_iter_cfg, injection_work.bh_target_pressure, lt.min_margin_MPa, lt.depth_min_m, ternary(lt.ok, 'OK', 'WARN'));
        end
        if ~lt.ok && isfinite(P_bh_limit) && injection_work.bh_target_pressure >= P_bh_limit - 1e-9
            warning('[液态下输] 已触及破裂上限 %.2f MPa，仍未满足全程液态。', P_bh_limit);
        end

        % 2) 若仍不满足，再尝试降低井口温度（平台冷却）
        Tmin = getfield_or_default(injection_work.liquid_transport, 'min_temp_inlet_C', -5.0);
        dT = max(getfield_or_default(injection_work.liquid_transport, 'temp_step_C', 1.0), 0.1);
        max_iter_need = ceil(max(injection_work.temp_inlet - Tmin, 0.0) / dT);
        max_iter_T = min(max_iter_cfg, max_iter_need);
        if ~lt.ok && max_iter_T > 0
            fprintf('[液态下输] 自动调整：在 P_bh_target=%.2f MPa 下，降低井口温度并重新边界匹配（最多%d次）...\n', ...
                injection_work.bh_target_pressure, max_iter_T);
            for it = 1:max_iter_T
                if lt.ok
                    break;
                end

                T_next = injection_work.temp_inlet - dT;
                if T_next < Tmin
                    break;
                end

                injection_work.temp_inlet = T_next;
                if use_boundary_solver
                    [injection_work, sim_solve] = solve_inlet_pressure_to_match_formation(...
                        depth_grid_solve, dz_solve, T_env_solve, injection_work, model, well, thermal, CO2, const, simopt_solver);
                else
                    sim_solve = simulate_olga_coupled(depth_grid_solve, dz_solve, T_env_solve, injection_work, model, well, thermal, CO2, const, simopt_solver);
                end

                lt = check_liquid_transport(depth_grid_solve, sim_solve, CO2, margin_req);
                fprintf('  [降温 %d/%d] T_inlet=%.1f°C -> min(P-Psat)=%.2f MPa @ %.0fm | %s\n', ...
                    it, max_iter_T, injection_work.temp_inlet, lt.min_margin_MPa, lt.depth_min_m, ternary(lt.ok, 'OK', 'WARN'));
            end
        elseif ~lt.ok
            fprintf('[液态下输] 井口温度已到下限 %.1f°C，跳过降温调整。\n', Tmin);
        end

        if ~lt.ok
            warning('[液态下输] 自动调整后仍未满足全程液态：min(P-Psat)=%.2f MPa @ %.0fm。建议提高井底超压/降低井口温度/提高流量或重新核对破裂上限。', ...
                lt.min_margin_MPa, lt.depth_min_m);
        end
    end
end


% 最终在主网格上运行一次主计算
injection = injection_work;
sim = simulate_olga_coupled(depth_grid, dz, T_env, injection, model, well, thermal, CO2, const, simopt);

fig_main = [];
if isfield(options, 'enable_plots') && options.enable_plots
    fprintf('\n[可视化] 生成井筒汇总图...\n');
    fig_main = plot_wellbore_summary(depth_grid, T_env, sim, well, injection, CO2);
    graphics_keepalive();
    maybe_export_figure(fig_main, 'Fig01_Wellbore_Summary');
end

% 最终液态校核（主网格）
if isfinite(margin_req)
    lt_final = check_liquid_transport(depth_grid, sim, CO2, margin_req);
    fprintf('[液态下输] 校核(主网格):  min(P-Psat)=%.2f MPa @ %.0fm | max(β_v)=%.3f | %s\n', ...
        lt_final.min_margin_MPa, lt_final.depth_min_m, lt_final.max_beta_vapor, ternary(lt_final.ok, 'OK', 'WARN'));
end

T_co2 = sim.T_co2;
P_co2 = sim.P_co2;
rho_co2 = sim.rho_co2;
mu_co2 = sim.mu_co2;
cp_co2 = sim.cp_co2;
k_co2 = sim.k_co2;
phase = sim.phase;
v_co2 = sim.v_co2;
Re_co2 = sim.Re_co2;
dpdz_fric_Pa_m = sim.dpdz_fric_Pa_m;
tau_w_Pa = sim.tau_w_Pa;
flow_regime = sim.flow_regime;
v_g_phase = sim.v_g_phase;
v_l_phase = sim.v_l_phase;
v_rel_gl = sim.v_rel_gl;
a_i_m2_m3 = sim.a_i_m2_m3;
tau_i_Pa = sim.tau_i_Pa;
K_i_N_s_m4 = sim.K_i_N_s_m4;
Q_heat = sim.Q_heat;
Q_hydrate_rxn = sim.Q_hydrate_rxn;
U_total = sim.U_total;
h_inner = sim.h_inner;
mu_JT = sim.mu_JT;
beta_vapor = sim.beta_vapor;
alpha_g = sim.alpha_g;
alpha_g_drift = getfield_or_default(sim, 'alpha_g_drift', NaN(size(alpha_g)));
alpha_g_interfacial = getfield_or_default(sim, 'alpha_g_interfacial', NaN(size(alpha_g)));
delta_dep = sim.delta_dep;
hydrate_risk = sim.hydrate_risk;
T_hydrate = sim.T_hydrate;
subcooling = sim.subcooling;
dep_flux = sim.dep_flux;
dep_thickness = sim.dep_thickness;
solid_phi = sim.solid_phi;
agg_d = sim.agg_d;
plugging_index = sim.plugging_index;
plugging_risk = sim.plugging_risk;

gas = sim.gas;
m_dot_gas = sim.m_dot_gas;
m_dot_water0 = sim.m_dot_water0;

m_dot = m_dot_gas + m_dot_water0; % kg/s（入口总质量流，未含生成固相的微小变化）
g = const.g;

fprintf('计算完成！\n');

%% 关键结果提取

fprintf(' 计算结果 \n');


% 关键位置索引
idx_mudline = find(depth_grid <= well.mudline_depth, 1, 'last');
if isempty(idx_mudline)
    idx_mudline = 1;
end
idx_packer = find(depth_grid >= well.packer_depth, 1);

ambient = struct();
ambient.rho_seawater_kg_m3 = 1028.0; % 与外对流换热段取值保持一致（数量级校核用）
ambient.P_seawater_mudline_MPa = 0.101325 + ambient.rho_seawater_kg_m3 * const.g * well.seawater_depth / 1e6;
ambient.dP_co2_minus_seawater_mudline_MPa = P_co2(idx_mudline) - ambient.P_seawater_mudline_MPa;
ambient.mudline_balance_flag = ternary(ambient.dP_co2_minus_seawater_mudline_MPa >= 0.0, 'OK', 'UNDERBALANCE');

% 水合物质量的理论上限及折减值。
% - 理论上限：假设 CO2 与海水充分接触并完全转化为结构I CO2 水合物（仅作上限，不代表可达量）。
% 按接触效率、传质效率和可用水量折减。
solidification = struct();
solidification.enabled = true;
solidification.is_theoretical_upper_bound = true;
solidopt = getfield_or_default(options, 'solidification', struct());
t_inj_s = getfield_or_default(getfield_or_default(options, 'storage', struct()), 'injection_duration_s', model.hydrate.simulation_time_s);
t_inj_s = max(t_inj_s, 0.0);
gas_in = create_gas_mixture(injection.gas_components, injection.gas_mole_fraction, injection.gas_kij);
Mw_mix = sum(gas_in.z .* gas_in.Mw_kgmol);
mfrac_co2 = gas_in.z(gas_in.idx_co2) * gas_in.Mw_kgmol(gas_in.idx_co2) / max(Mw_mix, 1e-12);
% 使用主计算得到的质量流量
m_dot_co2_in_kg_s = m_dot_gas * mfrac_co2;
frac_hyd = hydrate_mass_fractions(model.hydrate);

solidification.duration_s = t_inj_s;
solidification.m_dot_co2_in_kg_s = m_dot_co2_in_kg_s;
% 分阶段注入（ramp）下的累计注入：默认把 MRST 的软启动视为工程注入排量逐步扩大。
% - 对于线性 ramp（0->1）：等效注入时长 t_eff = t - 0.5*t_ramp（t>=t_ramp）；若 t<t_ramp，t_eff=0.5*t^2/t_ramp
t_eff_s = t_inj_s;
try
    ctrl = getfield_or_default(getfield_or_default(getfield_or_default(options, 'storage', struct()), 'mrst', struct()), 'injection_control', struct());
    ctrl_enabled = isstruct(ctrl) && getfield_or_default(ctrl, 'enabled', false);
    mode = '';
    if isfield(ctrl, 'mode'); mode = lower(strtrim(char(ctrl.mode))); end
    t_ramp_s = getfield_or_default(ctrl, 't_ramp_s', 0.0);
    if ctrl_enabled && strcmp(mode, 'flux') && isfinite(t_ramp_s) && t_ramp_s > 0 && isfinite(t_inj_s) && t_inj_s > 0
        if t_inj_s >= t_ramp_s
            t_eff_s = max(t_inj_s - 0.5*t_ramp_s, 0.0);
        else
            t_eff_s = max(0.5 * (t_inj_s^2) / t_ramp_s, 0.0);
        end
    end
catch
end
solidification.duration_effective_s = t_eff_s;
solidification.M_co2_injected_kg = m_dot_co2_in_kg_s * t_eff_s;
solidification.M_hydrate_max_kg = solidification.M_co2_injected_kg / max(frac_hyd.co2, 1e-12);
solidification.M_water_required_kg = solidification.M_hydrate_max_kg * frac_hyd.water;
solidification.ratio_hydrate_to_co2 = 1 / max(frac_hyd.co2, 1e-12);
solidification.ratio_water_to_co2 = frac_hyd.water / max(frac_hyd.co2, 1e-12);

% 受限修正：默认系数=1（即退化为理论上限），由 options.solidification 覆盖
solidification.contact_efficiency = min(max(getfield_or_default(solidopt, 'contact_efficiency', 1.0), 0.0), 1.0);
solidification.mass_transfer_efficiency = min(max(getfield_or_default(solidopt, 'mass_transfer_efficiency', 1.0), 0.0), 1.0);
solidification.available_water_kg = getfield_or_default(solidopt, 'available_water_kg', Inf); % Inf 表示不受供水限制
solidification.effective_fraction = min(max(solidification.contact_efficiency * solidification.mass_transfer_efficiency, 0.0), 1.0);

solidification.M_co2_effective_kg = solidification.M_co2_injected_kg * solidification.effective_fraction;
solidification.M_hydrate_est_kg = solidification.M_co2_effective_kg / max(frac_hyd.co2, 1e-12);
solidification.M_water_required_est_kg = solidification.M_hydrate_est_kg * frac_hyd.water;
solidification.M_hydrate_limited_kg = solidification.M_hydrate_est_kg;
solidification.water_limited = false;
if isfinite(solidification.available_water_kg) && solidification.available_water_kg >= 0
    M_hyd_from_water = solidification.available_water_kg / max(frac_hyd.water, 1e-12);
    if M_hyd_from_water < solidification.M_hydrate_limited_kg
        solidification.M_hydrate_limited_kg = M_hyd_from_water;
        solidification.water_limited = true;
    end
end
solidification.M_co2_trapped_est_kg = solidification.M_hydrate_limited_kg * frac_hyd.co2;

% 泥线界面温度取海底水温；井底另算地层条件下的平衡温度。
[solidification.Teq_mudline_C, solidification.risk_mudline, solidification.subcool_mudline_C] = hydrate_equilibrium_and_risk(...
    P_co2(idx_mudline), T_env(idx_mudline), 1.0, thermal.seawater_salinity, model.hydrate, []);
[solidification.Teq_bottom_C, solidification.risk_bottom, solidification.subcool_bottom_C] = hydrate_equilibrium_and_risk(...
    P_co2(end), T_env(end), 1.0, thermal.seawater_salinity, model.hydrate, []);

% 井口
fprintf('\n【井口】(%.0fm)\n', depth_grid(1));
fprintf('  温度: %.2f°C\n', T_co2(1));
fprintf('  压力: %.2f MPa (表压%.2f MPa)\n', P_co2(1), P_co2(1)-0.101325);
fprintf('  密度: %.1f kg/m3\n', rho_co2(1));
fprintf('  相态: %s\n', phase{1});

% 泥线
fprintf('\n【泥线】(%.0fm, 水深%.0fm)\n', well.mudline_depth, well.seawater_depth);
fprintf('  CO2温度: %.2f°C (环境%.2f°C, 温差%.2f°C)\n', ...
    T_co2(idx_mudline), T_env(idx_mudline), T_co2(idx_mudline)-T_env(idx_mudline));
fprintf('  压力: %.2f MPa\n', P_co2(idx_mudline));
fprintf('  环境静水压: %.2f MPa (ΔP=%.2f MPa) [%s]\n', ambient.P_seawater_mudline_MPa, ambient.dP_co2_minus_seawater_mudline_MPa, ambient.mudline_balance_flag);
fprintf('  相态: %s\n', phase{idx_mudline});
if hydrate_risk(idx_mudline) > 0
    fprintf('  [RISK] 水合物风险: 等级%d, 过冷度%.2f°C\n', hydrate_risk(idx_mudline), subcooling(idx_mudline));
end

% 封隔器
fprintf('\n【RTTS封隔器】(%.0fm)\n', well.packer_depth);
fprintf('  温度: %.2f°C\n', T_co2(idx_packer));
fprintf('  压力: %.2f MPa\n', P_co2(idx_packer));
fprintf('  相态: %s\n', phase{idx_packer});

% 井底
fprintf('\n【井底】(%.0fm)\n', well.total_depth);
fprintf('  CO2温度: %.2f°C (环境%.2f°C)\n', T_co2(end), T_env(end));
fprintf('  实际压力: %.2f MPa\n', P_co2(end));
target_bh_MPa = getfield_or_default(injection, 'bh_target_pressure', NaN);
if ~isfinite(target_bh_MPa)
    target_bh_MPa = getfield_or_default(injection, 'formation_pressure', NaN);
end
fprintf('  目标压力: %.2f MPa\n', target_bh_MPa);
fprintf('  地层孔隙压力: %.2f MPa\n', injection.formation_pressure);
pressure_err_pct = abs(P_co2(end)-target_bh_MPa)/max(target_bh_MPa, 1e-12)*100;
pressure_err_flag = ternary(pressure_err_pct <= 1.0, 'OK', 'WARN');
fprintf('  压力误差: %.2f%% [%s]\n', pressure_err_pct, pressure_err_flag);
if isfinite(P_bh_max_MPa)
    dP_to_frac = P_bh_max_MPa - P_co2(end);
    if dP_to_frac >= 0
        fprintf('  破裂上限: %.2f MPa (离破裂%.2f MPa) [OK]\n', P_bh_max_MPa, dP_to_frac);
    else
        fprintf('  破裂上限: %.2f MPa (超破裂%.2f MPa) [WARN]\n', P_bh_max_MPa, -dP_to_frac);
    end
end
fprintf('  密度: %.1f kg/m3\n', rho_co2(end));
fprintf('  相态: %s\n', phase{end});

fprintf('\n【深海水合物固化（理论上限/受限修正）】\n');
fprintf('  运行时间: %.1f 天 | 累计CO2≈%.1f t\n', solidification.duration_s/86400, solidification.M_co2_injected_kg/1000);
fprintf('  泥线: T_env=%.2f°C, P=%.2fMPa, T_eq=%.2f°C, ΔT=%.2f°C\n', ...
    T_env(idx_mudline), P_co2(idx_mudline), solidification.Teq_mudline_C, solidification.subcool_mudline_C);
fprintf('  理论上限: 水合物≈%.1f t (需水≈%.1f t, H2O/CO2≈%.2f)\n', ...
    solidification.M_hydrate_max_kg/1000, solidification.M_water_required_kg/1000, solidification.ratio_water_to_co2);
fprintf('  受限修正: η_contact=%.2f, η_mt=%.2f, 水供给=%s -> 水合物≈%.1f t, CO2固化≈%.1f t%s\n', ...
    solidification.contact_efficiency, solidification.mass_transfer_efficiency, ...
    ternary(isfinite(solidification.available_water_kg), sprintf('%.1f t', solidification.available_water_kg/1000), 'Inf'), ...
    solidification.M_hydrate_limited_kg/1000, solidification.M_co2_trapped_est_kg/1000, ...
    ternary(solidification.water_limited, ' [WATER-LIMIT]', ''));

% 温度变化统计
fprintf('\n【温度演化】\n');
fprintf('  井口->井底: %.2f°C -> %.2f°C (变化%.2f°C)\n', ...
    T_co2(1), T_co2(end), T_co2(end)-T_co2(1));
fprintf('  最低温度: %.2f°C @ %.0fm\n', min(T_co2), depth_grid(find(T_co2==min(T_co2),1)));
fprintf('  最大温度: %.2f°C @ %.0fm\n', max(T_co2), depth_grid(find(T_co2==max(T_co2),1)));

% 相态统计
phases_unique = unique(phase);
fprintf('\n【相态分布】\n');
for p = 1:length(phases_unique)
    phase_name = phases_unique{p};
    idx_phase = find(strcmp(phase, phase_name));
    if ~isempty(idx_phase)
        fprintf('  %s: %.0f-%.0fm (长度%.0fm, 占比%.1f%%)\n', ...
            phase_name, depth_grid(idx_phase(1)), depth_grid(idx_phase(end)), ...
            depth_grid(idx_phase(end))-depth_grid(idx_phase(1)), ...
            length(idx_phase)/n_nodes*100);
    end
end

% 水合物风险评估
hydrate_zones = find(hydrate_risk > 0);
if ~isempty(hydrate_zones)
    fprintf('\n【水合物风险区域】\n');
    fprintf('  风险深度: %.0f-%.0fm (长度%.0fm)\n', ...
        depth_grid(hydrate_zones(1)), depth_grid(hydrate_zones(end)), ...
        depth_grid(hydrate_zones(end))-depth_grid(hydrate_zones(1)));
    fprintf('  最大过冷度: %.2f°C\n', max(subcooling));
    fprintf('  建议: 优先控制含水量、提高注入温度或加强保温/伴热，降低水合物生成风险\n');
else
    fprintf('\n【水合物风险】安全\n');
end

% 流动特性
fprintf('\n【流动特性】\n');
fprintf('  体积流量: %.2f m³/h | 质量流量: %.3f kg/s (%.2f t/h)\n', injection.volume_flow_rate, m_dot, m_dot*3.6);
fprintf('  平均流速: %.2f m/s\n', mean(v_co2));
fprintf('  平均雷诺数: %.0f (%s)\n', mean(Re_co2), ...
    ternary(mean(Re_co2)>2300, '湍流', '层流'));
fprintf('  总压降: %.2f MPa (重力%.2f + 摩阻%.2f)\n', ...
    P_co2(end)-P_co2(1), sum(rho_co2(1:end-1)*g*dz)/1e6, ...
    P_co2(end)-P_co2(1)-sum(rho_co2(1:end-1)*g*dz)/1e6);

% 传热分析
fprintf('\n【传热分析】\n');
fprintf('  平均传热系数: %.2f W/(m2*K)\n', mean(U_total(1:end-1)));
fprintf('  总传热量: %.2f kW\n', sum(Q_heat)/1000);
fprintf('  平均传热功率: %.2f W/m\n', sum(Q_heat)/(well.total_depth));

%% 套管强度
fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                    安全性评估                              ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n');

% 9-5/8"套管抗内压强度（API 5CT标准）
% 假设N80钢级，壁厚约40.5 lb/ft
P_burst = 31.5;  % 套管额定爆破压（MPa）
P_max_casing = max(P_co2(depth_grid <= well.casing_depth));

fprintf('\n【套管强度】\n');
fprintf('  最大内压: %.2f MPa (深度%.0fm)\n', P_max_casing, ...
    depth_grid(find(P_co2==P_max_casing,1)));
fprintf('  额定抗内压: %.2f MPa\n', P_burst);
fprintf('  安全系数: %.2f %s\n', P_burst/P_max_casing, ...
    ternary(P_burst/P_max_casing>1.5, '[OK]', '[WARN]'));

% 温度应力
dT_max = max(abs(T_co2 - T_env));
fprintf('\n【温度应力】\n');
fprintf('  最大温差: %.2f°C\n', dT_max);
if dT_max > 50
    fprintf('  [WARN] 建议考虑热应力分析\n');
else
    fprintf('  [OK] 温差在安全范围内\n');
end

%% 数据导出
fprintf('\n[导出] 保存结果数据...\n');

results = struct();
results.depth = depth_grid;
results.T_co2 = T_co2;
results.T_env = T_env;
results.P_co2 = P_co2;
results.rho_co2 = rho_co2;
results.mu_co2 = mu_co2;
results.cp_co2 = cp_co2;
results.k_co2 = k_co2;
results.phase = phase;
results.v_co2 = v_co2;
results.Re_co2 = Re_co2;
results.dpdz_fric_Pa_m = dpdz_fric_Pa_m;
results.tau_w_Pa = tau_w_Pa;
results.flow_regime = flow_regime;
results.v_g_phase = v_g_phase;
results.v_l_phase = v_l_phase;
results.v_rel_gl = v_rel_gl;
results.a_i_m2_m3 = a_i_m2_m3;
results.tau_i_Pa = tau_i_Pa;
results.K_i_N_s_m4 = K_i_N_s_m4;
results.Q_heat = Q_heat;
results.Q_hydrate_rxn = Q_hydrate_rxn;
results.h_mix = sim.h_mix;
results.h_balance_err = sim.h_balance_err;
results.U_total = U_total;
results.h_inner = h_inner;
results.hydrate_risk = hydrate_risk;
results.T_hydrate = T_hydrate;
results.subcooling = subcooling;
results.mu_JT = mu_JT;
results.beta_vapor = beta_vapor;
results.alpha_g = alpha_g;
results.alpha_g_drift = alpha_g_drift;
results.alpha_g_interfacial = alpha_g_interfacial;
results.delta_dep = delta_dep;
results.dep_flux = dep_flux;
results.dep_flux_net = sim.dep_flux_net;
results.ero_flux = sim.ero_flux;
results.dep_thickness = dep_thickness;
results.solid_phi = solid_phi;
results.agg_d = agg_d;
results.plugging_index = plugging_index;
results.plugging_risk = plugging_risk;
results.m_dot_gas_free = sim.m_dot_gas_free;
results.z_co2 = sim.z_co2;
results.m_dot_water = sim.m_dot_water;
results.m_dot_hydrate_susp = sim.m_dot_hydrate_susp;
results.model = model;
results.environment = environment;
if isstruct(env_profile) && isfield(env_profile, 'ok') && env_profile.ok
    results.environment_profile = env_profile;
end
results.ambient = ambient;
results.solidification = solidification;

if options.save_mat
    save('CO2_injection_results_enhanced.mat', 'results', 'well', 'injection', 'thermal');
end

% 导出关键数据为Excel
if options.export_excel
    export_xlsx = 'CO2_injection_data.xlsx';
    T_export = table(depth_grid, T_co2, P_co2, rho_co2, mu_co2*1e6, Re_co2, v_co2, dpdz_fric_Pa_m, tau_w_Pa, v_rel_gl, tau_i_Pa, a_i_m2_m3, ...
        beta_vapor, alpha_g, alpha_g_drift, alpha_g_interfacial, results.z_co2, results.m_dot_gas_free, results.m_dot_water, results.m_dot_hydrate_susp, ...
        T_hydrate, subcooling, hydrate_risk, dep_thickness, dep_flux, results.dep_flux_net, results.ero_flux, ...
        solid_phi, agg_d, plugging_index, plugging_risk, h_inner, U_total, Q_heat, Q_hydrate_rxn, results.h_mix, results.h_balance_err, ...
        'VariableNames', {'Depth_m', 'Temperature_C', 'Pressure_MPa', 'Density_kg_m3', 'Viscosity_uPas', 'Re', 'Velocity_m_s', ...
        'dpdz_fric_Pa_m', 'tau_w_Pa', 'v_rel_gl_m_s', 'tau_i_Pa', 'a_i_m2_m3', 'BetaVapor', 'AlphaGas', 'AlphaGas_Drift', 'AlphaGas_Interfacial', 'Z_CO2', 'MdotGasFree_kg_s', 'MdotWater_kg_s', 'MdotHydrateSusp_kg_s', ...
        'T_HydrateEq_C', 'Subcool_C', 'HydrateRisk', 'HydrateDep_m', 'DepFlux_kg_m2_s', 'DepFluxNet_kg_m2_s', 'EroFlux_kg_m2_s', ...
        'SolidPhi', 'AggD_m', 'PlugIndex', 'PlugRisk', 'h_inner_W_m2_K', 'U_W_m2_K', 'QHeat_W', 'QHydrateRxn_W', 'h_mix_J_kg', 'h_balance_err_J_kg'});
    writetable(T_export, export_xlsx);

    % 同步导出的表格到 out
    try
        if ~isfield(out, 'export') || ~isstruct(out.export); out.export = struct(); end
        out.export.excel = struct();
        out.export.excel.file = export_xlsx;
        out.export.excel.tables = struct();
        out.export.excel.tables.Sheet1 = T_export;
        out.export.excel.sheets = {'Sheet1'};
    catch ME
        warning('[导出] Excel->工作区同步失败：%s', ME.message);
    end
end

% 导出所有可视化图的作图数据（一个工作簿，多 sheet）
if isfield(options, 'export_plot_excel') && logical(options.export_plot_excel)
    plot_xlsx = char(getfield_or_default(options, 'plot_excel_file', 'CO2_visualization_plotdata.xlsx'));
    try
        plot_export = export_visualization_plotdata_excel(plot_xlsx, options, depth_grid, sim, env_profile, well, injection, thermal, environment);
        try
            if ~isfield(out, 'export') || ~isstruct(out.export); out.export = struct(); end
            out.export.plot_excel = plot_export;
        catch ME2
            warning('[导出] plot_excel->工作区同步失败：%s', ME2.message);
        end
        fprintf('  [OK] 可视化作图数据已导出到 %s\n', plot_xlsx);
    catch ME
        warning('[导出] 可视化作图数据导出失败：%s', ME.message);
    end
end

if options.save_mat
    fprintf('  [OK] 数据已保存到 CO2_injection_results_enhanced.mat\n');
end
if options.export_excel
    fprintf('  [OK] Excel表格已导出到 CO2_injection_data.xlsx\n');
end

fprintf('计算节点: %d；计算时间: %.2f s\n', n_nodes, toc);
try
    out.options = options;
    out.model = model;
    out.well = well;
    out.injection = injection;
    out.thermal = thermal;
    if exist('environment','var'); out.environment = environment; end
    if exist('env_profile','var'); out.environment_profile = env_profile; end
    if exist('CO2','var'); out.CO2 = CO2; end
    if exist('const','var'); out.const = const; end
    if exist('sim','var'); out.sim = sim; end
    if exist('results','var'); out.results = results; end
catch ME
    warning('[输出] out 组装失败：%s', ME.message);
end

end

function sim = simulate_olga_coupled(depth_grid, dz, T_env, injection, model, well, thermal, CO2, const, simopt)
    if nargin < 10 || isempty(simopt)
        simopt = struct();
    end
    verbose = true;
    if isfield(simopt, 'verbose')
        verbose = logical(simopt.verbose);
    end
    progress_every = 500;
    if isfield(simopt, 'progress_every')
        progress_every = max(1, round(simopt.progress_every));
    end
    show_progress = verbose;
    if isfield(simopt, 'show_progress')
        show_progress = logical(simopt.show_progress);
    end
    progress_update_interval_s = 0.5;
    if isfield(simopt, 'progress_update_interval_s')
        progress_update_interval_s = max(0.05, getfield_or_default(simopt, 'progress_update_interval_s', progress_update_interval_s));
    end
    verbose_details = verbose && ~show_progress;

    n_nodes = numel(depth_grid);

    kij_in = [];
    if isfield(injection, 'gas_kij')
        kij_in = injection.gas_kij;
    end
    gas = create_gas_mixture(injection.gas_components, injection.gas_mole_fraction, kij_in);
    z0 = gas.z;

    % 体积流量 -> 质量流量换算：volume_flow_rate 被定义为入口总混相体积流量(m³/h)。
    % 先用相态(flash/forced_liquid)计算入口混相密度，再按 free_water_mass_ratio 加入自由水。
    T_inlet_K = injection.temp_inlet + 273.15;
    P_inlet_Pa = (injection.pressure_inlet + 0.101325) * 1e6;  % 表压转绝压

    Q_total_m3_s = max(getfield_or_default(injection, 'volume_flow_rate', 0.0), 0.0) / 3600.0;
    m_dot_override_kg_s = getfield_or_default(injection, 'mass_flow_rate_kg_s', NaN);
    free_water_mass_ratio = max(getfield_or_default(injection, 'free_water_mass_ratio', 0.0), 0.0); % kg_water/kg_gas

    phase_mode_in = "flash";
    if isfield(model, 'thermo') && isfield(model.thermo, 'phase_mode')
        phase_mode_in = lower(string(model.thermo.phase_mode));
    end

    flash_in = struct('beta', 0.0, 'x', z0, 'y', z0, 'K', ones(size(z0)), 'converged', true);
    use_sat_recondense_in = false;
    if isfield(model, 'thermo') && isfield(model.thermo, 'flash')
        use_sat_recondense_in = logical(getfield_or_default(model.thermo.flash, 'recondense_enabled', false));
    end
    if phase_mode_in == "liquid"
        flash_in.beta = 0.0;
    elseif use_sat_recondense_in && phase_mode_in == "flash"
        flash_in = co2_flash_saturation_recondense(T_inlet_K, P_inlet_Pa, gas, model);
    else
        flash_model = "wilson";
        if isfield(model, 'thermo') && isfield(model.thermo, 'flash_model')
            flash_model = lower(string(model.thermo.flash_model));
        end
        if flash_model == "phi_pr"
            flashopt = struct();
            if isfield(model, 'thermo') && isfield(model.thermo, 'flash')
                flashopt = model.thermo.flash;
            end
            flash_in = co2_flash_pr_phi(T_inlet_K, P_inlet_Pa, gas, [], flashopt);
        else
            flash_in = co2_flash_wilson(T_inlet_K, P_inlet_Pa, gas);
        end
    end

    % 入口气相摩尔分率
    % 可指定入口气相摩尔分率，
    % 模拟井口附近CO2部分气化的物理过程（WHP接近饱和压力时）。
    if isfield(injection, 'beta_inlet_override') && ~isempty(injection.beta_inlet_override)
        flash_in.beta = max(min(injection.beta_inlet_override, 1.0), 0.0);
    end

    beta_in = min(max(getfield_or_default(flash_in, 'beta', 0.0), 0.0), 1.0);
    x_in = z0;
    y_in = z0;
    if isfield(flash_in, 'x'); x_in = flash_in.x(:)'; end
    if isfield(flash_in, 'y'); y_in = flash_in.y(:)'; end
    x_in = x_in / max(sum(x_in), 1e-12);
    y_in = y_in / max(sum(y_in), 1e-12);

    props_l_in = co2_phase_props_pr(T_inlet_K, P_inlet_Pa, x_in, gas, 'liquid');
    props_v_in = co2_phase_props_pr(T_inlet_K, P_inlet_Pa, y_in, gas, 'vapor');
    rho_l_in = max(props_l_in.rho, 1e-9);
    rho_v_in = max(props_v_in.rho, 1e-9);

    Mw_v_in = sum(y_in .* gas.Mw_kgmol);
    Mw_l_in = sum(x_in .* gas.Mw_kgmol);
    m_v_unit = beta_in * Mw_v_in;
    m_l_unit = (1 - beta_in) * Mw_l_in;
    w_v_in = 0.0;
    if (m_v_unit + m_l_unit) > 0
        w_v_in = m_v_unit / (m_v_unit + m_l_unit);
    end
    w_v_in = min(max(w_v_in, 0.0), 1.0);

    aq_in = aqueous_props(T_inlet_K, P_inlet_Pa, thermal.seawater_salinity);
    rho_w_in = max(aq_in.rho, 1e-9);

    vol_per_kg_gas = w_v_in / rho_v_in + (1.0 - w_v_in) / rho_l_in + free_water_mass_ratio / rho_w_in; % m3/kg_gas
    vol_per_kg_gas = max(vol_per_kg_gas, 1e-12);

    m_dot_gas = Q_total_m3_s / vol_per_kg_gas; % 气体质量流量（kg/s）
    m_dot_water0 = m_dot_gas * free_water_mass_ratio;

    m_dot_override_kg_s = getfield_or_default(injection, 'mass_flow_rate_kg_s', NaN);
    if isfinite(m_dot_override_kg_s) && m_dot_override_kg_s > 0
        m_dot_gas = m_dot_override_kg_s / (1.0 + free_water_mass_ratio);
    end
    m_dot_water0 = m_dot_gas * free_water_mass_ratio;

    Mw_mix0 = sum(z0 .* gas.Mw_kgmol);
    n_dot_total0 = m_dot_gas / max(Mw_mix0, 1e-12); % mol/s
    n_dot_comp0 = n_dot_total0 * z0;                % mol/s（各组分）

    delta_dep_init = zeros(n_nodes, 1); % m（可用于循环注入的“续算初值”）
    if isfield(simopt, 'initial_delta_dep') && ~isempty(simopt.initial_delta_dep)
        dep0 = simopt.initial_delta_dep;
        if isscalar(dep0)
            delta_dep_init(:) = max(dep0, 0.0);
        else
            dep0 = dep0(:);
            if numel(dep0) ~= n_nodes
                error('simopt.initial_delta_dep 尺寸不匹配：期望 %d×1。', n_nodes);
            end
            delta_dep_init = max(dep0, 0.0);
        end
    end
    delta_dep = delta_dep_init;

    % 沉积-水力/传热几何回写开关：
    % - couple_to_flow=true ：沉积厚度进入有效内径 D_inner（可产生强非线性回耦，适用于堵塞/压降放大评估）。
    % - couple_to_flow=false：保持通径不变，仅输出沉积厚度。
    dep_couple_to_flow = true;
    try
        if ~(isfield(model, 'hydrate') && isstruct(model.hydrate) && isfield(model.hydrate, 'enabled') && model.hydrate.enabled)
            dep_couple_to_flow = false;
        elseif ~(isfield(model.hydrate, 'deposition') && isstruct(model.hydrate.deposition) ...
                && isfield(model.hydrate.deposition, 'enabled') && model.hydrate.deposition.enabled)
            dep_couple_to_flow = false;
        elseif isfield(model.hydrate.deposition, 'couple_to_flow')
            dep_couple_to_flow = logical(model.hydrate.deposition.couple_to_flow);
        end
    catch
    end
    if verbose && ~dep_couple_to_flow && n_nodes > 1
        fprintf('  [沉积] couple_to_flow=false：沉积不回写通径（避免浅层出现不现实的局部节流尖峰）。\n');
    end

    progress_line_active = false;
    t_progress_start = [];
    t_progress_last = [];
    if show_progress
        t_progress_start = tic;
        t_progress_last = tic;
    end

    n_seg = max(n_nodes - 1, 1);
    n_iter_max = max(1, round(model.iteration.max_iter));
    if ~dep_couple_to_flow
        n_iter_max = 1;
    end
    n_total = n_seg * n_iter_max;

    for iter = 1:n_iter_max
        if verbose_details
            fprintf('  [耦合迭代] %d/%d\n', iter, model.iteration.max_iter);
        end

        T_co2 = zeros(n_nodes, 1);
        P_co2 = zeros(n_nodes, 1);
        rho_co2 = zeros(n_nodes, 1);
        mu_co2 = zeros(n_nodes, 1);
        cp_co2 = zeros(n_nodes, 1);
        k_co2 = zeros(n_nodes, 1);

        phase = cell(n_nodes, 1);
        beta_vapor = zeros(n_nodes, 1);
        alpha_g = zeros(n_nodes, 1);
        alpha_g_drift = NaN(n_nodes, 1);
        alpha_g_interfacial = NaN(n_nodes, 1);
        holdup_model_used = cell(n_nodes, 1);
        holdup_model_used(:) = {''};

        v_co2 = zeros(n_nodes, 1);
        Re_co2 = zeros(n_nodes, 1);
        f_D = zeros(n_nodes, 1);
        dpdz_fric_Pa_m_profile = zeros(n_nodes, 1);
        tau_w_Pa_profile = zeros(n_nodes, 1);
        flow_regime = cell(n_nodes, 1);
        flow_regime(:) = {''};
        v_g_phase = zeros(n_nodes, 1);
        v_l_phase = zeros(n_nodes, 1);
        v_rel_gl = zeros(n_nodes, 1);
        a_i_m2_m3 = zeros(n_nodes, 1);
        tau_i_Pa = zeros(n_nodes, 1);
        K_i_N_s_m4 = zeros(n_nodes, 1);
        Re_i = zeros(n_nodes, 1);
        Cd_i = zeros(n_nodes, 1);
        f_i = zeros(n_nodes, 1);
        d_i_m = zeros(n_nodes, 1);

        Q_heat = zeros(n_nodes, 1);
        U_total = zeros(n_nodes, 1);
        h_inner = NaN(n_nodes, 1);          % 实际用于计算的内对流换热系数（W/m2-K）
        ht_regime = cell(n_nodes, 1);       % 两相内换热流型标识（用于校核）
        ht_regime(:) = {''};
        Q_hydrate_rxn = zeros(n_nodes, 1); % 反应放热功率（W，按dz段），正值表示向流体放热
        h_mix = zeros(n_nodes, 1);         % 混相比焓（J/kg），用于焓法能量闭合
        h_balance_err = zeros(n_nodes, 1); % 焓闭合残差（J/kg）：h_actual - h_target（仅enthalpy模式填充）

        hydrate_risk = zeros(n_nodes, 1);
        T_hydrate = zeros(n_nodes, 1);
        subcooling = zeros(n_nodes, 1);
        % 水合物平衡有效压力诊断：分压 vs PR 逸度（用于高压/混相一致性评估）
        hydrate_P_eff_pp_MPa = NaN(n_nodes, 1);     % P_eff= P*y（或P），单位MPa
        hydrate_P_eff_fug_MPa = NaN(n_nodes, 1);    % P_eff= f，单位MPa
        hydrate_P_eff_used_MPa = NaN(n_nodes, 1);   % 实际用于平衡线的 P_eff
        hydrate_Teq_pp_C = NaN(n_nodes, 1);
        hydrate_Teq_fug_C = NaN(n_nodes, 1);
        hydrate_use_fugacity = false(n_nodes, 1);
        dep_flux = zeros(n_nodes, 1);
        dep_flux_net = zeros(n_nodes, 1);
        ero_flux = zeros(n_nodes, 1);
        dep_thickness = delta_dep_init;
        solid_phi = zeros(n_nodes, 1);
        agg_d = zeros(n_nodes, 1);
        plugging_index = zeros(n_nodes, 1);
        plugging_risk = zeros(n_nodes, 1);

        m_dot_water = m_dot_water0;
        m_dot_hydrate_susp = 0.0;
        n_dot_comp = n_dot_comp0;
        gas.z = z0;

        m_dot_gas_free_profile = zeros(n_nodes, 1);
        z_co2_profile = zeros(n_nodes, 1);
        m_dot_water_profile = zeros(n_nodes, 1);
        m_dot_hydrate_susp_profile = zeros(n_nodes, 1);

        T_co2(1) = injection.temp_inlet;
        P_co2(1) = injection.pressure_inlet + 0.101325;
        agg_d(1) = model.hydrate.aggregation.d0_m;

        flash_guess = [];
        energy_model = "cp";
        if isfield(model, 'thermo') && isfield(model.thermo, 'energy_model')
            energy_model = lower(string(model.thermo.energy_model));
        end
        T_ref_K = 273.15;
        if isfield(model, 'thermo') && isfield(model.thermo, 'enthalpy_ref_C')
            T_ref_K = 273.15 + model.thermo.enthalpy_ref_C;
        end
        cp_hydrate = 2000.0;
        if isfield(model, 'thermo') && isfield(model.thermo, 'enthalpy')
            cp_hydrate = getfield_or_default(model.thermo.enthalpy, 'hydrate_cp_J_kg_K', 2000.0);
        end

        phase_mode = "flash";
        if isfield(model, 'thermo') && isfield(model.thermo, 'phase_mode')
            phase_mode = lower(string(model.thermo.phase_mode));
        end
        for i = 1:n_nodes
            z = depth_grid(i);

            D_inner = max(well.coil_tubing_id - 2*delta_dep(i), model.hydraulics.min_inner_diameter_m);
            A_flow = pi * (D_inner/2)^2;
            perimeter = pi * D_inner;

            T_K = T_co2(i) + 273.15;
            P_Pa = P_co2(i) * 1e6;

            % 组分随程更新：水合物生成会消耗CO2，导致混合气组分变化（对闪蒸/密度/水力有影响）
            n_dot_total = sum(n_dot_comp);
            if n_dot_total <= 0
                error('自由相气体摩尔流量耗尽，无法继续计算（检查注入量/水合物参数）。');
            end
            gas.z = n_dot_comp / n_dot_total;

            m_dot_gas_free = sum(n_dot_comp .* gas.Mw_kgmol);
            m_dot_gas_free_profile(i) = m_dot_gas_free;
            z_co2_profile(i) = gas.z(gas.idx_co2);
            m_dot_water_profile(i) = m_dot_water;
            m_dot_hydrate_susp_profile(i) = m_dot_hydrate_susp;

            use_sat_recondense = false;
            if isfield(model, 'thermo') && isfield(model.thermo, 'flash')
                use_sat_recondense = logical(getfield_or_default(model.thermo.flash, 'recondense_enabled', false));
            end

            if phase_mode == "liquid"
                flash = struct('beta', 0.0, 'x', gas.z, 'y', gas.z, 'ok', true, 'method', 'forced_liquid');
            elseif phase_mode == "match_bhp_twophase"
                % 匹配实测BHP模式：全程按液态计算flash（确保热力学性质正确），
                % 但在井筒上部（深度 < twophase_transition_depth）对密度施加折减因子，
                % 模拟近井口CO2部分闪蒸导致的有效密度降低。
                % 密度折减在后续 rho_co2 赋值处实施，此处仍用液态flash。
                flash = struct('beta', 0.0, 'x', gas.z, 'y', gas.z, 'ok', true, 'method', 'match_bhp_liquid');
            elseif use_sat_recondense && phase_mode == "flash"
                flash = co2_flash_saturation_recondense(T_K, P_Pa, gas, model);
                flash_guess = flash;
            else
                flash_model = 'wilson';
                if isfield(model, 'thermo') && isfield(model.thermo, 'flash_model')
                    flash_model = lower(string(model.thermo.flash_model));
                end
                if flash_model == "phi_pr"
                    flashopt = struct();
                    if isfield(model, 'thermo') && isfield(model.thermo, 'flash')
                        flashopt = model.thermo.flash;
                    end
                    flash = co2_flash_pr_phi(T_K, P_Pa, gas, flash_guess, flashopt);
                else
                    flash = co2_flash_wilson(T_K, P_Pa, gas);
                end
                flash_guess = flash;
            end
            beta_vapor(i) = flash.beta;
            phase{i} = phase_label(beta_vapor(i), P_co2(i), T_co2(i), CO2);

            props_l = co2_phase_props_pr(T_K, P_Pa, flash.x, gas, 'liquid');
            if phase_mode == "liquid"
                props_v = props_l;
            else
                props_v = co2_phase_props_pr(T_K, P_Pa, flash.y, gas, 'vapor');
            end

            aq = aqueous_props(T_K, P_Pa, thermal.seawater_salinity);

            % Wilson 闪蒸给出气相摩尔分率，按相分子量换算质量流率。
            Mw_v = sum(flash.y .* gas.Mw_kgmol);
            Mw_l = sum(flash.x .* gas.Mw_kgmol);
            n_dot_v = beta_vapor(i) * n_dot_total;
            n_dot_l = (1 - beta_vapor(i)) * n_dot_total;
            m_dot_co2_v = n_dot_v * Mw_v;
            m_dot_co2_l = n_dot_l * Mw_l;
            m_dot_co2 = m_dot_co2_v + m_dot_co2_l;

            Vdot_g = m_dot_co2_v / max(props_v.rho, 1e-9);
            Vdot_lc = m_dot_co2_l / max(props_l.rho, 1e-9);
            Vdot_aq = m_dot_water / max(aq.rho, 1e-9);
            Vdot_s = m_dot_hydrate_susp / max(model.hydrate.rho_hydrate, 1e-9);
            Vdot_liq = Vdot_lc + Vdot_aq + Vdot_s;
            Vdot_total = Vdot_g + Vdot_liq;

            j_g = Vdot_g / max(A_flow, 1e-12);
            j = Vdot_total / max(A_flow, 1e-12);
            j_l = max(j - j_g, 0.0);

            rho_liq_eff = props_l.rho;
            if Vdot_liq > 0
                m_dot_liq_eff = m_dot_co2_l + m_dot_water + m_dot_hydrate_susp;
                rho_liq_eff = m_dot_liq_eff / max(Vdot_liq, 1e-12);
            end

            alpha_g_noslip = ternary(Vdot_total > 0, Vdot_g / Vdot_total, 0.0);
            alpha_g(i) = alpha_g_noslip;
            holdup_model_used{i} = 'noslip';
            if isfield(model.hydraulics, 'two_fluid') && isfield(model.hydraulics.two_fluid, 'enabled') && model.hydraulics.two_fluid.enabled
                if Vdot_g > 0 && Vdot_liq > 0
                    two_fluid_cfg = model.hydraulics.two_fluid;
                    holdup_model = lower(string(getfield_or_default(two_fluid_cfg, 'holdup_model', 'drift_flux')));
                    if holdup_model == "momentum_balance" || holdup_model == "interfacial"
                        holdup_model = "interfacial_terminal";
                    end

                    reg_model = lower(string(getfield_or_default(two_fluid_cfg, 'regime_model', 'none')));
                    reg_here = "two_phase";
                    if reg_model == "simple"
                        [C0_eff, Vgj_override_m_s, reg_here] = drift_flux_params_vertical_simple(j_g, j, rho_liq_eff, props_v.rho, D_inner, const.g);
                        two_fluid_cfg.C0 = C0_eff;
                        two_fluid_cfg.Vgj_override_m_s = Vgj_override_m_s; % NaN 表示用默认Harmathy
                        flow_regime{i} = char(reg_here);
                    end

                    [alpha_df, ~, ~, ~] = void_fraction_drift_flux(j_g, j, rho_liq_eff, props_v.rho, two_fluid_cfg, const.g);
                    alpha_g_drift(i) = alpha_df;

                    if holdup_model == "drift_flux"
                        alpha_g(i) = alpha_df;
                        holdup_model_used{i} = 'drift_flux';
                    elseif holdup_model == "interfacial_terminal"
                        aq_frac_liq_mu = ternary(Vdot_liq > 0, Vdot_aq / Vdot_liq, 0.0);
                        mu_l_eff = max((1-aq_frac_liq_mu)*props_l.mu + aq_frac_liq_mu*aq.mu, 1e-12);
                        alpha_guess = ternary(isfinite(alpha_df), alpha_df, alpha_g_noslip);
                        [alpha_itf, holinfo] = void_fraction_interfacial_terminal(...
                            j_g, j_l, props_v.rho, rho_liq_eff, mu_l_eff, D_inner, const.g, two_fluid_cfg, alpha_guess, reg_here);
                        if holinfo.ok
                            alpha_g(i) = alpha_itf;
                            alpha_g_interfacial(i) = alpha_itf;
                            holdup_model_used{i} = holinfo.model;
                        else
                            alpha_g(i) = alpha_df;
                            holdup_model_used{i} = ternary(isfield(holinfo, 'model') && ~isempty(holinfo.model), holinfo.model, 'drift_flux_fallback');
                        end
                    else
                        alpha_g(i) = alpha_df;
                        holdup_model_used{i} = 'drift_flux_default';
                    end
                end
            end
            alpha_g(i) = min(max(alpha_g(i), 0.0), 1.0);
            if isempty(flow_regime{i})
                flow_regime{i} = ternary((Vdot_g > 0 && Vdot_liq > 0), 'two_phase', 'single');
            end

            alpha_l_total = min(max(1.0 - alpha_g(i), 0.0), 1.0);
            solid_phi(i) = ternary(Vdot_liq > 0, Vdot_s / Vdot_liq, 0.0); % 固相在液相中的体积分数（更贴近浆体流变）
            aq_frac_liq = ternary(Vdot_liq > 0, Vdot_aq / Vdot_liq, 0.0); % 水相在液相中的体积分数（用于液相混合物性）

            m_dot_total = m_dot_co2 + m_dot_water + m_dot_hydrate_susp;
            rho_co2(i) = ternary(Vdot_total > 0, m_dot_total / Vdot_total, props_l.rho);

            % match_bhp_twophase模式：在两相过渡段内对密度施加折减
            % 物理含义：井口附近CO2部分闪蒸，有效密度低于纯液态
            if phase_mode == "match_bhp_twophase"
                L_trans = injection.twophase_transition_depth;
                if z < L_trans && L_trans > 0
                    % 密度从井口的 rho_gas 线性过渡到 L_trans 处的 rho_liquid
                    rho_gas_ref = 60.0;  % 气态CO2参考密度 (kg/m³)
                    frac = z / L_trans;  % 0(井口) -> 1(过渡深度)
                    rho_eff = rho_gas_ref + (rho_co2(i) - rho_gas_ref) * frac;
                    rho_co2(i) = rho_eff;
                end
            end
            mu_co2(i) = mixture_viscosity(props_v.mu, props_l.mu, aq.mu, aq_frac_liq, solid_phi(i));
            cp_co2(i) = mixture_cp(props_v.cp, props_l.cp, aq.cp, m_dot_co2_v, m_dot_co2_l, m_dot_water, m_dot_hydrate_susp);
            k_co2(i) = mixture_k(props_v.k, props_l.k, aq.k, aq_frac_liq, alpha_g(i));

            % 相间拖曳/界面剪切（诊断量输出）：仅两相区定义相间相对速度与界面剪切
            v_rel_gl(i) = 0.0;
            if Vdot_g > 0 && Vdot_liq <= 0
                v_g_phase(i) = j;
                v_l_phase(i) = 0.0;
            elseif Vdot_liq > 0 && Vdot_g <= 0
                v_g_phase(i) = 0.0;
                v_l_phase(i) = j;
            else
                v_g_phase(i) = 0.0;
                v_l_phase(i) = 0.0;
            end

            if Vdot_g > 0 && Vdot_liq > 0 && alpha_g(i) > 0 && alpha_g(i) < 1
                v_g_phase(i) = j_g / max(alpha_g(i), 1e-12);
                alpha_l_here = min(max(1.0 - alpha_g(i), 0.0), 1.0);
                v_l_phase(i) = j_l / max(alpha_l_here, 1e-12);
                v_rel_gl(i) = v_g_phase(i) - v_l_phase(i);

                use_interfacial = false;
                if isfield(model.hydraulics, 'two_fluid') && isfield(model.hydraulics.two_fluid, 'interfacial') && isfield(model.hydraulics.two_fluid.interfacial, 'enabled')
                    use_interfacial = logical(model.hydraulics.two_fluid.interfacial.enabled);
                end
                if use_interfacial
                    reg = flow_regime{i};
                    if isempty(reg)
                        reg = 'two_phase';
                    end
                    mu_l_eff = max(mu_co2(i), 1e-12);
                    itf = interfacial_drag_shear_two_fluid(D_inner, j_g, j_l, alpha_g(i), ...
                        props_v.rho, props_v.mu, rho_liq_eff, mu_l_eff, reg, model.hydraulics.two_fluid);
                    a_i_m2_m3(i) = itf.a_i_m2_m3;
                    tau_i_Pa(i) = itf.tau_i_Pa;
                    K_i_N_s_m4(i) = itf.K_i_N_s_m4;
                    Re_i(i) = itf.Re_i;
                    Cd_i(i) = ternary(isfinite(itf.Cd), itf.Cd, 0.0);
                    f_i(i) = ternary(isfinite(itf.f_i), itf.f_i, 0.0);
                    d_i_m(i) = ternary(isfinite(itf.d_char_m), itf.d_char_m, 0.0);
                end
            end

            m_dot_in_seg = m_dot_total;
            Hdot_in_seg = NaN;
            if energy_model == "enthalpy"
                h_v = pr_eos_enthalpy_pr(T_K, P_Pa, flash.y, gas, 'vapor', T_ref_K);
                h_l = pr_eos_enthalpy_pr(T_K, P_Pa, flash.x, gas, 'liquid', T_ref_K);
                h_aq = aq.cp * (T_K - T_ref_K);
                h_s = cp_hydrate * (T_K - T_ref_K);
                h_mix(i) = (m_dot_co2_v*h_v + m_dot_co2_l*h_l + m_dot_water*h_aq + m_dot_hydrate_susp*h_s) / max(m_dot_in_seg, 1e-12);
                Hdot_in_seg = m_dot_in_seg * h_mix(i);
            end

            v_co2(i) = m_dot_total / max(rho_co2(i) * A_flow, 1e-12);
            Re_co2(i) = rho_co2(i) * v_co2(i) * D_inner / max(mu_co2(i), 1e-12);

            % 两流体摩阻闭合（OLGA风格：holdup+两相摩阻倍增）；输出等效Darcy系数与壁面剪切应力
            dpdz_fric_Pa_m = 0.0;
            tau_w_Pa = 0.0;
            if Vdot_total > 0
                j_l = max(j - j_g, 0.0);
                use_two_fluid_fric = false;
                if isfield(model.hydraulics, 'two_fluid') && isfield(model.hydraulics.two_fluid, 'enabled') && model.hydraulics.two_fluid.enabled
                    use_two_fluid_fric = (Vdot_g > 0 && Vdot_liq > 0);
                end

                if use_two_fluid_fric
                    mu_l_eff = max(mu_co2(i), 1e-12); % 液相（含水/固）有效黏度：用上面计算的浆体黏度
                    [dpdz_fric_Pa_m, ~] = two_phase_friction_gradient_lockhart_martinelli(...
                        j_g, j_l, props_v.rho, rho_liq_eff, props_v.mu, mu_l_eff, well.tubing_roughness, D_inner);
                else
                    f_sp = darcy_friction_factor(Re_co2(i), well.tubing_roughness, D_inner);
                    dpdz_fric_Pa_m = f_sp * (rho_co2(i) * v_co2(i)^2 / 2) / max(D_inner, 1e-12);
                end
                dpdz_fric_Pa_m = max(dpdz_fric_Pa_m, 0.0);
                tau_w_Pa = dpdz_fric_Pa_m * D_inner / 4;
            end
            dpdz_fric_Pa_m_profile(i) = dpdz_fric_Pa_m;
            tau_w_Pa_profile(i) = tau_w_Pa;

            denom_dyn = (rho_co2(i) * v_co2(i)^2 / 2);
            if denom_dyn > 0 && dpdz_fric_Pa_m > 0
                f_D(i) = (dpdz_fric_Pa_m * D_inner) / denom_dyn;
            else
                f_D(i) = 0.0;
            end

            % 本段净沉积质量流（kg/s），用于能量守恒（沉积质量从流动中移除）
            m_dot_dep_net_seg = 0.0;

            if model.hydrate.enabled
                y_co2_gas = max(flash.y(gas.idx_co2), 0.0);

                % 平衡有效压力：分压 vs PR 逸度
                %
                % - 平衡线拟合（hyd.equilibrium.lnp_poly_coeffs）默认是ln(P)相关式；若直接把逸度 f 当成 P 代入，会导致平衡线显著偏移。
                % - 因此：默认只用分压（P*y）做平衡线与风险判断；逸度仅用于诊断对比。
                % 用逸度驱动时，平衡线系数须按逸度标定。
                hyd = model.hydrate;
                use_pp = logical(getfield_or_default(hyd, 'use_partial_pressure', true));
                use_fugacity = logical(getfield_or_default(hyd, 'use_fugacity', false)); % 仅控制“是否用于平衡线/风险”，不影响诊断计算
                P_eff_pp_MPa = max(P_co2(i) * ternary(use_pp, y_co2_gas, 1.0), 0.0);

                fug_cfg = struct();
                if isfield(hyd, 'fugacity') && isstruct(hyd.fugacity)
                    fug_cfg = hyd.fugacity;
                end
                fug_mode = lower(string(getfield_or_default(fug_cfg, 'mode', ...
                    ternary(use_fugacity, 'on', 'off'))));
                P_thr = getfield_or_default(fug_cfg, 'P_threshold_MPa', 5.0);
                beta_min = getfield_or_default(fug_cfg, 'beta_mixed_min', 0.01);
                beta_max = getfield_or_default(fug_cfg, 'beta_mixed_max', 0.99);

                P_eff_fug_MPa = NaN;
                need_fug_calc = (fug_mode == "on") || (fug_mode == "auto" && (P_co2(i) >= P_thr || (beta_vapor(i) > beta_min && beta_vapor(i) < beta_max)));
                if need_fug_calc
                    y_eff = ternary(use_pp, y_co2_gas, 1.0);
                    x_eff = ternary(use_pp, max(flash.x(gas.idx_co2), 0.0), 1.0);
                    try
                        [lnphi_v, ~] = pr_eos_lnphi_pr(T_K, P_Pa, flash.y, gas, 'vapor');
                        [lnphi_l, ~] = pr_eos_lnphi_pr(T_K, P_Pa, flash.x, gas, 'liquid');
                        phi_v = exp(lnphi_v(gas.idx_co2));
                        phi_l = exp(lnphi_l(gas.idx_co2));
                        f_v = max(y_eff * phi_v * P_co2(i), 0.0); % MPa
                        f_l = max(x_eff * phi_l * P_co2(i), 0.0); % MPa
                        if beta_vapor(i) >= 0.99
                            P_eff_fug_MPa = f_v;
                        elseif beta_vapor(i) <= 0.01
                            P_eff_fug_MPa = f_l;
                        else
                            P_eff_fug_MPa = max(f_v, f_l);
                        end
                    catch
                        P_eff_fug_MPa = NaN;
                    end
                end

                use_fug_in_eq = (fug_mode == "on") || (fug_mode == "auto" && use_fugacity);
                use_fug_here = use_fug_in_eq && isfinite(P_eff_fug_MPa);
                P_eff_used_MPa = ternary(use_fug_here, P_eff_fug_MPa, P_eff_pp_MPa);

                hydrate_P_eff_pp_MPa(i) = P_eff_pp_MPa;
                hydrate_P_eff_fug_MPa(i) = P_eff_fug_MPa;
                hydrate_P_eff_used_MPa(i) = P_eff_used_MPa;
                hydrate_use_fugacity(i) = use_fug_here;

                % 以 P_eff_used 计算风险等级；并额外保存 pp/fug 下的 Teq 对比（误差诊断）
                [T_hydrate(i), risk_eq, subcooling(i)] = hydrate_equilibrium_and_risk(...
                    P_co2(i), T_co2(i), 1.0, thermal.seawater_salinity, hyd, P_eff_used_MPa);
                hydrate_Teq_pp_C(i) = hydrate_equilibrium_and_risk(...
                    P_co2(i), T_co2(i), 1.0, thermal.seawater_salinity, hyd, P_eff_pp_MPa);
                if isfinite(P_eff_fug_MPa)
                    hydrate_Teq_fug_C(i) = hydrate_equilibrium_and_risk(...
                        P_co2(i), T_co2(i), 1.0, thermal.seawater_salinity, hyd, P_eff_fug_MPa);
                end

                if m_dot_water > 0
                    hydrate_risk(i) = risk_eq;
                    if subcooling(i) > 0
                        tau_w = tau_w_Pa;
                        dep_frac = 0.0;
                        if isfield(model.hydrate, 'deposition') && isfield(model.hydrate.deposition, 'enabled') && model.hydrate.deposition.enabled
                            dep_frac = hydrate_deposition_fraction(tau_w, model.hydrate.deposition);
                        end

                        % 按水相体积分数估算水与 CO₂ 的接触面积。
                        alpha_aq_total = ternary(Vdot_total > 0, Vdot_aq / Vdot_total, 0.0);
                        a_int = hydrate_interfacial_area(alpha_aq_total, model.hydrate.drop_diameter_m);
                        k_mult_kin = hydrate_kinetics_multiplier(model.hydrate);
                        r_vol = model.hydrate.k_surf * k_mult_kin * a_int * (subcooling(i)^model.hydrate.n);

                        seg_vol = A_flow * dz;
                        m_dot_hyd_gen_seg = r_vol * seg_vol;
                        max_conv = 1.0;
                        if isfield(model.hydrate, 'max_water_conversion')
                            max_conv = min(max(model.hydrate.max_water_conversion, 0.0), 1.0);
                        end
                        m_dot_hyd_gen_seg = min(m_dot_hyd_gen_seg, max_hydrate_from_water(m_dot_water * max_conv, model.hydrate));
                        m_dot_co2_avail = n_dot_comp(gas.idx_co2) * gas.Mw_kgmol(gas.idx_co2);
                        m_dot_hyd_gen_seg = min(m_dot_hyd_gen_seg, max_hydrate_from_co2(m_dot_co2_avail, model.hydrate));

                        if isfield(model.hydrate, 'include_heat_of_formation') && model.hydrate.include_heat_of_formation
                            dH = 0.0;
                            if isfield(model.hydrate, 'deltaH_formation_J_per_kg')
                                dH = model.hydrate.deltaH_formation_J_per_kg;
                            end
                            Q_hydrate_rxn(i) = -dH * m_dot_hyd_gen_seg; % dH<0 放热 => Q>0
                        end

                        m_dot_water_consumed = hydrate_water_consumption(m_dot_hyd_gen_seg, model.hydrate);
                        m_dot_water = max(0.0, m_dot_water - m_dot_water_consumed);

                        % 沉积增长（可选：冲刷/再悬浮，使沉积趋于饱和而非无限增长）
                        t_op = max(getfield_or_default(simopt, 'injection_duration_s', model.hydrate.simulation_time_s), 0.0);
                        rho_h = max(model.hydrate.rho_hydrate, 1e-9);

                        m_dot_dep_attempt = dep_frac * m_dot_hyd_gen_seg;
                        dep_flux(i) = ternary(perimeter > 0, m_dot_dep_attempt / (perimeter * dz), 0.0); % kg/(m2·s)
                        dep_flux_net(i) = dep_flux(i);
                        ero_flux(i) = 0.0;

                        delta0 = delta_dep_init(i);
                        delta_end = delta0 + dep_flux(i) * t_op / rho_h;
                        m_dot_dep_net_seg = m_dot_dep_attempt;
                        m_dot_resusp_seg = 0.0;

                        if isfield(model.hydrate, 'deposition') && isfield(model.hydrate.deposition, 'erosion') && model.hydrate.deposition.erosion.enabled
                            ero = model.hydrate.deposition.erosion;
                            tau_c = max(getfield_or_default(ero, 'tau_crit_Pa', 80), 0.0);
                            k_ero = max(getfield_or_default(ero, 'k_ero_kg_m3_s_Pa', 0.0), 0.0);
                            dTau = max(tau_w - tau_c, 0.0);

                            if dTau > 0 && k_ero > 0 && dep_flux(i) > 0 && t_op > 0
                                % 线性冲刷：j_ero = k*(τ-τc)*δ，dδ/dt = (j_dep - j_ero)/ρ
                                delta_eq = dep_flux(i) / (k_ero * dTau);
                                lambda = (k_ero * dTau) / rho_h; % 1/s
                                delta_end = delta_eq + (delta0 - delta_eq) * exp(-lambda * t_op);
                                delta_end = max(delta_end, 0.0);

                                wall_area = perimeter * dz;
                                m_dot_net_seg = rho_h * (delta_end - delta0) * wall_area / max(t_op, 1e-12);
                                if m_dot_net_seg >= 0
                                    m_dot_dep_net_seg = min(m_dot_net_seg, m_dot_dep_attempt);
                                    m_dot_resusp_seg = max(m_dot_dep_attempt - m_dot_dep_net_seg, 0.0);
                                    dep_flux_net(i) = ternary(wall_area > 0, m_dot_dep_net_seg / wall_area, 0.0);
                                    ero_flux(i) = ternary(wall_area > 0, m_dot_resusp_seg / wall_area, 0.0);
                                else
                                    m_dot_dep_net_seg = 0.0;
                                    m_dot_resusp_seg = m_dot_dep_attempt - m_dot_net_seg; % net<0 => resusp > attempt
                                    dep_flux_net(i) = 0.0;
                                    ero_flux(i) = ternary(wall_area > 0, m_dot_resusp_seg / wall_area, 0.0);
                                end
                            end
                        end

                        dep_thickness(i) = delta_end;

                        % 悬浮相：未沉积部分 + 再悬浮部分
                        resusp_frac = 1.0;
                        if isfield(model.hydrate, 'deposition') && isfield(model.hydrate.deposition, 'erosion') && model.hydrate.deposition.erosion.enabled
                            resusp_frac = min(max(getfield_or_default(model.hydrate.deposition.erosion, 'resuspension_fraction', 1.0), 0.0), 1.0);
                        end
                        m_dot_hydrate_susp = m_dot_hydrate_susp + (1 - dep_frac) * m_dot_hyd_gen_seg + resusp_frac * m_dot_resusp_seg;

                        % 更新局部固相体积分数（用于聚集与堵塞评估）
                        Vdot_s_now = m_dot_hydrate_susp / rho_h;
                        Vdot_aq_now = m_dot_water / max(aq.rho, 1e-9);
                        Vdot_liq_now = Vdot_lc + Vdot_aq_now + Vdot_s_now;
                        solid_phi(i) = ternary(Vdot_liq_now > 0, Vdot_s_now / Vdot_liq_now, 0.0);
                        % CO2消耗：从自由相组分流量中扣除（沉积/悬浮都已计入水合物质量）
                        m_dot_co2_consumed = hydrate_co2_consumption(m_dot_hyd_gen_seg, model.hydrate);
                        n_dot_co2_consumed = m_dot_co2_consumed / max(gas.Mw_kgmol(gas.idx_co2), 1e-12);
                        n_dot_co2_consumed = min(n_dot_co2_consumed, n_dot_comp(gas.idx_co2));
                        n_dot_comp(gas.idx_co2) = max(0.0, n_dot_comp(gas.idx_co2) - n_dot_co2_consumed);

                        if model.hydrate.aggregation.enabled
                            dt = dz / max(v_co2(i), 1e-6);
                            agg_d(i) = max(agg_d(i), model.hydrate.aggregation.d0_m);
                            agg_d_next = hydrate_aggregation_step(agg_d(i), solid_phi(i), dt, model.hydrate.aggregation);
                            if i < n_nodes
                                agg_d(i+1) = max(agg_d(i+1), agg_d_next);
                            end
                        end

                        max_dep = 0.5 * (well.coil_tubing_id - model.hydraulics.min_inner_diameter_m);
                        idx_dep = dep_thickness(i) / max(max_dep, 1e-6);
                        idx_solid = solid_phi(i) / max(model.hydrate.aggregation.phi_crit, 1e-6);
                        idx_agg = agg_d(i) / max(model.hydrate.aggregation.d_crit_m, 1e-9);
                        plugging_index(i) = max([idx_dep, idx_solid, idx_agg]);

                        if plugging_index(i) < 0.5
                            plugging_risk(i) = 0;
                        elseif plugging_index(i) < 1.0
                            plugging_risk(i) = 1;
                        elseif plugging_index(i) < 2.0
                            plugging_risk(i) = 2;
                        else
                            plugging_risk(i) = 3;
                        end

                        hydrate_risk(i) = max(hydrate_risk(i), plugging_risk(i));
                    end
                end
            end

            if i < n_nodes
                T_env_avg = (T_env(i) + T_env(i+1)) / 2;

                % 内对流换热：两相时使用流型相关液膜/干涸闭合（OLGA风格的最小可用实现）
                h_inner_override = [];
                ht_regime{i} = '';
                if D_inner > 0 && Vdot_total > 0
                    inner_model = lower(string(getfield_or_default(model.heat_transfer, 'inner_model', 'mixture')));
                    use_tp_ht = (inner_model == "two_phase_regime");
                    if use_tp_ht
                        tp_enabled = true;
                        if isfield(model.heat_transfer, 'inner') && isfield(model.heat_transfer.inner, 'two_phase') && isfield(model.heat_transfer.inner.two_phase, 'enabled')
                            tp_enabled = logical(model.heat_transfer.inner.two_phase.enabled);
                        end
                        use_tp_ht = tp_enabled;
                    end

                    if use_tp_ht && Vdot_g > 0 && Vdot_liq > 0
                        j_l = max(j - j_g, 0.0);

                        k_liq_eff = (1-aq_frac_liq)*props_l.k + aq_frac_liq*aq.k;
                        mu_liq_eff = max(mu_co2(i), 1e-12);
                        m_dot_liq = m_dot_co2_l + m_dot_water + m_dot_hydrate_susp;
                        cp_liq_eff = props_l.cp;
                        if m_dot_liq > 0
                            cp_liq_eff = (m_dot_co2_l*props_l.cp + m_dot_water*aq.cp + m_dot_hydrate_susp*cp_hydrate) / max(m_dot_liq, 1e-12);
                        end

                        props_g_ht = struct('rho', props_v.rho, 'mu', props_v.mu, 'cp', props_v.cp, 'k', props_v.k);
                        props_l_ht = struct('rho', rho_liq_eff, 'mu', mu_liq_eff, 'cp', cp_liq_eff, 'k', k_liq_eff);

                        regime_hint = "unknown";
                        try
                            [~, ~, regime_hint] = drift_flux_params_vertical_simple(j_g, j, rho_liq_eff, props_v.rho, D_inner, const.g);
                        catch
                            regime_hint = "unknown";
                        end

                        [h_inner_override, htinfo] = inner_htc_two_phase_regime(D_inner, j_g, j_l, alpha_g(i), props_g_ht, props_l_ht, regime_hint, model.heat_transfer);
                        ht_regime{i} = char(htinfo.regime);
                        h_inner_override = max(h_inner_override, 1.0);
                    else
                        ht_regime{i} = ternary((Vdot_g > 0 && Vdot_liq > 0), 'two_phase_mixture', 'single');
                    end
                end

                [U_total(i), h_comp] = calculate_multilayer_heat_transfer(...
                    z, T_co2(i), T_env_avg, P_co2(i), v_co2(i), Re_co2(i), ...
                    well, CO2, thermal, k_co2(i), mu_co2(i), cp_co2(i), model, delta_dep(i), h_inner_override);
                h_inner(i) = h_comp.h_inner;

                A_heat = pi * D_inner * dz;
                Q = U_total(i) * A_heat * (T_env_avg - T_co2(i));
                Q_heat(i) = Q;

                dP_gravity = rho_co2(i) * const.g * dz / 1e6;
                dP_fric_loss = dpdz_fric_Pa_m * dz / 1e6;
                P_co2(i+1) = P_co2(i) + dP_gravity - dP_fric_loss;
                if isfield(model, 'hydraulics')
                    min_pressure_abs_MPa = getfield_or_default(model.hydraulics, 'min_pressure_abs_MPa', 0.12);
                    if isfinite(min_pressure_abs_MPa)
                        P_co2(i+1) = max(P_co2(i+1), min_pressure_abs_MPa);
                    end
                end

                dH_fric = (dP_fric_loss * 1e6) / max(rho_co2(i), 1e-9);
                m_dot_out_seg = total_massflow_from_components(n_dot_comp, gas, m_dot_water, m_dot_hydrate_susp);
                h_out_target = NaN;

                if energy_model == "enthalpy"
                    % 段能量守恒（忽略动能变化）：
                    % H_out + H_deposit = H_in + Q + Q_rxn
                    % 摩阻通过压降影响焓平衡，不另加摩擦热源项。
                    h_dep = cp_hydrate * (T_K - T_ref_K);
                    Hdot_dep = m_dot_dep_net_seg * h_dep;

                    % 井筒竖直向下流动：需要包含重力势能项（≈ṁ*g*dz），否则会在高压段产生系统性温度偏差。
                    Wdot_gravity = m_dot_in_seg * const.g * dz;
                    Hdot_out_target = Hdot_in_seg + Q + Q_hydrate_rxn(i) + Wdot_gravity - Hdot_dep;
                    h_out_target = Hdot_out_target / max(m_dot_out_seg, 1e-12);

                    dH_guess = (Q + Q_hydrate_rxn(i) + Wdot_gravity - Hdot_dep) / max(m_dot_in_seg, 1e-12);
                    dT_guess = dH_guess / max(cp_co2(i), 1e-9);
                    T_guess = T_co2(i) + dT_guess;

                    T_co2(i+1) = solve_temperature_from_enthalpy(T_guess, P_co2(i+1)*1e6, h_out_target, ...
                        n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_guess, T_ref_K);

                    [h_out_actual, ~] = mixture_specific_enthalpy(T_co2(i+1)+273.15, P_co2(i+1)*1e6, ...
                        n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_guess, T_ref_K);
                    h_balance_err(i) = h_out_actual - h_out_target;
                else
                    dH_heat = (Q + Q_hydrate_rxn(i)) / max(m_dot_in_seg, 1e-12);
                    dT = (dH_heat + dH_fric) / max(cp_co2(i), 1e-9);
                    T_co2(i+1) = T_co2(i) + dT;
                end

                % 可选：加速度项原型（默认关闭）。为避免重算闪蒸，暂仅在CO2单相区使用EOS估算ρ_out。
                if isfield(model, 'hydraulics') && isfield(model.hydraulics, 'include_acceleration') && model.hydraulics.include_acceleration
                    if (beta_vapor(i) <= 0.01 || beta_vapor(i) >= 0.99) && D_inner > 0
                        D_inner_next = max(well.coil_tubing_id - 2*delta_dep(i+1), model.hydraulics.min_inner_diameter_m);
                        A_flow_next = pi * (D_inner_next/2)^2;

                        n_dot_total_out = sum(n_dot_comp);
                        if n_dot_total_out > 0
                            gas_out = gas;
                            gas_out.z = n_dot_comp / n_dot_total_out;
                            m_dot_co2_out = sum(n_dot_comp .* gas_out.Mw_kgmol);

                            phase_out = ternary(beta_vapor(i) >= 0.99, 'vapor', 'liquid');
                            props_co2_out = co2_phase_props_pr(T_co2(i+1)+273.15, P_co2(i+1)*1e6, gas_out.z, gas_out, phase_out);
                            aq_out = aqueous_props(T_co2(i+1)+273.15, P_co2(i+1)*1e6, thermal.seawater_salinity);

                            Vdot_co2_out = m_dot_co2_out / max(props_co2_out.rho, 1e-9);
                            Vdot_aq_out = m_dot_water / max(aq_out.rho, 1e-9);
                            Vdot_s_out = m_dot_hydrate_susp / max(model.hydrate.rho_hydrate, 1e-9);
                            Vdot_total_out = Vdot_co2_out + Vdot_aq_out + Vdot_s_out;

                            rho_out = ternary(Vdot_total_out > 0, m_dot_out_seg / Vdot_total_out, rho_co2(i));
                            v_out = ternary(rho_out > 0 && A_flow_next > 0, m_dot_out_seg / (rho_out * A_flow_next), v_co2(i));

                            dP_accel_MPa = (rho_co2(i) * v_co2(i)^2 - rho_out * v_out^2) / 1e6;
                            if isfinite(dP_accel_MPa)
                                P_co2(i+1) = P_co2(i+1) + dP_accel_MPa;

                                if energy_model == "enthalpy" && isfinite(h_out_target)
                                    T_co2(i+1) = solve_temperature_from_enthalpy(T_co2(i+1), P_co2(i+1)*1e6, h_out_target, ...
                                        n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_guess, T_ref_K);
                                    [h_out_actual, ~] = mixture_specific_enthalpy(T_co2(i+1)+273.15, P_co2(i+1)*1e6, ...
                                        n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_guess, T_ref_K);
                                    h_balance_err(i) = h_out_actual - h_out_target;
                                end
                            end
                        end
                    end
                end
            end

            if show_progress
                seg_i = min(i, n_seg);
                n_done = (iter-1) * n_seg + seg_i;
                pct = 100.0 * (n_done / max(n_total, 1));

                do_print = false;
                if ~progress_line_active
                    do_print = true;
                elseif toc(t_progress_last) >= progress_update_interval_s
                    do_print = true;
                elseif i == n_nodes && iter == n_iter_max
                    do_print = true;
                end

                if do_print
                    elapsed_s = toc(t_progress_start);
                    eta_s = NaN;
                    if pct > 0
                        eta_s = elapsed_s * (100.0 / pct - 1.0);
                    end
                    if isfinite(eta_s) && eta_s >= 0
                        fprintf('\r[主计算] 总进度 %6.2f%% | 耦合迭代 %d/%d | 节点 %d/%d | 剩余 %.0f s', ...
                            pct, iter, n_iter_max, seg_i, n_seg, eta_s);
                    else
                        fprintf('\r[主计算] 总进度 %6.2f%% | 耦合迭代 %d/%d | 节点 %d/%d', ...
                            pct, iter, n_iter_max, seg_i, n_seg);
                    end
                    progress_line_active = true;
                    t_progress_last = tic;
                end
            end

            if verbose_details && (mod(i, progress_every) == 0 || i == n_nodes)
                fprintf('  %4.0fm: T=%6.2f°C, P=%6.2f MPa, ρ=%6.1f kg/m3, β_v=%.2f', ...
                    z, T_co2(i), P_co2(i), rho_co2(i), beta_vapor(i));
                if hydrate_risk(i) > 0
                    fprintf('  [水合物R=%d, ΔT=%.2f°C]', hydrate_risk(i), subcooling(i));
                end
                fprintf('\n');
                if usejava('desktop')
                    graphics_keepalive();
                end
            end
        end

        dep_thickness = min(dep_thickness, 0.49*well.coil_tubing_id);
        delta_err = 0.0;
        if dep_couple_to_flow
            omega = min(max(model.iteration.relaxation, 0.0), 1.0);
            delta_new = dep_thickness;
            delta_err = max(abs(delta_new - delta_dep));
            delta_dep = (1-omega)*delta_dep + omega*delta_new;
        else
            % 不做几何回写：保持通径为初始状态（通常为清管）。
            delta_dep = delta_dep_init;
        end

        if verbose
            if show_progress && progress_line_active
                fprintf('\n');
                progress_line_active = false;
            end
            fprintf('    收敛检查：max|Δδ| = %.3e m\n', delta_err);
        end
        if energy_model == "enthalpy" && verbose
            fprintf('    焓闭合残差: max|d_h| = %.3g J/kg\n', max(abs(h_balance_err)));
        end
        if delta_err < model.iteration.delta_tol_m
            break;
        end
    end

    if show_progress
        if progress_line_active
            fprintf('\r[主计算] 总进度 100.00%% | 完成\n');
        else
            fprintf('[主计算] 总进度 100.00%% | 完成\n');
        end
    end

    Q_heat(end) = ternary(n_nodes > 1, Q_heat(end-1), 0.0);
    U_total(end) = ternary(n_nodes > 1, U_total(end-1), 0.0);
    h_inner(end) = ternary(n_nodes > 1, h_inner(end-1), h_inner(end));
    ht_regime{end} = ternary(n_nodes > 1, ht_regime{end-1}, 'single');
    dpdz_fric_Pa_m_profile(end) = ternary(n_nodes > 1, dpdz_fric_Pa_m_profile(end-1), 0.0);
    tau_w_Pa_profile(end) = ternary(n_nodes > 1, tau_w_Pa_profile(end-1), 0.0);
    flow_regime{end} = ternary(n_nodes > 1, flow_regime{end-1}, 'single');
    mu_JT = zeros(n_nodes, 1);

    % 水合物平衡诊断汇总：分压 vs PR 逸度
    hydrate_diag = struct('ok', false);
    try
        hyd = model.hydrate;
        fug_cfg = struct();
        if isfield(hyd, 'fugacity') && isstruct(hyd.fugacity)
            fug_cfg = hyd.fugacity;
        end
        fug_mode_cfg = lower(string(getfield_or_default(fug_cfg, 'mode', ...
            ternary(getfield_or_default(hyd, 'use_fugacity', false), 'on', 'off'))));

        idx = find(isfinite(hydrate_Teq_fug_C) & isfinite(hydrate_Teq_pp_C));
        hydrate_diag.ok = true;
        hydrate_diag.fugacity_mode = char(fug_mode_cfg);
        hydrate_diag.n_fugacity_calc = sum(isfinite(hydrate_P_eff_fug_MPa));
        hydrate_diag.n_fugacity_used = sum(hydrate_use_fugacity);
        hydrate_diag.n_compare = numel(idx);
        hydrate_diag.max_abs_Teq_diff_C = NaN;
        hydrate_diag.mean_Teq_diff_C = NaN;
        hydrate_diag.depth_max_abs_diff_m = NaN;
        if ~isempty(idx)
            dTeq = hydrate_Teq_fug_C(idx) - hydrate_Teq_pp_C(idx);
            [mx, kmax] = max(abs(dTeq));
            idx_max = idx(kmax);
            hydrate_diag.max_abs_Teq_diff_C = mx;
            hydrate_diag.mean_Teq_diff_C = mean(dTeq);
            hydrate_diag.depth_max_abs_diff_m = depth_grid(idx_max);
        end

        report = logical(getfield_or_default(fug_cfg, 'report', false));
        if report && verbose && hydrate_diag.n_compare > 0
            fprintf('[水合物平衡] Teq: fugacity(%s) vs pp -> max|Δ|=%.2f°C @ %.0fm, meanΔ=%.2f°C (N=%d)\n', ...
                hydrate_diag.fugacity_mode, hydrate_diag.max_abs_Teq_diff_C, hydrate_diag.depth_max_abs_diff_m, ...
                hydrate_diag.mean_Teq_diff_C, hydrate_diag.n_compare);
        end
    catch
        hydrate_diag = struct('ok', false);
    end

    sim = struct();
    sim.T_co2 = T_co2;
    sim.P_co2 = P_co2;
    sim.rho_co2 = rho_co2;
    sim.mu_co2 = mu_co2;
    sim.cp_co2 = cp_co2;
    sim.k_co2 = k_co2;
    sim.phase = phase;
    sim.v_co2 = v_co2;
    sim.Re_co2 = Re_co2;
    sim.Q_heat = Q_heat;
    sim.Q_hydrate_rxn = Q_hydrate_rxn;
    sim.h_mix = h_mix;
    sim.h_balance_err = h_balance_err;
    sim.U_total = U_total;
    sim.h_inner = h_inner;
    sim.ht_regime = ht_regime;
    sim.dpdz_fric_Pa_m = dpdz_fric_Pa_m_profile;
    sim.tau_w_Pa = tau_w_Pa_profile;
    sim.flow_regime = flow_regime;
    sim.v_g_phase = v_g_phase;
    sim.v_l_phase = v_l_phase;
    sim.v_rel_gl = v_rel_gl;
    sim.a_i_m2_m3 = a_i_m2_m3;
    sim.tau_i_Pa = tau_i_Pa;
    sim.K_i_N_s_m4 = K_i_N_s_m4;
    sim.Re_i = Re_i;
    sim.Cd_i = Cd_i;
    sim.f_i = f_i;
    sim.d_i_m = d_i_m;
    sim.mu_JT = mu_JT;
    sim.beta_vapor = beta_vapor;
    sim.alpha_g = alpha_g;
    sim.alpha_g_drift = alpha_g_drift;
    sim.alpha_g_interfacial = alpha_g_interfacial;
    sim.holdup_model_used = holdup_model_used;
    sim.delta_dep = delta_dep;
    sim.hydrate_risk = hydrate_risk;
    sim.T_hydrate = T_hydrate;
    sim.subcooling = subcooling;
    sim.hydrate_P_eff_pp_MPa = hydrate_P_eff_pp_MPa;
    sim.hydrate_P_eff_fug_MPa = hydrate_P_eff_fug_MPa;
    sim.hydrate_P_eff_used_MPa = hydrate_P_eff_used_MPa;
    sim.hydrate_Teq_pp_C = hydrate_Teq_pp_C;
    sim.hydrate_Teq_fug_C = hydrate_Teq_fug_C;
    sim.hydrate_use_fugacity = hydrate_use_fugacity;
    sim.hydrate_diag = hydrate_diag;
    sim.dep_flux = dep_flux;
    sim.dep_flux_net = dep_flux_net;
    sim.ero_flux = ero_flux;
    sim.dep_thickness = dep_thickness;
    sim.solid_phi = solid_phi;
    sim.agg_d = agg_d;
    sim.plugging_index = plugging_index;
    sim.plugging_risk = plugging_risk;
    sim.T_env = T_env(:);
    sim.gas = gas;
    sim.m_dot_gas = m_dot_gas;
    sim.m_dot_water0 = m_dot_water0;
    sim.m_dot_gas_free = m_dot_gas_free_profile;
    sim.z_co2 = z_co2_profile;
    sim.m_dot_water = m_dot_water_profile;
    sim.m_dot_hydrate_susp = m_dot_hydrate_susp_profile;
end

% 井口边界反求（可选）：以目标井底压力作为约束，二分法求解井口表压。
% 适用前提：在给定流量/物性/传热条件下，井底压力随井口压力单调变化（对多数井筒稳态工况成立）。
function [injection_out, sim_out] = solve_inlet_pressure_to_match_formation(depth_grid, dz, T_env, injection_in, model, well, thermal, CO2, const, simopt)
    target = getfield_or_default(injection_in, 'bh_target_pressure', injection_in.formation_pressure); % MPa（绝压）
    tol = 0.01;
    if isfield(injection_in, 'solve_tol_MPa')
        tol = max(injection_in.solve_tol_MPa, 1e-6);
    end
    max_iter = 25;
    if isfield(injection_in, 'solve_max_iter')
        max_iter = max(1, round(injection_in.solve_max_iter));
    end

    simopt_eval = simopt;
    simopt_eval.verbose = false;
    simopt_eval.progress_every = 1000000000;
    simopt_eval.show_progress = false;

    show_progress = false;
    if isfield(simopt, 'show_progress')
        show_progress = logical(simopt.show_progress);
    end
    progress_update_interval_s = 0.5;
    if isfield(simopt, 'progress_update_interval_s')
        progress_update_interval_s = max(0.05, getfield_or_default(simopt, 'progress_update_interval_s', progress_update_interval_s));
    end

    % 速度优先：井底压力对井口压力的灵敏度通常≈1，因此用误差回代的固定点迭代即可快速收敛
    step_cap = 5.0;
    if isfield(injection_in, 'solve_step_cap_MPa')
        step_cap = max(injection_in.solve_step_cap_MPa, 0.1);
    end

    p = max(injection_in.pressure_inlet, 0.0);
    best_p = p;
    best_err = inf;
    best_sim = [];

    t_progress_start = [];
    t_progress_last = [];
    progress_line_active = false;
    if show_progress
        t_progress_start = tic;
        t_progress_last = tic;
        fprintf('[边界匹配] 总进度  0.00%%');
        progress_line_active = true;
    end

    for it = 1:max_iter
        inj = injection_in;
        inj.pressure_inlet = max(p, 0.0);

        sim_tmp = simulate_olga_coupled(depth_grid, dz, T_env, inj, model, well, thermal, CO2, const, simopt_eval);
        if ~isstruct(sim_tmp) || ~isfield(sim_tmp, 'P_co2') || isempty(sim_tmp.P_co2) || ~isfinite(sim_tmp.P_co2(end))
            best_sim = sim_tmp;
            best_p = p;
            best_err = NaN;
            break;
        end

        err = sim_tmp.P_co2(end) - target;
        if abs(err) < abs(best_err)
            best_err = err;
            best_p = p;
            best_sim = sim_tmp;
        end

        if show_progress && (toc(t_progress_last) >= progress_update_interval_s || it == 1 || abs(err) <= tol || it == max_iter)
            pct = 100.0 * (it / max(max_iter, 1));
            elapsed_s = toc(t_progress_start);
            eta_s = NaN;
            if it > 0
                eta_s = (max_iter - it) * (elapsed_s / it);
            end
            if isfinite(eta_s) && eta_s >= 0
                fprintf('\r[边界匹配] 总进度 %6.2f%% | iter %d/%d | err %.3f MPa | 剩余 %.0f s', pct, it, max_iter, err, eta_s);
            else
                fprintf('\r[边界匹配] 总进度 %6.2f%% | iter %d/%d | err %.3f MPa', pct, it, max_iter, err);
            end
            t_progress_last = tic;
            progress_line_active = true;
        end

        if abs(err) <= tol
            break;
        end

        step = -err; % MPa
        step = min(max(step, -step_cap), step_cap);
        p_next = max(p + step, 0.0);
        if abs(p_next - p) < 1e-6
            break;
        end
        p = p_next;
    end

    if show_progress && progress_line_active
        fprintf('\r[边界匹配] 总进度 100.00%% | 完成\n');
    end

    if ~isfinite(best_err)
        warning('入口压力反求失败（返回 NaN），保持原 pressure_inlet=%.3f MPa。', max(injection_in.pressure_inlet, 0.0));
        injection_out = injection_in;
        sim_out = best_sim;
        return;
    end

    if abs(best_err) > tol
        warning('入口压力反求未收敛：best_err=%.3f MPa（tol=%.3f, iter<=%d），将使用最优迭代点 pressure_inlet=%.3f MPa。', ...
            best_err, tol, max_iter, best_p);
    end

    injection_out = injection_in;
    injection_out.pressure_inlet = best_p;
    sim_out = best_sim;
    fprintf('  [边界匹配] 井口表压=%.3f MPa, 井底误差=%.3f MPa（tol=%.3f, iter<=%d）\n', best_p, best_err, tol, max_iter);
end

% 近井壳层回耦：外迭代将 ΔP_shell(t_end) 计入井口表压反求（保持流量不变）
function cfg = setup_graphics_runtime(options)
    % 图形后端与导出设置。
    cfg = struct();
    if nargin < 1 || ~isstruct(options)
        try
            setappdata(0, 'CO2_HYDRATE_GRAPHICS', cfg);
        catch
        end
        return;
    end

    if isfield(options, 'graphics') && isstruct(options.graphics)
        cfg = options.graphics;
    end

    if ~isfield(cfg, 'enabled'); cfg.enabled = true; end
    if ~isfield(cfg, 'opengl'); cfg.opengl = 'auto'; end
    if ~isfield(cfg, 'renderer'); cfg.renderer = 'opengl'; end
    if ~isfield(cfg, 'figure_pause_s'); cfg.figure_pause_s = 0.05; end
    if ~isfield(cfg, 'retry_pause_s'); cfg.retry_pause_s = 0.20; end
    if ~isfield(cfg, 'keepalive_pause_s'); cfg.keepalive_pause_s = 0.02; end
    if ~isfield(cfg, 'auto_export'); cfg.auto_export = false; end
    if ~isfield(cfg, 'export_dir'); cfg.export_dir = 'figures'; end
    if ~isfield(cfg, 'export_format'); cfg.export_format = 'png'; end
    if ~isfield(cfg, 'export_dpi'); cfg.export_dpi = 200; end
    if ~isfield(cfg, 'export_close'); cfg.export_close = false; end
    % 多子图中隐藏坐标区工具栏。
    if ~isfield(cfg, 'disable_axes_toolbar'); cfg.disable_axes_toolbar = true; end

    cfg.is_desktop = usejava('desktop');
    cfg.export_dir_abs = resolve_export_dir(cfg.export_dir, cfg.auto_export);

    try
        setappdata(0, 'CO2_HYDRATE_GRAPHICS', cfg);
        setappdata(0, 'CO2_HYDRATE_EXPORTS', {});
    catch
    end

    if ~cfg.enabled
        return;
    end

    need_graphics = false;
    if isfield(options, 'enable_plots') && options.enable_plots
        need_graphics = true;
    end
    if isfield(cfg, 'auto_export') && cfg.auto_export
        need_graphics = true;
    end
    if ~need_graphics
        return;
    end

    try
        set(groot, 'DefaultFigureVisible', 'on');
    catch
    end
    if cfg.is_desktop
        try
            set(groot, 'DefaultFigureWindowStyle', 'normal');
        catch
        end
    end

    % 中文显示与避免乱码（禁用 latex 解释器；统一字体）
    font_name = char(getfield_or_default(cfg, 'font_name', 'auto'));
    if strcmpi(strtrim(font_name), 'auto')
        if ispc
            font_name = 'Microsoft YaHei';
        else
            font_name = 'Helvetica';
        end
    end
    try; set(groot, 'DefaultAxesFontName', font_name); catch; end
    try; set(groot, 'DefaultTextFontName', font_name); catch; end
    try; set(groot, 'DefaultLegendFontName', font_name); catch; end
    try; set(groot, 'DefaultColorbarFontName', font_name); catch; end
    try; set(groot, 'DefaultTextInterpreter', 'none'); catch; end
    try; set(groot, 'DefaultLegendInterpreter', 'none'); catch; end
    try; set(groot, 'DefaultAxesTickLabelInterpreter', 'none'); catch; end

    try
        mode = lower(strtrim(string(cfg.opengl)));
        if mode == "software"
            ws = warning('query', 'MATLAB:opengl:switchRendererDeprecated');
            warning('off', 'MATLAB:opengl:switchRendererDeprecated');
            opengl('software');
            warning(ws.state, ws.identifier);
        elseif mode == "hardware"
            ws = warning('query', 'MATLAB:opengl:switchRendererDeprecated');
            warning('off', 'MATLAB:opengl:switchRendererDeprecated');
            opengl('hardware');
            warning(ws.state, ws.identifier);
        end
    catch
    end
end

function dir_abs = resolve_export_dir(export_dir, create_dir)
    if nargin < 2
        create_dir = false;
    end
    if nargin < 1 || isempty(export_dir)
        export_dir = 'figures';
    end
    if isstring(export_dir)
        export_dir = char(export_dir);
    end
    export_dir = strtrim(export_dir);
    if isempty(export_dir)
        export_dir = 'figures';
    end

    base_dir = fileparts(mfilename('fullpath'));
    if is_absolute_path(export_dir)
        dir_abs = export_dir;
    else
        dir_abs = fullfile(base_dir, export_dir);
    end

    if create_dir
        try
            if ~isfolder(dir_abs)
                mkdir(dir_abs);
            end
        catch ME
            warning('[导出] 无法创建导出目录：%s', ME.message);
        end
    end
end

function tf = is_absolute_path(p)
    tf = false;
    if nargin < 1 || isempty(p)
        return;
    end
    if isstring(p)
        p = char(p);
    end
    if ~ischar(p)
        return;
    end
    p = strtrim(p);
    if isempty(p)
        return;
    end

    if ispc
        if startsWith(p, '\\') || startsWith(p, '//')
            tf = true;
            return;
        end
        if numel(p) >= 3 && isletter(p(1)) && p(2) == ':' && (p(3) == '\' || p(3) == '/')
            tf = true;
        end
        return;
    end

    tf = startsWith(p, '/');
end

function maybe_export_figure(fig, basename)
    if nargin < 1 || isempty(fig) || ~ishandle(fig)
        return;
    end
    try
        if ~strcmp(get(fig, 'Type'), 'figure')
            return;
        end
    catch
        return;
    end

    cfg = [];
    try
        cfg = getappdata(0, 'CO2_HYDRATE_GRAPHICS');
    catch
    end
    if ~isstruct(cfg) || ~isfield(cfg, 'auto_export') || ~cfg.auto_export
        return;
    end

    export_dir = '';
    if isfield(cfg, 'export_dir_abs') && ~isempty(cfg.export_dir_abs)
        export_dir = cfg.export_dir_abs;
    end
    if isempty(export_dir)
        export_dir = resolve_export_dir(getfield_or_default(cfg, 'export_dir', 'figures'), true);
    end

    if nargin < 2 || isempty(basename)
        try
            basename = get(fig, 'Name');
        catch
            basename = '';
        end
    end
    basename = sanitize_filename(basename);
    if isempty(basename)
        try
            basename = sprintf('Figure%d', get(fig, 'Number'));
        catch
            basename = 'Figure';
        end
    end

    fmt = lower(strtrim(char(getfield_or_default(cfg, 'export_format', 'png'))));
    dpi = getfield_or_default(cfg, 'export_dpi', 200);

    out_path = fullfile(export_dir, sprintf('%s.%s', basename, fmt));
    try
        % 导出时隐藏坐标区工具栏。
        tb_state = {};
        try
            axs = findall(fig, 'Type', 'axes');
            for iax = 1:numel(axs)
                ax = axs(iax);
                if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar, 'Visible')
                    tb_state{end+1,1} = ax; %#ok<AGROW>
                    tb_state{end,2} = ax.Toolbar.Visible;
                    ax.Toolbar.Visible = 'off';
                end
            end
        catch
            tb_state = {};
        end

        if has_graphics_function('exportgraphics')
            if strcmp(fmt, 'png') || strcmp(fmt, 'jpg') || strcmp(fmt, 'jpeg') || strcmp(fmt, 'tif') || strcmp(fmt, 'tiff')
                exportgraphics(fig, out_path, 'Resolution', dpi);
            else
                exportgraphics(fig, out_path);
            end
        else
            print(fig, out_path, ['-d', fmt], sprintf('-r%d', dpi));
        end

        % 恢复坐标区工具栏可见性
        for k = 1:size(tb_state, 1)
            try
                ax = tb_state{k,1};
                vis = tb_state{k,2};
                if isvalid(ax) && isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar, 'Visible')
                    ax.Toolbar.Visible = vis;
                end
            catch
            end
        end

        exports = {};
        try
            exports = getappdata(0, 'CO2_HYDRATE_EXPORTS');
        catch
        end
        if isempty(exports)
            exports = {};
        end
        exports{end+1} = out_path; %#ok<AGROW>
        try
            setappdata(0, 'CO2_HYDRATE_EXPORTS', exports);
        catch
        end
    catch ME
        warning('[导出] 图像导出失败：%s', ME.message);
    end

    if isfield(cfg, 'export_close') && cfg.export_close
        try
            close(fig);
        catch
        end
    end
end

function tf = has_graphics_function(fname)
    % `exist(...,'file')` 在 p-code/内建函数场景下会返回 6/5，不应只判断 ==2。
    tf = false;
    if nargin < 1 || isempty(fname)
        return;
    end
    try
        tf = (exist(fname, 'file') ~= 0) || (exist(fname, 'builtin') ~= 0);
    catch
        tf = false;
    end
end

function tf = axes_contains_raster_content(ax)
    tf = false;
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end
    try
        tf = ~isempty(findall(ax, 'Type', 'image'));
    catch
        tf = false;
    end
end

function tf = axes_is_probably_3d(ax)
    tf = false;
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end

    try
        v = view(ax);
        if isnumeric(v) && numel(v) >= 2
            if abs(v(1)) > 1e-9 || abs(v(2) - 90) > 1e-9
                tf = true;
                return;
            end
        end
    catch
    end

    % 含非平面 ZData 的坐标区按三维处理。
    try
        ch = allchild(ax);
    catch
        ch = [];
    end
    for i = 1:numel(ch)
        h = ch(i);
        if isempty(h) || ~ishandle(h)
            continue;
        end
        z = [];
        try
            if isprop(h, 'ZData')
                z = get(h, 'ZData');
            end
        catch
            z = [];
        end
        if isempty(z)
            continue;
        end
        z = z(isfinite(z));
        if isempty(z)
            continue;
        end
        if any(abs(z - z(1)) > 1e-12)
            tf = true;
            return;
        end
    end
end

function install_axes_copy_vector_context_menu(fig)
    % 为多子图图窗安装仅复制当前子图(向量)菜单项，避免复制时把 sgtitle 等全局标题一并带走。
    % 用法：在完成绘图（axes 已创建）后调用一次即可。
    if nargin < 1 || isempty(fig) || ~ishandle(fig)
        return;
    end
    try
        if ~strcmp(get(fig, 'Type'), 'figure')
            return;
        end
    catch
        return;
    end

    axs = gobjects(0);
    try
        axs = findall(fig, 'Type', 'axes');
    catch
        axs = gobjects(0);
    end
    if isempty(axs)
        return;
    end

    cfg = [];
    disable_toolbar = false;
    try
        cfg = getappdata(0, 'CO2_HYDRATE_GRAPHICS');
    catch
    end
    try
        if isstruct(cfg) && isfield(cfg, 'disable_axes_toolbar')
            disable_toolbar = logical(cfg.disable_axes_toolbar);
        end
    catch
    end

    for iax = 1:numel(axs)
        ax = axs(iax);
        if ~ishandle(ax)
            continue;
        end
        try
            if ~strcmp(get(ax, 'Type'), 'axes')
                continue;
            end
        catch
            continue;
        end
        try
            ax_tag = get(ax, 'Tag');
            if strcmpi(ax_tag, 'legend') || strcmpi(ax_tag, 'Colorbar')
                continue;
            end
        catch
        end

        if disable_toolbar
            try
                if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar, 'Visible')
                    ax.Toolbar.Visible = 'off';
                end
            catch
            end
        end

        % 避免重复添加
        try
            if isappdata(ax, 'CO2_COPYVECTOR_MENUITEM')
                h = getappdata(ax, 'CO2_COPYVECTOR_MENUITEM');
                if ~isempty(h) && ishandle(h)
                    continue;
                end
            end
        catch
        end

        cm = [];
        try
            cm = ax.UIContextMenu;
        catch
            cm = [];
        end
        if isempty(cm) || ~ishandle(cm)
            try
                cm = uicontextmenu(fig);
                ax.UIContextMenu = cm;
            catch
                continue;
            end
        end

        try
            ax_cb = ax;
            mi = uimenu(cm, 'Label', '复制为向量图（仅当前子图，不含总标题）', ...
                'Callback', @(~,~)copy_axes_vector_to_clipboard(ax_cb));
            setappdata(ax, 'CO2_COPYVECTOR_MENUITEM', mi);
        catch
        end

        % 添加导出菜单项
        try
            uimenu(cm, 'Label', '导出子图为 SVG...', 'Separator', 'on', ...
                'Callback', @(~,~)export_single_axes(ax, 'svg'));
            uimenu(cm, 'Label', '导出子图为 PDF...', ...
                'Callback', @(~,~)export_single_axes(ax, 'pdf'));
            uimenu(cm, 'Label', '导出子图为 PNG (300dpi)...', ...
                'Callback', @(~,~)export_single_axes(ax, 'png'));
        catch
        end

        try
            if ~disable_toolbar && exist('axtoolbarbtn', 'file') == 2 && isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
                if isappdata(ax, 'CO2_COPYVECTOR_TOOLBTN')
                    hb = getappdata(ax, 'CO2_COPYVECTOR_TOOLBTN');
                    if ~isempty(hb) && ishandle(hb)
                        continue;
                    end
                end
                tb = ax.Toolbar;
                btn = axtoolbarbtn(tb, 'push');
                btn.Tooltip = '复制为向量图（仅当前子图，不含总标题）';
                ax_cb = ax;
                btn.ButtonPushedFcn = @(~,~)copy_axes_vector_to_clipboard(ax_cb);
                btn.Icon = local_copy_icon_16();
                setappdata(ax, 'CO2_COPYVECTOR_TOOLBTN', btn);
            end
        catch
        end
    end
end

function icon = local_copy_icon_16()
    % 生成一个简易 16×16 复制图标（避免依赖 matlabroot 图标路径）。
    icon = zeros(16, 16, 3);

    % 背景透明用白底近似（AxesToolbar 会自行融合主题色，这里保持中性）。
    bg = 1.0;
    icon(:,:,1) = bg; icon(:,:,2) = bg; icon(:,:,3) = bg;

    % 两个重叠纸张轮廓（深灰）
    col = [0.20 0.20 0.20];
    icon = draw_rect(icon, 5, 3, 9, 11, col, 1);
    icon = draw_rect(icon, 3, 5, 9, 11, col, 1);

    % 内部几条浅线（淡灰）
    col2 = [0.55 0.55 0.55];
    icon = draw_hline(icon, 5, 6, 10, col2);
    icon = draw_hline(icon, 7, 6, 10, col2);
    icon = draw_hline(icon, 9, 6, 10, col2);
end

function img = draw_rect(img, x0, y0, w, h, col, thickness)
    for t = 0:(thickness-1)
        img = draw_hline(img, y0+t, x0, x0+w-1, col);
        img = draw_hline(img, y0+h-1-t, x0, x0+w-1, col);
        img = draw_vline(img, x0+t, y0, y0+h-1, col);
        img = draw_vline(img, x0+w-1-t, y0, y0+h-1, col);
    end
end

function img = draw_hline(img, y, x1, x2, col)
    y = max(1, min(size(img, 1), y));
    x1 = max(1, min(size(img, 2), x1));
    x2 = max(1, min(size(img, 2), x2));
    if x2 < x1; return; end
    img(y, x1:x2, 1) = col(1);
    img(y, x1:x2, 2) = col(2);
    img(y, x1:x2, 3) = col(3);
end

function img = draw_vline(img, x, y1, y2, col)
    x = max(1, min(size(img, 2), x));
    y1 = max(1, min(size(img, 1), y1));
    y2 = max(1, min(size(img, 1), y2));
    if y2 < y1; return; end
    img(y1:y2, x, 1) = col(1);
    img(y1:y2, x, 2) = col(2);
    img(y1:y2, x, 3) = col(3);
end

function copy_axes_vector_to_clipboard(ax)
    % 将单个坐标区复制为"向量图"到剪贴板，完全不包含 sgtitle 等 figure 级元素。
    % 方法：创建合适尺寸的临时figure，复制axes并放大字体，确保输出美观。
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end
    try
        if ~strcmp(get(ax, 'Type'), 'axes')
            return;
        end
    catch
        return;
    end

    % 获取配置参数（可通过 options.graphics.copy_* 自定义）
    cfg = [];
    try
        cfg = getappdata(0, 'CO2_HYDRATE_GRAPHICS');
    catch
    end
    % 默认图幅与西文字体。
    fig_width_cm = 16;   % 图窗宽度 cm（单栏约8cm，双栏约16cm）
    fig_height_cm = 12;  % 图窗高度 cm
    font_size = 14;      % 坐标轴刻度字号（调大）
    label_size = 16;     % 轴标签字号（调大）
    title_size = 16;     % 标题字号（调大）
    font_name = 'Times New Roman'; % 新罗马字体
    if isstruct(cfg)
        if isfield(cfg, 'copy_fig_width_cm'); fig_width_cm = cfg.copy_fig_width_cm; end
        if isfield(cfg, 'copy_fig_height_cm'); fig_height_cm = cfg.copy_fig_height_cm; end
        if isfield(cfg, 'copy_font_size'); font_size = cfg.copy_font_size; end
        if isfield(cfg, 'copy_label_size'); label_size = cfg.copy_label_size; end
        if isfield(cfg, 'copy_title_size'); title_size = cfg.copy_title_size; end
        if isfield(cfg, 'copy_font_name'); font_name = cfg.copy_font_name; end
    end

    % 创建临时figure，使用合适的尺寸
    temp_fig = [];
    try
        % 转换 cm 到 pixels（假设96 DPI）
        px_per_cm = 96 / 2.54;
        fig_w_px = round(fig_width_cm * px_per_cm);
        fig_h_px = round(fig_height_cm * px_per_cm);

        temp_fig = figure('Visible', 'off', 'Color', 'w', ...
            'Units', 'pixels', 'Position', [100 100 fig_w_px fig_h_px], ...
            'PaperUnits', 'centimeters', 'PaperSize', [fig_width_cm fig_height_cm], ...
            'PaperPosition', [0 0 fig_width_cm fig_height_cm]);

        % 复制axes到临时figure；若遇到 yyaxis 等多坐标系限制，退化为直接复制当前坐标区
        try
            new_ax = copyobj(ax, temp_fig);
        catch MEcopy
            if ~isempty(temp_fig) && ishandle(temp_fig)
                try
                    close(temp_fig);
                catch
                end
            end
            if copy_axes_direct_fallback_with_temp_english(ax, font_name, font_size, label_size, title_size)
                mode_msg = 'vector';
                if axes_contains_raster_content(ax)
                    mode_msg = 'image(900dpi)';
                end
                fprintf('[Copy] Copied as %s (direct-axes fallback): %s\n', mode_msg, MEcopy.message);
                return;
            end

            % 双纵轴图可复制整图后再裁切。
            if ~axes_is_probably_3d(ax)
                if copy_axes_vector_by_cloned_figure(ax, fig_width_cm, fig_height_cm, font_name, font_size, label_size, title_size)
                    mode_msg = 'vector';
                    if axes_contains_raster_content(ax)
                        mode_msg = 'image(900dpi)';
                    end
                    fprintf('[Copy] Copied as %s (cloned-figure fallback): %s\n', mode_msg, MEcopy.message);
                    return;
                end
            end
            rethrow(MEcopy);
        end

        % 调整新axes的位置，留出足够边距给标签
        set(new_ax, 'Units', 'normalized', 'Position', [0.15 0.15 0.78 0.78]);

        % 放大字体并设置为 Times New Roman
        set(new_ax, 'FontSize', font_size, 'FontName', font_name);

        % 将所有文本转为英文并调整字体
        convert_axes_to_english(new_ax, font_name, font_size, label_size, title_size);

        % 调整 colorbar 字体（如果有）
        cb = findobj(new_ax.Parent, 'Type', 'colorbar');
        for icb = 1:numel(cb)
            try
                set(cb(icb), 'FontSize', font_size, 'FontName', font_name);
                if ~isempty(cb(icb).Label)
                    cb(icb).Label.String = translate_label(cb(icb).Label.String);
                    set(cb(icb).Label, 'FontSize', label_size, 'FontName', font_name);
                end
            catch
            end
        end

        % 从源坐标区复制图例。
        rebuild_legend_from_source(ax, new_ax, temp_fig, font_name, font_size);
        apply_axes_vector_copy_style(new_ax);

        % 确保临时figure没有任何sgtitle
        try
            delete(findall(temp_fig, 'Type', 'Text', 'Tag', 'sgtitle'));
        catch
        end

        % 使用copygraphics复制临时figure
        try
            if has_graphics_function('copygraphics')
                mode_msg = 'vector';
                if axes_contains_raster_content(new_ax)
                    % 云图/图像轴：高分辨率栅格复制，避免向量分块白线伪影
                    copygraphics(temp_fig, 'ContentType', 'image', 'Resolution', 900, 'BackgroundColor', 'white');
                    mode_msg = 'image(900dpi)';
                else
                    copygraphics(temp_fig, 'ContentType', 'vector', 'BackgroundColor', 'white');
                end
                close(temp_fig);
                fprintf('[Copy] Copied as %s (%.0fx%.0f cm, font %d pt, Times New Roman)\n', mode_msg, fig_width_cm, fig_height_cm, font_size);
                return;
            end
        catch
        end

        % 备用：使用print复制到剪贴板（Windows EMF）
        try
            if ispc
                print(temp_fig, '-dmeta');
                close(temp_fig);
                fprintf('[Copy] Copied as EMF (%.0fx%.0f cm)\n', fig_width_cm, fig_height_cm);
                return;
            end
        catch
        end

        % 如果都失败了，关闭临时figure
        close(temp_fig);

    catch ME
        % 出错时确保关闭临时figure
        if ~isempty(temp_fig) && ishandle(temp_fig)
            try
                close(temp_fig);
            catch
            end
        end
        warning('CO2:CopyFailed', '%s', ME.message);
    end

    warning('CO2:CopyUnavailable', 'copygraphics/print not available.');
end

function ok = copy_axes_direct_fallback_with_temp_english(ax, font_name, font_size, label_size, title_size)
    ok = false;
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end

    state = capture_axes_annotation_state(ax);
    cleanup_obj = onCleanup(@()restore_axes_annotation_state(state)); %#ok<NASGU>

    try
        set(ax, 'FontSize', font_size, 'FontName', font_name);
    catch
    end

    try
        convert_axes_to_english(ax, font_name, font_size, label_size, title_size);
    catch
    end

    cbs = find_colorbars_for_axes(ax);
    for i = 1:numel(cbs)
        try
            set(cbs(i), 'FontSize', font_size, 'FontName', font_name);
            if ~isempty(cbs(i).Label)
                cbs(i).Label.String = translate_label(cbs(i).Label.String);
                set(cbs(i).Label, 'FontSize', label_size, 'FontName', font_name);
            end
        catch
        end
    end

    legs = find_legends_for_axes(ax);
    for i = 1:numel(legs)
        try
            set(legs(i), 'FontSize', font_size-1, 'FontName', font_name);
            if isprop(legs(i), 'String')
                legs(i).String = cellfun(@translate_label, cellstr(legs(i).String), 'UniformOutput', false);
            end
            if ~isempty(legs(i).Title)
                legs(i).Title.String = translate_label(legs(i).Title.String);
                set(legs(i).Title, 'FontSize', font_size, 'FontName', font_name);
            end
        catch
        end
    end
    apply_axes_vector_copy_style(ax);

    try
        if has_graphics_function('copygraphics')
            if axes_contains_raster_content(ax)
                copygraphics(ax, 'ContentType', 'image', 'Resolution', 900, 'BackgroundColor', 'white');
            else
                copygraphics(ax, 'ContentType', 'vector', 'BackgroundColor', 'white');
            end
            ok = true;
            return;
        end
    catch
    end
end

function state = capture_axes_annotation_state(ax)
    state = struct();
    state.ax = [];
    state.ax_font_name = '';
    state.ax_font_size = [];
    state.ax_x_grid = '';
    state.ax_y_grid = '';
    state.ax_z_grid = '';
    state.ax_x_minor_grid = '';
    state.ax_y_minor_grid = '';
    state.ax_z_minor_grid = '';
    state.labels = struct('h', {}, 'String', {}, 'FontName', {}, 'FontSize', {});
    state.texts = struct('h', {}, 'String', {}, 'FontName', {}, 'FontSize', {});
    state.legends = struct('h', {}, 'String', {}, 'FontName', {}, 'FontSize', {}, ...
        'TitleString', {}, 'TitleFontName', {}, 'TitleFontSize', {}, 'Box', {});
    state.colorbars = struct('h', {}, 'FontName', {}, 'FontSize', {}, ...
        'LabelString', {}, 'LabelFontName', {}, 'LabelFontSize', {});
    state.constantlines = struct('h', {}, 'Label', {});

    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end
    state.ax = ax;
    try; state.ax_font_name = ax.FontName; catch; end
    try; state.ax_font_size = ax.FontSize; catch; end
    try; state.ax_x_grid = ax.XGrid; catch; end
    try; state.ax_y_grid = ax.YGrid; catch; end
    try; state.ax_z_grid = ax.ZGrid; catch; end
    try; state.ax_x_minor_grid = ax.XMinorGrid; catch; end
    try; state.ax_y_minor_grid = ax.YMinorGrid; catch; end
    try; state.ax_z_minor_grid = ax.ZMinorGrid; catch; end

    labels = {};
    try; labels{end+1} = ax.XLabel; catch; end %#ok<AGROW>
    try; labels{end+1} = ax.YLabel; catch; end %#ok<AGROW>
    try; labels{end+1} = ax.ZLabel; catch; end %#ok<AGROW>
    try; labels{end+1} = ax.Title;  catch; end %#ok<AGROW>
    il = 0;
    for i = 1:numel(labels)
        h = labels{i};
        if isempty(h) || ~ishandle(h)
            continue;
        end
        il = il + 1;
        state.labels(il).h = h;
        try; state.labels(il).String = h.String; catch; state.labels(il).String = ''; end
        try; state.labels(il).FontName = h.FontName; catch; state.labels(il).FontName = ''; end
        try; state.labels(il).FontSize = h.FontSize; catch; state.labels(il).FontSize = []; end
    end

    txt = [];
    try
        txt = findall(ax, 'Type', 'Text');
    catch
        txt = [];
    end
    it = 0;
    for i = 1:numel(txt)
        if isempty(txt(i)) || ~ishandle(txt(i))
            continue;
        end
        it = it + 1;
        state.texts(it).h = txt(i);
        try; state.texts(it).String = txt(i).String; catch; state.texts(it).String = ''; end
        try; state.texts(it).FontName = txt(i).FontName; catch; state.texts(it).FontName = ''; end
        try; state.texts(it).FontSize = txt(i).FontSize; catch; state.texts(it).FontSize = []; end
    end

    legs = find_legends_for_axes(ax);
    il = 0;
    for i = 1:numel(legs)
        if isempty(legs(i)) || ~ishandle(legs(i))
            continue;
        end
        il = il + 1;
        state.legends(il).h = legs(i);
        try; state.legends(il).String = legs(i).String; catch; state.legends(il).String = {}; end
        try; state.legends(il).FontName = legs(i).FontName; catch; state.legends(il).FontName = ''; end
        try; state.legends(il).FontSize = legs(i).FontSize; catch; state.legends(il).FontSize = []; end
        try; state.legends(il).TitleString = legs(i).Title.String; catch; state.legends(il).TitleString = ''; end
        try; state.legends(il).TitleFontName = legs(i).Title.FontName; catch; state.legends(il).TitleFontName = ''; end
        try; state.legends(il).TitleFontSize = legs(i).Title.FontSize; catch; state.legends(il).TitleFontSize = []; end
        try; state.legends(il).Box = legs(i).Box; catch; state.legends(il).Box = ''; end
    end

    cbs = find_colorbars_for_axes(ax);
    ic = 0;
    for i = 1:numel(cbs)
        if isempty(cbs(i)) || ~ishandle(cbs(i))
            continue;
        end
        ic = ic + 1;
        state.colorbars(ic).h = cbs(i);
        try; state.colorbars(ic).FontName = cbs(i).FontName; catch; state.colorbars(ic).FontName = ''; end
        try; state.colorbars(ic).FontSize = cbs(i).FontSize; catch; state.colorbars(ic).FontSize = []; end
        try; state.colorbars(ic).LabelString = cbs(i).Label.String; catch; state.colorbars(ic).LabelString = ''; end
        try; state.colorbars(ic).LabelFontName = cbs(i).Label.FontName; catch; state.colorbars(ic).LabelFontName = ''; end
        try; state.colorbars(ic).LabelFontSize = cbs(i).Label.FontSize; catch; state.colorbars(ic).LabelFontSize = []; end
    end

    cls = [];
    try
        cls = findall(ax, 'Type', 'ConstantLine');
    catch
        cls = [];
    end
    iln = 0;
    for i = 1:numel(cls)
        if isempty(cls(i)) || ~ishandle(cls(i))
            continue;
        end
        iln = iln + 1;
        state.constantlines(iln).h = cls(i);
        try; state.constantlines(iln).Label = cls(i).Label; catch; state.constantlines(iln).Label = ''; end
    end
end

function restore_axes_annotation_state(state)
    if ~isstruct(state)
        return;
    end

    try
        if isfield(state, 'ax') && ~isempty(state.ax) && ishandle(state.ax)
            if isfield(state, 'ax_font_name') && ~isempty(state.ax_font_name)
                set(state.ax, 'FontName', state.ax_font_name);
            end
            if isfield(state, 'ax_font_size') && ~isempty(state.ax_font_size)
                set(state.ax, 'FontSize', state.ax_font_size);
            end
            if isfield(state, 'ax_x_grid') && ~isempty(state.ax_x_grid)
                set(state.ax, 'XGrid', state.ax_x_grid);
            end
            if isfield(state, 'ax_y_grid') && ~isempty(state.ax_y_grid)
                set(state.ax, 'YGrid', state.ax_y_grid);
            end
            if isfield(state, 'ax_z_grid') && ~isempty(state.ax_z_grid)
                set(state.ax, 'ZGrid', state.ax_z_grid);
            end
            if isfield(state, 'ax_x_minor_grid') && ~isempty(state.ax_x_minor_grid)
                set(state.ax, 'XMinorGrid', state.ax_x_minor_grid);
            end
            if isfield(state, 'ax_y_minor_grid') && ~isempty(state.ax_y_minor_grid)
                set(state.ax, 'YMinorGrid', state.ax_y_minor_grid);
            end
            if isfield(state, 'ax_z_minor_grid') && ~isempty(state.ax_z_minor_grid)
                set(state.ax, 'ZMinorGrid', state.ax_z_minor_grid);
            end
        end
    catch
    end

    if isfield(state, 'labels')
        for i = 1:numel(state.labels)
            h = state.labels(i).h;
            if isempty(h) || ~ishandle(h)
                continue;
            end
            try; h.String = state.labels(i).String; catch; end
            try; h.FontName = state.labels(i).FontName; catch; end
            try; h.FontSize = state.labels(i).FontSize; catch; end
        end
    end

    if isfield(state, 'texts')
        for i = 1:numel(state.texts)
            h = state.texts(i).h;
            if isempty(h) || ~ishandle(h)
                continue;
            end
            try; h.String = state.texts(i).String; catch; end
            try; h.FontName = state.texts(i).FontName; catch; end
            try; h.FontSize = state.texts(i).FontSize; catch; end
        end
    end

    if isfield(state, 'legends')
        for i = 1:numel(state.legends)
            h = state.legends(i).h;
            if isempty(h) || ~ishandle(h)
                continue;
            end
            try; h.String = state.legends(i).String; catch; end
            try; h.FontName = state.legends(i).FontName; catch; end
            try; h.FontSize = state.legends(i).FontSize; catch; end
            try; h.Title.String = state.legends(i).TitleString; catch; end
            try; h.Title.FontName = state.legends(i).TitleFontName; catch; end
            try; h.Title.FontSize = state.legends(i).TitleFontSize; catch; end
            try; h.Box = state.legends(i).Box; catch; end
        end
    end

    if isfield(state, 'colorbars')
        for i = 1:numel(state.colorbars)
            h = state.colorbars(i).h;
            if isempty(h) || ~ishandle(h)
                continue;
            end
            try; h.FontName = state.colorbars(i).FontName; catch; end
            try; h.FontSize = state.colorbars(i).FontSize; catch; end
            try; h.Label.String = state.colorbars(i).LabelString; catch; end
            try; h.Label.FontName = state.colorbars(i).LabelFontName; catch; end
            try; h.Label.FontSize = state.colorbars(i).LabelFontSize; catch; end
        end
    end

    if isfield(state, 'constantlines')
        for i = 1:numel(state.constantlines)
            h = state.constantlines(i).h;
            if isempty(h) || ~ishandle(h)
                continue;
            end
            try; h.Label = state.constantlines(i).Label; catch; end
        end
    end
end

function legs = find_legends_for_axes(ax)
    legs = [];
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end
    try
        if isprop(ax, 'Legend') && ~isempty(ax.Legend) && ishandle(ax.Legend)
            legs = ax.Legend;
            return;
        end
    catch
    end
    fig = [];
    try
        fig = ancestor(ax, 'figure');
    catch
        fig = [];
    end
    if isempty(fig) || ~ishandle(fig)
        return;
    end

    cand = findobj(fig, 'Type', 'legend');
    for i = 1:numel(cand)
        if legend_belongs_to_axes(cand(i), ax)
            legs = [legs; cand(i)]; %#ok<AGROW>
        end
    end
end

function cbs = find_colorbars_for_axes(ax)
    cbs = [];
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end
    fig = [];
    try
        fig = ancestor(ax, 'figure');
    catch
        fig = [];
    end
    if isempty(fig) || ~ishandle(fig)
        return;
    end

    cand = findobj(fig, 'Type', 'colorbar');
    for i = 1:numel(cand)
        if colorbar_belongs_to_axes(cand(i), ax)
            cbs = [cbs; cand(i)]; %#ok<AGROW>
        end
    end
end

function tf = legend_belongs_to_axes(leg, ax)
    tf = false;
    if isempty(leg) || ~ishandle(leg) || isempty(ax) || ~ishandle(ax)
        return;
    end

    try
        if isprop(leg, 'Axes')
            ax_ref = leg.Axes;
            for i = 1:numel(ax_ref)
                if isequal(ax_ref(i), ax)
                    tf = true;
                    return;
                end
            end
        end
    catch
    end

    try
        if isprop(leg, 'PlotChildren')
            ph = leg.PlotChildren;
            for i = 1:numel(ph)
                try
                    if isequal(ancestor(ph(i), 'axes'), ax)
                        tf = true;
                        return;
                    end
                catch
                end
            end
        end
    catch
    end
end

function tf = colorbar_belongs_to_axes(cb, ax)
    tf = false;
    if isempty(cb) || ~ishandle(cb) || isempty(ax) || ~ishandle(ax)
        return;
    end

    try
        if isprop(cb, 'Axes')
            ax_ref = cb.Axes;
            for i = 1:numel(ax_ref)
                if isequal(ax_ref(i), ax)
                    tf = true;
                    return;
                end
            end
        end
    catch
    end

    try
        if isprop(cb, 'Peer') && isequal(cb.Peer, ax)
            tf = true;
            return;
        end
    catch
    end
end

function rebuild_legend_from_source(ax_src, ax_dst, fig_dst, font_name, font_size)
    if nargin < 5 || isempty(font_size)
        font_size = 14;
    end
    if nargin < 4 || isempty(font_name)
        font_name = 'Times New Roman';
    end
    if nargin < 3 || isempty(fig_dst) || ~ishandle(fig_dst)
        return;
    end
    if nargin < 2 || isempty(ax_dst) || ~ishandle(ax_dst)
        return;
    end
    if nargin < 1 || isempty(ax_src) || ~ishandle(ax_src)
        return;
    end

    src_leg = get_axes_primary_legend(ax_src);
    if isempty(src_leg) || ~ishandle(src_leg)
        return;
    end

    src_loc = '';
    src_orient = '';
    src_num_cols = [];
    src_box = '';
    src_pos_norm = [];
    src_ax_pos_norm = [];
    dst_ax_pos_norm = [];
    try; src_loc = char(src_leg.Location); catch; end
    try; src_orient = char(src_leg.Orientation); catch; end
    try; src_num_cols = src_leg.NumColumns; catch; end
    try; src_box = char(src_leg.Box); catch; end
    try
        old_units = src_leg.Units;
        src_leg.Units = 'normalized';
        src_pos_norm = src_leg.Position;
        src_leg.Units = old_units;
    catch
        src_pos_norm = [];
    end
    try
        old_units = ax_src.Units;
        ax_src.Units = 'normalized';
        src_ax_pos_norm = ax_src.Position;
        ax_src.Units = old_units;
    catch
        src_ax_pos_norm = [];
    end
    try
        old_units = ax_dst.Units;
        ax_dst.Units = 'normalized';
        dst_ax_pos_norm = ax_dst.Position;
        ax_dst.Units = old_units;
    catch
        dst_ax_pos_norm = [];
    end

    labels = {};
    try
        labels = cellstr(src_leg.String);
    catch
        try
            labels = cellstr(get(src_leg, 'String'));
        catch
            labels = {};
        end
    end
    if isempty(labels)
        return;
    end
    for i = 1:numel(labels)
        labels{i} = translate_label(labels{i});
    end

    src_handles = get_legend_plot_children(src_leg, ax_src);
    dst_handles = map_handles_from_source_axes(ax_src, ax_dst, src_handles);
    if isempty(dst_handles)
        dst_handles = collect_legendable_children(ax_dst);
    end
    if isempty(dst_handles)
        return;
    end

    n = min(numel(labels), numel(dst_handles));
    labels = labels(1:n);
    dst_handles = dst_handles(1:n);

    try
        old_legs = findobj(fig_dst, 'Type', 'legend');
        if ~isempty(old_legs)
            delete(old_legs);
        end
    catch
    end

    new_leg = [];
    try
        new_leg = legend(ax_dst, dst_handles, labels);
    catch
        try
            new_leg = legend(ax_dst, labels);
        catch
            new_leg = [];
        end
    end
    if isempty(new_leg) || ~ishandle(new_leg)
        return;
    end

    try; new_leg.FontName = font_name; catch; end
    try; new_leg.FontSize = font_size - 1; catch; end
    try; new_leg.Interpreter = 'none'; catch; end
    try; if isprop(new_leg, 'AutoUpdate'); new_leg.AutoUpdate = 'off'; end; catch; end
    try
        if ~isempty(src_loc)
            new_leg.Location = src_loc;
        end
    catch
    end
    try
        if ~isempty(src_orient)
            new_leg.Orientation = src_orient;
        end
    catch
    end
    try
        if isprop(new_leg, 'NumColumns') && ~isempty(src_num_cols)
            new_leg.NumColumns = src_num_cols;
        end
    catch
    end
    try
        if ~isempty(src_box)
            new_leg.Box = src_box;
        end
    catch
    end
    try; new_leg.Box = 'off'; catch; end

    % 关键：保留原图例的实际位置，避免 Location='best/none' 在复制图中重新布局。
    % 图例位置由图窗坐标换算到目标坐标区。
    try
        if ~isempty(src_pos_norm) && isnumeric(src_pos_norm) && numel(src_pos_norm) == 4 && all(isfinite(src_pos_norm))
            if strcmpi(src_loc, 'best') || strcmpi(src_loc, 'bestoutside') || strcmpi(src_loc, 'none')
                mapped_pos = map_legend_position_between_axes(src_pos_norm, src_ax_pos_norm, dst_ax_pos_norm);
                if ~isempty(mapped_pos)
                    new_leg.Units = 'normalized';
                    new_leg.Location = 'none';
                    new_leg.Position = mapped_pos;
                end
            end
        end
    catch
    end

    try
        if ~isempty(src_leg.Title)
            new_leg.Title.String = translate_label(src_leg.Title.String);
            new_leg.Title.FontName = font_name;
            new_leg.Title.FontSize = font_size;
        end
    catch
    end
end

function apply_axes_vector_copy_style(ax)
    % 复制为向量图时：去掉网格背景与图例边框，并移除顶部标题与井壁标注文字。
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end

    try; set(ax, 'XGrid', 'off', 'YGrid', 'off'); catch; end
    try; set(ax, 'ZGrid', 'off'); catch; end
    try; set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off'); catch; end
    try; set(ax, 'ZMinorGrid', 'off'); catch; end

    legs = find_legends_for_axes(ax);
    for i = 1:numel(legs)
        try
            legs(i).Box = 'off';
        catch
        end
    end

    try
        if ~isempty(ax.Title) && ishandle(ax.Title)
            ax.Title.String = '';
        end
    catch
    end

    % 处理常量线标签：去掉井壁 r_w 标签；其他标签翻译为英文
    cls = [];
    try
        cls = findall(ax, 'Type', 'ConstantLine');
    catch
        cls = [];
    end
    for i = 1:numel(cls)
        try
            lbl = cls(i).Label;
        catch
            lbl = '';
        end
        lbl_s = '';
        try
            lbl_s = char(lbl);
        catch
            try
                lbl_s = char(string(lbl));
            catch
                lbl_s = '';
            end
        end
        lbl_l = lower(strtrim(lbl_s));
        if contains(lbl_l, 'r_w') || contains(lbl_s, '井壁') || contains(lbl_l, 'well wall')
            try; cls(i).Label = ''; catch; end
        else
            try; cls(i).Label = translate_label(lbl_s); catch; end
        end
    end
end

function dst_pos = map_legend_position_between_axes(src_leg_pos, src_ax_pos, dst_ax_pos)
    dst_pos = [];
    if nargin < 3
        return;
    end
    if ~(isnumeric(src_leg_pos) && isnumeric(src_ax_pos) && isnumeric(dst_ax_pos))
        return;
    end
    if numel(src_leg_pos) ~= 4 || numel(src_ax_pos) ~= 4 || numel(dst_ax_pos) ~= 4
        return;
    end
    if any(~isfinite(src_leg_pos)) || any(~isfinite(src_ax_pos)) || any(~isfinite(dst_ax_pos))
        return;
    end
    if src_ax_pos(3) <= eps || src_ax_pos(4) <= eps || dst_ax_pos(3) <= eps || dst_ax_pos(4) <= eps
        return;
    end

    rel_x = (src_leg_pos(1) - src_ax_pos(1)) / src_ax_pos(3);
    rel_y = (src_leg_pos(2) - src_ax_pos(2)) / src_ax_pos(4);
    rel_w = src_leg_pos(3) / src_ax_pos(3);
    rel_h = src_leg_pos(4) / src_ax_pos(4);

    dst_pos = [ ...
        dst_ax_pos(1) + rel_x * dst_ax_pos(3), ...
        dst_ax_pos(2) + rel_y * dst_ax_pos(4), ...
        rel_w * dst_ax_pos(3), ...
        rel_h * dst_ax_pos(4)];
end

function leg = get_axes_primary_legend(ax)
    leg = [];
    legs = find_legends_for_axes(ax);
    if ~isempty(legs)
        leg = legs(1);
    end
end

function h = get_legend_plot_children(leg, ax_ref)
    h = [];
    if nargin < 1 || isempty(leg) || ~ishandle(leg)
        return;
    end
    try
        if isprop(leg, 'PlotChildren')
            h = leg.PlotChildren;
        end
    catch
        h = [];
    end

    if isempty(h) && nargin >= 2 && ~isempty(ax_ref) && ishandle(ax_ref)
        h = collect_legendable_children(ax_ref);
    end

    if isempty(h)
        return;
    end
    valid = false(numel(h), 1);
    for i = 1:numel(h)
        valid(i) = ~isempty(h(i)) && ishandle(h(i));
    end
    h = h(valid);
end

function dst_handles = map_handles_from_source_axes(ax_src, ax_dst, src_handles)
    dst_handles = [];
    if nargin < 3 || isempty(src_handles)
        return;
    end
    if isempty(ax_src) || ~ishandle(ax_src) || isempty(ax_dst) || ~ishandle(ax_dst)
        return;
    end

    src_children = [];
    dst_children = [];
    try; src_children = allchild(ax_src); catch; end
    try; dst_children = allchild(ax_dst); catch; end
    if isempty(src_children) || isempty(dst_children)
        return;
    end

    n = min(numel(src_children), numel(dst_children));
    src_children = src_children(1:n);
    dst_children = dst_children(1:n);

    for i = 1:numel(src_handles)
        idx = [];
        for j = 1:n
            if isequal(src_children(j), src_handles(i))
                idx = j;
                break;
            end
        end
        if ~isempty(idx) && ishandle(dst_children(idx))
            dst_handles = [dst_handles; dst_children(idx)]; %#ok<AGROW>
        end
    end
end

function h = collect_legendable_children(ax)
    h = [];
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end

    children = [];
    try
        children = allchild(ax);
    catch
        children = [];
    end
    if isempty(children)
        return;
    end

    for i = numel(children):-1:1
        child = children(i);
        if isempty(child) || ~ishandle(child)
            continue;
        end
        icon_style = '';
        try
            icon_style = child.Annotation.LegendInformation.IconDisplayStyle;
        catch
            icon_style = '';
        end
        if strcmpi(icon_style, 'off')
            continue;
        end
        if ~isempty(icon_style)
            h = [h; child]; %#ok<AGROW>
            continue;
        end

        tp = '';
        try
            tp = lower(get(child, 'Type'));
        catch
            tp = '';
        end
        if any(strcmp(tp, {'line', 'patch', 'surface', 'scatter', 'bar', 'area', 'histogram', 'contour', 'hggroup'}))
            h = [h; child]; %#ok<AGROW>
        end
    end
end

function ok = copy_axes_vector_by_cloned_figure(ax, fig_width_cm, fig_height_cm, font_name, font_size, label_size, title_size)
    ok = false;
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end

    fig_src = [];
    fig_clone = [];
    try
        fig_src = ancestor(ax, 'figure');
    catch
        fig_src = [];
    end
    if isempty(fig_src) || ~ishandle(fig_src)
        return;
    end

    try
        fig_clone = copyobj(fig_src, 0);
        px_per_cm = 96 / 2.54;
        fig_w_px = round(fig_width_cm * px_per_cm);
        fig_h_px = round(fig_height_cm * px_per_cm);
        set(fig_clone, 'Visible', 'off', 'Color', 'w', ...
            'Units', 'pixels', 'Position', [100 100 fig_w_px fig_h_px], ...
            'PaperUnits', 'centimeters', 'PaperSize', [fig_width_cm fig_height_cm], ...
            'PaperPosition', [0 0 fig_width_cm fig_height_cm]);

        try
            delete(findall(fig_clone, 'Type', 'Text', 'Tag', 'sgtitle'));
        catch
        end

        axs_clone = findall(fig_clone, 'Type', 'axes');
        if isempty(axs_clone)
            error('No axes in cloned figure');
        end

        pos_src = [];
        try
            pos_src = get(ax, 'Position');
        catch
            pos_src = [];
        end

        cand = gobjects(0);
        for i = 1:numel(axs_clone)
            tg = '';
            try; tg = get(axs_clone(i), 'Tag'); catch; end
            if strcmpi(tg, 'legend') || strcmpi(tg, 'Colorbar')
                continue;
            end
            cand(end+1) = axs_clone(i); %#ok<AGROW>
        end
        if isempty(cand)
            cand = axs_clone;
        end

        ax_clone = cand(1);
        if ~isempty(pos_src)
            best = inf;
            for i = 1:numel(cand)
                try
                    p = get(cand(i), 'Position');
                    d = sum(abs(p(:) - pos_src(:)));
                    if d < best
                        best = d;
                        ax_clone = cand(i);
                    end
                catch
                end
            end
        end

        for i = 1:numel(cand)
            if ~isequal(cand(i), ax_clone)
                try; delete(cand(i)); catch; end
            end
        end

        try
            set(ax_clone, 'Units', 'normalized', 'Position', [0.15 0.15 0.78 0.78]);
            set(ax_clone, 'FontSize', font_size, 'FontName', font_name);
        catch
        end

        convert_axes_to_english(ax_clone, font_name, font_size, label_size, title_size);

        cbs = findobj(fig_clone, 'Type', 'colorbar');
        for i = 1:numel(cbs)
            try
                set(cbs(i), 'FontSize', font_size, 'FontName', font_name);
                if ~isempty(cbs(i).Label)
                    cbs(i).Label.String = translate_label(cbs(i).Label.String);
                    set(cbs(i).Label, 'FontSize', label_size, 'FontName', font_name);
                end
            catch
            end
        end

        legs = findobj(fig_clone, 'Type', 'legend');
        for i = 1:numel(legs)
            try
                set(legs(i), 'FontSize', font_size-1, 'FontName', font_name);
                if isprop(legs(i), 'String')
                    legs(i).String = cellfun(@translate_label, cellstr(legs(i).String), 'UniformOutput', false);
                end
            catch
            end
        end
        apply_axes_vector_copy_style(ax_clone);

        if has_graphics_function('copygraphics')
            if axes_contains_raster_content(ax_clone)
                copygraphics(fig_clone, 'ContentType', 'image', 'Resolution', 900, 'BackgroundColor', 'white');
            else
                copygraphics(fig_clone, 'ContentType', 'vector', 'BackgroundColor', 'white');
            end
            ok = true;
        elseif ispc
            print(fig_clone, '-dmeta');
            ok = true;
        end
    catch
        ok = false;
    end

    if ~isempty(fig_clone) && ishandle(fig_clone)
        try; close(fig_clone); catch; end
    end
end

function export_single_axes(ax, fmt)
    % 将单个子图导出为文件（svg/pdf/png），自动弹出保存对话框
    % 所有文字转为英文，使用 Times New Roman 字体
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end
    if nargin < 2 || isempty(fmt)
        fmt = 'svg';
    end
    fmt = lower(strtrim(char(fmt)));

    % 获取配置参数（调大字号，使用 Times New Roman）
    cfg = [];
    try
        cfg = getappdata(0, 'CO2_HYDRATE_GRAPHICS');
    catch
    end
    fig_width_cm = 16;
    fig_height_cm = 12;
    font_size = 14;      % 调大
    label_size = 16;     % 调大
    title_size = 16;     % 调大
    font_name = 'Times New Roman';
    if isstruct(cfg)
        if isfield(cfg, 'copy_fig_width_cm'); fig_width_cm = cfg.copy_fig_width_cm; end
        if isfield(cfg, 'copy_fig_height_cm'); fig_height_cm = cfg.copy_fig_height_cm; end
        if isfield(cfg, 'copy_font_size'); font_size = cfg.copy_font_size; end
        if isfield(cfg, 'copy_label_size'); label_size = cfg.copy_label_size; end
        if isfield(cfg, 'copy_title_size'); title_size = cfg.copy_title_size; end
        if isfield(cfg, 'copy_font_name'); font_name = cfg.copy_font_name; end
    end

    % 弹出保存对话框
    switch fmt
        case 'svg'
            [fname, fpath] = uiputfile({'*.svg', 'SVG Vector (*.svg)'}, 'Export as SVG', 'subplot.svg');
        case 'pdf'
            [fname, fpath] = uiputfile({'*.pdf', 'PDF Document (*.pdf)'}, 'Export as PDF', 'subplot.pdf');
        case 'png'
            [fname, fpath] = uiputfile({'*.png', 'PNG Image (*.png)'}, 'Export as PNG', 'subplot.png');
        otherwise
            [fname, fpath] = uiputfile({'*.*', 'All Files (*.*)'}, 'Export', ['subplot.' fmt]);
    end
    if isequal(fname, 0) || isequal(fpath, 0)
        return;  % 用户取消
    end
    out_file = fullfile(fpath, fname);

    % 优先导出当前坐标区。
    if has_graphics_function('exportgraphics')
        try
            if strcmp(fmt, 'png')
                exportgraphics(ax, out_file, 'Resolution', 300, 'BackgroundColor', 'white');
            else
                exportgraphics(ax, out_file, 'ContentType', 'vector', 'BackgroundColor', 'white');
            end
            fprintf('[Export] Saved (preferred direct axes): %s\n', out_file);
            if ispc
                try
                    winopen(fpath);
                catch
                end
            end
            return;
        catch
            % 直导失败时，再退化到临时figure链路。
        end
    end

    % yyaxis 等多坐标系坐标区不再走 copyobj 临时图链路，直接导出当前坐标区。
    has_multi_coord = false;
    try
        if isprop(ax, 'YAxis') && numel(ax.YAxis) > 1
            has_multi_coord = true;
        end
        if isprop(ax, 'XAxis') && numel(ax.XAxis) > 1
            has_multi_coord = true;
        end
    catch
    end
    if has_multi_coord && has_graphics_function('exportgraphics')
        if strcmp(fmt, 'png')
            exportgraphics(ax, out_file, 'Resolution', 300, 'BackgroundColor', 'white');
        else
            exportgraphics(ax, out_file, 'ContentType', 'vector', 'BackgroundColor', 'white');
        end
        fprintf('[Export] Saved (multi-coordinate direct): %s\n', out_file);
        if ispc
            try
                winopen(fpath);
            catch
            end
        end
        return;
    end

    % 创建临时 figure
    temp_fig = [];
    try
        px_per_cm = 96 / 2.54;
        fig_w_px = round(fig_width_cm * px_per_cm);
        fig_h_px = round(fig_height_cm * px_per_cm);

        temp_fig = figure('Visible', 'off', 'Color', 'w', ...
            'Units', 'pixels', 'Position', [100 100 fig_w_px fig_h_px], ...
            'PaperUnits', 'centimeters', 'PaperSize', [fig_width_cm fig_height_cm], ...
            'PaperPosition', [0 0 fig_width_cm fig_height_cm]);

        try
            new_ax = copyobj(ax, temp_fig);
        catch MEcopy
            if ~isempty(temp_fig) && ishandle(temp_fig)
                try
                    close(temp_fig);
                catch
                end
                temp_fig = [];
            end
            if has_graphics_function('exportgraphics')
                if strcmp(fmt, 'png')
                    exportgraphics(ax, out_file, 'Resolution', 300, 'BackgroundColor', 'white');
                else
                    exportgraphics(ax, out_file, 'ContentType', 'vector', 'BackgroundColor', 'white');
                end
                fprintf('[Export] Saved (direct axes): %s\n', out_file);
                if ispc
                    try
                        winopen(fpath);
                    catch
                    end
                end
                return;
            end
            rethrow(MEcopy);
        end
        set(new_ax, 'Units', 'normalized', 'Position', [0.15 0.15 0.78 0.78]);

        % 设置字体为 Times New Roman
        set(new_ax, 'FontSize', font_size, 'FontName', font_name);

        % 将所有文本转为英文
        convert_axes_to_english(new_ax, font_name, font_size, label_size, title_size);

        % 调整 colorbar
        cb = findobj(temp_fig, 'Type', 'colorbar');
        for icb = 1:numel(cb)
            try
                set(cb(icb), 'FontSize', font_size, 'FontName', font_name);
                if ~isempty(cb(icb).Label)
                    cb(icb).Label.String = translate_label(cb(icb).Label.String);
                    set(cb(icb).Label, 'FontSize', label_size, 'FontName', font_name);
                end
            catch
            end
        end

        % 调整图例
        leg = findobj(temp_fig, 'Type', 'legend');
        for ileg = 1:numel(leg)
            try
                set(leg(ileg), 'FontSize', font_size-1, 'FontName', font_name);
                if isprop(leg(ileg), 'String')
                    leg(ileg).String = cellfun(@translate_label, cellstr(leg(ileg).String), 'UniformOutput', false);
                end
            catch
            end
        end

        % 导出
        if has_graphics_function('exportgraphics')
            if strcmp(fmt, 'png')
                exportgraphics(temp_fig, out_file, 'Resolution', 300, 'BackgroundColor', 'white');
            else
                exportgraphics(temp_fig, out_file, 'ContentType', 'vector', 'BackgroundColor', 'white');
            end
        else
            switch fmt
                case 'svg'
                    print(temp_fig, out_file, '-dsvg');
                case 'pdf'
                    print(temp_fig, out_file, '-dpdf', '-bestfit');
                case 'png'
                    print(temp_fig, out_file, '-dpng', '-r300');
                otherwise
                    print(temp_fig, out_file, ['-d' fmt]);
            end
        end

        close(temp_fig);
        fprintf('[Export] Saved: %s (%.0fx%.0f cm, Times New Roman)\n', out_file, fig_width_cm, fig_height_cm);

        % 打开所在文件夹（可选）
        if ispc
            try
                winopen(fpath);
            catch
            end
        end

    catch ME
        if ~isempty(temp_fig) && ishandle(temp_fig)
            try; close(temp_fig); catch; end
        end
        % 临时图窗导出失败时直接导出坐标区。
        try
            if has_graphics_function('exportgraphics')
                if strcmp(fmt, 'png')
                    exportgraphics(ax, out_file, 'Resolution', 300, 'BackgroundColor', 'white');
                else
                    exportgraphics(ax, out_file, 'ContentType', 'vector', 'BackgroundColor', 'white');
                end
                fprintf('[Export] Saved (final direct fallback): %s\n', out_file);
                if ispc
                    try
                        winopen(fpath);
                    catch
                    end
                end
                return;
            end
        catch ME2
            errordlg(sprintf('Export failed: %s | fallback failed: %s', ME.message, ME2.message), 'Export Error');
            return;
        end
        errordlg(sprintf('Export failed: %s', ME.message), 'Export Error');
    end
end

function convert_axes_to_english(ax, font_name, font_size, label_size, title_size)
    % 将坐标轴的所有文本转换为英文，并设置字体
    if nargin < 2; font_name = 'Times New Roman'; end
    if nargin < 3; font_size = 14; end
    if nargin < 4; label_size = 16; end
    if nargin < 5; title_size = 16; end

    % X轴标签
    if ~isempty(ax.XLabel) && ishandle(ax.XLabel)
        ax.XLabel.String = translate_label(ax.XLabel.String);
        set(ax.XLabel, 'FontSize', label_size, 'FontName', font_name);
    end
    % Y轴标签
    if ~isempty(ax.YLabel) && ishandle(ax.YLabel)
        ax.YLabel.String = translate_label(ax.YLabel.String);
        set(ax.YLabel, 'FontSize', label_size, 'FontName', font_name);
    end
    % Z轴标签
    if ~isempty(ax.ZLabel) && ishandle(ax.ZLabel)
        ax.ZLabel.String = translate_label(ax.ZLabel.String);
        set(ax.ZLabel, 'FontSize', label_size, 'FontName', font_name);
    end
    % 标题
    if ~isempty(ax.Title) && ishandle(ax.Title)
        ax.Title.String = translate_label(ax.Title.String);
        set(ax.Title, 'FontSize', title_size, 'FontName', font_name);
    end

    % 兼容 yyaxis：分别处理全部 X/Y 轴标注（左/右轴）
    try
        if isprop(ax, 'XAxis') && ~isempty(ax.XAxis)
            for ix = 1:numel(ax.XAxis)
                lbl = ax.XAxis(ix).Label;
                if ~isempty(lbl) && ishandle(lbl)
                    lbl.String = translate_label(lbl.String);
                    set(lbl, 'FontSize', label_size, 'FontName', font_name);
                end
            end
        end
    catch
    end
    try
        if isprop(ax, 'YAxis') && ~isempty(ax.YAxis)
            for iy = 1:numel(ax.YAxis)
                lbl = ax.YAxis(iy).Label;
                if ~isempty(lbl) && ishandle(lbl)
                    lbl.String = translate_label(lbl.String);
                    set(lbl, 'FontSize', label_size, 'FontName', font_name);
                end
            end
        end
    catch
    end

    % 处理所有 Text 对象
    txt_objs = findall(ax, 'Type', 'Text');
    for it = 1:numel(txt_objs)
        try
            txt_objs(it).String = translate_label(txt_objs(it).String);
            set(txt_objs(it), 'FontName', font_name, 'FontSize', font_size);
        catch
        end
    end
end

function out = translate_label(in)
    % 将中文标签翻译为英文
    if isempty(in)
        out = in;
        return;
    end

    % 处理 cell 数组
    if iscell(in)
        out = cellfun(@translate_label, in, 'UniformOutput', false);
        return;
    end

    % 转为字符串
    s = char(in);

    % 中英文对照表（常用标签）
    dict = {
        % 复合标签（优先完整短语，避免部分替换后中英混杂）
        'z_{com} (m, 海平面基准深度)', 'z_{com} (m, sea-level referenced depth)';
        '固化封存占注入量（%）', 'Trapped CO2 / injected CO2 (%)';
        '最大 Sh (-)', 'Max Sh (-)';

        % 近井壳层图（优先完整短语，避免部分替换残留中文）
        '壳层水合物饱和度 Sh(r,t)', 'Shell hydrate saturation Sh(r,t)';
        '壳层温度场 T(r,t)', 'Shell temperature field T(r,t)';
        '壳层生长强度 / 固化占比', 'Shell growth intensity / trapping ratio';
        '注入能力衰减（按径向压降线性缩放）', 'Injectivity decline (scaled by radial pressure drop)';
        '径向剖面：Sh(r)（选时刻）', 'Radial profile: Sh(r) (selected times)';
        '井壁 r_w=', 'Well wall r_w=';
        '无温度场结果', 'No temperature-field result';

        % 温度耦合图例（优先完整短语，避免部分替换残留中文）
        '近井单元 T(r_{well})', 'Near-well cell T(r_{well})';
        '前缘 T(r_{hyd,front})', 'Hydrate front T(r_{hyd,front})';
        '远场 T(r_{max})', 'Far field T(r_{max})';
        '羽流区：CO2加权均温', 'Plume zone: CO2 weighted mean temperature';
        '固化区：水合物加权均温', 'Solidification zone: Hydrate weighted mean temperature';

        % 坐标轴标签
        '时间', 'Time';
        '天', 'day';
        '小时', 'hour';
        '深度', 'Depth';
        '半径', 'Radius';
        '距离', 'Distance';
        '温度', 'Temperature';
        '压力', 'Pressure';
        '饱和度', 'Saturation';
        '密度', 'Density';
        '渗透率', 'Permeability';
        '孔隙度', 'Porosity';
        '质量', 'Mass';
        '流量', 'Flow rate';
        '速度', 'Velocity';
        '粘度', 'Viscosity';

        % 物质名称
        'CO2羽流', 'CO2 Plume';
        'CO2饱和度', 'CO2 Saturation';
        '水合物', 'Hydrate';
        '固化水合物', 'Solidified Hydrate';
        '气相', 'Gas Phase';
        '液相', 'Liquid Phase';
        '超临界', 'Supercritical';

        % 图标题常用词
        '云图', 'Contour';
        '剖面', 'Profile';
        '分布', 'Distribution';
        '演化', 'Evolution';
        '变化', 'Variation';
        '顶视', 'Top View';
        '等值面', 'Isosurface';
        '末期', 'Final';
        '初始', 'Initial';

        % 单位相关
        '米', 'm';
        '千克', 'kg';
        '吨', 't';
        '秒', 's';
        '年', 'year';
        '摄氏度', 'degC';

        % MRST相关
        'MRST储层', 'MRST Reservoir';
        '3D羽流', '3D Plume';
        '径向', 'Radial';
        '注入', 'Injection';
        '储层', 'Reservoir';
        '封存', 'Storage';
        '井筒', 'Wellbore';
        '井底', 'Bottomhole';
        '井口', 'Wellhead';
        '泥线', 'Mudline';
        '海底', 'Seabed';
        '海平面基准深度', 'sea-level referenced depth';
        '海平面基准', 'sea-level referenced';
        '海平面', 'sea level';
        '最大', 'Max';
        '壳层', 'Shell';
        '井壁', 'Well wall';
        '近井单元', 'Near-well cell';
        '前缘', 'Front';
        '远场', 'Far field';
        '羽流区', 'Plume zone';
        '固化区', 'Solidification zone';
        '加权均温', 'weighted mean temperature';

        % 符号替换
        '（', '(';
        '）', ')';
        '，', ', ';
        '：', ': ';
        '℃', 'degC';
    };

    % 执行替换
    out = s;
    for i = 1:size(dict, 1)
        out = strrep(out, dict{i, 1}, dict{i, 2});
    end

    % 处理常见模式：Sg(r,t) 等保持不变

    % 如果仍有中文字符，尝试移除或替换
    % 检测中文字符（Unicode范围）
    has_chinese = any(out >= char(hex2dec('4E00')) & out <= char(hex2dec('9FFF')));
    if has_chinese
        % 对于无法翻译的中文，保留原样但添加警告
        % 或者可以选择移除
        warning('Label contains untranslated Chinese: %s', out);
    end
end

function out = sanitize_filename(name)
    if nargin < 1 || isempty(name)
        out = '';
        return;
    end
    if isstring(name)
        name = char(name);
    end
    if ~ischar(name)
        out = '';
        return;
    end
    out = strtrim(name);
    out = regexprep(out, '[\r\n\t]', '_');
    out = regexprep(out, '[<>:"/\\|?*]', '_');
    out = regexprep(out, '\s+', '_');
    out = regexprep(out, '_+', '_');
    out = regexprep(out, '^_+|_+$', '');
    if numel(out) > 160
        out = out(1:160);
    end
end

function fig = safe_figure(varargin)
    max_try = 3;
    fig = [];
    cfg = [];
    try
        cfg = getappdata(0, 'CO2_HYDRATE_GRAPHICS');
    catch
    end
    pause_s = 0.05;
    retry_pause_s = 0.20;
    if isstruct(cfg)
        pause_s = getfield_or_default(cfg, 'figure_pause_s', pause_s);
        retry_pause_s = getfield_or_default(cfg, 'retry_pause_s', retry_pause_s);
    end

    fig_num = [];
    prop_args = varargin;
    if ~isempty(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1}) && isfinite(varargin{1})
        fig_num = varargin{1};
        prop_args = varargin(2:end);
    end
    for k = 1:max_try
        try
            if isempty(fig_num)
                fig = figure(prop_args{:});
            else
                fig = figure(fig_num);
                if ~isempty(prop_args)
                    set(fig, prop_args{:});
                end
            end
            try
                set(fig, 'Visible', 'on', 'WindowStyle', 'normal');
            catch
            end
            if isstruct(cfg) && isfield(cfg, 'renderer') && ~isempty(cfg.renderer)
                try
                    set(fig, 'Renderer', char(cfg.renderer));
                catch
                end
            end
            try
                ensure_figure_on_screen(fig);
            catch
            end
            if usejava('desktop')
                drawnow;
                pause(pause_s);
            end
            return;
        catch ME
            if k >= max_try
                warning('[绘图] 图窗创建失败：%s', ME.message);
                return;
            end
            pause(retry_pause_s);
        end
    end
end

function ensure_figure_on_screen(fig)
    if isempty(fig) || ~ishandle(fig)
        return;
    end

    try
        if isprop(fig, 'WindowState')
            ws = '';
            try
                ws = get(fig, 'WindowState');
            catch
            end
            if ischar(ws) && strcmpi(ws, 'maximized')
                return;
            end
        end
    catch
    end

    old_units = '';
    try
        old_units = get(fig, 'Units');
        set(fig, 'Units', 'pixels');
    catch
    end

    pos = [];
    try
        pos = get(fig, 'Position'); % [x y w h]
    catch
    end
    if isempty(pos) || numel(pos) ~= 4 || any(~isfinite(pos))
        try
            if ~isempty(old_units)
                set(fig, 'Units', old_units);
            end
        catch
        end
        return;
    end

    mp = get(0, 'MonitorPositions'); % [left bottom width height]
    if isempty(mp)
        try
            if ~isempty(old_units)
                set(fig, 'Units', old_units);
            end
        catch
        end
        return;
    end

    fig_center = [pos(1) + 0.5*pos(3), pos(2) + 0.5*pos(4)];
    best = 1;
    for i = 1:size(mp, 1)
        left = mp(i, 1); bottom = mp(i, 2); w = mp(i, 3); h = mp(i, 4);
        if fig_center(1) >= left && fig_center(1) <= (left + w) && fig_center(2) >= bottom && fig_center(2) <= (bottom + h)
            best = i;
            break;
        end
    end
    mon = mp(best, :);
    mon_left = mon(1);
    mon_bottom = mon(2);
    mon_w = mon(3);
    mon_h = mon(4);
    mon_right = mon_left + mon_w;
    mon_top = mon_bottom + mon_h;

    margin_l = 10;
    margin_b = 10;
    margin_r = 10;
    margin_top = 60; % 预留标题栏/任务栏，避免最大化/最小化按钮出屏

    max_w = max(200, mon_w - (margin_l + margin_r));
    max_h = max(200, mon_h - (margin_b + margin_top));
    pos(3) = min(pos(3), max_w);
    pos(4) = min(pos(4), max_h);

    pos(1) = min(max(pos(1), mon_left + margin_l), mon_right - margin_r - pos(3));
    pos(2) = min(max(pos(2), mon_bottom + margin_b), mon_top - margin_top - pos(4));

    try
        set(fig, 'Position', pos);
    catch
    end

    try
        if ~isempty(old_units)
            set(fig, 'Units', old_units);
        end
    catch
    end
end

function graphics_keepalive()
    % 在长循环/批量扫描中调用，避免 GUI 图形超时/握手失败
    if ~usejava('desktop')
        return;
    end
    cfg = [];
    try
        cfg = getappdata(0, 'CO2_HYDRATE_GRAPHICS');
    catch
    end
    pause_s = 0.02;
    if isstruct(cfg)
        pause_s = getfield_or_default(cfg, 'keepalive_pause_s', pause_s);
    end
    try
        drawnow limitrate
    catch
        drawnow
    end
    pause(pause_s);
end

function gas = create_gas_mixture(component_names, z, kij)
    component_names = cellstr(component_names);
    z = z(:)';
    if numel(component_names) ~= numel(z)
        error('gas_components 与 gas_mole_fraction 长度不一致');
    end
    if any(z < 0)
        error('gas_mole_fraction 不能为负');
    end
    if sum(z) <= 0
        error('gas_mole_fraction 之和必须>0');
    end
    z = z / sum(z);

    n = numel(component_names);
    Tc = zeros(1, n);
    Pc = zeros(1, n);
    omega = zeros(1, n);
    Mw = zeros(1, n);
    cp_ig = zeros(1, n);

    for i = 1:n
        p = component_props(component_names{i});
        Tc(i) = p.Tc_K;
        Pc(i) = p.Pc_Pa;
        omega(i) = p.omega;
        Mw(i) = p.Mw_kgmol;
        cp_ig(i) = p.cp_ig_J_mol_K;
    end

    kij_mat = zeros(n, n);
    if nargin >= 3 && ~isempty(kij)
        if ~isnumeric(kij) || any(~isfinite(kij(:)))
            error('gas_kij 必须为有限数值矩阵（n×n）');
        end
        if ~isequal(size(kij), [n, n])
            error('gas_kij 维度不匹配：期望 %dx%d（与 gas_components 顺序一致）', n, n);
        end
        kij_mat = 0.5 * (kij + kij'); % 强制对称
        kij_mat(1:n+1:end) = 0.0;     % 对角线必须为0
    end

    idx_co2 = find(strcmpi(component_names, 'CO2'), 1);
    if isempty(idx_co2)
        error('gas_components 必须包含 CO2（用于水合物平衡/驱动力计算）');
    end

    gas = struct();
    gas.names = component_names(:)';
    gas.z = z;
    gas.Tc_K = Tc;
    gas.Pc_Pa = Pc;
    gas.omega = omega;
    gas.Mw_kgmol = Mw;
    gas.Mw_mix_kgmol = sum(z .* Mw);
    gas.cp_ig_J_mol_K = cp_ig;
    gas.kij = kij_mat;
    gas.idx_co2 = idx_co2;
end

function p = component_props(name)
    key = upper(strtrim(name));
    switch key
        case 'CO2'
            p.Tc_K = 304.1282;
            p.Pc_Pa = 7.3773e6;
            p.omega = 0.22394;
            p.Mw_kgmol = 44.0095e-3;
            p.cp_ig_J_mol_K = 37.1; % 近似常数（用于焓闭合的理想气部分）
        case 'N2'
            p.Tc_K = 126.192;
            p.Pc_Pa = 3.3958e6;
            p.omega = 0.0372;
            p.Mw_kgmol = 28.0134e-3;
            p.cp_ig_J_mol_K = 29.1;
        case 'CH4'
            p.Tc_K = 190.564;
            p.Pc_Pa = 4.5992e6;
            p.omega = 0.011;
            p.Mw_kgmol = 16.043e-3;
            p.cp_ig_J_mol_K = 35.7;
        case 'H2S'
            p.Tc_K = 373.2;
            p.Pc_Pa = 9.0e6;
            p.omega = 0.10;
            p.Mw_kgmol = 34.0809e-3;
            p.cp_ig_J_mol_K = 34.6;
        otherwise
            error('未支持的组分：%s（请在 component_props 中补充）', name);
    end
end

function flash = co2_flash_wilson(T_K, P_Pa, gas)
    % Wilson K值闪蒸（用于是否两相的工程近似；密度仍由PR-EOS计算）
    K = (gas.Pc_Pa ./ P_Pa) .* exp(5.37 * (1 + gas.omega) .* (1 - gas.Tc_K ./ T_K));

    z = gas.z;
    f0 = sum(z .* (K - 1));
    f1 = sum(z .* (K - 1) ./ K);

    if f0 < 0
        beta = 0.0; % 全液
    elseif f1 > 0
        beta = 1.0; % 全气
    else
        lo = 0.0; hi = 1.0;
        for it = 1:80
            mid = 0.5*(lo+hi);
            denom = 1 + mid*(K-1);
            f = sum(z .* (K - 1) ./ denom);
            if abs(f) < 1e-10
                break;
            end
            if f > 0
                lo = mid;
            else
                hi = mid;
            end
        end
        beta = 0.5*(lo+hi);
    end

    denom = 1 + beta*(K-1);
    x = z ./ denom;
    y = K .* x;
    x = x / sum(x);
    y = y / sum(y);

    flash = struct('beta', beta, 'K', K, 'x', x, 'y', y);
end

function flash = co2_flash_pr_phi(T_K, P_Pa, gas, flash_guess, flashopt)
    % PR-EOS (phi-phi) 闪蒸：K_i = phi_L/phi_V
    if nargin < 4
        flash_guess = [];
    end
    if nargin < 5 || isempty(flashopt)
        flashopt = struct();
    end
    max_iter = getfield_or_default(flashopt, 'max_iter', 30);
    tol_lnK = getfield_or_default(flashopt, 'tol_lnK', 1e-8);
    damping = min(max(getfield_or_default(flashopt, 'damping', 0.5), 0.0), 1.0);

    z = gas.z(:)' / sum(gas.z);

    % 先做相态可行性判断：若PR在整体组分下只有单实根，则认为单相，避免 K≈1 时 beta=0.5 的数值假解
    rZ = pr_eos_real_roots_pr(T_K, P_Pa, z, gas);
    if numel(rZ) <= 1
        flash_w = co2_flash_wilson(T_K, P_Pa, gas);
        beta_sp = 0.0;
        if flash_w.beta > 0.5
            beta_sp = 1.0;
        elseif flash_w.beta == 0.5
            beta_sp = ternary(rZ(1) > 0.5, 1.0, 0.0);
        end
        x = z;
        y = z;
        K = ones(size(z));
        flash = struct('beta', beta_sp, 'K', K, 'x', x, 'y', y, 'converged', true);
        return;
    end

    % 初值：优先用上一点结果，否则用Wilson
    if isstruct(flash_guess) && isfield(flash_guess, 'K') && numel(flash_guess.K) == numel(z)
        K = max(flash_guess.K(:)', 1e-12);
    else
        flash0 = co2_flash_wilson(T_K, P_Pa, gas);
        K = max(flash0.K(:)', 1e-12);
    end
    K = min(max(K, 1e-8), 1e8);

    lnK = log(K);
    converged = false;
    beta = 0.0;
    x = z;
    y = z;

    for it = 1:max_iter
        beta = rachford_rice_beta(z, exp(lnK));
        denom = 1 + beta*(exp(lnK) - 1);
        denom = max(denom, 1e-15);
        x = z ./ denom;
        y = exp(lnK) .* x;
        x = x / sum(x);
        y = y / sum(y);

        [lnphi_L, ~] = pr_eos_lnphi_pr(T_K, P_Pa, x, gas, 'liquid');
        [lnphi_V, ~] = pr_eos_lnphi_pr(T_K, P_Pa, y, gas, 'vapor');

        lnK_new = lnphi_L - lnphi_V;
        if all(isfinite(lnK_new))
            if max(abs(lnK_new - lnK)) < tol_lnK
                lnK = lnK_new;
                converged = true;
                break;
            end
            lnK = (1 - damping) * lnK + damping * lnK_new;
            lnK = min(max(lnK, log(1e-12)), log(1e12));
        else
            break;
        end
    end

    K = exp(lnK);
    if max(abs(K - 1)) < 1e-6
        % 近似单相：用Wilson给出偏气/偏液方向，避免 beta=0.5
        flash_w = co2_flash_wilson(T_K, P_Pa, gas);
        beta = ternary(flash_w.beta >= 0.5, 1.0, 0.0);
        x = z;
        y = z;
        flash = struct('beta', beta, 'K', K, 'x', x, 'y', y, 'converged', converged);
        return;
    end
    beta = rachford_rice_beta(z, K);
    denom = 1 + beta*(K - 1);
    denom = max(denom, 1e-15);
    x = z ./ denom;
    y = K .* x;
    x = x / sum(x);
    y = y / sum(y);

    flash = struct('beta', beta, 'K', K, 'x', x, 'y', y, 'converged', converged);
end

function beta = rachford_rice_beta(z, K)
    % 求解 Rachford-Rice 方程的摩尔汽相分率 beta
    z = z(:)';
    K = K(:)';
    K = max(K, 1e-12);

    f0 = sum(z .* (K - 1));
    f1 = sum(z .* (K - 1) ./ K);

    if f0 < 0
        beta = 0.0;
        return;
    end
    if f1 > 0
        beta = 1.0;
        return;
    end

    lo = 0.0;
    hi = 1.0;
    for it = 1:80
        mid = 0.5 * (lo + hi);
        denom = 1 + mid*(K - 1);
        denom = max(denom, 1e-15);
        f = sum(z .* (K - 1) ./ denom);
        if abs(f) < 1e-12
            break;
        end
        if f > 0
            lo = mid;
        else
            hi = mid;
        end
    end
    beta = 0.5 * (lo + hi);
end

function [lnphi, Z] = pr_eos_lnphi_pr(T_K, P_Pa, comp, gas, phase)
    % Peng-Robinson EOS fugacity coefficient for mixture (supports kij via gas.kij)
    R = 8.314462618;
    comp = comp(:)' / sum(comp);

    Tc = gas.Tc_K(:)';
    Pc = gas.Pc_Pa(:)';
    omega = gas.omega(:)';

    kappa = 0.37464 + 1.54226*omega - 0.26992*omega.^2;
    Tr = T_K ./ Tc;
    sqrtTr = sqrt(max(Tr, 1e-12));
    alpha = (1 + kappa .* (1 - sqrtTr)).^2;

    ai0 = 0.45724 * (R^2 .* Tc.^2 ./ Pc);
    ai = ai0 .* alpha;
    bi = 0.07780 * (R .* Tc ./ Pc);

    n = numel(comp);
    kij = zeros(n, n);
    if isfield(gas, 'kij') && isnumeric(gas.kij) && isequal(size(gas.kij), [n, n])
        kij = 0.5 * (gas.kij + gas.kij');
        kij(1:n+1:end) = 0.0;
    end

    b_m = sum(comp .* bi);
    sqrt_ai = sqrt(max(ai, 0.0));
    a_ij = (sqrt_ai(:) * sqrt_ai(:).') .* (1.0 - kij);
    a_m = comp * a_ij * comp(:);
    a_m = max(a_m, 1e-30);
    a_i_mix = (a_ij * comp(:))';

    A = a_m * P_Pa / (R^2 * T_K^2);
    B = b_m * P_Pa / (R * T_K);

    coeffs = [1, -(1-B), (A - 3*B^2 - 2*B), -(A*B - B^2 - B^3)];
    r = roots(coeffs);
    r = r(abs(imag(r)) < 1e-9);
    if isempty(r)
        error('PR-EOS fugacity 失败：无实根（T=%.2fK, P=%.2fMPa）', T_K, P_Pa/1e6);
    end
    r = sort(real(r));
    if strcmpi(phase, 'vapor')
        Z = r(end);
    else
        Z = r(1);
    end
    Z = max(Z, 1e-9);

    Bi = bi / max(b_m, 1e-30);
    term1 = Bi .* (Z - 1) - log(max(Z - B, 1e-12));

    sqrt2 = sqrt(2);
    logterm = log(max((Z + (1+sqrt2)*B) / max(Z + (1-sqrt2)*B, 1e-12), 1e-12));
    term2 = (A / max(2*sqrt2*B, 1e-30)) .* (2*a_i_mix/max(a_m, 1e-30) - Bi) .* logterm;

    lnphi = term1 - term2;
end

function r = pr_eos_real_roots_pr(T_K, P_Pa, comp, gas)
    % 返回 PR EOS 的实根（升序）
    R = 8.314462618;
    comp = comp(:)' / sum(comp);

    Tc = gas.Tc_K(:)';
    Pc = gas.Pc_Pa(:)';
    omega = gas.omega(:)';

    kappa = 0.37464 + 1.54226*omega - 0.26992*omega.^2;
    Tr = T_K ./ Tc;
    sqrtTr = sqrt(max(Tr, 1e-12));
    alpha = (1 + kappa .* (1 - sqrtTr)).^2;

    ai0 = 0.45724 * (R^2 .* Tc.^2 ./ Pc);
    ai = ai0 .* alpha;
    bi = 0.07780 * (R .* Tc ./ Pc);

    n = numel(comp);
    kij = zeros(n, n);
    if isfield(gas, 'kij') && isnumeric(gas.kij) && isequal(size(gas.kij), [n, n])
        kij = 0.5 * (gas.kij + gas.kij');
        kij(1:n+1:end) = 0.0;
    end

    b_m = sum(comp .* bi);
    sqrt_ai = sqrt(max(ai, 0.0));
    a_ij = (sqrt_ai(:) * sqrt_ai(:).') .* (1.0 - kij);
    a_m = comp * a_ij * comp(:);
    a_m = max(a_m, 1e-30);

    A = a_m * P_Pa / (R^2 * T_K^2);
    B = b_m * P_Pa / (R * T_K);

    coeffs = [1, -(1-B), (A - 3*B^2 - 2*B), -(A*B - B^2 - B^3)];
    rr = roots(coeffs);
    rr = rr(abs(imag(rr)) < 1e-9);
    if isempty(rr)
        r = [];
        return;
    end
    r = sort(real(rr));
end

function m_dot_total = total_massflow_from_components(n_dot_comp, gas, m_dot_water, m_dot_hydrate_susp)
    n_dot_comp = n_dot_comp(:)';
    m_dot_gas = sum(n_dot_comp .* gas.Mw_kgmol);
    m_dot_total = m_dot_gas + m_dot_water + m_dot_hydrate_susp;
end

function h_mass = pr_eos_enthalpy_pr(T_K, P_Pa, comp, gas, phase, T_ref_K)
    % PR EOS: h = h_ig(T) + h_res(T,P)
    if nargin < 6 || isempty(T_ref_K)
        T_ref_K = 273.15;
    end
    comp = comp(:)' / sum(comp);
    Mw_mix = sum(comp .* gas.Mw_kgmol);

    h_ig_molar = ideal_gas_enthalpy_molar(T_K, comp, gas, T_ref_K);
    h_res_molar = pr_eos_residual_enthalpy_molar_pr(T_K, P_Pa, comp, gas, phase);

    h_mass = (h_ig_molar + h_res_molar) / max(Mw_mix, 1e-12);
end

function h_ig_molar = ideal_gas_enthalpy_molar(T_K, comp, gas, T_ref_K)
    if nargin < 4 || isempty(T_ref_K)
        T_ref_K = 273.15;
    end
    comp = comp(:)' / sum(comp);
    cp_mix = sum(comp .* gas.cp_ig_J_mol_K);
    h_ig_molar = cp_mix * (T_K - T_ref_K);
end

function [h_res_molar, Z] = pr_eos_residual_enthalpy_molar_pr(T_K, P_Pa, comp, gas, phase)
    % Peng-Robinson residual enthalpy (molar) for mixture (supports kij via gas.kij)
    R = 8.314462618;
    comp = comp(:)' / sum(comp);

    Tc = gas.Tc_K(:)';
    Pc = gas.Pc_Pa(:)';
    omega = gas.omega(:)';

    kappa = 0.37464 + 1.54226*omega - 0.26992*omega.^2;
    Tr = T_K ./ Tc;
    sqrtTr = sqrt(max(Tr, 1e-12));
    m = 1 + kappa .* (1 - sqrtTr);
    alpha = m.^2;

    ai0 = 0.45724 * (R^2 .* Tc.^2 ./ Pc);
    ai = ai0 .* alpha;
    bi = 0.07780 * (R .* Tc ./ Pc);

    d_alpha_dT = -(m .* kappa) ./ max(sqrtTr .* Tc, 1e-12);
    dai_dT = ai0 .* d_alpha_dT;

    n = numel(comp);
    kij = zeros(n, n);
    if isfield(gas, 'kij') && isnumeric(gas.kij) && isequal(size(gas.kij), [n, n])
        kij = 0.5 * (gas.kij + gas.kij');
        kij(1:n+1:end) = 0.0;
    end

    b_m = sum(comp .* bi);
    sqrt_ai = sqrt(max(ai, 0.0));
    d_sqrt_ai_dT = 0.5 * dai_dT ./ max(sqrt_ai, 1e-30);

    a_ij = (sqrt_ai(:) * sqrt_ai(:).') .* (1.0 - kij);
    a_m = comp * a_ij * comp(:);
    a_m = max(a_m, 1e-30);

    da_ij_dT = (d_sqrt_ai_dT(:) * sqrt_ai(:).' + sqrt_ai(:) * d_sqrt_ai_dT(:).') .* (1.0 - kij);
    da_m_dT = comp * da_ij_dT * comp(:);

    A = a_m * P_Pa / (R^2 * T_K^2);
    B = b_m * P_Pa / (R * T_K);

    coeffs = [1, -(1-B), (A - 3*B^2 - 2*B), -(A*B - B^2 - B^3)];
    r = roots(coeffs);
    r = r(abs(imag(r)) < 1e-9);
    if isempty(r)
        error('PR-EOS 焓计算失败：无实根（T=%.2fK, P=%.2fMPa）', T_K, P_Pa/1e6);
    end
    r = sort(real(r));
    if strcmpi(phase, 'vapor')
        Z = r(end);
    else
        Z = r(1);
    end
    Z = max(Z, 1e-9);

    sqrt2 = sqrt(2);
    logterm = log(max((Z + (1+sqrt2)*B) / max(Z + (1-sqrt2)*B, 1e-12), 1e-12));
    h_res_molar = R*T_K*(Z - 1) + ((T_K*da_m_dT - a_m) / max(2*sqrt2*b_m, 1e-30)) * logterm;
end

function [h_mix, flash] = mixture_specific_enthalpy(T_K, P_Pa, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_guess, T_ref_K)
    if nargin < 9
        flash_guess = [];
    end
    if nargin < 10 || isempty(T_ref_K)
        T_ref_K = 273.15;
    end

    n_dot_comp = n_dot_comp(:)';
    n_dot_total = sum(n_dot_comp);
    if n_dot_total <= 0
        error('气相组分摩尔流量<=0，无法计算焓。');
    end
    gas_local = gas;
    gas_local.z = n_dot_comp / n_dot_total;

    flash_model = "wilson";
    if isfield(model, 'thermo') && isfield(model.thermo, 'flash_model')
        flash_model = lower(string(model.thermo.flash_model));
    end
    if flash_model == "phi_pr"
        flashopt = struct();
        if isfield(model, 'thermo') && isfield(model.thermo, 'flash')
            flashopt = model.thermo.flash;
        end
        flash = co2_flash_pr_phi(T_K, P_Pa, gas_local, flash_guess, flashopt);
    else
        flash = co2_flash_wilson(T_K, P_Pa, gas_local);
    end

    beta = flash.beta;
    n_dot_v = beta * n_dot_total;
    n_dot_l = (1 - beta) * n_dot_total;
    Mw_v = sum(flash.y .* gas_local.Mw_kgmol);
    Mw_l = sum(flash.x .* gas_local.Mw_kgmol);
    m_dot_co2_v = n_dot_v * Mw_v;
    m_dot_co2_l = n_dot_l * Mw_l;

    h_v = pr_eos_enthalpy_pr(T_K, P_Pa, flash.y, gas_local, 'vapor', T_ref_K);
    h_l = pr_eos_enthalpy_pr(T_K, P_Pa, flash.x, gas_local, 'liquid', T_ref_K);

    aq = aqueous_props(T_K, P_Pa, thermal.seawater_salinity);
    h_aq = aq.cp * (T_K - T_ref_K);

    cp_hyd = 2000.0;
    if isfield(model, 'thermo') && isfield(model.thermo, 'enthalpy')
        cp_hyd = getfield_or_default(model.thermo.enthalpy, 'hydrate_cp_J_kg_K', 2000.0);
    end
    h_s = cp_hyd * (T_K - T_ref_K);

    m_dot_total = m_dot_co2_v + m_dot_co2_l + m_dot_water + m_dot_hydrate_susp;
    h_mix = (m_dot_co2_v*h_v + m_dot_co2_l*h_l + m_dot_water*h_aq + m_dot_hydrate_susp*h_s) / max(m_dot_total, 1e-12);
end

function T_out_C = solve_temperature_from_enthalpy(T_init_C, P_Pa, h_target, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_guess, T_ref_K)
    if nargin < 10
        flash_guess = [];
    end
    if nargin < 11 || isempty(T_ref_K)
        T_ref_K = 273.15;
    end

    cfg = struct();
    if isfield(model, 'thermo') && isfield(model.thermo, 'enthalpy')
        cfg = model.thermo.enthalpy;
    end
    max_iter = max(4, round(getfield_or_default(cfg, 'max_iter', 10)));
    tol = max(getfield_or_default(cfg, 'tol_J_kg', 200.0), 1.0);
    max_step = max(getfield_or_default(cfg, 'max_step_C', 5.0), 0.1);

    T_min = -80.0;
    T_max = 200.0;

    T0 = min(max(T_init_C, T_min), T_max);
    flash0 = flash_guess;
    [h0, flash0] = mixture_specific_enthalpy(T0+273.15, P_Pa, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash0, T_ref_K);
    e0 = h0 - h_target;

    best_T = T0;
    best_err = abs(e0);
    if best_err <= tol
        T_out_C = best_T;
        return;
    end

    % 以 T0 为中心构造括区，再用二分法求温度。
    step = max_step;
    T_low = max(T0 - step, T_min);
    T_high = min(T0 + step, T_max);

    flash_low = flash0;
    [h_low, flash_low] = mixture_specific_enthalpy(T_low+273.15, P_Pa, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_low, T_ref_K);
    e_low = h_low - h_target;

    flash_high = flash0;
    [h_high, flash_high] = mixture_specific_enthalpy(T_high+273.15, P_Pa, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_high, T_ref_K);
    e_high = h_high - h_target;

    if abs(e_low) < best_err; best_err = abs(e_low); best_T = T_low; end
    if abs(e_high) < best_err; best_err = abs(e_high); best_T = T_high; end
    if best_err <= tol
        T_out_C = best_T;
        return;
    end

    max_expand = 12;
    n_expand = 0;
    while (e_low * e_high > 0) && (n_expand < max_expand) && (T_low > T_min || T_high < T_max)
        n_expand = n_expand + 1;
        step = step * 2;
        new_low = max(T0 - step, T_min);
        new_high = min(T0 + step, T_max);

        if new_low < T_low - 1e-12
            T_low = new_low;
            flash_low = flash0;
            [h_low, flash_low] = mixture_specific_enthalpy(T_low+273.15, P_Pa, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_low, T_ref_K);
            e_low = h_low - h_target;
            if abs(e_low) < best_err; best_err = abs(e_low); best_T = T_low; end
        end
        if new_high > T_high + 1e-12
            T_high = new_high;
            flash_high = flash0;
            [h_high, flash_high] = mixture_specific_enthalpy(T_high+273.15, P_Pa, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_high, T_ref_K);
            e_high = h_high - h_target;
            if abs(e_high) < best_err; best_err = abs(e_high); best_T = T_high; end
        end
        if best_err <= tol
            T_out_C = best_T;
            return;
        end
    end

    if e_low * e_high > 0
        T_out_C = best_T;
        return;
    end

    if e_low > 0
        [T_low, T_high] = deal(T_high, T_low);
        [e_low, e_high] = deal(e_high, e_low);
        [flash_low, flash_high] = deal(flash_high, flash_low);
    end

    for it = 1:max_iter
        T_mid = 0.5 * (T_low + T_high);
        flash_mid = flash0;
        if abs(T_mid - T_low) <= abs(T_mid - T_high)
            flash_mid = flash_low;
        else
            flash_mid = flash_high;
        end
        [h_mid, flash_mid] = mixture_specific_enthalpy(T_mid+273.15, P_Pa, n_dot_comp, m_dot_water, m_dot_hydrate_susp, gas, thermal, model, flash_mid, T_ref_K);
        e_mid = h_mid - h_target;

        if abs(e_mid) < best_err
            best_err = abs(e_mid);
            best_T = T_mid;
        end
        if best_err <= tol
            T_out_C = best_T;
            return;
        end

        if e_low * e_mid <= 0
            T_high = T_mid;
            e_high = e_mid;
            flash_high = flash_mid;
        else
            T_low = T_mid;
            e_low = e_mid;
            flash_low = flash_mid;
        end
    end

    T_out_C = best_T;
end

function props = co2_phase_props_pr(T_K, P_Pa, comp, gas, phase)
    R = 8.314462618;
    comp = comp(:)' / sum(comp);

    [Z, rho] = pr_eos_Z_rho(T_K, P_Pa, comp, gas, phase);
    Mw_mix = sum(comp .* gas.Mw_kgmol);

    % 物性（简化：用于水力/传热的工程近似；如需精确请替换为专用物性库）
    T_C = T_K - 273.15;
    is_vapor = strcmpi(phase, 'vapor');
    if is_vapor
        mu = 1.4e-5 * (T_K/300)^0.70 * sqrt(max(Mw_mix,1e-9)/0.044);
        cp = 850 + 200*(T_K/273.15);
        k = 0.0145 * (T_K/273.15)^0.85;
    else
        mu = 1.0e-4 * exp(300/max(T_K,1e-6)) * sqrt(max(Mw_mix,1e-9)/0.044);
        cp = 2100 + 3*T_C;
        k = 0.105 - 0.00025*T_C + 0.00001*(P_Pa/1e6);
    end
    mu = max(mu, 5e-6);
    cp = max(cp, 500);
    k = max(k, 0.01);

    props = struct('Z', Z, 'rho', rho, 'mu', mu, 'cp', cp, 'k', k, 'Mw_kgmol', Mw_mix, 'R', R);
end

function [Z, rho] = pr_eos_Z_rho(T_K, P_Pa, comp, gas, phase)
    R = 8.314462618;
    comp = comp(:)' / sum(comp);

    Tc = gas.Tc_K;
    Pc = gas.Pc_Pa;
    omega = gas.omega;

    kappa = 0.37464 + 1.54226*omega - 0.26992*omega.^2;
    Tr = T_K ./ Tc;
    alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;

    ai = 0.45724 * (R^2 .* Tc.^2 ./ Pc) .* alpha;
    bi = 0.07780 * (R .* Tc ./ Pc);

    n = numel(comp);
    kij = zeros(n, n);
    if isfield(gas, 'kij') && isnumeric(gas.kij) && isequal(size(gas.kij), [n, n])
        kij = 0.5 * (gas.kij + gas.kij');
        kij(1:n+1:end) = 0.0;
    end

    sqrt_ai = sqrt(max(ai, 0.0));
    a_ij = (sqrt_ai(:) * sqrt_ai(:).') .* (1.0 - kij);
    a_m = comp * a_ij * comp(:);
    b_m = sum(comp .* bi);

    A = a_m * P_Pa / (R^2 * T_K^2);
    B = b_m * P_Pa / (R * T_K);

    coeffs = [1, -(1-B), (A - 3*B^2 - 2*B), -(A*B - B^2 - B^3)];
    r = roots(coeffs);
    r = r(abs(imag(r)) < 1e-8);
    if isempty(r)
        error('PR-EOS 求解失败：无实根（T=%.2fK, P=%.2fMPa）', T_K, P_Pa/1e6);
    end
    r = real(r);

    if strcmpi(phase, 'vapor')
        Z = max(r);
    else
        Z = min(r);
    end
    Z = max(Z, 1e-6);

    Mw_mix = sum(comp .* gas.Mw_kgmol);
    rho = (P_Pa * Mw_mix) / (Z * R * T_K);
    rho = max(rho, 0.1);
end

function aq = aqueous_props(T_K, P_Pa, salinity_wt)
    T_C = T_K - 273.15;
    salinity_wt = max(salinity_wt, 0.0);

    rho = (1000 + 8.0*salinity_wt) - 0.3*(T_C - 4);
    rho = max(rho, 950);

    mu = 1e-3 * exp(-0.033*(T_C - 20));
    mu = mu * (1 + 0.02*salinity_wt);

    cp = 4180 * (1 - 0.01*salinity_wt);
    k = 0.60 * (1 - 0.003*(T_C - 20));
    cp = max(cp, 3000);
    k = max(k, 0.4);

    aq = struct('rho', rho, 'mu', mu, 'cp', cp, 'k', k);
end

function mu_mix = mixture_viscosity(mu_v, mu_l, mu_aq, aq_frac_liq, solid_phi)
    % 液相（CO2液+水相）有效黏度：按体积分数做对数混合，再用Kriger-Dougherty描述固相增黏
    aq_frac_liq = min(max(aq_frac_liq, 0.0), 1.0);
    mu_l = max(mu_l, 1e-12);
    mu_aq = max(mu_aq, 1e-12);

    ln_mu_cont = (1-aq_frac_liq)*log(mu_l) + aq_frac_liq*log(mu_aq);
    mu_cont = exp(ln_mu_cont);

    phi = min(max(solid_phi, 0.0), 0.59);
    phi_max = 0.60;
    intrinsic = 2.5;
    mu_slurry = mu_cont * (1 - phi/phi_max)^(-intrinsic*phi_max);

    % 混相：两流体框架下黏度由连续液相/浆体控制，给一个小下限避免数值问题
    mu_mix = max(mu_slurry, min([mu_v, mu_l, mu_aq, 1e-6]));
end

function cp_mix = mixture_cp(cp_v, cp_l, cp_aq, m_v, m_l, m_w, m_s)
    cp_h = 2100; % 水合物固体比热（近似）
    m_total = m_v + m_l + m_w + m_s;
    if m_total <= 0
        cp_mix = cp_l;
        return;
    end
    cp_mix = (m_v*cp_v + m_l*cp_l + m_w*cp_aq + m_s*cp_h) / m_total;
    cp_mix = max(cp_mix, 500);
end

function k_mix = mixture_k(k_v, k_l, k_aq, aq_frac_liq, alpha_g)
    aq_frac_liq = min(max(aq_frac_liq, 0.0), 1.0);
    k_liq = (1-aq_frac_liq)*k_l + aq_frac_liq*k_aq;
    alpha_g = min(max(alpha_g, 0.0), 1.0);
    k_mix = (1-alpha_g)*k_liq + alpha_g*k_v;
    k_mix = max(k_mix, 0.01);
end

function f = darcy_friction_factor(Re, eps_m, D_m)
    if Re <= 0
        f = 0;
        return;
    end
    if Re < 2300
        f = 64 / Re;
        return;
    end
    rel = eps_m / max(D_m, 1e-12);
    f = 1 / (-1.8*log10((rel/3.7)^1.11 + 6.9/Re))^2;
end

function [dpdz_fric_Pa_m, info] = two_phase_friction_gradient_lockhart_martinelli(j_g, j_l, rho_g, rho_l, mu_g, mu_l, eps_m, D_m)
    % 两相摩阻压降（Lockhart-Martinelli，采用Darcy摩阻因子）
    % 输入为各相表观速度j，密度/黏度为相/液相有效物性。
    % 该相关是工程闭合，用于替代HEM混相摩阻；与drift-flux holdup一起构成最小两流体水力内核。
    if nargin < 8
        error('two_phase_friction_gradient_lockhart_martinelli: 参数不足');
    end

    j_g = max(j_g, 0.0);
    j_l = max(j_l, 0.0);
    rho_g = max(rho_g, 1e-9);
    rho_l = max(rho_l, 1e-9);
    mu_g = max(mu_g, 1e-12);
    mu_l = max(mu_l, 1e-12);
    D_m = max(D_m, 1e-12);
    eps_m = max(eps_m, 0.0);

    info = struct();

    if j_g <= 0 && j_l <= 0
        dpdz_fric_Pa_m = 0.0;
        return;
    end

    % 单相退化
    if j_g <= 0
        Re_l = rho_l * j_l * D_m / mu_l;
        f_l = darcy_friction_factor(Re_l, eps_m, D_m);
        dpdz_fric_Pa_m = f_l * (rho_l * j_l^2 / 2) / D_m;
        info.Re_l = Re_l; info.f_l = f_l; info.Re_g = 0; info.f_g = 0;
        return;
    end
    if j_l <= 0
        Re_g = rho_g * j_g * D_m / mu_g;
        f_g = darcy_friction_factor(Re_g, eps_m, D_m);
        dpdz_fric_Pa_m = f_g * (rho_g * j_g^2 / 2) / D_m;
        info.Re_g = Re_g; info.f_g = f_g; info.Re_l = 0; info.f_l = 0;
        return;
    end

    Re_g = rho_g * j_g * D_m / mu_g;
    Re_l = rho_l * j_l * D_m / mu_l;
    f_g = darcy_friction_factor(Re_g, eps_m, D_m);
    f_l = darcy_friction_factor(Re_l, eps_m, D_m);

    dpdz_g0 = f_g * (rho_g * j_g^2 / 2) / D_m;
    dpdz_l0 = f_l * (rho_l * j_l^2 / 2) / D_m;

    dpdz_g0 = max(dpdz_g0, 1e-12);
    dpdz_l0 = max(dpdz_l0, 1e-12);

    X = sqrt(dpdz_l0 / dpdz_g0);
    X = max(X, 1e-12);

    is_g_lam = Re_g < 2300;
    is_l_lam = Re_l < 2300;
    if ~is_g_lam && ~is_l_lam
        C = 20;
    elseif is_g_lam && is_l_lam
        C = 5;
    elseif ~is_g_lam && is_l_lam
        C = 12;
    else
        C = 10;
    end

    phi_l2 = 1 + C / X + 1 / (X^2);
    dpdz_fric_Pa_m = phi_l2 * dpdz_l0;
    dpdz_fric_Pa_m = max(dpdz_fric_Pa_m, 0.0);

    info.Re_g = Re_g;
    info.Re_l = Re_l;
    info.f_g = f_g;
    info.f_l = f_l;
    info.C = C;
    info.X = X;
    info.phi_l2 = phi_l2;
    info.dpdz_g0 = dpdz_g0;
    info.dpdz_l0 = dpdz_l0;
end

function info = interfacial_drag_shear_two_fluid(D_m, j_g, j_l, alpha_g, rho_g, mu_g, rho_l, mu_l, regime, two_fluid)
    % 相间拖曳/界面剪切（诊断量）：给出界面面积浓度a_i、界面剪切τ_i与动量交换系数K_i
    %
    %
    % - 当前版本仍以 drift-flux 给定 holdup（α_g）与相速，界面项先作为输出；
    if nargin < 10 || isempty(two_fluid)
        two_fluid = struct();
    end

    cfg = struct();
    if isfield(two_fluid, 'interfacial')
        cfg = two_fluid.interfacial;
    end
    enabled = true;
    if isfield(cfg, 'enabled')
        enabled = logical(cfg.enabled);
    end

    info = struct('enabled', enabled, 'regime', '', 'model', '', 'v_g_m_s', 0.0, 'v_l_m_s', 0.0, ...
                  'v_rel_m_s', 0.0, 'a_i_m2_m3', 0.0, 'tau_i_Pa', 0.0, 'F_i_N_m3', 0.0, ...
                  'K_i_N_s_m4', 0.0, 'Re_i', 0.0, 'Cd', NaN, 'f_i', NaN, 'd_char_m', NaN);
    if ~enabled
        return;
    end

    D_m = max(D_m, 1e-12);
    j_g = max(j_g, 0.0);
    j_l = max(j_l, 0.0);
    rho_g = max(rho_g, 1e-9);
    rho_l = max(rho_l, 1e-9);
    mu_g = max(mu_g, 1e-12);
    mu_l = max(mu_l, 1e-12);
    alpha_g = min(max(alpha_g, 0.0), 1.0);
    alpha_l = min(max(1.0 - alpha_g, 0.0), 1.0);

    if j_g <= 0 || j_l <= 0 || alpha_g <= 0 || alpha_l <= 0
        return;
    end

    v_g = j_g / max(alpha_g, 1e-12);
    v_l = j_l / max(alpha_l, 1e-12);
    v_rel = v_g - v_l;
    v_rel_mag = abs(v_rel);
    if v_rel_mag <= 1e-9
        info.v_g_m_s = v_g;
        info.v_l_m_s = v_l;
        info.v_rel_m_s = v_rel;
        return;
    end

    hindered_exp = getfield_or_default(cfg, 'hindered_exp', 0.0);
    hindered = 1.0;
    if hindered_exp > 0
        hindered = (max(alpha_l, 1e-6))^(-hindered_exp);
    end
    tau_max = getfield_or_default(cfg, 'tau_i_max_Pa', 5e4);
    tau_max = max(tau_max, 1.0);

    reg = lower(string(regime));
    if strlength(reg) == 0
        reg = "two_phase";
    end
    info.regime = char(reg);

    if reg == "annular"
        % 环状流：界面视为气核-液膜边界，采用类壁面摩阻闭合
        eps_i = getfield_or_default(cfg, 'interface_roughness_m', 0.0);
        eps_i = max(eps_i, 0.0);

        D_core = D_m * sqrt(max(alpha_g, 0.0));
        D_core = min(max(D_core, 1e-6), D_m);
        a_i = 4.0 * D_core / max(D_m^2, 1e-12); % ≈ 4*sqrt(alpha_g)/D
        Re_i = rho_g * v_rel_mag * D_core / mu_g;
        f_i = darcy_friction_factor(Re_i, eps_i, D_core);
        tau_i = hindered * f_i * (rho_g * v_rel_mag^2) / 8.0;
        tau_i = min(max(tau_i, 0.0), tau_max);
        F_i = tau_i * a_i;
        K_i = F_i / max(v_rel_mag, 1e-12);

        info.model = 'annular_friction';
        info.v_g_m_s = v_g;
        info.v_l_m_s = v_l;
        info.v_rel_m_s = v_rel;
        info.a_i_m2_m3 = a_i;
        info.tau_i_Pa = tau_i;
        info.F_i_N_m3 = F_i;
        info.K_i_N_s_m4 = K_i;
        info.Re_i = Re_i;
        info.f_i = f_i;
        info.d_char_m = D_core;
        return;
    end

    % 分散流（bubbly/slug近似）：球形分散相拖曳（Schiller-Naumann）
    d_b = getfield_or_default(cfg, 'dispersed_bubble_diameter_m', 3e-3);
    d_b = min(max(d_b, 1e-5), 0.5*D_m);

    Re_b = rho_l * v_rel_mag * d_b / mu_l;
    Cd = drag_coefficient_schiller_naumann(Re_b);

    a_i = 6.0 * alpha_g / max(d_b, 1e-12); % m2/m3
    tau_i = hindered * Cd * rho_l * v_rel_mag^2 / 8.0; % Pa（与分散相数密度无关）
    tau_i = min(max(tau_i, 0.0), tau_max);
    F_i = tau_i * a_i; % N/m3
    K_i = F_i / max(v_rel_mag, 1e-12);

    info.model = 'dispersed_drag';
    info.v_g_m_s = v_g;
    info.v_l_m_s = v_l;
    info.v_rel_m_s = v_rel;
    info.a_i_m2_m3 = a_i;
    info.tau_i_Pa = tau_i;
    info.F_i_N_m3 = F_i;
    info.K_i_N_s_m4 = K_i;
    info.Re_i = Re_b;
    info.Cd = Cd;
    info.d_char_m = d_b;
end

function Cd = drag_coefficient_schiller_naumann(Re)
    % Schiller-Naumann 拖曳系数（球形分散相的工程相关）
    Re = max(Re, 0.0);
    if Re <= 1e-12
        Cd = 0.0;
        return;
    end
    if Re < 1000
        Cd = 24.0 / Re * (1.0 + 0.15 * Re^0.687);
    else
        Cd = 0.44;
    end
    Cd = min(max(Cd, 0.0), 10.0);
end

function [C0_eff, Vgj_override_m_s, regime] = drift_flux_params_vertical_simple(j_g, j, rho_l, rho_g, D_m, g)
    % 简化垂直两相流型：用无滑移空隙率alpha0分段（bubbly/slug-annular）来给出drift-flux参数
    % 输出：
    % - C0_eff：分布参数
    % - Vgj_override_m_s：漂移速度幅值（m/s）；NaN表示用默认Harmathy
    % - regime：流型标签（仅用于诊断）
    j_g = max(j_g, 0.0);
    j = max(j, 0.0);
    rho_l = max(rho_l, 1e-9);
    rho_g = max(rho_g, 1e-9);
    D_m = max(D_m, 1e-12);
    g = max(g, 0.0);

    alpha0 = ternary(j > 0, j_g / j, 0.0);
    alpha0 = min(max(alpha0, 0.0), 1.0);

    alpha_bubble_max = 0.25;
    alpha_annular_min = 0.75;

    if alpha0 < alpha_bubble_max
        regime = "bubbly";
        C0_eff = 1.2 - 0.2*sqrt(rho_g / rho_l);
        C0_eff = min(max(C0_eff, 1.0), 1.2);
        Vgj_override_m_s = NaN; % 用Harmathy
        return;
    end

    if alpha0 < alpha_annular_min
        regime = "slug";
        C0_eff = 1.0;
        % Taylor bubble 数量级：0.35*sqrt(g*D*Δρ/ρ_l)
        Vgj_override_m_s = 0.35 * sqrt(g * D_m * max(rho_l - rho_g, 0.0) / rho_l);
        Vgj_override_m_s = max(Vgj_override_m_s, 0.0);
        return;
    end

    regime = "annular";
    C0_eff = 1.0;
    Vgj_override_m_s = 0.0;
end

function [alpha_g, v_g, v_l, Vgj] = void_fraction_drift_flux(j_g, j, rho_l, rho_g, two_fluid, g)
    % Drift-flux holdup（Zuber-Findlay）：v_g = C0*j + Vgj
    % 本模型坐标约定下行注入为正；气相浮力导致的漂移速度取反号（Vgj<0）。
    j_g = max(j_g, 0.0);
    j = max(j, 0.0);

    if j <= 0 || j_g <= 0
        alpha_g = 0.0;
        v_g = max(getfield_or_default(two_fluid, 'min_vg_m_s', 1e-3), 1e-6);
        v_l = j;
        Vgj = 0.0;
        return;
    end
    if j_g >= j
        alpha_g = 1.0;
        v_g = max(j_g, getfield_or_default(two_fluid, 'min_vg_m_s', 1e-3));
        v_l = 0.0;
        Vgj = 0.0;
        return;
    end

    C0 = getfield_or_default(two_fluid, 'C0', 1.2);
    sigma = getfield_or_default(two_fluid, 'sigma_N_m', 0.02);
    direction = string(getfield_or_default(two_fluid, 'flow_direction', 'down'));

    rho_l = max(rho_l, 1e-6);
    rho_g = max(rho_g, 1e-6);
    g = max(g, 0.0);

    % 漂移速度：默认用Harmathy（bubbly数量级），也允许外部覆盖（用于slug/annular等简化流型闭合）
    Vgj_override = getfield_or_default(two_fluid, 'Vgj_override_m_s', NaN);
    if isfinite(Vgj_override)
        Vgj_mag = max(Vgj_override, 0.0);
    else
        Vgj_mag = 1.53 * (g * sigma * max(rho_l - rho_g, 0.0) / (rho_l^2))^(0.25);
    end
    if strcmpi(direction, "down")
        Vgj = -Vgj_mag;
    else
        Vgj = Vgj_mag;
    end

    v_g = C0*j + Vgj;
    v_g = max(v_g, getfield_or_default(two_fluid, 'min_vg_m_s', 1e-3));

    alpha_g = j_g / max(v_g, 1e-12);
    alpha_g = min(max(alpha_g, 0.0), 1.0);

    v_l = (j - alpha_g*v_g) / max(1.0 - alpha_g, 1e-12);
    v_l = max(v_l, 0.0);
end

function [alpha_g, info] = void_fraction_interfacial_terminal(j_g, j_l, rho_g, rho_l, mu_l, D_m, g, two_fluid, alpha_guess, regime)
    % 终端滑移（分散流）反算holdup：以浮力=拖曳(Schiller-Naumann)闭合|v_rel|，
    % 再由质量守恒解二次方程得到alpha_g；环状流等复杂流型由上层回退处理。
    if nargin < 10 || isempty(regime)
        regime = "two_phase";
    end
    if nargin < 9 || isempty(alpha_guess) || ~isfinite(alpha_guess)
        alpha_guess = 0.5;
    end
    if nargin < 8 || isempty(two_fluid)
        two_fluid = struct();
    end

    alpha_g = NaN;
    info = struct('ok', false, 'model', 'interfacial_terminal', 'v_rel_m_s', 0.0, 'v_term_m_s', 0.0, ...
                  'Re', 0.0, 'Cd', 0.0, 'iterations', 0);

    reg = lower(string(regime));
    if reg == "annular"
        info.model = 'annular_fallback';
        return;
    end

    j_g = max(j_g, 0.0);
    j_l = max(j_l, 0.0);
    rho_g = max(rho_g, 1e-9);
    rho_l = max(rho_l, 1e-9);
    mu_l = max(mu_l, 1e-12);
    D_m = max(D_m, 1e-12);
    g = max(g, 0.0);

    if j_g <= 0 && j_l <= 0
        alpha_g = 0.0;
        info.ok = true;
        return;
    end
    if j_g <= 0
        alpha_g = 0.0;
        info.ok = true;
        return;
    end
    if j_l <= 0
        alpha_g = 1.0;
        info.ok = true;
        return;
    end
    if g <= 0
        alpha_g = j_g / max(j_g + j_l, 1e-12);
        info.ok = true;
        return;
    end

    cfg = struct();
    if isfield(two_fluid, 'interfacial')
        cfg = two_fluid.interfacial;
    end
    d_b = getfield_or_default(cfg, 'dispersed_bubble_diameter_m', 3e-3);
    d_b = min(max(d_b, 1e-5), 0.5*D_m);
    hindered_exp = getfield_or_default(cfg, 'hindered_exp', 0.0);

    direction = lower(string(getfield_or_default(two_fluid, 'flow_direction', 'down')));
    sgn = ternary(direction == "down", -1.0, 1.0);

    alpha = min(max(alpha_guess, 1e-6), 1.0 - 1e-6);
    max_iter = 20;
    tol = 1e-4;
    for k = 1:max_iter
        alpha_l = min(max(1.0 - alpha, 1e-6), 1.0);
        hindered = 1.0;
        if hindered_exp > 0
            hindered = (alpha_l)^(-hindered_exp);
        end

        [v_term, Re_b, Cd] = terminal_velocity_sphere_schiller_naumann(rho_l, rho_g, mu_l, d_b, g, hindered);
        v_rel = sgn * v_term;

        [alpha_new, ~, ~, ok] = void_fraction_from_slip_velocity(j_g, j_l, v_rel, alpha_guess);
        if ~ok
            info.iterations = k;
            return;
        end

        info.v_rel_m_s = v_rel;
        info.v_term_m_s = v_term;
        info.Re = Re_b;
        info.Cd = Cd;
        info.iterations = k;

        if abs(alpha_new - alpha) < tol
            alpha = alpha_new;
            info.ok = true;
            alpha_g = alpha;
            return;
        end
        alpha = 0.5*alpha + 0.5*alpha_new;
        alpha = min(max(alpha, 1e-6), 1.0 - 1e-6);
    end

    alpha_g = alpha;
    info.ok = true;
end

function [alpha_g, v_g, v_l, ok] = void_fraction_from_slip_velocity(j_g, j_l, v_rel, alpha_pref)
    j_g = max(j_g, 0.0);
    j_l = max(j_l, 0.0);
    j = j_g + j_l;

    alpha_g = 0.0;
    v_g = 0.0;
    v_l = 0.0;
    ok = false;

    if j <= 0
        ok = true;
        return;
    end
    if j_g <= 0
        alpha_g = 0.0;
        v_l = j;
        ok = true;
        return;
    end
    if j_l <= 0
        alpha_g = 1.0;
        v_g = j;
        ok = true;
        return;
    end

    if nargin < 4 || isempty(alpha_pref) || ~isfinite(alpha_pref)
        alpha_pref = j_g / max(j, 1e-12);
    end
    alpha_pref = min(max(alpha_pref, 0.0), 1.0);

    if abs(v_rel) <= 1e-9
        alpha_g = j_g / max(j, 1e-12);
        v_g = j;
        v_l = j;
        ok = true;
        return;
    end

    a = v_rel;
    b = -(j + v_rel);
    c = j_g;
    disc = b^2 - 4*a*c;
    if disc < 0 || ~isfinite(disc)
        return;
    end
    sd = sqrt(disc);

    alpha1 = (-b + sd) / (2*a);
    alpha2 = (-b - sd) / (2*a);
    cand = [alpha1, alpha2];
    cand = cand(isfinite(cand) & cand > 0.0 & cand < 1.0);
    if isempty(cand)
        return;
    end
    [~, idx] = min(abs(cand - alpha_pref));
    alpha_g = cand(idx);

    v_g = j_g / max(alpha_g, 1e-12);
    v_l = j_l / max(1.0 - alpha_g, 1e-12);
    ok = isfinite(v_g) && isfinite(v_l) && v_g >= 0 && v_l >= 0;
end

function [v_term, Re_b, Cd] = terminal_velocity_sphere_schiller_naumann(rho_l, rho_g, mu_l, d_b, g, hindered)
    rho_l = max(rho_l, 1e-9);
    rho_g = max(rho_g, 1e-9);
    mu_l = max(mu_l, 1e-12);
    d_b = max(d_b, 1e-12);
    g = max(g, 0.0);
    hindered = max(hindered, 1e-9);

    d_rho = max(rho_l - rho_g, 0.0);
    if d_rho <= 0 || g <= 0
        v_term = 0.0;
        Re_b = 0.0;
        Cd = 0.0;
        return;
    end

    Cd = 0.44;
    v_term = 0.0;
    Re_b = 0.0;
    for k = 1:30
        v_new = sqrt((4.0/3.0) * (d_rho/rho_l) * g * d_b / max(Cd*hindered, 1e-12));
        Re_b = rho_l * v_new * d_b / mu_l;
        Cd_new = drag_coefficient_schiller_naumann(Re_b);
        if Cd_new <= 0
            v_term = 0.0;
            Cd = 0.0;
            Re_b = 0.0;
            return;
        end
        if k > 1 && abs(v_new - v_term) / max(v_new, 1e-12) < 1e-3
            v_term = v_new;
            Cd = Cd_new;
            return;
        end
        v_term = v_new;
        Cd = 0.5*Cd + 0.5*Cd_new;
    end
end

function prof = woa18_profile_ncss(lat_deg, lon_deg, month, resolution, timeout_s, cache_mat, enable_cache, verbose)
    % 通过 THREDDS NetCDF Subset Service (NCSS) 获取 WOA18 某点的温度/盐度垂向剖面（CSV）
    %
    % WOA18 采用多年平均气候态（decav）。
    % - 该函数只取单点（lat, lon）所有标准深度层，避免 OPeNDAP 读取较慢的问题。
    prof = struct();
    prof.ok = false;
    prof.source = 'woa18_ncss';
    prof.lat_req = lat_deg;
    prof.lon_req = lon_deg;
    prof.month = month;
    prof.resolution = resolution;
    prof.depth_m = [];
    prof.T_C = [];
    prof.S_psu = [];

    month = min(max(round(month), 1), 12);
    lon_wrapped = mod(lon_deg, 360.0);
    if lon_wrapped < 0
        lon_wrapped = lon_wrapped + 360.0;
    end
    prof.lon_req = lon_wrapped;

    use_cache = enable_cache && ~isempty(cache_mat) && isfile(cache_mat);
    if use_cache
        try
            s = load(cache_mat, 'woa18_cache');
            if isfield(s, 'woa18_cache') && isstruct(s.woa18_cache)
                c = s.woa18_cache;
                same = isfield(c, 'ok') && c.ok ...
                    && isfield(c, 'month') && c.month == month ...
                    && isfield(c, 'resolution') && strcmp(char(c.resolution), char(resolution)) ...
                    && isfield(c, 'lat_req') && abs(c.lat_req - lat_deg) < 1e-6 ...
                    && isfield(c, 'lon_req') && abs(c.lon_req - lon_wrapped) < 1e-6;
                if same
                    prof = c;
                    if verbose
                        fprintf('  [WOA18] 使用缓存: %s\n', cache_mat);
                    end
                    return;
                end
            end
        catch
        end
    end

    opt = weboptions('Timeout', timeout_s);

    try
        url_t = sprintf(['https://www.ncei.noaa.gov/thredds-ocean/ncss/grid/ncei/woa/temperature/decav/%s/' ...
            'woa18_decav_t%02d_01.nc?var=t_an&latitude=%.4f&longitude=%.4f&accept=csv'], ...
            char(resolution), month, lat_deg, lon_wrapped);
        txt_t = webread(url_t, opt);
        [z_t, t_val] = ncss_parse_point_profile_csv(txt_t, "t_an");

        url_s = sprintf(['https://www.ncei.noaa.gov/thredds-ocean/ncss/grid/ncei/woa/salinity/decav/%s/' ...
            'woa18_decav_s%02d_01.nc?var=s_an&latitude=%.4f&longitude=%.4f&accept=csv'], ...
            char(resolution), month, lat_deg, lon_wrapped);
        txt_s = webread(url_s, opt);
        [z_s, s_val] = ncss_parse_point_profile_csv(txt_s, "s_an");
    catch ex
        if verbose
            fprintf('  [WOA18] NCSS 获取失败：%s\n', ex.message);
        end
        return;
    end

    if isempty(z_t) || isempty(z_s)
        if verbose
            fprintf('  [WOA18] NCSS 返回为空（可能为陆点/浅水缺测）。\n');
        end
        return;
    end

    if numel(z_t) ~= numel(z_s) || any(abs(z_t - z_s) > 1e-6)
        z = unique([z_t; z_s]);
        t_val = interp1(z_t, t_val, z, 'linear', NaN);
        s_val = interp1(z_s, s_val, z, 'linear', NaN);
    else
        z = z_t;
    end

    [z, idx] = sort(z(:));
    t_val = t_val(idx);
    s_val = s_val(idx);

    prof.depth_m = z;
    prof.T_C = t_val;
    prof.S_psu = s_val;
    prof.ok = any(isfinite(prof.T_C)) && any(isfinite(prof.S_psu));

    if enable_cache && prof.ok && ~isempty(cache_mat)
        try
            woa18_cache = prof;
            save(cache_mat, 'woa18_cache');
        catch
        end
    end
end

function [z_m, v] = ncss_parse_point_profile_csv(csv_text, var_prefix)
    z_m = [];
    v = [];

    lines = splitlines(string(csv_text));
    lines = lines(strlength(lines) > 0);
    if numel(lines) < 2
        return;
    end

    header = split(lines(1), ',');
    idx_alt = find(startsWith(header, "alt"), 1);
    idx_var = find(startsWith(header, var_prefix), 1);
    if isempty(idx_alt) || isempty(idx_var)
        return;
    end

    n = numel(lines) - 1;
    z_tmp = NaN(n, 1);
    v_tmp = NaN(n, 1);
    for i = 1:n
        parts = split(lines(i+1), ',');
        if numel(parts) >= max(idx_alt, idx_var)
            z_tmp(i) = str2double(parts(idx_alt));
            v_tmp(i) = str2double(parts(idx_var));
        end
    end

    ok = isfinite(z_tmp) & isfinite(v_tmp);
    z_m = z_tmp(ok);
    v = v_tmp(ok);
end

function val = getfield_or_default(s, field, default_val)
    field = char(field);
    val = default_val;
    if isstruct(s) && isfield(s, field)
        val = s.(field);
    end
end

function [T_eq_C, risk_level, subcool, diag] = hydrate_equilibrium_and_risk(P_MPa, T_C, y_co2_gas, salinity_wt, hyd, P_eff_override_MPa)
    % CO2 水合物平衡温度（可配置/可标定）：
    % - 有效压力：P_eff = P*y（分压）或 PR 逸度（由上层传入 P_eff_override_MPa）
    % - 平衡线模型：默认 ln(P_eff) 多项式，可用 table 插值替代
    % - 盐度：线性降温修正 + 裕度 margin

    if nargin < 5 || isempty(hyd) || ~isstruct(hyd)
        hyd = struct();
    end

    % ---- 有效压力（MPa）----
    use_partial_pressure = logical(getfield_or_default(hyd, 'use_partial_pressure', true));
    if nargin >= 6 && ~isempty(P_eff_override_MPa) && isfinite(P_eff_override_MPa)
        P_eff = max(P_eff_override_MPa, 0.0); % 直接使用外部给定的有效压力（如逸度）
    elseif use_partial_pressure
        P_eff = max(P_MPa * max(y_co2_gas, 0.0), 0.0);
    else
        P_eff = max(P_MPa, 0.0); % 更保守：近似为纯CO2分压
    end

    % ---- 平衡线配置（可标定）----
    eqcfg = struct();
    if isfield(hyd, 'equilibrium') && isstruct(hyd.equilibrium)
        eqcfg = hyd.equilibrium;
    end

    P_min_MPa = getfield_or_default(eqcfg, 'P_min_MPa', 1.3);
    if ~(isfinite(P_eff) && P_eff >= max(P_min_MPa, 0.0))
        T_eq_C = NaN;  % 修正：使用NaN代替-100，避免绘图时出现奇异点
        risk_level = 0;
        subcool = -1;
        if nargout >= 4
            diag = struct('ok', false, 'P_eff_MPa', P_eff, 'P_min_MPa', P_min_MPa);
        end
        return;
    end

    model_name = lower(string(getfield_or_default(eqcfg, 'model', 'lnp_poly')));

    % ln(P_eff) 多项式系数：T_pure = c0 + c1*lnP + c2*lnP^2 + ...
    coeffs = getfield_or_default(eqcfg, 'lnp_poly_coeffs', []);
    if isempty(coeffs)
        coeffs = getfield_or_default(eqcfg, 'coeffs', []);
    end
    if isempty(coeffs)
        coeffs = [-3.9477 7.8436 -0.4887]; % 默认：现有拟合
    end
    coeffs = coeffs(:)';

    lnP = log(max(P_eff, 1e-12));
    T_pure = NaN;
    if model_name == "lnp_poly" || model_name == "poly_lnp" || model_name == "poly_lnp_power"
        T_pure = 0.0;
        for k = 1:numel(coeffs)
            T_pure = T_pure + coeffs(k) * lnP^(k-1);
        end
    elseif model_name == "lnp_table" || model_name == "table"
        P_tab = getfield_or_default(eqcfg, 'table_P_MPa', []);
        T_tab = getfield_or_default(eqcfg, 'table_Teq_C', []);
        if isempty(P_tab) || isempty(T_tab) || numel(P_tab) ~= numel(T_tab)
            error('hydrate_equilibrium_and_risk: table 模式需要 equilibrium.table_P_MPa 与 equilibrium.table_Teq_C（等长向量）。');
        end
        P_tab = P_tab(:);
        T_tab = T_tab(:);
        ok = isfinite(P_tab) & isfinite(T_tab) & (P_tab > 0);
        P_tab = P_tab(ok);
        T_tab = T_tab(ok);
        [P_tab, idx_sort] = sort(P_tab);
        T_tab = T_tab(idx_sort);
        lnP_tab = log(P_tab);
        T_pure = interp1(lnP_tab, T_tab, lnP, 'linear', 'extrap');
    else
        error('hydrate_equilibrium_and_risk: 不支持的 equilibrium.model="%s"（可选：lnp_poly / lnp_table）。', char(model_name));
    end

    % ---- 盐度修正 + 裕度 ----
    eq_margin_C = getfield_or_default(eqcfg, 'margin_C', getfield_or_default(hyd, 'equilibrium_margin_C', 0.0));
    salt_k = getfield_or_default(eqcfg, 'salt_deltaT_C_per_wt', getfield_or_default(hyd, 'salt_deltaT_C_per_wt', -0.6));
    delta_salt = salt_k * max(salinity_wt, 0.0);

    T_eq_C = T_pure + delta_salt + eq_margin_C;
    subcool = T_eq_C - T_C;

    % 风险分级
    if subcool < 0
        risk_level = 0;
    elseif subcool < 2
        risk_level = 1;
    elseif subcool < 5
        risk_level = 2;
    else
        risk_level = 3;
    end

    if nargout >= 4
        diag = struct();
        diag.ok = true;
        diag.model = char(model_name);
        diag.coeffs = coeffs;
        diag.use_partial_pressure = use_partial_pressure;
        diag.P_MPa = P_MPa;
        diag.y_co2 = y_co2_gas;
        diag.P_eff_MPa = P_eff;
        diag.P_min_MPa = P_min_MPa;
        diag.T_pure_C = T_pure;
        diag.delta_salt_C = delta_salt;
        diag.margin_C = eq_margin_C;
        diag.T_eq_C = T_eq_C;
        diag.subcool_C = subcool;
        diag.risk_level = risk_level;
    end
end

function dep_frac = hydrate_deposition_fraction(tau_w, dep_model)
    tau_crit = max(dep_model.tau_crit_Pa, 1e-6);
    m = max(dep_model.m, 0.1);
    base = 1 / (1 + (tau_w/tau_crit)^m);
    dep_frac = min(max(base, 0.0), 1.0);
end

function a_int = hydrate_interfacial_area(alpha_l, d_drop)
    alpha_l = min(max(alpha_l, 0.0), 1.0);
    d_drop = max(d_drop, 1e-12);
    a_int = 6 * alpha_l / d_drop; % 1/m
end

function mult = hydrate_kinetics_multiplier(hyd)
    mult = 1.0;
    if ~isstruct(hyd)
        return;
    end

    if isfield(hyd, 'catalyst') && isstruct(hyd.catalyst) ...
            && isfield(hyd.catalyst, 'enabled') && hyd.catalyst.enabled
        k_mult = getfield_or_default(hyd.catalyst, 'k_mult', 1.0);
        mult = mult * max(k_mult, 0.0);
    end
end

function frac = hydrate_mass_fractions(hyd)
    mw_h2o = 18.01528e-3;
    mw_co2 = 44.0095e-3;
    h2o_per = hyd.h2o_per_guest;
    m_hyd_per_mol = mw_co2 + h2o_per*mw_h2o;
    frac = struct();
    frac.water = (h2o_per*mw_h2o) / m_hyd_per_mol;
    frac.co2 = mw_co2 / m_hyd_per_mol;
end

function m_dot_hyd_max = max_hydrate_from_water(m_dot_water, hyd)
    frac = hydrate_mass_fractions(hyd);
    m_dot_hyd_max = m_dot_water / max(frac.water, 1e-9);
end

function m_dot_hyd_max = max_hydrate_from_co2(m_dot_co2, hyd)
    frac = hydrate_mass_fractions(hyd);
    m_dot_hyd_max = m_dot_co2 / max(frac.co2, 1e-9);
end

function m_dot_water = hydrate_water_consumption(m_dot_hydrate, hyd)
    frac = hydrate_mass_fractions(hyd);
    m_dot_water = m_dot_hydrate * frac.water;
end

function m_dot_co2 = hydrate_co2_consumption(m_dot_hydrate, hyd)
    frac = hydrate_mass_fractions(hyd);
    m_dot_co2 = m_dot_hydrate * frac.co2;
end

function d_next = hydrate_aggregation_step(d_curr, solid_phi, dt, agg)
    solid_phi = min(max(solid_phi, 0.0), 0.6);
    dt = max(dt, 0.0);

    k_eff = agg.k_agg;

    d_next = d_curr + k_eff * solid_phi^2 * dt;
    d_next = max(d_next, agg.d0_m);
    d_next = min(d_next, agg.d_crit_m*10);
end

function label = phase_label(beta_vapor, P_MPa, T_C, CO2)
    if beta_vapor > 0.01 && beta_vapor < 0.99
        label = '两相';
        return;
    end
    if T_C > CO2.Tc && P_MPa > CO2.Pc
        label = '超临界';
    elseif beta_vapor >= 0.99
        label = '气态';
    else
        label = '液态';
    end
end

%% 多层传热模型
%
% - 采用径向热阻串联的准稳态近似，返回以内壁面积为基准的 U_total（W/m2-K）。
% - 上层（海水段）外侧可选横掠圆柱对流；地层段外侧用半无限体导热的等效换热系数近似（含时间尺度假设）。
% - 不求解井筒/环空流体的瞬态对流与环空循环，仅提供工程量级的温度剖面闭合。
function [U_total, h_struct] = calculate_multilayer_heat_transfer(...
    depth, T_co2, T_env_local, P_co2, v_co2, Re_co2, well, CO2, thermal, k_co2, mu_co2, cp_co2, model, delta_dep, h_inner_override)

    % 层次（简化，贴近OLGA的热阻串联）：流体 -> (可选水合物沉积) -> 油管钢 -> 环空（仅导热） -> (可选套管/水泥) -> 海水/地层

    % 0. 有效内径（考虑沉积层）
    r_ti = well.coil_tubing_id / 2;
    r_f = max(r_ti - max(delta_dep, 0), model.hydraulics.min_inner_diameter_m/2);
    D_inner = 2 * r_f;

    % 1. 内对流传热系数
    use_override = false;
    if nargin >= 15 && ~isempty(h_inner_override) && isfinite(h_inner_override) && h_inner_override > 0
        use_override = true;
    end
    if use_override
        h_inner = h_inner_override;
        Nu = NaN;
        Pr = NaN;
    else
        [h_inner, Nu, Pr] = internal_htc_db(Re_co2, cp_co2, mu_co2, k_co2, D_inner);
    end
    R_inner_conv = 1 / max(2*pi*r_f*h_inner, 1e-12);

    % 2. 水合物沉积导热（可选）
    R_deposit = 0;
    if max(delta_dep,0) > 0
        R_deposit = log(r_ti/r_f) / (2*pi*max(model.hydrate.k_hydrate,1e-12));
    end

    % 3. 油管壁导热
    k_steel = 45;  % W/(m·K)
    r_to = well.coil_tubing_od / 2;
    R_tubing_wall = log(r_to/r_ti) / (2*pi*k_steel);

    % 环空几何（海水段：5-7/8"隔水通道/钻杆；地层段：9-5/8"套管；裸眼段：8-1/2"井眼）
    is_seawater = depth <= well.mudline_depth;
    use_riser = false;
    if is_seawater ...
            && isfield(well, 'riser_id') && isfield(well, 'riser_od') ...
            && isfinite(well.riser_id) && isfinite(well.riser_od) ...
            && well.riser_id > well.coil_tubing_od * 1.001 ...
            && well.riser_od > well.riser_id * 1.001
        use_riser = true;
    end

    if use_riser
        D_annulus_outer = well.riser_id;
    elseif depth <= well.casing_depth
        % 套管段（含地层套管段）：环空外径为套管内径
        D_annulus_outer = well.casing_id;
    else
        % 裸眼段：环空外径为井眼直径
        D_annulus_outer = well.openhole_diameter;
    end
    D_annulus_outer = max(D_annulus_outer, well.coil_tubing_od * 1.001);
    L_char = max((D_annulus_outer - well.coil_tubing_od) / 2, 1e-6);  % 环空特征长度（m）

    % 3. 环空流体传热（基于物理模型计算）
    if depth < well.platform_elevation
        % 平台段：环空填充海水（或完井液），但外部环境是空气
        % 环空传热特性与海水段相同
        k_seawater = 0.6;  % W/(m·K)
        beta_seawater = 2e-4;  % 体积膨胀系数 (1/K)
        nu_seawater = 1e-6;  % 运动粘度 (m2/s)
        alpha_seawater = 1.4e-7;  % 热扩散系数 (m2/s)

        dT_annulus = abs(T_co2 - T_env_local);  % 温差

        % Rayleigh数（自然对流）
        Ra = 9.81 * beta_seawater * dT_annulus * L_char^3 / (nu_seawater * alpha_seawater);

        % 环空Nusselt数（Raithby-Hollands关联式）
        if Ra > 1e3
            Nu_annulus = 0.18 * (Ra * (well.coil_tubing_od/D_annulus_outer)^(3/5))^0.25;
        else
            Nu_annulus = 1.0;  % 纯导热
        end

        h_annulus = Nu_annulus * k_seawater / L_char;
    elseif depth <= well.mudline_depth
        % 海水段：海水环空（自然对流）
        k_seawater = 0.6;  % W/(m·K)
        beta_seawater = 2e-4;  % 体积膨胀系数 (1/K)
        nu_seawater = 1e-6;  % 运动粘度 (m2/s)
        alpha_seawater = 1.4e-7;  % 热扩散系数 (m2/s)

        dT_annulus = abs(T_co2 - T_env_local);  % 温差

        % Rayleigh数（自然对流）
        Ra = 9.81 * beta_seawater * dT_annulus * L_char^3 / (nu_seawater * alpha_seawater);

        % 环空Nusselt数（Raithby-Hollands关联式）
        if Ra > 1e3
            Nu_annulus = 0.18 * (Ra * (well.coil_tubing_od/D_annulus_outer)^(3/5))^0.25;
        else
            Nu_annulus = 1.0;  % 纯导热
        end

        h_annulus = Nu_annulus * k_seawater / L_char;
    else
        % 地层段
        if depth > well.packer_depth
            % 封隔器以下：静止流体（钻井液残留）
            k_mud = 0.8;  % 钻井液导热系数
            h_annulus = k_mud / L_char;  % 纯导热
        else
            % 钻井液环空（微循环）
            k_mud = 0.8;
            h_annulus = 2 * k_mud / L_char;  % 轻微对流
        end
    end

    % 覆盖：环空按导热处理（不考虑对流分布）
    if strcmpi(model.heat_transfer.annulus_model, 'conduction')
        if depth <= well.mudline_depth
            % 平台段和海水段：环空都填充海水
            k_annulus = 0.6; % 海水
        else
            % 地层段：钻井液环空
            k_annulus = 0.8; % 钻井液/完井液（简化）
        end
        r_annulus_outer = max(D_annulus_outer / 2, r_to * 1.001);
        R_annulus = log(r_annulus_outer/r_to) / (2*pi*max(k_annulus,1e-12));
    else
        R_annulus = 1/(2*pi*r_to*h_annulus);
    end

    % 4. 外管壁导热（海水段：隔水通道/钻杆；地层段：套管；裸眼段：无外管）
    k_cement = 0.8;  % W/(m·K)
    r_cement = well.cement_outer_radius;
    if use_riser
        r_ci = well.riser_id / 2;
        r_co = well.riser_od / 2;
        R_casing_wall = log(r_co/r_ci) / (2*pi*k_steel);
        R_cement = 0.0;     % 海水段无水泥环
        r_outer = r_co;     % 海水段外边界为隔水通道外径
    else
        r_ci = well.casing_id / 2;
        r_co = well.casing_od / 2;
        if depth <= well.casing_depth
            % 套管段（地层内）
            R_casing_wall = log(r_co/r_ci) / (2*pi*k_steel);
            % 5. 水泥环导热
            R_cement = log(r_cement/r_co) / (2*pi*k_cement);
            % 外边界半径（用于外部传热）：水泥环外径
            r_outer = r_cement;
        else
            % 裸眼段
            R_casing_wall = 0.0;
            R_cement = 0.0;
            r_outer = well.openhole_diameter / 2;
        end
    end

    % 6. 地层/海水/空气传热（基于实际物理过程）
    if depth < well.platform_elevation
        % 平台段：空气自然对流
        k_air = 0.026;  % W/(m·K)
        D_outer = 2 * r_outer;  % 外径

        % 空气物性（25°C）
        beta_air = 3.4e-3;  % 体积膨胀系数 (1/K)
        nu_air = 1.5e-5;  % 运动粘度 (m2/s)
        alpha_air = 2.2e-5;  % 热扩散系数 (m2/s)
        Pr_air = 0.71;  % Prandtl数

        dT_air = abs(T_co2 - T_env_local);  % 温差

        % Rayleigh数（自然对流，横向圆柱）
        Ra_D = 9.81 * beta_air * dT_air * D_outer^3 / (nu_air * alpha_air);

        % Churchill-Chu关联式（横向圆柱自然对流）
        if Ra_D < 1e12
            Nu_air = (0.6 + 0.387 * Ra_D^(1/6) / (1 + (0.559/Pr_air)^(9/16))^(8/27))^2;
        else
            Nu_air = 0.1 * Ra_D^(1/3);  % 高Ra数近似
        end
        Nu_air = max(Nu_air, 1.0);  % 最小为纯导热

        h_formation = Nu_air * k_air / D_outer;
        R_formation = 1 / (2*pi*r_outer*h_formation);
    elseif depth <= well.mudline_depth
        % 海水段：对流传热（考虑海流）
        k_seawater = 0.6;  % W/(m·K)
        v_current = 0.05;  % 海流速度 m/s（典型值）
        D_outer = 2 * r_outer;  % 外径

        % 海水物性（3°C）
        rho_seawater = 1028;  % kg/m3
        mu_seawater = 1.6e-3;  % Pa·s
        Pr_seawater = 11.5;  % Prandtl数

        % Reynolds数（横掠圆柱）
        Re_seawater = rho_seawater * v_current * D_outer / mu_seawater;

        % Churchill-Bernstein关联式（横掠圆柱）
        if Re_seawater < 1
            Nu_seawater = 0.3 + 0.62 * Re_seawater^0.5 * Pr_seawater^(1/3);
        else
            Nu_seawater = 0.3 + 0.62 * Re_seawater^0.5 * Pr_seawater^(1/3) * ...
                         (1 + (Re_seawater/282000)^(5/8))^(4/5);
        end

        h_formation = Nu_seawater * k_seawater / D_outer;
        R_formation = 1 / (2*pi*r_outer*h_formation);
    else
        % 地层：稳态导热（考虑时间效应的近似）
        k_formation = well.formation_thermal_conductivity;

        % 半无限大介质的有效传热系数（考虑时间）
        % 假设注入时间 t_injection = 1天 = 86400秒
        t_injection = 86400;  % 秒
        alpha_formation = 1e-6;  % 地层热扩散系数 m2/s

        % 渗透深度
        delta_thermal = sqrt(alpha_formation * t_injection);  % 约0.3m

        h_formation = k_formation / delta_thermal;
        R_formation = 1 / (2*pi*r_outer*h_formation);
    end
    if depth <= well.mudline_depth
        sw_mult = getfield_or_default(model.heat_transfer, 'seawater_external_multiplier', 1.0);
        sw_hmax = getfield_or_default(model.heat_transfer, 'seawater_external_h_max_W_m2K', Inf);
        h_formation = max(h_formation * max(sw_mult, 0.0), 1e-9);
        if isfinite(sw_hmax)
            h_formation = min(h_formation, max(sw_hmax, 1e-9));
        end
        R_formation = 1 / (2*pi*r_outer*h_formation);
    end

    % 总热阻（per unit length）
    R_total = R_inner_conv + R_deposit + R_tubing_wall + R_annulus + ...
              R_casing_wall + R_cement + R_formation;

    % 以内壁（流体侧）面积为基准的总传热系数（A_i = 2*pi*r_f）
    U_total = 1 / max(R_total * 2*pi*r_f, 1e-12);  % W/(m2·K)

    % 返回关键中间量（用于分析/校核）
    h_struct = struct('h_inner', h_inner, 'h_annulus', h_annulus, 'h_external', h_formation, ...
                      'R_deposit', R_deposit, 'R_annulus', R_annulus);
end

function [h, Nu, Pr] = internal_htc_db(Re, cp, mu, k, D)
    % 管内强制对流：Dittus–Boelter（湍流）+ 常温壁面层流Nu=4.36
    Re = max(Re, 0.0);
    cp = max(cp, 0.0);
    mu = max(mu, 1e-12);
    k = max(k, 1e-12);
    D = max(D, 1e-12);
    Pr = max(cp * mu / k, 0.0);
    if Re > 2300
        Nu = 0.023 * Re^0.8 * Pr^0.4;
    else
        Nu = 4.36;
    end
    h = Nu * k / D;
end

function [h_inner, info] = inner_htc_two_phase_regime(D_m, j_g, j_l, alpha_g, props_g, props_l, regime_hint, heat_transfer)
    % 两相内对流换热：按流型处理液相连续/环状液膜/干涸三类主导机制（用于注入保守包络）
    if nargin < 8 || isempty(heat_transfer)
        heat_transfer = struct();
    end

    tp = struct();
    if isfield(heat_transfer, 'inner')
        tp = heat_transfer.inner;
    end
    tp2 = struct();
    if isfield(tp, 'two_phase')
        tp2 = tp.two_phase;
    end
    h_min = getfield_or_default(tp2, 'h_min_W_m2K', 1.0);
    h_max = getfield_or_default(tp2, 'h_max_W_m2K', 1e5);
    dryout_alpha_l = getfield_or_default(tp2, 'dryout_alpha_l', 1e-3);
    dryout_j_l_m_s = getfield_or_default(tp2, 'dryout_j_l_m_s', 1e-4);
    film_min_thickness_m = getfield_or_default(tp2, 'film_min_thickness_m', 5e-6);

    D_m = max(D_m, 1e-12);
    j_g = max(j_g, 0.0);
    j_l = max(j_l, 0.0);
    alpha_g = min(max(alpha_g, 0.0), 1.0);
    alpha_l = min(max(1.0 - alpha_g, 0.0), 1.0);

    v_g = ternary(alpha_g > 0, j_g / max(alpha_g, 1e-12), 0.0);
    v_l = ternary(alpha_l > 0, j_l / max(alpha_l, 1e-12), 0.0);

    Re_g = props_g.rho * max(v_g, 1e-6) * D_m / max(props_g.mu, 1e-12);
    [h_g, ~, ~] = internal_htc_db(Re_g, props_g.cp, props_g.mu, props_g.k, D_m);

    Re_l = props_l.rho * max(v_l, 1e-6) * D_m / max(props_l.mu, 1e-12);
    [h_l_bulk, ~, ~] = internal_htc_db(Re_l, props_l.cp, props_l.mu, props_l.k, D_m);

    regime = lower(string(regime_hint));
    if strlength(regime) == 0 || regime == "unknown"
        alpha0 = ternary((j_g + j_l) > 0, j_g / max(j_g + j_l, 1e-12), alpha_g);
        if alpha0 < 0.25
            regime = "bubbly";
        elseif alpha0 < 0.85
            regime = "slug";
        else
            regime = "annular";
        end
    end

    dry = (j_l <= 0) || (alpha_l <= 0) || ((alpha_l < dryout_alpha_l) && (j_l < dryout_j_l_m_s));

    delta_film = NaN;
    h_l_film = NaN;
    if ~dry && regime == "annular"
        r = D_m / 2;
        alpha_l_eff = min(max(alpha_l, 0.0), 0.99);
        r_core = r * sqrt(max(1.0 - alpha_l_eff, 0.0));
        delta_film = max(r - r_core, film_min_thickness_m);
        delta_film = min(delta_film, 0.49 * r);
        D_h_film = max(4.0 * delta_film, 1e-12);

        Re_film = props_l.rho * max(v_l, 1e-6) * D_h_film / max(props_l.mu, 1e-12);
        [h_l_film, ~, ~] = internal_htc_db(Re_film, props_l.cp, props_l.mu, props_l.k, D_h_film);
    end

    if dry
        h_inner = h_g;
    else
        h_inner = ternary(regime == "annular" && isfinite(h_l_film), h_l_film, h_l_bulk);
    end
    if ~isfinite(h_inner)
        h_inner = h_min;
    end
    h_inner = min(max(h_inner, h_min), h_max);

    info = struct('regime', regime, 'dry', dry, 'alpha_g', alpha_g, 'alpha_l', alpha_l, ...
                  'h_inner', h_inner, 'h_g', h_g, 'h_l_bulk', h_l_bulk, 'h_l_film', h_l_film, ...
                  'delta_film_m', delta_film);
end

%% 饱和压力计算（与原代码一致）
function P_sat = saturation_pressure_co2(T)
    Tc = 31.1;
    Pc = 7.38;

    if T < -56.6
        P_sat = 0.1;
    elseif T >= Tc
        P_sat = Pc;
    else
        T_K = T + 273.15;
        Tc_K = Tc + 273.15;
        tau = 1 - T_K/Tc_K;

        a1 = -7.0602;
        a2 = 1.9391;
        a3 = -1.6463;
        a4 = -3.2995;

        ln_Pr = (a1*tau + a2*tau^1.5 + a3*tau^2 + a4*tau^4)*Tc_K/T_K;
        P_sat = Pc * exp(ln_Pr);
    end

    P_sat = max(0.1, min(Pc, P_sat));
end

function flash = co2_flash_saturation_recondense(T_K, P_Pa, gas, model)
    % Saturation-constrained flash for CO2-rich injection streams.
    % It stabilizes the shallow gas/two-phase interval and lets the stream
    % recondense automatically when hydrostatic pressure exceeds Psat(T).
    T_C = T_K - 273.15;
    P_MPa = P_Pa / 1e6;
    z = gas.z(:)';
    z = z / max(sum(z), 1e-12);

    cfg = struct();
    if isfield(model, 'thermo') && isfield(model.thermo, 'flash')
        cfg = model.thermo.flash;
    end

    min_z_co2 = min(max(getfield_or_default(cfg, 'recondense_min_z_co2', 0.70), 0.0), 1.0);
    beta_min = min(max(getfield_or_default(cfg, 'recondense_beta_min', 0.0), 0.0), 1.0);
    beta_max = min(max(getfield_or_default(cfg, 'recondense_beta_max', 1.0), 0.0), 1.0);
    if beta_max < beta_min
        tmp = beta_min;
        beta_min = beta_max;
        beta_max = tmp;
    end

    flash = struct('beta', 0.0, 'x', z, 'y', z, 'K', ones(size(z)), ...
        'converged', true, 'ok', true, 'method', 'sat_recondense');

    if gas.z(gas.idx_co2) < min_z_co2 || T_C >= 31.1
        flash.method = 'sat_recondense_liquid_or_supercritical';
        return;
    end

    P_sat = saturation_pressure_co2(T_C);
    liq_margin = max(getfield_or_default(cfg, 'recondense_liquid_margin_MPa', 0.03), 0.0);
    vap_margin = max(getfield_or_default(cfg, 'recondense_vapor_margin_MPa', 0.03), 0.0);

    if P_MPa >= P_sat + liq_margin
        beta = 0.0;
        method = 'sat_recondense_liquid';
    elseif P_MPa <= max(P_sat - vap_margin, 0.0)
        beta = 1.0;
        method = 'sat_recondense_vapor';
    else
        width = max(liq_margin + vap_margin, 1e-6);
        beta = (P_sat + liq_margin - P_MPa) / width;
        beta = min(max(beta, 0.0), 1.0);
        method = 'sat_recondense_transition';
    end

    flash.beta = min(max(beta, beta_min), beta_max);
    flash.method = method;
    flash.P_sat_MPa = P_sat;
    flash.margin_MPa = P_MPa - P_sat;
end

%% 破裂压力估算（量级）
function [P_frac_MPa, info] = estimate_fracture_pressure_MPa(well, injection, const)
    P_frac_MPa = NaN;
    info = struct('method', 'unavailable', 'method_str', '', 'depth_ref_m', NaN);
    if nargin < 2 || ~isstruct(injection)
        return;
    end

    if isfield(injection, 'fracture_pressure_MPa') && isfinite(injection.fracture_pressure_MPa)
        P_frac_MPa = injection.fracture_pressure_MPa;
        info.method = 'direct';
        info.method_str = sprintf('指定P_{frac}=%.2f', P_frac_MPa);
        return;
    end

    z_td = getfield_or_default(well, 'total_depth', NaN);
    z_ref0 = getfield_or_default(well, 'platform_elevation', 0.0);
    depth_ref_m = z_td - z_ref0; % 近似：以海平面为参考的水深+井深
    if ~isfinite(depth_ref_m) || depth_ref_m <= 0
        depth_ref_m = z_td;
    end
    depth_ref_m = max(depth_ref_m, 0.0);
    info.depth_ref_m = depth_ref_m;

    g0 = 9.80665;
    if nargin >= 3 && isstruct(const) && isfield(const, 'g') && isfinite(const.g)
        g0 = const.g;
    end

    if isfield(injection, 'fracture_gradient_MPa_per_100m') && isfinite(injection.fracture_gradient_MPa_per_100m) ...
            && injection.fracture_gradient_MPa_per_100m > 0
        grad = injection.fracture_gradient_MPa_per_100m; % MPa/100m
        P_frac_MPa = 0.101325 + grad * depth_ref_m / 100.0;
        info.method = 'gradient_MPa_per_100m';
        info.method_str = sprintf('梯度%.2fMPa/100m@%.0fm', grad, depth_ref_m);
        return;
    end

    if isfield(injection, 'fracture_gradient_g_cm3') && isfinite(injection.fracture_gradient_g_cm3) ...
            && injection.fracture_gradient_g_cm3 > 0
        rho = 1000.0 * injection.fracture_gradient_g_cm3; % kg/m3
        P_frac_MPa = 0.101325 + rho * g0 * depth_ref_m / 1e6;
        info.method = 'gradient_g_cm3';
        info.method_str = sprintf('梯度%.2fg/cm³@%.0fm', injection.fracture_gradient_g_cm3, depth_ref_m);
        return;
    end

    % 破裂压力系数：Kf = P_frac / P_pore（常见现场给法）
    if isfield(injection, 'fracture_factor') && isfinite(injection.fracture_factor) && injection.fracture_factor > 0
        P_pore = getfield_or_default(injection, 'formation_pressure', NaN);
        if ~isfinite(P_pore)
            P_pore = getfield_or_default(injection, 'bh_target_pressure', NaN);
        end
        if isfinite(P_pore) && P_pore > 0
            P_frac_MPa = injection.fracture_factor * P_pore;
            info.method = 'factor';
            info.method_str = sprintf('Kf=%.3f×P_{pore}=%.2f', injection.fracture_factor, P_pore);
            return;
        end
    end
end

%% 液态下输校核（量级：以 CO2 饱和压力曲线判别）
function out = check_liquid_transport(depth_grid, sim, CO2, margin_req_MPa)
    if nargin < 4 || isempty(margin_req_MPa)
        margin_req_MPa = 0.0;
    end
    if nargin < 3 || isempty(CO2)
        CO2 = struct();
    end

    out = struct('ok', false, ...
        'min_margin_MPa', NaN, 'depth_min_m', NaN, 'T_min_C', NaN, 'P_min_MPa', NaN, ...
        'P_sat_min_MPa', NaN, 'max_beta_vapor', NaN);

    if nargin < 2 || ~isstruct(sim) || ~isfield(sim, 'T_co2') || ~isfield(sim, 'P_co2')
        return;
    end

    z = depth_grid(:);
    T = sim.T_co2(:);
    P = sim.P_co2(:);
    n = min([numel(z), numel(T), numel(P)]);
    if n <= 0
        return;
    end
    z = z(1:n);
    T = T(1:n);
    P = P(1:n);

    beta = NaN(n, 1);
    if isfield(sim, 'beta_vapor') && ~isempty(sim.beta_vapor)
        beta = sim.beta_vapor(:);
        beta = beta(1:n);
    end

    Tc = getfield_or_default(CO2, 'Tc', 31.1);

    P_sat = zeros(n, 1);
    for i = 1:n
        P_sat(i) = saturation_pressure_co2(T(i));
    end
    margin = P - P_sat;
    [min_margin, idx] = min(margin);

    out.min_margin_MPa = min_margin;
    out.depth_min_m = z(idx);
    out.T_min_C = T(idx);
    out.P_min_MPa = P(idx);
    out.P_sat_min_MPa = P_sat(idx);

    beta_valid = beta(isfinite(beta));
    if ~isempty(beta_valid)
        out.max_beta_vapor = max(beta_valid);
    end

    ok_margin = all(margin >= margin_req_MPa);
    ok_T = all(T < Tc - 1e-9);
    ok_beta = true;
    if ~isempty(beta_valid)
        ok_beta = all(beta_valid <= 0.01);
    end
    out.ok = ok_margin && ok_T && ok_beta;
end

%% 井口压力扫掠：同时检查全程液态 + 井底不超破裂上限
function sweep = sweep_inlet_pressure_constraints(depth_grid_base, dz_base, T_env_base, injection_base, model_base, well, thermal, CO2, const, P_bh_limit_MPa, cfg)
    if nargin < 11 || isempty(cfg)
        cfg = struct();
    end

    sweep = struct();
    sweep.ok = false;

    % 固定井口温度
    T_inlet_C = getfield_or_default(cfg, 'temp_inlet_C', NaN);
    if ~isfinite(T_inlet_C)
        T_inlet_C = getfield_or_default(injection_base, 'temp_inlet', NaN);
    end
    if ~isfinite(T_inlet_C)
        warning('[调参] 未提供有效 temp_inlet_C，跳过扫掠。');
        return;
    end

    % 扫描点
    P_list = getfield_or_default(cfg, 'pressure_list_MPa', []);
    if isempty(P_list)
        P_rng = getfield_or_default(cfg, 'pressure_range_MPa', [0.0 6.0]);
        if numel(P_rng) < 2 || ~all(isfinite(P_rng))
            P_rng = [0.0 6.0];
        end
        P_lo = min(P_rng(1), P_rng(2));
        P_hi = max(P_rng(1), P_rng(2));
        n_points = max(3, round(getfield_or_default(cfg, 'n_points', 41)));
        P_list = linspace(P_lo, P_hi, n_points);
    end
    P_list = unique(max(P_list(:), 0.0));
    P_list = sort(P_list);

    % 约束：液态裕度
    margin_req = getfield_or_default(cfg, 'liquid_margin_MPa', NaN);
    if ~isfinite(margin_req)
        lt0 = getfield_or_default(injection_base, 'liquid_transport', struct());
        margin_req = getfield_or_default(lt0, 'margin_MPa', 0.0);
    end
    margin_req = max(margin_req, 0.0);

    % 压力扫描使用较粗的深度网格。
    dz_sweep = max(getfield_or_default(cfg, 'dz_m', 5.0), dz_base);
    n_seg = max(1, ceil(well.total_depth / dz_sweep));
    dz_sweep = well.total_depth / n_seg;
    z_sweep = (0:dz_sweep:well.total_depth)';
    T_env_sweep = interp1(depth_grid_base(:), T_env_base(:), z_sweep, 'linear', 'extrap');

    % 扫描时关闭水合物子模块（仅为液态/破裂约束，不需要固相/沉积闭合）
    model_sweep = model_base;
    if getfield_or_default(cfg, 'disable_hydrate', true)
        model_sweep.hydrate.enabled = false;
    end
    if isfield(model_sweep, 'iteration') && isstruct(model_sweep.iteration)
        model_sweep.iteration.max_iter = 1;
    end

    simopt = struct();
    simopt.verbose = false;
    simopt.progress_every = 1000000000;
    simopt.show_progress = false;

    nP = numel(P_list);
    P_bh = NaN(nP, 1);
    min_margin = NaN(nP, 1);
    z_min = NaN(nP, 1);
    T_min = NaN(nP, 1);
    P_min = NaN(nP, 1);
    P_sat_min = NaN(nP, 1);
    ok_liquid = false(nP, 1);
    ok_fracture = true(nP, 1);

    fprintf('\n[调参] 井口T/P扫掠：T_inlet=%.1f°C | 扫描P_inlet_gauge=%.2f~%.2f MPa | margin=%.2f MPa | dz=%.1f m\n', ...
        T_inlet_C, P_list(1), P_list(end), margin_req, dz_sweep);
    if isfinite(P_bh_limit_MPa)
        fprintf('[调参] 破裂上限（井底）=%.2f MPa\n', P_bh_limit_MPa);
    else
        fprintf('[调参] 破裂上限（井底）=NaN（未提供/不可用）\n');
    end

    for i = 1:nP
        inj = injection_base;
        inj.temp_inlet = T_inlet_C;
        inj.pressure_inlet = P_list(i);
        inj.solve_inlet_pressure = false;
        if isfield(inj, 'liquid_transport') && isstruct(inj.liquid_transport)
            inj.liquid_transport.enforce = false;
        end

        sim = simulate_olga_coupled(z_sweep, dz_sweep, T_env_sweep, inj, model_sweep, well, thermal, CO2, const, simopt);
        P_bh(i) = sim.P_co2(end);

        lt = check_liquid_transport(z_sweep, sim, CO2, margin_req);
        min_margin(i) = lt.min_margin_MPa;
        z_min(i) = lt.depth_min_m;
        T_min(i) = lt.T_min_C;
        P_min(i) = lt.P_min_MPa;
        P_sat_min(i) = lt.P_sat_min_MPa;
        ok_liquid(i) = lt.ok;

        if isfinite(P_bh_limit_MPa) && isfinite(P_bh(i))
            ok_fracture(i) = (P_bh(i) <= P_bh_limit_MPa + 1e-9);
        end
    end

    ok_both = ok_liquid & ok_fracture;

    T = table(P_list, repmat(T_inlet_C, nP, 1), P_bh, min_margin, z_min, T_min, P_min, P_sat_min, ok_liquid, ok_fracture, ok_both, ...
        'VariableNames', {'P_inlet_gauge_MPa', 'T_inlet_C', 'P_bh_MPa', 'MinMargin_MPa', 'DepthMin_m', ...
        'T_at_MinMargin_C', 'P_at_MinMargin_MPa', 'P_sat_at_MinMargin_MPa', 'OkLiquid', 'OkFracture', 'OkBoth'});

    sweep.ok = true;
    sweep.table = T;
    sweep.P_bh_limit_MPa = P_bh_limit_MPa;
    sweep.liquid_margin_req_MPa = margin_req;

    idx_liq = find(ok_liquid, 1, 'first');
    if ~isempty(idx_liq)
        sweep.P_liquid_min_gauge_MPa = P_list(idx_liq);
    else
        sweep.P_liquid_min_gauge_MPa = NaN;
    end
    idx_frac = find(ok_fracture, 1, 'last');
    if ~isempty(idx_frac)
        sweep.P_fracture_max_gauge_MPa = P_list(idx_frac);
    else
        sweep.P_fracture_max_gauge_MPa = NaN;
    end
    idx_both = find(ok_both, 1, 'first');
    if ~isempty(idx_both)
        sweep.P_recommended_gauge_MPa = P_list(idx_both);
    else
        sweep.P_recommended_gauge_MPa = NaN;
    end

    if isfinite(sweep.P_recommended_gauge_MPa)
        fprintf('[调参] 可行：推荐 P_inlet_gauge≈%.2f MPa（同时满足液态+破裂）\n', sweep.P_recommended_gauge_MPa);
    else
        fprintf('[调参] 无同时满足解：液态所需 P_inlet_gauge_min≈%.2f MPa；破裂允许 P_inlet_gauge_max≈%.2f MPa\n', ...
            sweep.P_liquid_min_gauge_MPa, sweep.P_fracture_max_gauge_MPa);
    end

    if getfield_or_default(cfg, 'print_table', true)
        disp(T(:, {'P_inlet_gauge_MPa','P_bh_MPa','MinMargin_MPa','OkLiquid','OkFracture','OkBoth'}));
    end
end

%% 三元运算符辅助函数
function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

%% 注入流中CO2质量流量拆分
function plot_export = export_visualization_plotdata_excel(xlsx_path, options, depth_grid, sim, env_profile, well, injection, thermal, environment)
    if nargin < 1 || isempty(xlsx_path)
        xlsx_path = 'CO2_visualization_plotdata.xlsx';
    end
    if isstring(xlsx_path); xlsx_path = char(xlsx_path); end

    depth_grid = depth_grid(:);

    plot_export = struct();
    plot_export.file = xlsx_path;
    plot_export.tables = struct();
    plot_export.sheets = {};

    sheets = {};
    write_sheet = @(name) set_sheet_written(name);

    % Run info / manifest
    run_ts = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    model_ver = '1.0.1';
    plot_level = '';
    if nargin >= 2 && isstruct(options) && isfield(options, 'plot_level')
        try; plot_level = char(options.plot_level); catch; end
    end
    T_info = table({run_ts}, {model_ver}, {plot_level}, ...
        'VariableNames', {'RunTime', 'ModelVersion', 'PlotLevel'});
    writetable(T_info, xlsx_path, 'Sheet', 'RunInfo');
    plot_export.tables.RunInfo = T_info;
    write_sheet('RunInfo');

    % Fig01：井筒主结果（汇总）
    T_wb = table();
    try; T_wb.Depth_m = depth_grid; catch; end
    try; T_wb.T_CO2_C = sim.T_co2(:); catch; end
    try; T_wb.T_env_C = sim.T_env(:); catch; end
    try; T_wb.P_CO2_MPa = sim.P_co2(:); catch; end
    try; T_wb.Rho_CO2_kg_m3 = sim.rho_co2(:); catch; end
    try; T_wb.Mu_CO2_uPas = sim.mu_co2(:) * 1e6; catch; end
    try; T_wb.k_CO2_W_m_K = sim.k_co2(:); catch; end
    try; T_wb.V_CO2_m_s = sim.v_co2(:); catch; end
    try; T_wb.Re = sim.Re_co2(:); catch; end
    try; T_wb.Phase = sim.phase(:); catch; end
    try; T_wb.Subcool_C = sim.subcooling(:); catch; end
    % 液态裕度（P - Psat）：用于全程液态约束对照（Fig01 缺失时可由 T/P 回算）
    try
        if ismember('T_CO2_C', T_wb.Properties.VariableNames) && ismember('P_CO2_MPa', T_wb.Properties.VariableNames)
            T_tmp = T_wb.T_CO2_C(:);
            P_tmp = T_wb.P_CO2_MPa(:);
            n_tmp = min(numel(T_tmp), numel(P_tmp));
            T_tmp = T_tmp(1:n_tmp);
            P_tmp = P_tmp(1:n_tmp);

            P_sat_tmp = NaN(n_tmp, 1);
            for i = 1:n_tmp
                if isfinite(T_tmp(i))
                    P_sat_tmp(i) = saturation_pressure_co2(T_tmp(i));
                end
            end
            T_wb.P_sat_MPa = P_sat_tmp;
            T_wb.LiquidMargin_MPa = P_tmp - P_sat_tmp;
        end
    catch
    end
    % 传热相关（Fig01 中常用的导热/换热系数等）
    try; T_wb.h_inner_W_m2_K = sim.h_inner(:); catch; end
    try; T_wb.U_W_m2_K = sim.U_total(:); catch; end
    try; T_wb.Q_heat_W = sim.Q_heat(:); catch; end
    try; T_wb.Q_hydrate_rxn_W = sim.Q_hydrate_rxn(:); catch; end
    % 环境压力（与 Fig01 图中的环境静水压基线一致）
    try
        rho_env = 1028.0;
        g_env = 9.80665;
        z_ref = getfield_or_default(well, 'platform_elevation', 0.0);
        depth_below_sea_m = max(depth_grid - z_ref, 0.0);
        T_wb.P_env_MPa = 0.101325 + rho_env * g_env * depth_below_sea_m / 1e6;
    catch
    end
    % 地层基线：地温梯度 + 孔隙压力（用于曲线对照；压力为泥线静水压->井底孔压的线性插值）
    try
        z_mud = getfield_or_default(well, 'mudline_depth', NaN);
        z_td = getfield_or_default(well, 'total_depth', NaN);
        if ~isfinite(z_td); z_td = depth_grid(end); end

        if isfinite(z_mud) && isfinite(z_td) && z_td > z_mud
            msk_form = depth_grid > z_mud;
            z_form_m = depth_grid(msk_form) - z_mud;

            T_form = thermal.seabed_temp + thermal.geothermal_gradient * z_form_m;
            T_form_col = NaN(size(depth_grid));
            T_form_col(msk_form) = T_form;
            T_wb.T_formation_C = T_form_col;

            P_pore_col = NaN(size(depth_grid));
            P_bh_pore_MPa = getfield_or_default(injection, 'formation_pressure', NaN);
            if isfinite(P_bh_pore_MPa)
                rho_env = 1028.0;
                g_env = 9.80665;
                z_sea = getfield_or_default(well, 'seawater_depth', NaN);
                if ~isfinite(z_sea)
                    z_ref = getfield_or_default(well, 'platform_elevation', 0.0);
                    if isfinite(z_ref)
                        z_sea = z_mud - z_ref;
                    end
                end

                P_mud_hydro_MPa = NaN;
                if isfinite(z_sea) && z_sea >= 0
                    P_mud_hydro_MPa = 0.101325 + rho_env * g_env * z_sea / 1e6;
                end
                z_form_max_m = max(z_td - z_mud, 1e-12);
                if isfinite(P_mud_hydro_MPa)
                    P_pore_col(msk_form) = P_mud_hydro_MPa + (P_bh_pore_MPa - P_mud_hydro_MPa) * (z_form_m / z_form_max_m);
                end
            end
            T_wb.P_pore_linear_MPa = P_pore_col;
        end
    catch
    end
    writetable(T_wb, xlsx_path, 'Sheet', 'Fig01_Wellbore');
    plot_export.tables.Fig01_Wellbore = T_wb;
    write_sheet('Fig01_Wellbore');

    % Fig02b：地层剖面（地温梯度 + 孔隙压力基线）
    try
        z_mud = getfield_or_default(well, 'mudline_depth', NaN);
        z_td = getfield_or_default(well, 'total_depth', NaN);
        if ~isfinite(z_td); z_td = depth_grid(end); end

        if isfinite(z_mud) && isfinite(z_td) && z_td > z_mud
            msk_form = depth_grid > z_mud;
            z = depth_grid(msk_form);
            z_form_m = z - z_mud;

            T_geo_C = thermal.seabed_temp + thermal.geothermal_gradient * z_form_m;

            rho_env = 1028.0;
            g_env = 9.80665;
            z_sea = getfield_or_default(well, 'seawater_depth', NaN);
            if ~isfinite(z_sea)
                z_ref = getfield_or_default(well, 'platform_elevation', 0.0);
                if isfinite(z_ref)
                    z_sea = z_mud - z_ref;
                end
            end

            P_mud_hydro_MPa = NaN;
            if isfinite(z_sea) && z_sea >= 0
                P_mud_hydro_MPa = 0.101325 + rho_env * g_env * z_sea / 1e6;
            end
            P_bh_pore_MPa = getfield_or_default(injection, 'formation_pressure', NaN);
            P_pore_linear_MPa = NaN(size(z));
            if isfinite(P_mud_hydro_MPa) && isfinite(P_bh_pore_MPa)
                z_form_max_m = max(z_td - z_mud, 1e-12);
                P_pore_linear_MPa = P_mud_hydro_MPa + (P_bh_pore_MPa - P_mud_hydro_MPa) * (z_form_m / z_form_max_m);
            end

            T_form = table(z, z_form_m, T_geo_C, P_pore_linear_MPa, ...
                'VariableNames', {'Depth_mudlineRef_m', 'Depth_below_mudline_m', 'T_geo_C', 'P_pore_linear_MPa'});
            writetable(T_form, xlsx_path, 'Sheet', 'Fig02_Formation');
            plot_export.tables.Fig02_Formation = T_form;
            write_sheet('Fig02_Formation');
        end
    catch
    end

    % 工作表索引
    if ~isempty(sheets)
        C = [{'SheetName'}; sheets(:)];
        writecell(C, xlsx_path, 'Sheet', 'Manifest');
    else
        C = {'SheetName'; 'RunInfo'};
        writecell(C, xlsx_path, 'Sheet', 'Manifest');
    end

    plot_export.sheets = sheets;
    plot_export.manifest = C;

    function set_sheet_written(name)
        sheets{end+1,1} = name; %#ok<AGROW>
    end
end

function fig = plot_wellbore_summary(depth_grid, T_env, sim, well, injection, CO2)
    if nargin < 6 || isempty(CO2); CO2 = struct(); end
    if nargin < 5 || isempty(injection); injection = struct(); end
    if nargin < 4 || isempty(well); well = struct(); end
    if nargin < 3 || ~isstruct(sim)
        fig = [];
        return;
    end

    depth_grid = depth_grid(:);
    T_env = T_env(:);
    T_co2 = sim.T_co2(:);
    P_co2 = sim.P_co2(:);

    Tc = getfield_or_default(CO2, 'Tc', 31.1);
    Pc = getfield_or_default(CO2, 'Pc', 7.38);

    z_mud = getfield_or_default(well, 'mudline_depth', NaN);
    z_packer = getfield_or_default(well, 'packer_depth', NaN);
    z_td = getfield_or_default(well, 'total_depth', max(depth_grid));

    color_temp = [0.85 0.33 0.10];
    color_pressure = [0.00 0.45 0.74];
    color_env = [0.47 0.67 0.19];

    fig = safe_figure(1, 'Position', [120 80 1500 820], 'Color', 'w', 'Name', '井筒主结果（汇总）', 'NumberTitle', 'off');
    % 避免重复运行脚本时图窗复用导致对象叠加（尤其是双 X 轴叠加轴）。
    try; clf(fig); catch; end
    tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

    % 温度
    nexttile(tl, 1);
    h_env = plot(T_env, depth_grid, '--', 'Color', color_env, 'LineWidth', 2.0); hold on;
    h_co2 = plot(T_co2, depth_grid, '-', 'Color', color_temp, 'LineWidth', 2.6);
    h_eq = [];
    if isfield(sim, 'hydrate_Teq_pp_C') && ~isempty(sim.hydrate_Teq_pp_C)
        h_eq = plot(sim.hydrate_Teq_pp_C(:), depth_grid, ':', 'Color', [0.45 0.45 0.95], 'LineWidth', 1.8);
        lgd = legend([h_env, h_co2, h_eq], {'环境', 'CO2', '水合物平衡线(分压)'}, 'Location', 'best');
    elseif isfield(sim, 'T_hydrate') && ~isempty(sim.T_hydrate)
        h_eq = plot(sim.T_hydrate(:), depth_grid, ':', 'Color', [0.45 0.45 0.95], 'LineWidth', 1.8);
        lgd = legend([h_env, h_co2, h_eq], {'环境', 'CO2', '水合物平衡线'}, 'Location', 'best');
    else
        lgd = legend([h_env, h_co2], {'环境', 'CO2'}, 'Location', 'best');
    end
    % 防止后续 xline/yline 被自动加入图例（data1/data2/...）
    try
        set(lgd, 'AutoUpdate', 'off');
    catch
    end

    h_tc = xline(Tc, 'k:', 'LineWidth', 1.2, 'Label', 'Tc');
    try; set(h_tc, 'HandleVisibility', 'off'); catch; end
    add_depth_markers(gca, z_mud, z_packer, z_td);
    set(gca, 'YDir', 'reverse');
    xlabel('温度（°C）'); ylabel('深度（m）');
    title('温度', 'FontWeight', 'bold');
    grid on; box on;

    % 压力
    nexttile(tl, 2);
    h_p_co2 = plot(P_co2, depth_grid, '-', 'Color', color_pressure, 'LineWidth', 2.6); hold on;

    % 环境压力：以海水静水压近似（上部为海水段；下部作为静水压基线延拓，用于对照压差/欠平衡）
    rho_env = 1028.0;
    g_env = 9.80665;
    z_ref = getfield_or_default(well, 'platform_elevation', 0.0);
    depth_below_sea_m = max(depth_grid - z_ref, 0.0);
    P_env_MPa = 0.101325 + rho_env * g_env * depth_below_sea_m / 1e6;
    h_p_env = plot(P_env_MPa, depth_grid, '--', 'Color', color_env, 'LineWidth', 2.0);

    h_pc = xline(Pc, 'k:', 'LineWidth', 1.2, 'Label', 'Pc');
    try; set(h_pc, 'HandleVisibility', 'off'); catch; end
    add_depth_markers(gca, z_mud, z_packer, z_td);
    target_bh_MPa = getfield_or_default(injection, 'bh_target_pressure', NaN);
    if ~isfinite(target_bh_MPa)
        target_bh_MPa = getfield_or_default(injection, 'formation_pressure', NaN);
    end
    if isfinite(target_bh_MPa)
        h_target = plot(target_bh_MPa, z_td, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        try; set(h_target, 'HandleVisibility', 'off'); catch; end
        try
            text(target_bh_MPa, z_td, sprintf('  目标井底 %.2f MPa', target_bh_MPa), ...
                'Color', [0.75 0.0 0.0], 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        catch
        end
    end

    % 实际井底压力（曲线终点；用于和目标/约束区分，避免误读红点）
    try
        P_bh_actual_MPa = P_co2(end);
        h_bh = plot(P_bh_actual_MPa, z_td, 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        try; set(h_bh, 'HandleVisibility', 'off'); catch; end
        text(P_bh_actual_MPa, z_td, sprintf('  实际井底 %.2f MPa', P_bh_actual_MPa), ...
            'Color', [0.10 0.10 0.10], 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    catch
    end

    % 破裂上限（井底约束线，基于 estimate_fracture_pressure_MPa；用于直观看是否超压）
    try
        P_frac_MPa = estimate_fracture_pressure_MPa(well, injection, struct('g', 9.80665));
        P_bh_limit_MPa = NaN;
        if isfinite(P_frac_MPa)
            P_bh_limit_MPa = P_frac_MPa - max(getfield_or_default(injection, 'fracture_margin_MPa', 0.0), 0.0);
        end
        if isfinite(P_bh_limit_MPa)
            h_frac = xline(P_bh_limit_MPa, 'r--', 'LineWidth', 1.3, ...
                'Label', sprintf('破裂上限 %.2f', P_bh_limit_MPa));
            try; set(h_frac, 'HandleVisibility', 'off'); catch; end
        end
    catch
    end

    lgd_p = legend([h_p_co2, h_p_env], {'CO2', '环境静水压(近似)'}, 'Location', 'best');
    try; set(lgd_p, 'AutoUpdate', 'off'); catch; end
    set(gca, 'YDir', 'reverse');
    xlabel('压力（MPa，绝压）'); ylabel('深度（m）');
    title('压力', 'FontWeight', 'bold');
    grid on; box on;

    % CO2 相图（P–T）+ 井筒轨迹
    nexttile(tl, 3);
    try
        T_traj = T_co2(:);
        P_traj = P_co2(:); % MPa（绝压）

        T_pad = 2.0;
        P_pad = 0.5;
        T_min = floor(min(T_traj) - T_pad);
        T_max = ceil(max(T_traj) + T_pad);
        P_min = max(0.1, floor((min(P_traj) - P_pad) * 10) / 10);
        P_max = ceil((max(P_traj) + P_pad) * 10) / 10;

        % 坐标范围覆盖临界点 Tc、Pc。
        T_min = min(T_min, floor(Tc - 5));
        T_max = max(T_max, ceil(Tc + 5));
        P_min = min(P_min, floor((Pc - 1) * 10) / 10);
        P_max = max(P_max, ceil((Pc + 1) * 10) / 10);

        nT = 240;
        nP = 220;
        T_grid = linspace(T_min, T_max, nT);
        P_grid = linspace(P_min, P_max, nP);
        P_sat_T = arrayfun(@saturation_pressure_co2, T_grid);

        [TT, PP] = meshgrid(T_grid, P_grid);
        P_sat_mat = repmat(P_sat_T, nP, 1);

        % 相态分区（纯CO2工程近似）：气态/液态/超临界（按需移除两相图例与分区）
        phase_code = ones(size(TT)); % 1=气态
        phase_code(TT < Tc & PP > P_sat_mat) = 2; % 2=液态
        phase_code(TT > Tc & PP > Pc) = 3; % 3=超临界

        imagesc(T_grid, P_grid, phase_code);
        set(gca, 'YDir', 'normal');
        colormap(gca, [1.00 0.55 0.55; 0.50 0.70 1.00; 0.70 0.55 0.95]);
        caxis([0.5 3.5]);
        cb = colorbar;
        cb.Ticks = 1:3;
        cb.TickLabels = {'Gas', 'Liquid', 'Supercritical'};
        cb.TickDirection = 'out';
        hold on;

        % 饱和曲线 + 液态裕度线（便于对照全程液态约束）
        idx_sat = T_grid <= Tc;
        plot(T_grid(idx_sat), P_sat_T(idx_sat), 'k-', 'LineWidth', 2.4);
        if isfield(injection, 'liquid_transport') && isstruct(injection.liquid_transport) ...
                && isfield(injection.liquid_transport, 'enabled') && injection.liquid_transport.enabled ...
                && isfield(injection.liquid_transport, 'margin_MPa') && isfinite(injection.liquid_transport.margin_MPa) ...
                && injection.liquid_transport.margin_MPa > 0
            plot(T_grid(idx_sat), P_sat_T(idx_sat) + injection.liquid_transport.margin_MPa, 'k--', 'LineWidth', 1.0);
        end

        % 绘图采用分压平衡线，逸度用于对照。
        if isfield(sim, 'hydrate_Teq_pp_C') && ~isempty(sim.hydrate_Teq_pp_C)
            plot(sim.hydrate_Teq_pp_C(:), P_traj, ':', 'Color', [0.45 0.45 0.95], 'LineWidth', 1.8);
        elseif isfield(sim, 'T_hydrate') && ~isempty(sim.T_hydrate)
            plot(sim.T_hydrate(:), P_traj, ':', 'Color', [0.45 0.45 0.95], 'LineWidth', 1.8);
        end

        % 井筒轨迹（用双层线条保证在任意底色上可见）
        plot(T_traj, P_traj, 'w-', 'LineWidth', 3.0);
        plot(T_traj, P_traj, 'k-', 'LineWidth', 1.6);
        scatter(T_traj(1), P_traj(1), 36, 'k', 'filled');
        scatter(T_traj(end), P_traj(end), 52, 'k', 'filled', 'Marker', '^');

        plot(Tc, Pc, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'w');
        xlabel('Temperature (°C)');
        ylabel('Pressure (MPa, abs)');
        grid on; box on;
    catch
        text(0.5, 0.5, 'P–T 相图绘制失败', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        axis off;
    end

    % 密度 + 流速
    nexttile(tl, 4);
    if isfield(sim, 'rho_co2') && isfield(sim, 'v_co2') && ~isempty(sim.rho_co2) && ~isempty(sim.v_co2)
        ax = gca;
        plot_dual_xaxis(ax, depth_grid, sim.rho_co2(:), '密度（kg/m³）', [0.10 0.10 0.10], ...
            sim.v_co2(:), '流速（m/s）', [0.00 0.45 0.74]);
        add_depth_markers(ax, z_mud, z_packer, z_td);
        % 仅保留双 X 轴标签（顶轴：流速），避免标题与顶轴文字重叠。
    else
        text(0.5, 0.5, '无密度/流速结果', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        axis off;
    end

    % 黏度 + 导热系数
    nexttile(tl, 5);
    if isfield(sim, 'mu_co2') && isfield(sim, 'k_co2') && ~isempty(sim.mu_co2) && ~isempty(sim.k_co2)
        mu_mPa_s = 1e3 * sim.mu_co2(:);
        ax = gca;
        plot_dual_xaxis(ax, depth_grid, mu_mPa_s, '黏度（mPa·s）', [0.85 0.33 0.10], ...
            sim.k_co2(:), '导热系数（W/m·K）', [0.30 0.65 0.30]);
        add_depth_markers(ax, z_mud, z_packer, z_td);
        % 仅保留双 X 轴标签（顶轴：导热系数），避免标题与顶轴文字重叠。
    else
        text(0.5, 0.5, '无黏度/导热结果', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        axis off;
    end

    % 过冷度
    nexttile(tl, 6);
    if isfield(sim, 'subcooling') && ~isempty(sim.subcooling)
        subcool = sim.subcooling(:);
        subcool_valid = subcool(isfinite(subcool));
        x_max = max([6.0; subcool_valid; 0.0]);

        hold on;
        patch([0 2 2 0], [0 0 z_td z_td], [0.85 1.00 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
        patch([2 5 5 2], [0 0 z_td z_td], [1.00 1.00 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
        patch([5 x_max x_max 5], [0 0 z_td z_td], [1.00 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.45);

        plot(subcool, depth_grid, 'b-', 'LineWidth', 2.4);
        xline(0, 'k:', 'LineWidth', 1.2);
        xline(2, 'g-', 'LineWidth', 1.4);
        xline(5, 'r-', 'LineWidth', 1.4);
        add_depth_markers(gca, z_mud, z_packer, z_td);
        set(gca, 'YDir', 'reverse');

        x_min = min([0.0; subcool_valid]);
        xlim([x_min x_max]);
        xlabel('过冷度（°C）'); ylabel('深度（m）');
        title('过冷度', 'FontWeight', 'bold');
        grid on; box on;
    else
        text(0.5, 0.5, '未启用水合物模块', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        axis off;
    end

    Q_m3_h = getfield_or_default(injection, 'volume_flow_rate', NaN);
    P_in_g = P_co2(1) - 0.101325;
    sgtitle(sprintf('深水 CO2 注入 | 流量=%.2f m³/h | 井深=%.0f m | 井口 %.1f°C/%.2fMPa(g) [%.2f abs] -> 井底 %.1f°C/%.2fMPa(abs)', ...
        Q_m3_h, z_td, T_co2(1), P_in_g, P_co2(1), T_co2(end), P_co2(end)), 'FontWeight', 'bold');
    install_axes_copy_vector_context_menu(fig);

    % tiledlayout 在首次 drawnow 后可能微调子图位置；同步双 X 轴顶轴位置，保证两张双轴图顶轴高度一致。
    try
        drawnow;
        sync_dual_xaxis_positions(fig);
    catch
    end
end

function hs = add_depth_markers(ax, z_mud, z_packer, z_td)
    hs = gobjects(0);
    if nargin < 1 || isempty(ax) || ~ishandle(ax)
        return;
    end

    if nargin >= 2 && isfinite(z_mud)
        h = yline(ax, z_mud, 'b--', 'LineWidth', 1.0);
        try; set(h, 'HandleVisibility', 'off'); catch; end
        hs(end+1, 1) = h; %#ok<AGROW>
    end

    if nargin >= 3 && isfinite(z_packer)
        h = yline(ax, z_packer, 'm--', 'LineWidth', 1.0);
        try; set(h, 'HandleVisibility', 'off'); catch; end
        hs(end+1, 1) = h; %#ok<AGROW>
    end

    if nargin >= 4 && isfinite(z_td)
        h = yline(ax, z_td, 'k-.', 'LineWidth', 0.9, ...
            'Label', sprintf('井底 %.0f m', z_td), ...
            'LabelHorizontalAlignment', 'left', ...
            'LabelVerticalAlignment', 'top');
        try; set(h, 'HandleVisibility', 'off'); catch; end
        hs(end+1, 1) = h; %#ok<AGROW>
    end
end

function ax_top = plot_dual_xaxis(ax_bottom, depth_m, x_bottom, x_bottom_label, x_bottom_color, x_top, x_top_label, x_top_color)
    depth_m = depth_m(:);
    x_bottom = x_bottom(:);
    x_top = x_top(:);

    % 重绘前移除同一位置的顶轴。
    if isstruct(ax_bottom.UserData) && isfield(ax_bottom.UserData, 'dual_xaxis') ...
            && isfield(ax_bottom.UserData.dual_xaxis, 'ax_top') && isvalid(ax_bottom.UserData.dual_xaxis.ax_top)
        try; delete(ax_bottom.UserData.dual_xaxis.ax_top); catch; end
    end

    h_bottom = plot(ax_bottom, x_bottom, depth_m, '-', 'Color', x_bottom_color, 'LineWidth', 2.2);
    hold(ax_bottom, 'on');
    set(ax_bottom, 'YDir', 'reverse');
    xlabel(ax_bottom, x_bottom_label, 'Color', x_bottom_color, 'Interpreter', 'none');
    ylabel(ax_bottom, '深度（m）', 'Interpreter', 'none');
    grid(ax_bottom, 'on');
    box(ax_bottom, 'on');
    ax_bottom.TickDir = 'out';
    ax_bottom.XColor = x_bottom_color;

    fig = ancestor(ax_bottom, 'figure');
    ax_top = axes('Parent', fig, 'Units', ax_bottom.Units, 'Position', ax_bottom.Position, 'Color', 'none', ...
        'XAxisLocation', 'top', 'YAxisLocation', 'right');
    % 关键：避免 plot/newplot 将 Color 从 'none' 重置为白底，从而遮挡底轴曲线。
    ax_top.NextPlot = 'add';
    ax_top.FontName = ax_bottom.FontName;
    ax_top.FontSize = ax_bottom.FontSize;
    ax_top.LineWidth = ax_bottom.LineWidth;

    % 顶轴曲线用虚线。
    h_top = plot(ax_top, x_top, depth_m, '--', 'Color', x_top_color, 'LineWidth', 1.8);
    ax_top.Color = 'none';
    set(ax_top, 'YDir', 'reverse', 'YTick', [], 'Box', 'off');
    xlabel(ax_top, x_top_label, 'Color', x_top_color, 'Interpreter', 'none');
    ax_top.TickDir = 'out';
    ax_top.XColor = x_top_color;
    ax_top.YColor = 'none';
    ax_top.HitTest = 'off';
    ax_top.XAxisLocation = 'top';
    linkaxes([ax_bottom ax_top], 'y');

    ax_bottom.UserData.dual_xaxis = struct('ax_top', ax_top, 'h_bottom', h_bottom, 'h_top', h_top);
    axes(ax_bottom);
end

function sync_dual_xaxis_positions(fig)
    if nargin < 1 || isempty(fig) || ~ishandle(fig)
        return;
    end
    try
        if ~strcmp(get(fig, 'Type'), 'figure')
            return;
        end
    catch
        return;
    end

    axs = findall(fig, 'Type', 'axes');
    for i = 1:numel(axs)
        ax = axs(i);
        try
            ud = ax.UserData;
        catch
            continue;
        end
        if ~isstruct(ud) || ~isfield(ud, 'dual_xaxis')
            continue;
        end
        dx = ud.dual_xaxis;
        if ~isstruct(dx) || ~isfield(dx, 'ax_top') || ~isvalid(dx.ax_top)
            continue;
        end
        ax_top = dx.ax_top;
        try
            ax_top.Units = ax.Units;
            ax_top.Position = ax.Position;
        catch
        end
    end
end
