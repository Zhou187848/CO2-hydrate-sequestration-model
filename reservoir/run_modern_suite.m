function run_modern_suite(group,partition,parts,overrides)
% 按分区运行地层参数扫描或数值对照。
if nargin<2, partition=1; end
if nargin<3, parts=1; end
if nargin<4, overrides=struct(); end
root=fileparts(mfilename('fullpath')); addpath(root);
cases={};
switch group
    case 'grid'
        for phi=[.25,.30,.35,.40]
            for k=[50,80,100,140,200,300,350,400,500]
                cases{end+1}=struct('name',sprintf('grid_phi%02d_k%03d',round(phi*100),k),...
                    'sand_phi',phi,'sand_k',k*9.869233e-16); %#ok<AGROW>
            end
        end
    case 'checks'
        cases={struct('name','matched_1d','dimension',1),...
            struct('name','farfield_1000','radius',1000,'nr',58),...
            struct('name','domain_deep','bottom',400),...
            struct('name','salt_DB','salt_equilibrium','DB'),...
            struct('name','kinetic_0p1','kinetic_multiplier',.1),...
            struct('name','kinetic_10','kinetic_multiplier',10),...
            struct('name','lambda_half','lambda_sand',1.25,'lambda_seal',1.25),...
            struct('name','lambda_double','lambda_sand',5,'lambda_seal',5),...
            struct('name','rockcp_half','rock_cp',400),...
            struct('name','rockcp_double','rock_cp',1600),...
            struct('name','seal_low','seal_k',.13*9.869233e-16),...
            struct('name','seal_high','seal_k',10*9.869233e-16),...
            struct('name','isotropic','anisotropy',1),...
            struct('name','lower_completion','inject_top',140,'inject_bottom',150),...
            struct('name','shutin_30d','duration',30*86400,'inject_end_s',10*86400,'snapshot_days',[1,5,10,11,15,20,30])};
    otherwise
        error('group 必须是 grid / checks');
end
for j=partition:parts:numel(cases)
    base=overrides;
    names=fieldnames(cases{j});
    for k=1:numel(names), base.(names{k})=cases{j}.(names{k}); end
    c=modern_configuration(base);
    file=fullfile(root,c.results_folder,[c.name '.mat']);
    if isfile(file)
        s=load(file,'out'); o=s.out;
        fields=fieldnames(c);
        for q=1:numel(fields), assert(isequaln(c.(fields{q}),o.cfg.(fields{q})),'缓存参数不匹配: %s/%s',c.name,fields{q}); end
        assert(abs(o.history(end,1)-c.duration)<1e-6,'不完整结果不能复用。');
        fprintf('Verified completed %s\n',c.name);
    else
        o=axisymmetric_co2(c);
    end
    assert(abs(o.summary.CO2_defect_kg)/o.summary.injected_kg<1e-6,'CO2 守恒不满足容差。');
    fprintf('COMPLETED %s partition %d/%d case %d/%d\n',group,partition,parts,j,numel(cases));
end
end
