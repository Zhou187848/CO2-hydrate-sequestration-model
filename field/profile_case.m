function out = profile_case(user)
% 二维轴对称 CO2–盐水–水合物全隐式有限体积模型。
% 饱和度以初始孔隙体积为基准。
% 两移动相为不互溶近似；未计 CO2 溶解、盐析和地质力学。
if nargin == 0, user = struct(); end
root = fileparts(mfilename('fullpath'));
if isempty(which('initVariablesADI'))
    error('Start MRST 2025b before running this solver.');
end
user = modern_configuration(user);
assert(strcmp(user.closure,'li2024_moridis2025_v1'),'Unsupported hydrate closure.');
cfg = defaults();
names = fieldnames(user);
for j=1:numel(names)
    assert(isfield(cfg,names{j}), '未知配置项: %s',names{j});
    cfg.(names{j})=user.(names{j});
end
% 可从本地 MAT 文件读取层顶和层底深度（m）。
if any(~isfinite([cfg.sand_top,cfg.sand_bottom]))
    input_file=getenv('CO2_MODEL_GEOMETRY');
    assert(isfile(input_file),'Supply sand_top and sand_bottom in metres, or set CO2_MODEL_GEOMETRY to a local MAT file containing both scalars.');
    geom=load(input_file,'sand_top','sand_bottom');
    if ~isfinite(cfg.sand_top), cfg.sand_top=geom.sand_top; end
    if ~isfinite(cfg.sand_bottom), cfg.sand_bottom=geom.sand_bottom; end
end
assert(isfinite(cfg.sand_top) && cfg.sand_top>=0 && cfg.sand_bottom>cfg.sand_top && cfg.bottom>=cfg.sand_bottom);
raw = load(fullfile(fileparts(root),'reservoir','properties.mat'));
prop = make_properties(raw);
m = make_model(cfg,prop);
x = [m.p0/1e6, zeros(m.n,2), m.t0, cfg.salt*ones(m.n,1)];
if cfg.initial_gas>0, x(:,2)=cfg.initial_gas; end
if cfg.initial_hydrate>0, x(:,3)=cfg.initial_hydrate; end
old = inventory(x,m,prop,cfg);
initial = sum(old(:,[1,2,3,5]),1);
t=0; dt=cfg.dt_initial; accepted=0; rejected=0;
budget=zeros(1,4); injection=0; hist=[]; snapshots={}; diagnostics=[];
clockstart=tic;
while t < cfg.duration-1e-8
    dt=min([dt,cfg.dt_max,cfg.duration-t]);
    future=cfg.snapshot_days*86400; future=future(future>t+1e-6);
    if ~isempty(future), dt=min(dt,min(future)-t); end
    if t<cfg.inject_end_s, dt=min(dt,cfg.inject_end_s-t); end
    if ~isempty(cfg.field_schedule)
        sc=cfg.field_schedule; cuts=sc(:,1:2); cuts=cuts(cuts>t+1e-7);
        if ~isempty(cuts), dt=min(dt,min(cuts)-t); end
        idx=find(t+dt/2>=sc(:,1) & t+dt/2<sc(:,2)); rate=0;
        if ~isempty(idx), rate=sc(idx,3); end
        amount=rate*dt;
        cfg.inject_temperature=interp1(cfg.field_temp(:,1),cfg.field_temp(:,2),t+dt/2,'linear');
    else
        [rate,amount]=injection_rate(t,dt,cfg);
    end
    [xn,ok,it,res]=newton(x,old,dt,rate,m,prop,cfg);
    if ~ok
        rejected=rejected+1;
        dt=dt/2;
        if cfg.verbose, fprintf('%s reject t=%.5fd dt=%.3fs res=%.3g\n',cfg.name,t/86400,dt,res); end
        assert(dt>=cfg.dt_min, 'Newton 未收敛，禁止输出伪完成结果。t=%.6g d, residual=%.4g', t/86400,res);
        continue
    end
    [~,flux,source,wellinfo] = equations(xn,old,dt,rate,m,prop,cfg,false);
    x=xn; t=t+dt; accepted=accepted+1;
    budget=budget+dt*(source-sum(flux(m.boundary,:),1));
    injection=injection+amount;
    now=inventory(x,m,prop,cfg);
    totals=sum(now(:,[1,2,3,5]),1);
    defects=totals-initial-budget;
    mh=sum(cfg.wc*cfg.rho_h*m.pv.*x(:,3));
    mg=sum(now(:,1))-mh;
    vgas=m.pv.*x(:,2);
    if sum(vgas)>0, zg=sum(vgas.*m.z)/sum(vgas); else, zg=NaN; end
    dp=max((x(m.inj,1)*1e6-m.p0(m.inj)))/1e6;
    hist(end+1,:)=[t,dt,it,res,injection,mh,mg,max(x(:,3)),max(x(:,2)),dp,zg,defects]; %#ok<AGROW>
    rf_g=max([0;m.r(x(:,2)>=.05)]); rf_h=max([0;m.r(x(:,3)>=.05)]);
    diagnostics(end+1,:)=[t,wellinfo,rf_g,rf_h,max(x(:,5)),...
        min((1-x(:,3)).^cfg.permeability_exponent)]; %#ok<AGROW>
    if any(abs(t-cfg.snapshot_days*86400)<dt/2) || t>=cfg.duration
        snapshots{end+1}=struct('t',t,'x',x); %#ok<AGROW>
    end
    if cfg.verbose && (mod(accepted,20)==0 || t>=cfg.duration)
        fprintf('%s t=%.4fd steps=%d dt=%.0fs hydrate=%.4f%% Sh=%.4f dp=%.4fMPa Cdef=%.3gkg\n',...
            cfg.name,t/86400,accepted,dt,100*mh/max(injection,eps),max(x(:,3)),dp,defects(1));
    end
    old=now;
    if it<=6, dt=dt*1.5; elseif it>12, dt=dt*.75; else, dt=dt*1.15; end
end
out=struct('cfg',cfg,'x',x,'initial',initial,'budget',budget,'history',hist,...
    'snapshots',{snapshots},'accepted',accepted,'rejected',rejected,'elapsed',toc(clockstart));
out.grid=rmfield(m,{'D','a','b','rh1','rh2','heat','hydraulic','ghost_x','face_z','boundary','inj_weight'});
out.history_columns={'time_s','dt_s','Newton_iterations','scaled_residual','injected_CO2_kg',...
    'hydrate_CO2_kg','mobile_CO2_kg','peak_Sh','peak_Sg','max_cell_overpressure_MPa',...
    'gas_volume_centroid_bml_m','CO2_defect_kg','water_defect_kg','salt_defect_kg','energy_defect_J'};
out.final_inventory=sum(old(:,[1,2,3,5]),1);
out.diagnostics=diagnostics;
out.diagnostic_columns={'time_s','reference_well_pressure_MPa','injection_mass_sum_kg_s',...
    'front_CO2_S005_m','front_hydrate_S005_m','max_salt_mass_fraction','min_k_over_k0'};
out.operator_checks=check_geometry(m,cfg);
out.summary=struct('name',cfg.name,'nr',m.nr,'nz',m.nz,'cells',m.n,...
    'injected_kg',injection,'hydrate_fraction_pct',100*mh/max(injection,eps),...
    'peak_Sh',max(x(:,3)),'peak_Sg',max(x(:,2)),'overpressure_MPa',dp,'centroid_bml_m',zg,...
    'CO2_defect_kg',defects(1),'water_defect_kg',defects(2),'salt_defect_kg',defects(3),...
    'energy_defect_J',defects(4),'min_T_C',min(x(:,4)),'max_T_C',max(x(:,4)),...
    'max_salt',max(x(:,5)),'steps',accepted,'rejected',rejected,'elapsed_s',out.elapsed);
if cfg.save
    folder=fullfile(root,cfg.results_folder); if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,[cfg.name '.mat']),'out','-v7');
    writetable(array2table(hist,'VariableNames',out.history_columns),fullfile(folder,[cfg.name '_history.csv']));
    writetable(array2table(diagnostics,'VariableNames',out.diagnostic_columns),fullfile(folder,[cfg.name '_diagnostics.csv']));
    fid=fopen(fullfile(folder,[cfg.name '_summary.json']),'w'); fprintf(fid,'%s',jsonencode(out.summary,PrettyPrint=true)); fclose(fid);
end
end

function v=check_geometry(m,c)
v.volume_relative=abs(sum(m.vol)-pi*(c.radius^2-c.rw^2)*(m.z_edges(end)-m.z_edges(1)))/sum(m.vol);
ng=size(m.ghost_x,1);
re=[m.r;c.radius*ones(m.nz,1);m.r(m.a(m.boundary(m.nz+1:end)))];
assert(numel(re)==m.n+ng);
radial=m.face_z==0;
q=m.heat(radial).*(log(re(m.a(radial)))-log(re(m.b(radial))));
j=ceil(m.a(radial)/m.nr); h=diff(m.z_edges);
exact=-2*pi*m.lambda(m.a(radial)).*h(j);
% 封闭热边界不参与定值导热的解析比较。
active=m.heat(radial)>0;
v.radial_log_flux_relative=max(abs(q(active)-exact(active))./abs(exact(active)));
if isempty(v.radial_log_flux_relative), v.radial_log_flux_relative=0; end
vertical=m.face_z~=0 & m.b<=m.n;
same=abs(m.lambda(m.a(vertical))-m.lambda(m.b(vertical)))<1e-12;
ids=find(vertical); ids=ids(same);
if isempty(ids), v.vertical_linear_flux_relative=0;
else
    i=mod(m.a(ids)-1,m.nr)+1; A=pi*diff(m.r_edges.^2);
    q=-m.heat(ids).*m.face_z(ids); exact=-m.lambda(m.a(ids)).*A(i);
    v.vertical_linear_flux_relative=max(abs(q-exact)./abs(exact));
end
end

function c=defaults()
c=struct('name','baseline_2d','nr',40,'z_refine',1,'radius',60,'bottom',250,...
    'dimension',2,'duration',864000,'dt_initial',30,'dt_max',3600,'dt_min',.01,...
    'rate_kg_s',10000/3600,'ramp_s',3600,'inject_end_s',Inf,'inject_temperature',5.7,...
    'inject_top',118,'inject_bottom',128,'rw',.10795,'salt',.0346,...
    'sand_top',NaN,'sand_bottom',NaN,'sand_phi',.35,'seal_phi',.40,...
    'sand_k',1.4e-13,'seal_k',1.3*9.869233e-16,'anisotropy',.1,...
    'lambda_sand',2.5,'lambda_seal',2.5,'rock_density',2650,'rock_cp',800,...
    'rho_h',1100,'cp_h',2000,'wc',.298,'reaction_heat',5e5,'k1',3e-6,...
    'nearwell_refine',1,'nearwell_radius',2,'permeability_exponent',4,...
    'sh_cap',.7,'swc',.2,'pc_sand',.2e6,'pc_seal',.5e6,'eq_shift',0,...
    'seafloor_temperature',3,'geothermal_gradient',.061,'reference_pressure',18.4e6,...
    'gravity',9.80665,'reaction',true,'closed',false,'initial_gas',0,'initial_hydrate',0,...
    'uniform_temperature',NaN,'uniform_pressure',NaN,'vertical_flow',true,...
    'tol',1e-9,'max_newton',24,'save',true,'verbose',true,'snapshot_days',[1,5,10]);
c.closure='li2024_moridis2025_v1'; c.results_folder='results';
c.salt_equilibrium='DL'; c.grain_radius=1.645e-4;
c.activation_energy=1.0288e5; c.formation_prefactor=8.4e11;
c.dissociation_prefactor=2.783e7; c.area_beta=5.3;
c.kinetic_multiplier=1; c.area_regularization=1e-10;
c.cp_h_slope=0; c.hydrate_u0=NaN;
c.well_control='distributed_flux'; c.band_only=false; c.field_schedule=[]; c.field_temp=[]; c.layer_profile=[];
end

function p=make_properties(raw)
for name={'gc_rho','gc_mu','gc_h'}
    p.(name{1})=griddedInterpolant({raw.p_axis(:),raw.t_axis(:)},raw.(name{1}),'linear','none');
end
for name={'aq_rho','aq_mu','aq_h'}
    p.(name{1})=griddedInterpolant({raw.p_axis(:),raw.t_axis(:),raw.s_axis(:)},raw.(name{1}),'linear','none');
end
end

function y=lookup(F,p,T,s)
% 插值梯度传递至 MRST ADI；导数在表格小间隔内中心差分。
pv=value(p); tv=value(T)+273.15;
if nargin==4
    sv=value(s); y0=F(pv,tv,sv);
else
    y0=F(pv,tv);
end
if ~isa(p,'ADI'), y=y0; return; end
ep=10; et=1e-4; es=1e-7;
if nargin==4
    dp=(F(pv+ep,tv,sv)-F(pv-ep,tv,sv))/(2*ep);
    dT=(F(pv,tv+et,sv)-F(pv,tv-et,sv))/(2*et);
    ds=(F(pv,tv,sv+es)-F(pv,tv,sv-es))/(2*es);
else
    dp=(F(pv+ep,tv)-F(pv-ep,tv))/(2*ep);
    dT=(F(pv,tv+et)-F(pv,tv-et))/(2*et);
end
y=y0+dp.*(p-pv)+dT.*(T-value(T));
if nargin==4, y=y+ds.*(s-sv); end
end

function m=make_model(c,prop)
r=logspace(log10(c.rw),log10(c.radius),c.nr+1)';
% 仅细分近井环形单元，保留原网格边界与远场离散。
assert(c.nearwell_refine>=1 && c.nearwell_refine==floor(c.nearwell_refine));
assert(c.permeability_exponent>=0);
if c.nearwell_refine>1
    extra=[];
    for ir=1:numel(r)-1
        if r(ir)<c.nearwell_radius
            sub=logspace(log10(r(ir)),log10(r(ir+1)),c.nearwell_refine+1)';
            extra=[extra;sub(2:end-1)]; %#ok<AGROW>
        end
    end
    r=sort([r;extra]);
end
if c.dimension==1
    ze=[c.sand_top;c.sand_bottom];
else
    % 显式保留模型分层边界与注入窗口，向上下逐渐变粗。
    breaks=[0,50,100,c.sand_top,c.inject_top,c.inject_bottom,c.sand_bottom,175,c.bottom];
    breaks=unique(sort(breaks(breaks>=0 & breaks<=c.bottom)));
    ze=breaks(1);
    for j=1:numel(breaks)-1
        dz=breaks(j+1)-breaks(j);
        if breaks(j)>=100 && breaks(j)<175, step=2/c.z_refine; else, step=12.5/c.z_refine; end
        q=linspace(breaks(j),breaks(j+1),max(1,ceil(dz/step))+1)';
        ze=[ze;q(2:end)]; %#ok<AGROW>
    end
end
nr=numel(r)-1; nz=numel(ze)-1;
% MRST 网格用于拓扑核对；物理算子使用精确环形体积和面半传导系数。
G=computeGeometry(tensorGrid(r,ze));
rc=sqrt(r(1:end-1).*r(2:end)); zc=(ze(1:end-1)+ze(2:end))/2;
[rr,zz]=ndgrid(rc,zc); [area,dz]=ndgrid(pi*diff(r.^2),diff(ze));
m=struct('nr',nr,'nz',nz,'n',nr*nz,'r_edges',r,'z_edges',ze,'r',rr(:),'z',zz(:),'vol',area(:).*dz(:));
assert(G.cells.num==m.n);
isSand=m.z>=c.sand_top & m.z<=c.sand_bottom;
m.phi=c.seal_phi*ones(m.n,1); m.phi(isSand)=c.sand_phi;
m.k=c.seal_k*ones(m.n,1); m.k(isSand)=c.sand_k;
% 各区间依次给出顶深、底深（m）和渗透率（mD）。
bands=c.layer_profile;
assert(~isempty(bands) && size(bands,2)==3 && all(isfinite(bands),'all') && all(bands(:,2)>bands(:,1)) && all(bands(:,3)>0),'Supply top, base and permeability for each interval.');
for row=bands'
    m.k(m.z>=row(1) & m.z<row(2))=row(3)*9.869233e-16;
end
if c.band_only, m.k(m.z<c.inject_top | m.z>c.inject_bottom)=c.seal_k; end
m.lambda=c.lambda_seal*ones(m.n,1); m.lambda(isSand)=c.lambda_sand;
m.pc=c.pc_seal*ones(m.n,1); m.pc(isSand)=c.pc_sand;
m.pv=m.vol.*m.phi;
% 相同连续初始剖面映射到各个网格；水压通过面平均密度迭代至静水平衡。
m.t0=c.seafloor_temperature+c.geothermal_gradient*m.z;
if isfinite(c.uniform_temperature), m.t0(:)=c.uniform_temperature; end
pz=c.reference_pressure+1025*c.gravity*(zc-mean([c.sand_top,c.sand_bottom]));
tz=m.t0(1:nr:end);
for it=1:15
    rho=lookup(prop.aq_rho,pz,tz,c.salt*ones(nz,1));
    dp=.5*(rho(1:end-1)+rho(2:end))*c.gravity.*diff(zc);
    temp=[0;cumsum(dp)];
    if nz>1, off=interp1(zc,temp,mean([c.sand_top,c.sand_bottom]),'linear','extrap'); else, off=0; end
    pz=c.reference_pressure+temp-off;
end
m.p0=repelem(pz(:),nr,1);
if isfinite(c.uniform_pressure), m.p0(:)=c.uniform_pressure; end
injmask=m.r==rc(1) & m.z>=c.inject_top & m.z<=c.inject_bottom;
if c.dimension==1, injmask(1)=true; end
m.inj=find(injmask); assert(~isempty(m.inj),'注入窗口未覆盖任何单元。');
m.inj_weight=dz(m.inj)/sum(dz(m.inj));
m.well_index=2*pi*dz(m.inj).*m.k(m.inj)./log(m.r(m.inj)/c.rw);
% 面方向统一由 a 单元指向 b 单元/外部定值状态。
a=[]; b=[]; rh1=[]; rh2=[]; ht=[]; hyd=[]; fz=[]; ghost=[];
for j=1:nz
    for i=1:nr-1
        ia=i+(j-1)*nr; ib=ia+1; rf=r(i+1); h=ze(j+1)-ze(j);
        addface(ia,ib,log(rf/rc(i))/(2*pi*h),log(rc(i+1)/rf)/(2*pi*h),1,true);
    end
end
for j=1:nz-1
    for i=1:nr
        ia=i+(j-1)*nr; ib=ia+nr; A=pi*(r(i+1)^2-r(i)^2);
        addface(ia,ib,(ze(j+1)-zc(j))/A,(zc(j+1)-ze(j+1))/A,c.anisotropy,c.vertical_flow);
    end
end
ni=numel(a);
for j=1:nz
    ia=j*nr; h=ze(j+1)-ze(j);
    addboundary(ia,log(r(end)/rc(end))/(2*pi*h),m.z(ia),1,~c.closed,m.p0(ia),m.t0(ia));
end
if c.dimension==2
    for i=1:nr
        A=pi*(r(i+1)^2-r(i)^2);
        ia=i;
        tb=c.seafloor_temperature;
        pb=m.p0(ia)-1025*c.gravity*m.z(ia);
        for it=1:10
            rho0=lookup(prop.aq_rho,m.p0(ia),m.t0(ia),c.salt);
            rhob=lookup(prop.aq_rho,pb,tb,c.salt);
            pb=m.p0(ia)-.5*(rho0+rhob)*c.gravity*m.z(ia);
        end
        if isfinite(c.uniform_temperature), tb=c.uniform_temperature; end
        if isfinite(c.uniform_pressure), pb=c.uniform_pressure; end
        addboundary(ia,(zc(1)-ze(1))/A,ze(1),c.anisotropy,~c.closed && c.vertical_flow,pb,tb);
        ia=i+(nz-1)*nr;
        tb=c.seafloor_temperature+c.geothermal_gradient*c.bottom;
        if isfinite(c.uniform_temperature), tb=c.uniform_temperature; end
        pb=m.p0(ia)+1025*c.gravity*(c.bottom-m.z(ia));
        addboundary(ia,(ze(end)-zc(end))/A,ze(end),c.anisotropy,false,pb,tb);
    end
end
nf=numel(a); internal=b<=m.n;
m.D=sparse([a;b(internal)],[(1:nf)';find(internal)],[ones(nf,1);-ones(sum(internal),1)],m.n,nf);
m.a=a; m.b=b; m.rh1=rh1; m.rh2=rh2; m.heat=ht; m.hydraulic=hyd;
m.face_z=fz; m.ghost_x=ghost; m.boundary=(ni+1:nf)';
    function addface(ia,ib,g1,g2,an,flow)
        a(end+1,1)=ia; b(end+1,1)=ib;
        rh1(end+1,1)=g1/(m.k(ia)*an); rh2(end+1,1)=g2/(m.k(ib)*an);
        ht(end+1,1)=1/(g1/m.lambda(ia)+g2/m.lambda(ib));
        hyd(end+1,1)=double(flow); fz(end+1,1)=m.z(ib)-m.z(ia);
    end
    function addboundary(ia,g1,zb,an,flow,pb,tb)
        ghost(end+1,:)=[pb/1e6,0,0,tb,c.salt];
        a(end+1,1)=ia; b(end+1,1)=m.n+size(ghost,1);
        rh1(end+1,1)=g1/(m.k(ia)*an); rh2(end+1,1)=0;
        ht(end+1,1)=double(~c.closed)*m.lambda(ia)/g1;
        hyd(end+1,1)=double(flow); fz(end+1,1)=zb-m.z(ia);
    end
end

function q=phase_properties(p,sg,sh,T,s,pcscale,prop,c)
sf=sg./(1-sh);
se=max(0,min(1,(1-sf-c.swc)/(1-c.swc)));
q.pc=pcscale.*(1-se).^2;
q.rhog=lookup(prop.gc_rho,p+q.pc,T);
q.mug=lookup(prop.gc_mu,p+q.pc,T);
q.hg=lookup(prop.gc_h,p+q.pc,T);
q.rhow=lookup(prop.aq_rho,p,T,s);
q.muw=lookup(prop.aq_mu,p,T,s);
q.hw=lookup(prop.aq_h,p,T,s);
q.krw=se.^2; q.krg=(1-se).^2;
end

function inv=inventory(x,m,prop,c)
[p,sg,sh,T,s]=split(x,false);
q=phase_properties(p,sg,sh,T,s,m.pc,prop,c);
mg=m.pv.*sg.*q.rhog; mw=m.pv.*(1-sg-sh).*q.rhow; mh=m.pv.*sh*c.rho_h;
zpot=-c.gravity*m.z;
uh=hydrate_internal_energy(T,c);
energy=(1-m.phi).*m.vol*c.rock_density*c.rock_cp.*T+...
    mg.*(q.hg-(p+q.pc)./q.rhog+zpot)+mw.*(q.hw-p./q.rhow+zpot)+mh.*(uh+zpot);
inv=[mg+c.wc*mh,(1-s).*mw+(1-c.wc)*mh,s.*mw,c.wc*mh,energy];
end

function [eq,flux,source,wellinfo]=equations(x,old,dt,rate,m,prop,c,usead)
[p,sg,sh,T,s]=split(x,usead);
ng=size(m.ghost_x,1);
pe=[p;m.ghost_x(:,1)*1e6]; ge=[sg;zeros(ng,1)]; he=[sh;zeros(ng,1)];
te=[T;m.ghost_x(:,4)]; se=[s;m.ghost_x(:,5)];
pcs=[m.pc;zeros(ng,1)];
q=phase_properties(pe,ge,he,te,se,pcs,prop,c);
a=m.a; b=m.b; zpot=-c.gravity*m.z;
dam=(1-he).^c.permeability_exponent;
trans=m.hydraulic./(m.rh1./dam(a)+m.rh2./dam(b));
dw=pe(a)-pe(b)+.5*(q.rhow(a)+q.rhow(b))*c.gravity.*m.face_z;
dg=pe(a)+q.pc(a)-pe(b)-q.pc(b)+.5*(q.rhog(a)+q.rhog(b))*c.gravity.*m.face_z;
uw=a; ug=a; uw(value(dw)<0)=b(value(dw)<0); ug(value(dg)<0)=b(value(dg)<0);
fw=trans.*dw.*q.krw(uw)./q.muw(uw).*q.rhow(uw);
fg=trans.*dg.*q.krg(ug)./q.mug(ug).*q.rhog(ug);
% 重力势能与焓使用同一上游状态和单元高程。
ze=[m.z; m.z(a(m.boundary))+m.face_z(m.boundary)];
ef=fw.*(q.hw(uw)-c.gravity*ze(uw))+fg.*(q.hg(ug)-c.gravity*ze(ug))+m.heat.*(te(a)-te(b));
fwater=(1-se(uw)).*fw; fsalt=se(uw).*fw;
mg=m.pv.*sg.*q.rhog(1:m.n); mw=m.pv.*(1-sg-sh).*q.rhow(1:m.n); mh=m.pv.*sh*c.rho_h;
uh=hydrate_internal_energy(T,c);
E=(1-m.phi).*m.vol*c.rock_density*c.rock_cp.*T+mg.*(q.hg(1:m.n)-(p+q.pc(1:m.n))./q.rhog(1:m.n)+zpot)+...
    mw.*(q.hw(1:m.n)-p./q.rhow(1:m.n)+zpot)+mh.*(uh+zpot);
reaction=hydrate_modern_closure(p+q.pc(1:m.n),T,s,sg,sh,m.phi,m.vol,c);
injT=c.inject_temperature+0*p;
wellinfo=[NaN,rate];
if strcmp(c.well_control,'common_bhp') && rate>0
    % 共用参考井压；井壁至首单元的圆柱半传导系数，纯CO2注入上游流度。
    ids=m.inj; pg=p(ids)+q.pc(ids);
    rhoi=lookup(prop.gc_rho,pg,injT(ids));
    mui=lookup(prop.gc_mu,pg,injT(ids));
    conduct=m.well_index.*(1-sh(ids)).^c.permeability_exponent.*rhoi./mui;
    hydro=(sum(rhoi)/numel(ids))*c.gravity.*(m.z(ids)-mean([c.inject_top,c.inject_bottom]));
    threshold=pg-hydro;
    active=true(numel(ids),1);
    for wi=1:numel(ids)+1
        pw=(rate+sum(conduct(active).*threshold(active)))./sum(conduct(active));
        next=value(pw)>value(threshold);
        if isequal(next,active), break; end
        active=next;
    end
    assert(any(active),'注入井无可用连接。');
    local=conduct.*max(pw-threshold,0);
    W=sparse(ids,1:numel(ids),1,m.n,numel(ids));
    inj=W*local;
    wellinfo=[value(pw)/1e6,sum(value(local))];
    % 注入焓以井壁CO2压力计算，与所需供给压力一致。
    hp=p+W*(pw+hydro-p(ids));
    hinj=lookup(prop.gc_h,hp,injT);
else
    inj=zeros(m.n,1); inj(m.inj)=rate*m.inj_weight;
    hinj=lookup(prop.gc_h,p,injT);
end
qin=inj.*(hinj+zpot);
eq={(mg+c.wc*mh-old(:,1)+dt*(m.D*fg-inj))./(m.pv*1000),...
    ((1-s).*mw+(1-c.wc)*mh-old(:,2)+dt*m.D*fwater)./(m.pv*1000),...
    (s.*mw-old(:,3)+dt*m.D*fsalt)./(m.pv*30),...
    (c.wc*mh-old(:,4)-dt*reaction)./(m.pv*1000),...
    (E-old(:,5)+dt*(m.D*ef-qin))./(m.vol*3e7)};
if nargout>1
    flux=[value(fg),value(fwater),value(fsalt),value(ef)];
    source=[sum(inj),0,0,sum(value(qin))];
end
end

function uh=hydrate_internal_energy(T,c)
% 不可压缩水合物的内能与压力无关；相焓为 u_h+p/rho_h。
if isfinite(c.hydrate_u0)
    uh=c.cp_h*T+.5*c.cp_h_slope*T.^2+c.hydrate_u0;
else
    uh=c.cp_h*T-c.wc*c.reaction_heat-c.reference_pressure/c.rho_h;
end
end

function [p,sg,sh,T,s]=split(x,ad)
if ad
    [p,sg,sh,T,s]=initVariablesADI(x(:,1),x(:,2),x(:,3),x(:,4),x(:,5));
else
    p=x(:,1); sg=x(:,2); sh=x(:,3); T=x(:,4); s=x(:,5);
end
p=p*1e6;
end

function [x,ok,it,res]=newton(x,old,dt,rate,m,prop,c)
ok=false; res=Inf;
for it=1:c.max_newton
    eq=equations(x,old,dt,rate,m,prop,c,true);
    R=vertcat(eq{:});
    assert(isreal(R.val) && all(isfinite(R.val)),'方程残差包含复数或非有限数值。');
    res=norm(R.val,inf);
    if res<c.tol, ok=true; return; end
    blocks=R.jac; J=horzcat(blocks{:});
    assert(size(J,2)==5*m.n,'Jacobian size %d x %d; expected %d unknowns',size(J,1),size(J,2),5*m.n);
    dx=reshape(-(J\R.val),m.n,5);
    if any(~isfinite(dx(:))), return; end
    alpha=1;
    for ls=1:24
        y=x+alpha*dx;
        feasible=all(y(:,1)>8.001 & y(:,1)<34 & y(:,2)>-1e-10 & y(:,3)>-1e-10 & ...
            y(:,2)+y(:,3)<.999999 & y(:,4)>.001 & y(:,4)<34.9 & y(:,5)>.000001 & y(:,5)<.19999);
        if feasible
            ey=equations(y,old,dt,rate,m,prop,c,false);
            rn=norm(vertcat(ey{:}),inf);
            if isfinite(rn) && (rn<(1-1e-4*alpha)*res || rn<c.tol), x=y; break; end
        end
        alpha=alpha/2;
    end
    if ls==24, return; end
end
end

function [rate,amount]=injection_rate(t,dt,c)
t0=min(t,c.inject_end_s); t1=min(t+dt,c.inject_end_s);
if c.ramp_s>0
    integral=@(x) min(x,c.ramp_s).^2/(2*c.ramp_s)+max(x-c.ramp_s,0);
    amount=c.rate_kg_s*(integral(t1)-integral(t0));
else
    amount=c.rate_kg_s*(t1-t0);
end
rate=amount/dt;
end
