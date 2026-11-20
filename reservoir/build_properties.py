"""生成低温高压 CO2 / NaCl 盐水物性表，记录软件版本与适用域。"""
from pathlib import Path
import json
ROOT = Path(__file__).resolve().parent
import numpy as np
from scipy.io import savemat
import CoolProp
from CoolProp.CoolProp import PropsSI

p = np.linspace(8e6, 35e6, 55)
t = np.linspace(273.15, 308.15, 71)
s = np.linspace(0, .20, 21)
# CO2 的同一参考焓用于库存与边界，避免任意参考态造成假热源。
P, T = np.meshgrid(p, t, indexing='ij')
data = {'p_axis':p, 't_axis':t, 's_axis':s}
for name, key in [('gc_rho','Dmass'), ('gc_mu','V'), ('gc_h','Hmass')]:
    data[name] = PropsSI(key, 'P', P.ravel(), 'T', T.ravel(), 'CO2').reshape(P.shape)
data['gc_h'] -= PropsSI('Hmass','P',18.4e6,'T',273.15,'CO2')
# Melinder NaCl 溶液: 密度、黏度和热容随温度和盐度变化；盐度取 NaCl 质量分数。
# 在 MNA 基准上加入恒定等温压缩系数，焓由同一 Gibbs 压力积分构建。
# 每个盐度的焓参考态统一为 0 °C、18.4 MPa，明确忽略混盐热。
shape=(len(p),len(t),len(s))
for name in ['aq_rho','aq_mu','aq_h']:
    data[name]=np.empty(shape)
for k, sal in enumerate(s):
    fluid=f'INCOMP::MNA[{sal:.8f}]'
    href=PropsSI('Hmass','P',18.4e6,'T',273.15,fluid)
    pref=18.4e6
    comp=4.5e-10
    rho0=PropsSI('Dmass','P',pref,'T',T.ravel(),fluid).reshape(P.shape)
    v0=1/rho0
    vplus=1/PropsSI('Dmass','P',pref,'T',(T+.001).ravel(),fluid).reshape(P.shape)
    vminus=1/PropsSI('Dmass','P',pref,'T',(T-.001).ravel(),fluid).reshape(P.shape)
    dvdt=(vplus-vminus)/.002
    h0=PropsSI('Hmass','P',pref,'T',T.ravel(),fluid).reshape(P.shape)-href
    data['aq_rho'][:,:,k]=rho0*np.exp(comp*(P-pref))
    data['aq_mu'][:,:,k]=PropsSI('V','P',P.ravel(),'T',T.ravel(),fluid).reshape(P.shape)
    data['aq_h'][:,:,k]=h0+(v0-T*dvdt)*(-np.expm1(-comp*(P-pref)))/comp
data['property_note']='CoolProp HEOS CO2; Melinder MNA brine with thermodynamically consistent constant-compressibility pressure extension. No CO2 dissolution or heat of salt mixing.'
savemat(ROOT/'properties.mat',data)
meta={'CoolProp':CoolProp.__version__,'pressure_MPa':[8,35], 'temperature_C':[0,35],
      'NaCl_mass_fraction':[0,.20], 'CO2':'HEOS; Span-Wagner EOS and bundled viscosity correlation',
      'brine':'CoolProp MNA (Melinder) NaCl-water surrogate; c=4.5e-10 /Pa; density and enthalpy from one Gibbs pressure integral',
      'enthalpy_reference':'0 C, 18.4 MPa separately at each salinity; no heat of mixing',
      'checks':{}}
for fluid in ['CO2','INCOMP::MNA[0.0346]']:
    meta['checks'][fluid]={key:PropsSI(key,'P',18.4e6,'T',282.65,fluid) for key in ['Dmass','V','Cpmass']}
(ROOT/'properties_provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
print(json.dumps(meta,indent=2))
