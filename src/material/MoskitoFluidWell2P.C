/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*                                                                        */
/*  This program is free software: you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation, either version 3 of the License, or     */
/*  (at your option) any later version.                                   */
/*                                                                        */
/*  This program is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          */
/*  GNU General Public License for more details.                          */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>  */
/**************************************************************************/

#include "MoskitoFluidWell2P.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell2P);

template <>
InputParameters
validParams<MoskitoFluidWell2P>()
{
  InputParameters params = validParams<MoskitoFluidWellGeneral>();
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for 2 phase EOS");
  params.addRequiredParam<UserObjectName>("viscosity_uo",
        "The name of the userobject for 2 phase viscosity Eq");
  params.addRequiredParam<UserObjectName>("drift_flux_uo",
        "The name of the userobject for drift flux model");
  return params;
}

MoskitoFluidWell2P::MoskitoFluidWell2P(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    eos_uo(getUserObject<MoskitoEOS2P>("eos_uo")),
    viscosity_uo(getUserObject<MoskitoViscosity2P>("viscosity_uo")),
    dfm_uo(getUserObject<MoskitoDriftFlux>("drift_flux_uo")),
    _cp_m(declareProperty<Real>("specific_heat")),
    _rho_g(declareProperty<Real>("gas_density")),
    _rho_l(declareProperty<Real>("liquid_density")),
    _rho_m(declareProperty<Real>("density")),
    _rho_pam(declareProperty<Real>("profile_mixture_density")),
    _drho_m_dp(declareProperty<Real>("drho_dp")),
    _drho_m_dp_2(declareProperty<Real>("drho_dp_2")),
    _drho_m_dT(declareProperty<Real>("drho_dT")),
    _drho_m_dh(declareProperty<Real>("drho_dh")),
    _drho_m_dh_2(declareProperty<Real>("drho_dh_2")),
    _vmfrac(declareProperty<Real>("mass_fraction")),
    _u_g(declareProperty<Real>("gas_velocity")),
    _u_l(declareProperty<Real>("liquid_velocity")),
    _vfrac(declareProperty<Real>("void_fraction")),
    _phase(declareProperty<Real>("current_phase")),
    _u_d(declareProperty<Real>("drift_velocity")),
    _c0(declareProperty<Real>("flow_type_c0")),
    _flow_pat(declareProperty<Real>("flow_pattern")),
    _v_sg(declareProperty<Real>("superficial_gas_velocity")),
    _v_sl(declareProperty<Real>("superficial_liquid_velocity")),
    _dgamma_dz(declareProperty<Real>("dgamma_dz")),
    _dgamma_dz_uj_gphi(declareProperty<Real>("dgamma_dz_uj_gphi")),
    _dgamma_dz_uj_phi(declareProperty<Real>("dgamma_dz_uj_phi")),
    _dgamma_dz_pj_gphi(declareProperty<Real>("dgamma_dz_pj_gphi")),
    _dgamma_dz_pj_phi(declareProperty<Real>("dgamma_dz_pj_phi")),
    _dgamma_dz_hj_gphi(declareProperty<Real>("dgamma_dz_hj_gphi")),
    _dgamma_dz_hj_phi(declareProperty<Real>("dgamma_dz_hj_phi")),
    _grad_flow(coupledGradient("flowrate")),
    _grad_h(coupledGradient("enthalpy")),
    _grad_p(coupledGradient("pressure"))
{
}

void
MoskitoFluidWell2P::computeQpProperties()
{
  // calculate required properties based on the given EOS
  eos_uo.VMFrac_T_from_p_h(fabs(_P[_qp]), _h[_qp], _vmfrac[_qp], _T[_qp], _phase[_qp]);
  eos_uo.rho_m_by_p(fabs(_P[_qp]), _h[_qp], _rho_m[_qp], _drho_m_dp[_qp], _drho_m_dp_2[_qp]);
  eos_uo.rho_m_by_h(fabs(_P[_qp]), _h[_qp], _rho_m[_qp], _drho_m_dh[_qp], _drho_m_dh_2[_qp]);
  _rho_l[_qp] = eos_uo.rho_l_from_p_T(fabs(_P[_qp]), _T[_qp], _phase[_qp]);
  _rho_g[_qp] = eos_uo.rho_g_from_p_T(fabs(_P[_qp]), _T[_qp], _phase[_qp]);
  _vfrac[_qp]  = (_rho_m[_qp] - _rho_l[_qp]) / (_rho_g[_qp] - _rho_l[_qp]);
  _rho_pam[_qp] = _rho_g[_qp] * _c0[_qp]  * _vfrac[_qp] + (1.0 - _vfrac[_qp] * _c0[_qp]) * _rho_l[_qp];
  _cp_m[_qp]  = eos_uo.cp_m_from_p_T(fabs(_P[_qp]), _T[_qp], _vmfrac[_qp], _phase[_qp]);

  _dia[_qp] = _d;
  _area[_qp] = PI * _d * _d / 4.0;
  _u[_qp] = _flow[_qp] / _area[_qp];
  _Re[_qp] = _rho_m[_qp] * _dia[_qp] * fabs(_u[_qp]) / viscosity_uo.mixture_mu(_P[_qp], _T[_qp], _vmfrac[_qp]);

  MoskitoFluidWellGeneral::computeQpProperties();

  // drift-flux calculator section
  if (_phase[_qp] == 2.0)
  {
    MoskitoDFGVar DFinp(_u[_qp], _rho_g[_qp], _rho_l[_qp], _vmfrac[_qp],
      _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
    dfm_uo.DFMCalculator(DFinp);
      DFinp.DFMOutput(_flow_pat[_qp], _v_sg[_qp], _v_sl[_qp], _c0[_qp], _u_d[_qp]);
  }
  else
  {
    _flow_pat[_qp] = 0;
    _c0[_qp] = 1.0;
    _u_d[_qp] = 0.0;
  }

  // based on mass weighted flow rate
  // momentum eq is valid only by mass mixing flow rate
  if (_phase[_qp] == 2.0)
  {
    _u_g[_qp]  = _c0[_qp] * _rho_m[_qp] * _u[_qp] + _rho_l[_qp] * _u_d[_qp];
    _u_g[_qp] /= _rho_pam[_qp];
    _u_l[_qp]  = (1.0 - _vfrac[_qp] * _c0[_qp]) * _rho_m[_qp]  * _u[_qp] - _rho_g[_qp] * _vfrac[_qp] * _u_d[_qp];
    _u_l[_qp] /= (1.0 - _vfrac[_qp]) * _rho_pam[_qp];
  }
  else
  {
    if (_vmfrac[_qp] == 0.0)
    {
      _u_g[_qp] = 0.0;
      _u_l[_qp] = _u[_qp];
    }
    else
    {
      _u_l[_qp] = 0.0;
      _u_g[_qp] = _u[_qp];
    }
  }

  if (_phase[_qp] == 2.0)
    DriftFluxMomentumEq();
  else
  {
    _dgamma_dz[_qp]         = 0.0;
    _dgamma_dz_uj_gphi[_qp] = 0.0;
    _dgamma_dz_uj_phi[_qp]  = 0.0;
    _dgamma_dz_pj_gphi[_qp] = 0.0;
    _dgamma_dz_pj_phi[_qp]  = 0.0;
    _dgamma_dz_hj_gphi[_qp] = 0.0;
    _dgamma_dz_hj_phi[_qp]  = 0.0;
  }

  // _lambda[_qp]  = (1.0 - (_d * _d) / std::pow(_d + _thickness , 2.0)) * _lambda0;
  // _lambda[_qp] += (_d * _d) / std::pow(_d + _thickness , 2.0) * eos_uo._lambda;
}

void
MoskitoFluidWell2P::DriftFluxMomentumEq()
{
  Real _drho_g_dp, _drho_g_dT, _drho_l_dp, _drho_l_dT, temp;

  eos_uo.rho_l_from_p_T(fabs(_P[_qp]), _T[_qp], temp, _drho_l_dp, _drho_l_dT, _phase[_qp]);
  eos_uo.rho_g_from_p_T(fabs(_P[_qp]), _T[_qp], temp, _drho_g_dp, _drho_g_dT, _phase[_qp]);

  /*
  All required coupled coeffient are derived below for momentum conservation
  In 2 phase momentum equation, the following rules should applied:
  1. gphi should be multiplied by (_grad_phi[j][_qp] * _well_dir[_qp])
  2. phi should be multiplied by (_phi[j][_qp])
  */

  // defined and used local variables
  Real _drho_g_dz, _drho_l_dz, _drho_m_dz, _drho_pam_dz, _dc_dz;
  Real _dc_dp, _drho_pam_dp, _dc_dh, _drho_pam_dh;
  Real part1, part2, temp1, temp2;

  // parameterizing gamma equation for simpilification
  Real a, b, c, d;
  a = _vfrac[_qp] / (1.0 - _vfrac[_qp]);
  b = _c0[_qp] - 1.0;
  c = _rho_g[_qp] * _rho_l[_qp] * _rho_m[_qp] / (_rho_pam[_qp] * _rho_pam[_qp]);
  d = b * _u[_qp] + _u_d[_qp];

  // parameterizing dgamma_dz equation for simpilification
  _drho_g_dz  = (_drho_g_dp * _grad_p[_qp] + _drho_g_dT * _grad_h[_qp] / _cp_m[_qp]) * _well_dir[_qp];
  _drho_l_dz  = (_drho_l_dp * _grad_p[_qp] + _drho_l_dT * _grad_h[_qp] / _cp_m[_qp]) * _well_dir[_qp];
  _drho_m_dz  =  _vfrac[_qp] * _drho_g_dz + (1.0 - _vfrac[_qp]) * _drho_l_dz;
  _drho_pam_dz = _c0[_qp] * _vfrac[_qp] * _drho_g_dz + (1.0 - _vfrac[_qp] * _c0[_qp]) * _drho_l_dz;
  part1  = _drho_g_dz  * _rho_l[_qp]     * _rho_m[_qp];
  part1 += _rho_g[_qp] * _drho_l_dz      * _rho_m[_qp];
  part1 += _rho_g[_qp] * _rho_l[_qp]     * _drho_m_dz;
  part2  = -2.0 * _drho_pam_dz * _rho_g[_qp] * _rho_l[_qp] * _rho_m[_qp];
  _dc_dz = (part1 * _rho_pam[_qp] + part2) / std::pow(_rho_pam[_qp] , 3.0);

  // residual for dgamma_dz in the momentum conservation
  _dgamma_dz[_qp]  = a * _dc_dz * std::pow(d, 2.0);
  _dgamma_dz[_qp] += 2.0 * a * b * c * d * _grad_flow[_qp] / _area[_qp] * _well_dir[_qp];

  // diagonal jacobian of the residual wrt uj for dgamma_dz in the momentum conservation
  _dgamma_dz_uj_gphi[_qp]  = 2.0 * a * b * c * d / _area[_qp];
  _dgamma_dz_uj_phi[_qp]  = _dc_dz * d + c * b * _grad_flow[_qp] / _area[_qp] * _well_dir[_qp];
  _dgamma_dz_uj_phi[_qp] *= 2.0 * a * b;

  _drho_pam_dp = _c0[_qp] * _drho_g_dp * _vfrac[_qp] + (1.0 - _vfrac[_qp] * _c0[_qp]) * _drho_l_dp;
  temp1  = _drho_g_dp  * _rho_l[_qp]     * _rho_m[_qp];
  temp1 += _rho_g[_qp] * _drho_l_dp      * _rho_m[_qp];
  temp1 += _rho_g[_qp] * _rho_l[_qp]     * _drho_m_dp[_qp];
  temp2  = -2.0 * _drho_pam_dp * _rho_g[_qp] * _rho_l[_qp] * _rho_m[_qp];
  _dc_dp  = (temp1 * _rho_pam[_qp] + temp2) / std::pow(_rho_pam[_qp] , 3.0);

  // off diagonal jacobian of the residual wrt pj for dgamma_dz in the momentum conservation
  _dgamma_dz_pj_gphi[_qp] = a * d * d * _dc_dp;
  _dgamma_dz_pj_phi[_qp]  = _drho_g_dz  * _drho_l_dp  * _rho_m[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_g_dz  * _rho_l[_qp] * _drho_m_dp[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_g_dp  * _drho_l_dz  * _rho_m[_qp];
  _dgamma_dz_pj_phi[_qp] += _rho_g[_qp] * _drho_l_dz  * _drho_m_dp[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_g_dp  * _rho_l[_qp] * _drho_m_dz;
  _dgamma_dz_pj_phi[_qp] += _rho_g[_qp] * _drho_l_dp  * _drho_m_dz;
  _dgamma_dz_pj_phi[_qp] *= _rho_pam[_qp];
  _dgamma_dz_pj_phi[_qp] += _drho_pam_dp * part1;
  _dgamma_dz_pj_phi[_qp] -= 2.0 * _drho_m_dz * temp1;
  _dgamma_dz_pj_phi[_qp] *= _rho_pam[_qp];
  _dgamma_dz_pj_phi[_qp] -= 3.0 * _drho_pam_dp * (part1 * _rho_pam[_qp] + part2);
  _dgamma_dz_pj_phi[_qp] /= std::pow(_rho_pam[_qp] , 4.0);
  _dgamma_dz_pj_phi[_qp] *= a * d * d;
  _dgamma_dz_pj_phi[_qp] += 2.0 * a * b * _dc_dp * d * _grad_flow[_qp] / _area[_qp] * _well_dir[_qp];

  _drho_pam_dh  = _c0[_qp] * _drho_g_dT * _vfrac[_qp] + (1.0 - _vfrac[_qp] * _c0[_qp]) * _drho_l_dT;
  _drho_pam_dh /= _cp_m[_qp];
  temp1  = _drho_g_dT  * _rho_l[_qp]     * _rho_m[_qp];
  temp1 += _rho_g[_qp] * _drho_l_dT      * _rho_m[_qp];
  temp1 += _rho_g[_qp] * _rho_l[_qp]     * _drho_m_dT[_qp];
  temp1 /= _cp_m[_qp];
  temp2  = -2.0 * _drho_pam_dh * _rho_g[_qp] * _rho_l[_qp] * _rho_m[_qp];
  _dc_dh  = (temp1 * _rho_pam[_qp] + temp2) / std::pow(_rho_pam[_qp] , 3.0);

  // off diagonal jacobian of the residual wrt hj for dgamma_dz in the momentum conservation
  _dgamma_dz_hj_gphi[_qp] = a * d * d * _dc_dh;
  _dgamma_dz_hj_phi[_qp]  = _drho_g_dz  * _drho_l_dT  * _rho_m[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_g_dz  * _rho_l[_qp] * _drho_m_dT[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_g_dT  * _drho_l_dz  * _rho_m[_qp];
  _dgamma_dz_hj_phi[_qp] += _rho_g[_qp] * _drho_l_dz  * _drho_m_dT[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_g_dT  * _rho_l[_qp] * _drho_m_dz;
  _dgamma_dz_hj_phi[_qp] += _rho_g[_qp] * _drho_l_dT  * _drho_m_dz;
  _dgamma_dz_hj_phi[_qp] *= _rho_pam[_qp] / _cp_m[_qp];
  _dgamma_dz_hj_phi[_qp] += _drho_pam_dh * part1;
  _dgamma_dz_hj_phi[_qp] -= 2.0 * _drho_m_dz * temp1;
  _dgamma_dz_hj_phi[_qp] *= _rho_pam[_qp];
  _dgamma_dz_hj_phi[_qp] -= 3.0 * _drho_pam_dh * (part1 * _rho_pam[_qp] + part2);
  _dgamma_dz_hj_phi[_qp] /= std::pow(_rho_pam[_qp] , 4.0);
  _dgamma_dz_hj_phi[_qp] *= a * d * d;
  _dgamma_dz_hj_phi[_qp] += 2.0 * a * b * _dc_dh * d * _grad_flow[_qp] / _area[_qp] * _well_dir[_qp];
}
