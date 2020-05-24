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

#include "MoskitoFluidWell_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell_2p1c);

template <>
InputParameters
validParams<MoskitoFluidWell_2p1c>()
{
  InputParameters params = validParams<MoskitoFluidWellGeneral>();
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable (J/kg)");
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for 2 phase EOS");
  params.addRequiredParam<UserObjectName>("viscosity_uo",
        "The name of the userobject for 2 phase viscosity Eq");
  params.addRequiredParam<UserObjectName>("drift_flux_uo",
        "The name of the userobject for drift flux model");
  return params;
}

MoskitoFluidWell_2p1c::MoskitoFluidWell_2p1c(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    eos_uo(getUserObject<MoskitoEOS2P>("eos_uo")),
    viscosity_uo(getUserObject<MoskitoViscosity2P>("viscosity_uo")),
    dfm_uo(getUserObject<MoskitoDriftFlux>("drift_flux_uo")),
    _T(declareProperty<Real>("temperature")),
    _cp_m(declareProperty<Real>("specific_heat")),
    _rho_g(declareProperty<Real>("gas_density")),
    _rho_l(declareProperty<Real>("liquid_density")),
    _rho_m(declareProperty<Real>("density")),
    _rho_pam(declareProperty<Real>("profile_mixture_density")),
    _drho_m_dp(declareProperty<Real>("drho_dp")),
    _drho_m_dp_2(declareProperty<Real>("drho_dp_2")),
    _drho_m_dh(declareProperty<Real>("drho_dh")),
    _drho_m_dh_2(declareProperty<Real>("drho_dh_2")),
    _drho_m_dph(declareProperty<Real>("drho_dph")),
    _vmfrac(declareProperty<Real>("mass_fraction")),
    _u_g(declareProperty<Real>("gas_velocity")),
    _u_l(declareProperty<Real>("liquid_velocity")),
    _vfrac(declareProperty<Real>("void_fraction")),
    _phase(declareProperty<Real>("current_phase")),
    _u_d(declareProperty<Real>("drift_velocity")),
    _c0(declareProperty<Real>("flow_type_c0")),
    _flow_pat(declareProperty<Real>("flow_pattern")),
    _dgamma_dh(declareProperty<Real>("dgamma_dh")),
    _dgamma_dp(declareProperty<Real>("dgamma_dp")),
    _dgamma_dq(declareProperty<Real>("dgamma_dq")),
    _dgamma2_dhq(declareProperty<Real>("dgamma2_dhq")),
    _dgamma2_dpq(declareProperty<Real>("dgamma2_dpq")),
    _dgamma2_dq2(declareProperty<Real>("dgamma2_dq2")),
    _dkappa_dh(declareProperty<Real>("dkappa_dh")),
    _dkappa_dp(declareProperty<Real>("dkappa_dp")),
    _dkappa_dq(declareProperty<Real>("dkappa_dq")),
    _domega_dh(declareProperty<Real>("domega_dh")),
    _domega_dp(declareProperty<Real>("domega_dp")),
    _domega_dq(declareProperty<Real>("domega_dq")),
    _h(coupledValue("enthalpy")),
    _grad_flow(coupledGradient("flowrate")),
    _grad_h(coupledGradient("enthalpy")),
    _grad_p(coupledGradient("pressure"))
{
}

void
MoskitoFluidWell_2p1c::computeQpProperties()
{
  MoskitoFluidWellGeneral::computeQpProperties();

  // To calculate required properties based on the given EOS
  eos_uo.VMFrac_T_from_p_h(_P[_qp], _h[_qp], _vmfrac[_qp], _T[_qp], _phase[_qp]);
  eos_uo.rho_m_by_p(_P[_qp], _h[_qp], _rho_m[_qp], _drho_m_dp[_qp], _drho_m_dp_2[_qp]);
  eos_uo.rho_m_by_h(_P[_qp], _h[_qp], _rho_m[_qp], _drho_m_dh[_qp], _drho_m_dh_2[_qp]);
  eos_uo.rho_m_by_ph(_P[_qp], _h[_qp], _drho_m_dph[_qp]);
  _rho_l[_qp] = eos_uo.rho_l_from_p_T(_P[_qp], _T[_qp], _phase[_qp]);
  _rho_g[_qp] = eos_uo.rho_g_from_p_T(_P[_qp], _T[_qp], _phase[_qp]);
  _vfrac[_qp]  = (_rho_m[_qp] - _rho_l[_qp]) / (_rho_g[_qp] - _rho_l[_qp]);
  _cp_m[_qp]  = eos_uo.cp_m_from_p_T(_P[_qp], _T[_qp], _vmfrac[_qp], _phase[_qp]);

  // To calculate the friction factor and Re No
  _u[_qp] = _flow[_qp] / _area[_qp];
  _Re[_qp] = _rho_m[_qp] * _dia[_qp] * _u[_qp] / viscosity_uo.mixture_mu(_P[_qp], _T[_qp], _vmfrac[_qp]);
  if (_f_defined)
    _friction[_qp] = _u_f;
  else
    MoodyFrictionFactor(_friction[_qp], _rel_roughness, _Re[_qp], _roughness_type);

  // drift-flux calculator section
    MoskitoDFGVar DFinp(_u[_qp], _rho_g[_qp], _rho_l[_qp], _vmfrac[_qp], _vfrac[_qp],
      _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
    dfm_uo.DFMCalculator(DFinp);
    DFinp.DFMOutput(_flow_pat[_qp], _c0[_qp], _u_d[_qp]);

  _rho_pam[_qp] = _rho_g[_qp] * _c0[_qp]  * _vfrac[_qp] + (1.0 - _vfrac[_qp] * _c0[_qp]) * _rho_l[_qp];

  PhaseVelocities();
  GammaDerivatives();
  KappaDerivatives();
  OmegaDerivatives();

  // _lambda[_qp]  = (1.0 - (_d * _d) / std::pow(_d + _thickness , 2.0)) * _lambda0;
  // _lambda[_qp] += (_d * _d) / std::pow(_d + _thickness , 2.0) * eos_uo._lambda;
}
void
MoskitoFluidWell_2p1c::PhaseVelocities()
{
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
}

Real
MoskitoFluidWell_2p1c::gamma(const Real & h, const Real & p, const Real & q)
{
  Real vmfrac, vfrac, T, phase, rho_l, rho_g, rho_m, rho_pam, dummy, c0, u_d;
  eos_uo.VMFrac_T_from_p_h(p, h, vmfrac, T, phase);
  rho_l = eos_uo.rho_l_from_p_T(p, T, phase);
  rho_g = eos_uo.rho_g_from_p_T(p, T, phase);
  rho_m = rho_l * rho_g / (vmfrac * (rho_l - rho_g) + rho_g);
  vfrac = (rho_m - rho_l) / (rho_g - rho_l);

  MoskitoDFGVar DFinp(q / _area[_qp], rho_g, rho_l, vmfrac, vfrac,
      _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
  dfm_uo.DFMCalculator(DFinp);
  DFinp.DFMOutput(dummy, c0, u_d);

  rho_pam = rho_g * c0  * vfrac + (1.0 - vfrac * c0) * rho_l;

  Real gamma = 0.0;
  gamma  = vfrac / (1.0 - vfrac);
  gamma *= rho_g * rho_l * rho_m / (rho_pam * rho_pam);
  gamma *= std::pow((c0 - 1.0) * q / _area[_qp] + u_d , 2.0);

  return gamma;
}

Real
MoskitoFluidWell_2p1c::kappa(const Real & h, const Real & p, const Real & q)
{
  Real vmfrac, vfrac, T, phase, rho_l, rho_g, rho_m, rho_pam, h_g, h_l, dummy, c0, u_d;
  eos_uo.VMFrac_T_from_p_h(p, h, vmfrac, T, phase);
  rho_l = eos_uo.rho_l_from_p_T(p, T, phase);
  rho_g = eos_uo.rho_g_from_p_T(p, T, phase);
  rho_m = rho_l * rho_g / (vmfrac * (rho_l - rho_g) + rho_g);
  vfrac = (rho_m - rho_l) / (rho_g - rho_l);

  MoskitoDFGVar DFinp(q / _area[_qp], rho_g, rho_l, vmfrac, vfrac,
      _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
  dfm_uo.DFMCalculator(DFinp);
  DFinp.DFMOutput(dummy, c0, u_d);

  rho_pam = rho_g * c0  * vfrac + (1.0 - vfrac * c0) * rho_l;
  eos_uo.h_lat(p, dummy, h_l, h_g);

  Real kappa = 0.0;
  kappa  = vfrac * rho_g * rho_l / rho_pam * (h_g - h_l);
  kappa *= std::pow((c0 - 1.0) * q / _area[_qp] + u_d , 2.0);

  return kappa;
}

Real
MoskitoFluidWell_2p1c::omega(const Real & h, const Real & p, const Real & q)
{
  Real vmfrac, vfrac, T, phase, rho_l, rho_g, rho_m, u_g, u_l, rho_pam, dummy, c0, u_d;
  eos_uo.VMFrac_T_from_p_h(p, h, vmfrac, T, phase);
  rho_l = eos_uo.rho_l_from_p_T(p, T, phase);
  rho_g = eos_uo.rho_g_from_p_T(p, T, phase);
  rho_m = rho_l * rho_g / (vmfrac * (rho_l - rho_g) + rho_g);
  vfrac = (rho_m - rho_l) / (rho_g - rho_l);

  MoskitoDFGVar DFinp(q / _area[_qp], rho_g, rho_l, vmfrac, vfrac,
      _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
  dfm_uo.DFMCalculator(DFinp);
  DFinp.DFMOutput(dummy, c0, u_d);

  rho_pam = rho_g * c0  * vfrac + (1.0 - vfrac * c0) * rho_l;

  u_g  = (c0 * rho_m * q / _area[_qp] + rho_l * u_d) / rho_pam;
  u_l  = (1.0 - vfrac * c0) * rho_m * q / _area[_qp] - rho_g * vfrac * u_d;
  u_l /= (1.0 - vfrac) * rho_pam;

  Real omega = 0.0;
  omega  = 3.0 * u_g * u_l * q / _area[_qp];
  omega -= (rho_g * vfrac * std::pow(u_l,3.0) + rho_l * (1.0 - vfrac) * std::pow(u_g,3.0)) / rho_m;
  omega *= 0.5 * vfrac * (1.0 - vfrac) * rho_g * rho_l / rho_m;

  return omega;
}

void
MoskitoFluidWell_2p1c::GammaDerivatives()
{
  _dgamma_dh[_qp] = 0.0; _dgamma_dp[_qp] = 0.0; _dgamma_dq[_qp] = 0.0;
  _dgamma2_dhq[_qp] = 0.0; _dgamma2_dpq[_qp] = 0.0; _dgamma2_dq2[_qp] = 0.0;

  if (_phase[_qp] == 2.0)
  {
    Real dh, dq, dp;
    Real tol = 1.0e-5;
    dh = tol * _h[_qp]; dp = tol * _P[_qp]; dq = tol * _flow[_qp];

    if (dh != 0.0)
    {
      _dgamma_dh[_qp]  = gamma(_h[_qp] + dh, _P[_qp], _flow[_qp]) - gamma(_h[_qp] - dh, _P[_qp], _flow[_qp]);
      _dgamma_dh[_qp] /= 2.0 * dh;
    }

    if (dp != 0.0)
    {
    _dgamma_dp[_qp]  = gamma(_h[_qp], _P[_qp] + dp, _flow[_qp]) - gamma(_h[_qp], _P[_qp] - dp, _flow[_qp]);
    _dgamma_dp[_qp] /= 2.0 * dp;
    }

    if (dq != 0.0)
    {
    _dgamma_dq[_qp]  = gamma(_h[_qp], _P[_qp], _flow[_qp] + dq) - gamma(_h[_qp], _P[_qp], _flow[_qp] - dq);
    _dgamma_dq[_qp] /= 2.0 * dq;
    }

    if (dh * dq != 0.0)
    {
    _dgamma2_dhq[_qp]  = gamma(_h[_qp] + dh, _P[_qp], _flow[_qp] + dq) + gamma(_h[_qp] - dh, _P[_qp], _flow[_qp] - dq);
    _dgamma2_dhq[_qp] -= gamma(_h[_qp] + dh, _P[_qp], _flow[_qp] - dq) + gamma(_h[_qp] - dh, _P[_qp], _flow[_qp] + dq);
    _dgamma2_dhq[_qp] /= 4.0 * dh * dq;
    }

    if (dp * dq != 0.0)
    {
    _dgamma2_dpq[_qp]  = gamma(_h[_qp], _P[_qp] + dp, _flow[_qp] + dq) + gamma(_h[_qp], _P[_qp] - dp, _flow[_qp] - dq);
    _dgamma2_dpq[_qp] -= gamma(_h[_qp], _P[_qp] + dp, _flow[_qp] - dq) + gamma(_h[_qp], _P[_qp] - dp, _flow[_qp] + dq);
    _dgamma2_dpq[_qp] /= 4.0 * dp * dq;
    }

    if (dq != 0.0)
    {
    _dgamma2_dq2[_qp]  = gamma(_h[_qp], _P[_qp], _flow[_qp] + dq) + gamma(_h[_qp], _P[_qp], _flow[_qp] - dq);
    _dgamma2_dq2[_qp] -= 2.0 * gamma(_h[_qp], _P[_qp], _flow[_qp]);
    _dgamma2_dq2[_qp] /=  dq * dq;
    }
  }
}

void
MoskitoFluidWell_2p1c::KappaDerivatives()
{
  _dkappa_dh[_qp] = 0.0; _dkappa_dp[_qp] = 0.0; _dkappa_dq[_qp] = 0.0;

  if (_phase[_qp] == 2.0)
  {
    Real dh, dq, dp;
    Real tol = 1.0e-5;
    dh = tol * _h[_qp]; dp = tol * _P[_qp]; dq = tol * _flow[_qp];

    if (dh != 0.0)
    {
      _dkappa_dh[_qp]  = kappa(_h[_qp] + dh, _P[_qp], _flow[_qp]) - kappa(_h[_qp] - dh, _P[_qp], _flow[_qp]);
      _dkappa_dh[_qp] /= 2.0 * dh;
    }

    if (dp != 0.0)
    {
    _dkappa_dp[_qp]  = kappa(_h[_qp], _P[_qp] + dp, _flow[_qp]) - kappa(_h[_qp], _P[_qp] - dp, _flow[_qp]);
    _dkappa_dp[_qp] /= 2.0 * dp;
    }

    if (dq != 0.0)
    {
    _dkappa_dq[_qp]  = kappa(_h[_qp], _P[_qp], _flow[_qp] + dq) - kappa(_h[_qp], _P[_qp], _flow[_qp] - dq);
    _dkappa_dq[_qp] /= 2.0 * dq;
    }
  }
}

void
MoskitoFluidWell_2p1c::OmegaDerivatives()
{
  _domega_dh[_qp] = 0.0; _domega_dp[_qp] = 0.0; _domega_dq[_qp] = 0.0;

  if (_phase[_qp] == 2.0)
  {
    Real dh, dq, dp;
    Real tol = 1.0e-5;
    dh = tol * _h[_qp]; dp = tol * _P[_qp]; dq = tol * _flow[_qp];

    if (dh != 0.0)
    {
      _domega_dh[_qp]  = omega(_h[_qp] + dh, _P[_qp], _flow[_qp]) - omega(_h[_qp] - dh, _P[_qp], _flow[_qp]);
      _domega_dh[_qp] /= 2.0 * dh;
    }

    if (dp != 0.0)
    {
    _domega_dp[_qp]  = omega(_h[_qp], _P[_qp] + dp, _flow[_qp]) - omega(_h[_qp], _P[_qp] - dp, _flow[_qp]);
    _domega_dp[_qp] /= 2.0 * dp;
    }

    if (dq != 0.0)
    {
    _domega_dq[_qp]  = omega(_h[_qp], _P[_qp], _flow[_qp] + dq) - omega(_h[_qp], _P[_qp], _flow[_qp] - dq);
    _domega_dq[_qp] /= 2.0 * dq;
    }
  }
}
