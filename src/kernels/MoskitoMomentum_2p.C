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

#include "MoskitoMomentum_2p.h"

registerMooseObject("MoskitoApp", MoskitoMomentum_2p);

template <>
InputParameters
validParams<MoskitoMomentum_2p>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable");
  params.addClassDescription("Momentum conservation equation for 2 phase "
        "(liquid and gas) pipe flow and it returns flowrate. A positive "
        "value in BCs is the production scenario and the flow is in the opposite"
        " way of the well direction as defined in the material, and vice versa");
  return params;
}

MoskitoMomentum_2p::MoskitoMomentum_2p(const InputParameters & parameters)
  : Kernel(parameters),
    _grad_p(coupledGradient("pressure")),
    _grad_h(coupledGradient("enthalpy")),
    _p_var_number(coupled("pressure")),
    _h_var_number(coupled("enthalpy")),
    _rho(getMaterialProperty<Real>("density")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dp_2(getMaterialProperty<Real>("drho_dp_2")),
    _drho_dh(getMaterialProperty<Real>("drho_dh")),
    _drho_dh_2(getMaterialProperty<Real>("drho_dh_2")),
    _d(getMaterialProperty<Real>("well_diameter")),
    _f(getMaterialProperty<Real>("well_moody_friction")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    _area(getMaterialProperty<Real>("well_area")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
    _flow_dir(getMaterialProperty<Real>("flow_direction_sign"))
{
  if (hasMaterialProperty<Real>("dgamma_dz"))
  {
    _dgamma_dz = &getMaterialProperty<Real>("dgamma_dz");
    _dgamma_dz_uj_gphi = &getMaterialProperty<Real>("dgamma_dz_uj_gphi");
    _dgamma_dz_uj_phi = &getMaterialProperty<Real>("dgamma_dz_uj_phi");
    _dgamma_dz_pj_gphi = &getMaterialProperty<Real>("dgamma_dz_pj_gphi");
    _dgamma_dz_pj_phi = &getMaterialProperty<Real>("dgamma_dz_pj_phi");
    _dgamma_dz_hj_gphi = &getMaterialProperty<Real>("dgamma_dz_hj_gphi");
    _dgamma_dz_hj_phi = &getMaterialProperty<Real>("dgamma_dz_hj_phi");
    _2p_ind = true;
  }
}

Real
MoskitoMomentum_2p::computeQpResidual()
{
  RealVectorValue r = 0.0;

  r += _drho_dp[_qp] * _grad_p[_qp];
  r += _drho_dh[_qp] * _grad_h[_qp];
  r *= _flow_dir[_qp] * _u[_qp] * _u[_qp];
  r += _flow_dir[_qp] * 2.0 * _rho[_qp] * _u[_qp] * _grad_u[_qp];
  r += _flow_dir[_qp] * _f[_qp] * _rho[_qp] * _u[_qp] * _u[_qp]
        * _well_dir[_qp] / (2.0 * _d[_qp]);
  r /= (_area[_qp] * _area[_qp]);
  if (_2p_ind)
    r += (*_dgamma_dz)[_qp] * _well_dir[_qp];
  r -= _grad_p[_qp];
  r += _rho[_qp] * _gravity[_qp];
  r *= _test[_i][_qp];

  return r * _well_dir[_qp];
}

Real
MoskitoMomentum_2p::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  j += _drho_dp[_qp] * _grad_p[_qp];
  j += _drho_dh[_qp] * _grad_h[_qp];
  j *= _flow_dir[_qp] * 2.0 * _phi[_j][_qp] * _u[_qp];
  j += _flow_dir[_qp] * 2.0 * _rho[_qp] * (_phi[_j][_qp]  * _grad_u[_qp]
        + _u[_qp] * _grad_phi[_j][_qp]);
  j += _flow_dir[_qp] * _f[_qp] * _rho[_qp] * _phi[_j][_qp]
        * _u[_qp]  * _well_dir[_qp] / _d[_qp];
  j /= (_area[_qp] * _area[_qp]);
  if (_2p_ind)
  {
    j += (*_dgamma_dz_uj_phi)[_qp] * _phi[_j][_qp] * _well_dir[_qp];
    j += (*_dgamma_dz_uj_gphi)[_qp] * _grad_phi[_j][_qp];
  }
  j *= _test[_i][_qp];

  return j * _well_dir[_qp];
}

Real
MoskitoMomentum_2p::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;

  if (jvar == _p_var_number)
  {
    j += _drho_dp_2[_qp] * _phi[_j][_qp] * _grad_p[_qp];
    j += _drho_dp[_qp] * _grad_phi[_j][_qp];
    j *= _flow_dir[_qp] * _u[_qp] * _u[_qp];
    j += _flow_dir[_qp] * 2.0 * _drho_dp[_qp] * _phi[_j][_qp] * _u[_qp] * _grad_u[_qp];
    j += _flow_dir[_qp] * _f[_qp] * _drho_dp[_qp] * _phi[_j][_qp] * _u[_qp]
          * _u[_qp] * _well_dir[_qp] / (2.0 * _d[_qp]);
    j /= (_area[_qp] * _area[_qp]);
    j -= _grad_phi[_j][_qp];
    j += _drho_dp[_qp] * _phi[_j][_qp] * _gravity[_qp];
    if (_2p_ind)
    {
      j += (*_dgamma_dz_pj_phi)[_qp] * _phi[_j][_qp] * _well_dir[_qp];
      j += (*_dgamma_dz_pj_gphi)[_qp] * _grad_phi[_j][_qp];
    }
    j *= _test[_i][_qp];
  }

  if (jvar == _h_var_number)
  {
    j += _drho_dh_2[_qp] * _phi[_j][_qp] * _grad_h[_qp];
    j += _drho_dh[_qp] * _grad_phi[_j][_qp];
    j *= _flow_dir[_qp] * _u[_qp] * _u[_qp];
    j += _flow_dir[_qp] * 2.0 * _drho_dh[_qp] * _phi[_j][_qp] * _u[_qp] * _grad_u[_qp];
    j += _flow_dir[_qp] * _f[_qp] * _drho_dh[_qp] * _phi[_j][_qp] * _u[_qp]
          * _u[_qp] * _well_dir[_qp] / (2.0 * _d[_qp]);
    j /= (_area[_qp] * _area[_qp]);
    j += _drho_dh[_qp] * _phi[_j][_qp] * _gravity[_qp];
    if (_2p_ind)
    {
      j += (*_dgamma_dz_hj_phi)[_qp] * _phi[_j][_qp] * _well_dir[_qp];
      j += (*_dgamma_dz_hj_gphi)[_qp] * _grad_phi[_j][_qp];
    }
    j *= _test[_i][_qp];
  }

  return j * _well_dir[_qp];
}
