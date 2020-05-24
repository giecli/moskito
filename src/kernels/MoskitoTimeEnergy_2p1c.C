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

#include "MoskitoTimeEnergy_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoTimeEnergy_2p1c);

template <>
InputParameters
validParams<MoskitoTimeEnergy_2p1c>()
{
  InputParameters params = validParams<TimeKernel>();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addClassDescription("Time derivative part of energy conservation equation for "
                  "1 phase (either liquid or gas) pipe flow and it returns enthalpy");

  return params;
}

MoskitoTimeEnergy_2p1c::MoskitoTimeEnergy_2p1c(const InputParameters & parameters)
  : TimeKernel(parameters),
    _q(coupledValue("flowrate")),
    _p_dot(coupledDot("pressure")),
    _q_dot(coupledDot("flowrate")),
    _dp_dot(coupledDotDu("pressure")),
    _dq_dot(coupledDotDu("flowrate")),
    _p_var_number(coupled("pressure")),
    _q_var_number(coupled("flowrate")),
    _area(getMaterialProperty<Real>("well_area")),
    _rho(getMaterialProperty<Real>("density")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dh(getMaterialProperty<Real>("drho_dh")),
    _drho_dp_2(getMaterialProperty<Real>("drho_dp_2")),
    _drho_dh_2(getMaterialProperty<Real>("drho_dh_2")),
    _drho_dph(getMaterialProperty<Real>("drho_dph")),
    _dgamma_dh(getMaterialProperty<Real>("dgamma_dh")),
    _dgamma_dp(getMaterialProperty<Real>("dgamma_dp")),
    _dgamma_dq(getMaterialProperty<Real>("dgamma_dq")),
    _dgamma2_dhq(getMaterialProperty<Real>("dgamma2_dhq")),
    _dgamma2_dpq(getMaterialProperty<Real>("dgamma2_dpq")),
    _dgamma2_dq2(getMaterialProperty<Real>("dgamma2_dq2"))
{
}

Real
MoskitoTimeEnergy_2p1c::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_dp[_qp] * _p_dot[_qp];
  r += _drho_dh[_qp] * _u_dot[_qp];
  r *= _u[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]);
  r += _rho[_qp] * (_u_dot[_qp] + _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]));
  r -= _p_dot[_qp];
  r += _dgamma_dh[_qp] * _u_dot[_qp];
  r += _dgamma_dp[_qp] * _p_dot[_qp];
  r += _dgamma_dq[_qp] * _q_dot[_qp];

  return r * _test[_i][_qp];
}

Real
MoskitoTimeEnergy_2p1c::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_dph[_qp] * _p_dot[_qp];
  j += _drho_dh_2[_qp] * _u_dot[_qp];
  j += _drho_dh[_qp] * _du_dot_du[_qp];
  j *= _u[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]);
  j += (_drho_dp[_qp] * _p_dot[_qp] + _drho_dh[_qp] * _u_dot[_qp]);
  j += _drho_dh[_qp] * (_u_dot[_qp] + _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]));
  j += _rho[_qp] * _du_dot_du[_qp];
  j += _dgamma_dh[_qp] * _du_dot_du[_qp];
  j += _dgamma2_dhq[_qp] * _q_dot[_qp];

  return j * _test[_i][_qp] * _phi[_j][_qp];
}

Real
MoskitoTimeEnergy_2p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _q_var_number)
  {
    j += (_drho_dp[_qp] * _p_dot[_qp] + _drho_dh[_qp] * _u_dot[_qp]) * _q[_qp];
    j += _rho[_qp] * _q_dot[_qp];
    j += _rho[_qp] * _q[_qp] * _dq_dot[_qp];
    j /= _area[_qp] * _area[_qp];
    j += _dgamma2_dhq[_qp] * _u_dot[_qp];
    j += _dgamma2_dpq[_qp] * _p_dot[_qp];
    j += _dgamma2_dq2[_qp] * _q_dot[_qp];
    j += _dgamma_dq[_qp] * _dq_dot[_qp];
  }

  if (jvar == _p_var_number)
  {
    j += _drho_dp_2[_qp] * _p_dot[_qp];
    j += _drho_dph[_qp] * _u_dot[_qp];
    j += _drho_dp[_qp] * _dp_dot[_qp];
    j *= (_u[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]));
    j += _drho_dp[_qp] * (_u_dot[_qp] + _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]));
    j -= _dp_dot[_qp];
    j += _dgamma_dp[_qp] * _dp_dot[_qp];
    j += _dgamma2_dpq[_qp] * _q_dot[_qp];
  }

  return j * _test[_i][_qp] * _phi[_j][_qp];
}
