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

#include "MoskitoEnergy_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoEnergy_2p1c);

template <>
InputParameters
validParams<MoskitoEnergy_2p1c>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addClassDescription("Energy conservation equation for 2 phase "
                                      "pipe flow and it returns enthalpy");

  return params;
}

MoskitoEnergy_2p1c::MoskitoEnergy_2p1c(const InputParameters & parameters)
  : Kernel(parameters),
  _q(coupledValue("flowrate")),
  _grad_q(coupledGradient("flowrate")),
  _grad_p(coupledGradient("pressure")),
  _q_var_number(coupled("flowrate")),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
  _lambda(getMaterialProperty<Real>("thermal_conductivity")),
  _cp(getMaterialProperty<Real>("specific_heat")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dh(getMaterialProperty<Real>("drho_dh")),
  _gravity(getMaterialProperty<RealVectorValue>("gravity")),
  _dkappa_dh(getMaterialProperty<Real>("dkappa_dh")),
  _dkappa_dp(getMaterialProperty<Real>("dkappa_dp")),
  _dkappa_dq(getMaterialProperty<Real>("dkappa_dq")),
  _domega_dh(getMaterialProperty<Real>("domega_dh")),
  _domega_dp(getMaterialProperty<Real>("domega_dp")),
  _domega_dq(getMaterialProperty<Real>("domega_dq"))
{
}

Real
MoskitoEnergy_2p1c::computeQpResidual()
{
  // r += _grad_test[_i][_qp] * _lambda[_qp] * _grad_u[_qp] / _cp[_qp];

  RealVectorValue r = 0.0;

  r += (_drho_dp[_qp] * _grad_p[_qp] + _drho_dh[_qp] * _grad_u[_qp]) * _q[_qp] * _u[_qp];
  r += _grad_q[_qp] * _rho[_qp] * _u[_qp];
  r += _rho[_qp] * _q[_qp] * (_grad_u[_qp] - _gravity[_qp]);
  r += _dkappa_dh[_qp] * _grad_u[_qp] + _dkappa_dp[_qp] * _grad_p[_qp] + _dkappa_dq[_qp] * _grad_q[_qp];
  r += _domega_dh[_qp] * _grad_u[_qp] + _domega_dp[_qp] * _grad_p[_qp] + _domega_dq[_qp] * _grad_q[_qp];
  r *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  return r * _well_dir[_qp];
}

Real
MoskitoEnergy_2p1c::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  j += _drho_dh[_qp] * _grad_phi[_j][_qp] * _q[_qp] * _u[_qp];
  j += (_drho_dp[_qp] * _grad_p[_qp] + _drho_dh[_qp] * _grad_u[_qp]) * _q[_qp] * _phi[_j][_qp];
  j += _grad_q[_qp] * _drho_dh[_qp] * _phi[_j][_qp] * _u[_qp];
  j += _grad_q[_qp] * _rho[_qp] * _phi[_j][_qp];
  j += _drho_dh[_qp] * _phi[_j][_qp] * _q[_qp] * (_grad_u[_qp] - _gravity[_qp]);
  j += _rho[_qp] * _q[_qp] * _grad_phi[_j][_qp];
  j += _dkappa_dh[_qp] * _grad_phi[_j][_qp];
  j += _domega_dh[_qp] * _grad_phi[_j][_qp];
  j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  return j * _well_dir[_qp];
}

Real
MoskitoEnergy_2p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;

  if (jvar == _q_var_number)
  {
    j += _drho_dp[_qp] * _grad_p[_qp] + _drho_dh[_qp] * _grad_u[_qp] * _phi[_j][_qp] * _u[_qp];
    j += _grad_phi[_j][_qp] * _rho[_qp] * _u[_qp];
    j += _rho[_qp] * _phi[_j][_qp] * (_grad_u[_qp] - _gravity[_qp]);
    j += _dkappa_dq[_qp] * _grad_phi[_j][_qp];
    j += _domega_dq[_qp] * _grad_phi[_j][_qp];
    j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];
  }

  if (jvar == _p_var_number)
  {
    j += _drho_dp[_qp] * _grad_phi[_j][_qp] * _q[_qp] * _u[_qp];
    j += _grad_q[_qp] * _drho_dp[_qp] * _phi[_j][_qp] * _u[_qp];
    j += _drho_dp[_qp] * _phi[_j][_qp] * _q[_qp] * (_grad_u[_qp] - _gravity[_qp]);
    j += _dkappa_dp[_qp] * _grad_phi[_j][_qp];
    j += _domega_dp[_qp] * _grad_phi[_j][_qp];
    j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];
  }

  return j * _well_dir[_qp];
}
