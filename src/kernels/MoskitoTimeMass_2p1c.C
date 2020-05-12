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

#include "MoskitoTimeMass_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoTimeMass_2p1c);

template <>
InputParameters
validParams<MoskitoTimeMass_2p1c>()
{
  InputParameters params = validParams<TimeKernel>();

  params.addRequiredCoupledVar("enthalpy", "Enthalpy nonlinear variable");
  params.addClassDescription("Time derivative part of mass conservation equation for "
                  "2 phase pipe flow and it returns pressure");

  return params;
}

MoskitoTimeMass_2p1c::MoskitoTimeMass_2p1c(const InputParameters & parameters)
  : TimeKernel(parameters),
    _h_dot(coupledDot("enthalpy")),
    _dh_dot(coupledDotDu("enthalpy")),
    _h_var_number(coupled("enthalpy")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dh(getMaterialProperty<Real>("drho_dh"))
{
}

Real
MoskitoTimeMass_2p1c::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_dp[_qp] * _u_dot[_qp];
  r += _drho_dh[_qp] * _h_dot[_qp];
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoTimeMass_2p1c::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_dp[_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
  j *= _test[_i][_qp];

  return j;
}

Real
MoskitoTimeMass_2p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _h_var_number)
  {
    j += _drho_dh[_qp] * _phi[_j][_qp] * _dh_dot[_qp];
    j *= _test[_i][_qp];
  }

  return j;
}
