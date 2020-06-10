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

#include "MoskitoFluidWell_1p1c.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell_1p1c);

template <>
InputParameters
validParams<MoskitoFluidWell_1p1c>()
{
  InputParameters params = validParams<MoskitoFluidWellGeneral>();
  params.addRequiredCoupledVar("temperature", "Temperature nonlinear variable (K)");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable (Pa)");
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for EOS");
  params.addRequiredParam<UserObjectName>("viscosity_uo",
        "The name of the userobject for viscosity Eq");

  return params;
}

MoskitoFluidWell_1p1c::MoskitoFluidWell_1p1c(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    eos_uo(getUserObject<MoskitoEOS1P>("eos_uo")),
    viscosity_uo(getUserObject<MoskitoViscosity1P>("viscosity_uo")),
    _cp(declareProperty<Real>("specific_heat")),
    _rho(declareProperty<Real>("density")),
    _drho_dp(declareProperty<Real>("drho_dp")),
    _drho_dT(declareProperty<Real>("drho_dT")),
    _h(declareProperty<Real>("h_from_p_T")),
    _P(coupledValue("pressure")),
    _T(coupledValue("temperature"))
{
}

void
MoskitoFluidWell_1p1c::computeQpProperties()
{
  MoskitoFluidWellGeneral::computeQpProperties();

  _cp[_qp] = eos_uo.cp(_P[_qp], _T[_qp]);
  _h[_qp] = eos_uo.h_from_p_T(_P[_qp], _T[_qp]);
  eos_uo.rho_from_p_T(_P[_qp], _T[_qp], _rho[_qp], _drho_dp[_qp], _drho_dT[_qp]);

  _u[_qp] = _flow[_qp] / _area[_qp];
  _Re[_qp] = _rho[_qp] * _dia[_qp] * fabs(_u[_qp]) / viscosity_uo.mu(_P[_qp], _T[_qp]);
  if (_f_defined)
    _friction[_qp] = _u_f;
  else
    MoodyFrictionFactor(_friction[_qp], _rel_roughness, _Re[_qp], _roughness_type);

  _lambda[_qp]  = (1.0 - (_d * _d) / std::pow(_d + _thickness , 2.0)) * _lambda0;
  _lambda[_qp] += (_d * _d) / std::pow(_d + _thickness , 2.0) * eos_uo.lambda(_P[_qp], _T[_qp]);
}
