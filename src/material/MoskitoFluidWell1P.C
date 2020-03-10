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

#include "MoskitoFluidWell1P.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell1P);

template <>
InputParameters
validParams<MoskitoFluidWell1P>()
{
  InputParameters params = validParams<MoskitoFluidWellGeneral>();
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for EOS");
  params.addRequiredParam<UserObjectName>("viscosity_uo",
        "The name of the userobject for viscosity Eq");

  return params;
}

MoskitoFluidWell1P::MoskitoFluidWell1P(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    eos_uo(getUserObject<MoskitoEOS1P>("eos_uo")),
    viscosity_uo(getUserObject<MoskitoViscosity1P>("viscosity_uo")),
    _cp(declareProperty<Real>("specific_heat")),
    _rho(declareProperty<Real>("density")),
    _drho_dp(declareProperty<Real>("drho_dp")),
    _drho_dT(declareProperty<Real>("drho_dT")),
    _drho_dh(declareProperty<Real>("drho_dh"))
{
}

void
MoskitoFluidWell1P::computeQpProperties()
{
  _T[_qp] = eos_uo.h_to_T(_P[_qp], _h[_qp]);
  _cp[_qp] = eos_uo.cp(_P[_qp], _T[_qp]);

  eos_uo.rho_from_p_T(_P[_qp], _T[_qp], _rho[_qp], _drho_dp[_qp], _drho_dT[_qp]);
  _drho_dh[_qp] = _drho_dT[_qp] / _cp[_qp];

  _dia[_qp] = _d;
  if (!(_area_defined*_perimeter_defined))
  {
    _area[_qp] = PI * _dia[_qp] * _dia[_qp] / 4.0;
    _perimeter[_qp] = PI * _dia[_qp];
  }
  else
  {
    _area[_qp] = _u_area;
    _perimeter[_qp] = _u_perimeter;
  }

  _u[_qp] = _flow[_qp] / _area[_qp];
  _Re[_qp] = _rho[_qp] * _dia[_qp] * fabs(_u[_qp]) / viscosity_uo.mu(_P[_qp], _T[_qp]);

  _lambda[_qp]  = (1.0 - (_d * _d) / std::pow(_d + _thickness , 2.0)) * _lambda0;
  _lambda[_qp] += (_d * _d) / std::pow(_d + _thickness , 2.0) * eos_uo.lambda(_P[_qp], _T[_qp]);

  MoskitoFluidWellGeneral::computeQpProperties();
}
