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

#include "MoskitoLateralHeat_1p_matrix.h"

registerMooseObject("MoskitoApp", MoskitoLateralHeat_1p_matrix);

template <>
InputParameters
validParams<MoskitoLateralHeat_1p_matrix>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Lateral heat exchange between wellbore "
        "and formation including tubing (mandatory), insulation, liquid filled "
        "annulus and cementation");
  return params;
}

MoskitoLateralHeat_1p_matrix::MoskitoLateralHeat_1p_matrix(const InputParameters & parameters)
  : Kernel(parameters),
  _rto(getMaterialProperty<Real>("radius_tubbing_outer")),
  _Uto(getMaterialProperty<Real>("thermal_resistivity_well")),
  _Twb(getMaterialProperty<Real>("temperature_well_formation_interface")),
  _diameter_liquid(getMaterialProperty<Real>("well_diameter"))
  {
  }

Real
MoskitoLateralHeat_1p_matrix::computeQpResidual()
{
  Real r = 0.0;
  r =  2.0 * PI * _rto[_qp] * _Uto[_qp] * (_u[_qp] - _Twb[_qp]);
  r /=  PI * _diameter_liquid[_qp] * _diameter_liquid[_qp] / 4.0;

  return  r * _test[_i][_qp];
}

Real
MoskitoLateralHeat_1p_matrix::computeQpJacobian()
{
  Real j = 0.0;
  j =  2.0 * PI * _rto[_qp] * _Uto[_qp] * _phi[_j][_qp];
  j /=  PI * _diameter_liquid[_qp] * _diameter_liquid[_qp] / 4.0;

  return  j * _test[_i][_qp];
}
