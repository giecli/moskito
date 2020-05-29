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

#pragma once

#include "Kernel.h"

class MoskitoLateralHeat_1p_matrix;

template <>
InputParameters validParams<MoskitoLateralHeat_1p_matrix>();

class MoskitoLateralHeat_1p_matrix : public Kernel
{
public:
  MoskitoLateralHeat_1p_matrix(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // Radius tubing outer
  const MaterialProperty<Real> & _rto;
  // Thermal wellbore resistivity
  const MaterialProperty<Real> & _Uto;
  // Temperature at formation - cement boundary
  const MaterialProperty<Real> & _Twb;
  // Diameter filled with liquid = _rti
  const MaterialProperty<Real> & _diameter_liquid;
  const Real gradC_to_gradR = 1.8;
  const Real Watt_to_Btu_per_h = 3.412141633;
  const Real m_to_ft    = 3.280839895;
  const Real Rankine_absol = 491.67;
  const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;
};
