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

#ifndef MOSKITOFLUIDWELL1P_H
#define MOSKITOFLUIDWELL1P_H

#include "MoskitoFluidWellGeneral.h"
#include "MoskitoEOS1P.h"
#include "MoskitoViscosity1P.h"

class MoskitoFluidWell1P;

template <>
InputParameters validParams<MoskitoFluidWell1P>();

class MoskitoFluidWell1P : public MoskitoFluidWellGeneral
{
public:
  MoskitoFluidWell1P(const InputParameters & parameters);
  virtual void computeQpProperties() override;

protected:
  // Userobject to equation of state
  const MoskitoEOS1P & eos_uo;
  // Userobject to Viscosity Eq
  const MoskitoViscosity1P & viscosity_uo;

  // The specific heat at constant pressure
  MaterialProperty<Real> & _cp;
  // The density
  MaterialProperty<Real> & _rho;
  // The first derivative of density wrt pressure
  MaterialProperty<Real> & _drho_dp;
  // The first derivative of density wrt temperature
  MaterialProperty<Real> & _drho_dT;
  // The first derivative of density wrt enthalpy
  MaterialProperty<Real> & _drho_dh;

};

#endif /* MOSKITOFLUIDWELL1P_H */
