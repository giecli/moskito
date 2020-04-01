
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

class MoskitoMomentum_2p;

template <>
InputParameters validParams<MoskitoMomentum_2p>();

class MoskitoMomentum_2p : public Kernel
{
public:
  MoskitoMomentum_2p(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // The gradient of the coupled pressure
  const VariableGradient & _grad_p;
  // The gradient of the coupled specific enthalpy
  const VariableGradient & _grad_h;

  // Variable numberings
  unsigned _p_var_number;
  unsigned _h_var_number;

  // The mixture density
  const MaterialProperty<Real> & _rho;
  // The first derivative of mixture density wrt pressure
  const MaterialProperty<Real> & _drho_dp;
  // The second derivative of mixture density wrt pressure
  const MaterialProperty<Real> & _drho_dp_2;
  // The first derivative of mixture density wrt enthalpy
  const MaterialProperty<Real> & _drho_dh;
  // The second derivative of mixture density wrt enthalpy
  const MaterialProperty<Real> & _drho_dh_2;
  // The pipe diameter
  const MaterialProperty<Real> & _d;
  // The pipe Moody friction factor
  const MaterialProperty<Real> & _f;
  // The gravity acceleration as a vector
  const MaterialProperty<RealVectorValue> & _gravity;
  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
  // The flow direction
  const MaterialProperty<Real> & _flow_dir;

  // residual for dgamma_dz in the momentum conservation
  const MaterialProperty<Real> * _dgamma_dz = NULL;
  // diagonal jacobian of the residual wrt uj for dgamma_dz in the momentum conservation
  const MaterialProperty<Real> * _dgamma_dz_uj_gphi = NULL;
  const MaterialProperty<Real> * _dgamma_dz_uj_phi = NULL;
  // diagonal jacobian of the residual wrt pj for dgamma_dz in the momentum conservation
  const MaterialProperty<Real> * _dgamma_dz_pj_gphi = NULL;
  const MaterialProperty<Real> * _dgamma_dz_pj_phi = NULL;
  // diagonal jacobian of the residual wrt hj for dgamma_dz in the momentum conservation
  const MaterialProperty<Real> * _dgamma_dz_hj_gphi = NULL;
  const MaterialProperty<Real> * _dgamma_dz_hj_phi = NULL;
private:
  // an indicator for 2 phase flow model
  bool _2p_ind = false;
};
