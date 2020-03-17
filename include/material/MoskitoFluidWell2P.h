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

#include "MoskitoFluidWellGeneral.h"
#include "MoskitoEOS2P.h"
#include "MoskitoViscosity2P.h"
#include "MoskitoDriftFlux.h"

class MoskitoFluidWell2P;

template <>
InputParameters validParams<MoskitoFluidWell2P>();

class MoskitoFluidWell2P : public MoskitoFluidWellGeneral
{
public:
  MoskitoFluidWell2P(const InputParameters & parameters);
  virtual void computeQpProperties() override;
  void DriftFluxMomentumEq();

protected:
  // Userobject to equation of state
  const MoskitoEOS2P & eos_uo;
  // Userobject to Viscosity Eq
  const MoskitoViscosity2P & viscosity_uo;
  // Userobject to Drift Fluc model
  const MoskitoDriftFlux & dfm_uo;

  // temperature
  MaterialProperty<Real> & _T;
  // The specific heat of mixture at constant pressure
  MaterialProperty<Real> & _cp_m;
  // Density of gas
  MaterialProperty<Real> & _rho_g;
  // Density of liquid
  MaterialProperty<Real> & _rho_l;
  // Density of mixture
  MaterialProperty<Real> & _rho_m;
  // Profile-adjusted density of mixture
  MaterialProperty<Real> & _rho_pam;
  // The first derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_m_dp;
  // The second derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_m_dp_2;
  // The first derivative of mixture density wrt temperature
  MaterialProperty<Real> & _drho_m_dT;
  // The first derivative of mixture density wrt enthalpy
  MaterialProperty<Real> & _drho_m_dh;
  // The second derivative of mixture density wrt enthalpy
  MaterialProperty<Real> & _drho_m_dh_2;
  // mass_fraction
  MaterialProperty<Real> & _vmfrac;
  // Gas velocity
  MaterialProperty<Real> & _u_g;
  // Liquid velocity
  MaterialProperty<Real> & _u_l;

  // void_fraction
  MaterialProperty<Real> & _vfrac;
  // current phase
  MaterialProperty<Real> & _phase;
  // drift velocity
  MaterialProperty<Real> & _u_d;
  // flow type parameter
  MaterialProperty<Real> & _c0;
  // flow pattern
  MaterialProperty<Real> & _flow_pat;
  // superficial velocities calculated by HK
  MaterialProperty<Real> &_v_sg;
  MaterialProperty<Real> &_v_sl;

  //refer to DriftFluxMomentumEq function
  // residual for dgamma_dz in the momentum conservation
  MaterialProperty<Real> & _dgamma_dz;
  // diagonal jacobian of the residual wrt uj for dgamma_dz in the momentum conservation
  MaterialProperty<Real> & _dgamma_dz_uj_gphi;
  MaterialProperty<Real> & _dgamma_dz_uj_phi;
  // diagonal jacobian of the residual wrt pj for dgamma_dz in the momentum conservation
  MaterialProperty<Real> & _dgamma_dz_pj_gphi;
  MaterialProperty<Real> & _dgamma_dz_pj_phi;
  // diagonal jacobian of the residual wrt hj for dgamma_dz in the momentum conservation
  MaterialProperty<Real> & _dgamma_dz_hj_gphi;
  MaterialProperty<Real> & _dgamma_dz_hj_phi;

  // The coupled enthalpy
  const VariableValue & _h;

  // The gradient of the coupled variables
  const VariableGradient & _grad_flow;
  const VariableGradient & _grad_h;
  const VariableGradient & _grad_p;
};
