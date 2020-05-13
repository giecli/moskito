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

class MoskitoFluidWell_2p1c;

template <>
InputParameters validParams<MoskitoFluidWell_2p1c>();

class MoskitoFluidWell_2p1c : public MoskitoFluidWellGeneral
{
public:
  MoskitoFluidWell_2p1c(const InputParameters & parameters);
  virtual void computeQpProperties() override;
  void PhaseVelocities();
  void GammaDerivatives();
  void KappaDerivatives();
  void OmegaDerivatives();
  Real gamma(const Real & h, const Real & p, const Real & q);
  Real kappa(const Real & h, const Real & p, const Real & q);
  Real omega(const Real & h, const Real & p, const Real & q);

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

  // The gamma first derivatives
  MaterialProperty<Real> & _dgamma_dh;
  // The gamma first derivatives
  MaterialProperty<Real> & _dgamma_dp;
  // The gamma first derivatives
  MaterialProperty<Real> & _dgamma_dq;
  // The gamma second derivatives
  MaterialProperty<Real> & _dgamma2_dhq;
  // The gamma second derivatives
  MaterialProperty<Real> & _dgamma2_dpq;
  // The gamma second derivatives
  MaterialProperty<Real> & _dgamma2_dq2;

  // The kappa first derivatives
  MaterialProperty<Real> & _dkappa_dh;
  // The kappa first derivatives
  MaterialProperty<Real> & _dkappa_dp;
  // The kappa first derivatives
  MaterialProperty<Real> & _dkappa_dq;

  // The omega first derivatives
  MaterialProperty<Real> & _domega_dh;
  // The omega first derivatives
  MaterialProperty<Real> & _domega_dp;
  // The omega first derivatives
  MaterialProperty<Real> & _domega_dq;

  // The coupled enthalpy
  const VariableValue & _h;

  // The gradient of the coupled variables
  const VariableGradient & _grad_flow;
  const VariableGradient & _grad_h;
  const VariableGradient & _grad_p;
};
