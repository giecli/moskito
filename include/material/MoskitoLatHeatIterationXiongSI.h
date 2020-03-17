/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2019 by Maziar Gholami Korzani                          */
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

#include "Material.h"
#include "Function.h"
#include "NewtonIteration.h"

class MoskitoLatHeatIterationXiongSI;

template <>
InputParameters validParams<MoskitoLatHeatIterationXiong>();

class MoskitoLatHeatIterationXiongSI : public Material, public NewtonIteration
{
public:
  MoskitoLatHeatIterationXiongSI(const InputParameters & parameters);
  virtual void computeQpProperties() override;

  virtual Real computeReferenceResidual(const Real trail_value, const Real scalar) override;
  virtual Real computeResidual(const Real trail_value, const Real scalar) override;
  virtual Real computeDerivative(const Real trailn_value, const Real scalar) override;
  virtual Real initialGuess(const Real trail_value) override;
  virtual Real minimumPermissibleValue(const Real trail_value) const override
  {
    return 0.000000001;
  }
  virtual Real maximumPermissibleValue(const Real trail_value) const override
  {
    return 100000.0;
  }

protected:
  // Calculate transient time function according to Ramey
  Real transienttimefunction(Real Uto);
  // Calculate deoth dependent formation temperature
  Real TemperatureFormation(Real _Tsurf);
  // Calculate temperature at formation wellbore
  Real TemperatureWBinterface(Real Uto, Real Temp);
  // Calculate temperature at annulus casing boundary
  Real TemperatureCasingAnnulusInterface(Real Uto, Real Temp);
  // Calculate tubing external temperature
  Real TemperatureTubingOuter(Real Uto, Real Temp);
  // Calculation of Annuls effect - radial heat transfer coefficient hr
  Real RadialHeatTransferCoefficient(Real Uto, Real T);
  // Calculation of Annuls effect - convective heat transfer coefficient hc
  Real ConvectiveHeatTransferCoefficient(Real Uto, Real T, Real grav);
  // Calculate Grashof Number needed for hc calculation
  Real Grashof(Real Uto, Real Temp, Real grav);
  // Calculate Rayleigh Number needed for hc calculation
  Real Rayleigh(Real grav,Real Uto, Real Temp);
  // temperature
  const MaterialProperty<Real> & _T;
  // Radius wellbore
  Real _rwb;
  // Radius tubing outer
  MaterialProperty<Real> & _RadTubout;
  // Variable to output heat loss
  MaterialProperty<Real> & _may;
  // Variable to output formation temperature
  MaterialProperty<Real> & _TRock;
  // Local parameter Radius tubing outer
  Real _rto;
  // Radius insulation
  Real _rins;
  // Radius casing inner
  Real _rci;
  // Radius casing outer respectivly cement
  Real _rcem;
  // Emissivity of inside casing surface
  Real _eao;
  // Emissivity of outside tubin/insulation surface
  Real _eai;
  // Density annulus fluid
  Real _rhoAnnulus;
  // Dyn Viscosity annulus fluid
  Real _nuAnnulus;
  // Thermal conductivity annulus fluid
  Real _lambdaAnnulus;
  // Heat capacity annulus fluid
  Real _cpAnnulus;
  //  Thermal volumetric expansion coefficient annulus fluid
  Real _betaAnnulus;
  // Annulus fluid thermal diffusivity = Temperaturleitfähigkeit
  Real _alphaAnnulus;
  // Thermal conductivoty formation (earth)
  Real _lambdaRock;
  // Formation thermal diffusivity
  Real _alphaRock;
  // Thermal wellbore resistivity
  MaterialProperty<Real> & _Uto;
  //  Temperature at wellbore formation interface
  MaterialProperty<Real> & _Twb;
  // Surface temperature
  Real _Tsurf;
  // Diameter of the pipe
  const MaterialProperty<Real> & _diameter;
  // Radius tubing inner
  Real _rti;
  // Thermal conductivity of the cement
  Real _lambdaCem;
  // Thermal conductivity of the casing
  Real _lambdaCas;
  // Thermal conductivity of the insulation
  Real _lambdaIns;
  // Thermal conductivity of the tubing
  Real _lambdaTub;
  // The gravity acceleration as a vector
  const MaterialProperty<RealVectorValue> & _gravity;
  // Function for geothermal gradient, has to be defined in the Input
  const Function & gradT;
  // mixing approach
    MooseEnum _hc;
  enum HC_case {Dropkin_Sommerscales, Raithby_Hollands, Churchill};
  // Ramey function parameter
  Real ft;
  // Definition of user time for steady state simulation or transient simulation time
  MooseEnum _td;
  enum Zeit {user_time, simulation_time};
  Real Timing;
  Real _ut;
  // Annulus outer and inner Rasius, is determined within this material
  Real _rai, _rao;
  // Tolerance of finite difference derivation
  const Real _tol;
  // Independent gravity for calculation of Raleigh and Grashof numers in case gravity is set to "0"
  RealVectorValue _independ_gravity;
  // Well direction for correction of gravity vector in terms of deviated well
  const  MaterialProperty<RealVectorValue> & _well_dir;
  // Stefan Bolzmann Konstante in SI units (W/(m² * K⁴))
  Real Boltz = 0.00000005670367;
  const Real s_to_h = 3600;
  // Convert from Si units to American system
  const Real m_to_ft    = 3.280839895;
  const Real m2_to_ft2  = 10.7639079;
  const Real m3_to_ft3  = 35.3146667;
  const Real kg_to_lbm  = 2.2046225;
  const Real gradC_to_gradR = 1.8;
  const Real Watt_to_Btu_per_h = 3.412141633;
  const Real J_to_Btu = 0.00094781712;
  const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;
  const Real Rankine_absol = 491.67;

};
