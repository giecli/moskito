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

class MoskitoLatHeatIterationXiong;

template <>
InputParameters validParams<MoskitoLatHeatIterationXiong>();

class MoskitoLatHeatIterationXiong : public Material, public NewtonIteration
{
public:
  MoskitoLatHeatIterationXiong(const InputParameters & parameters);
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
  Real TemperatureFormation();
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
  // Radius tubing outer
  MaterialProperty<Real> & _RadTubout;
  // Variable to output formation temperature
  MaterialProperty<Real> & _TRock;
    // Emissivity of inside casing surface
  Real _Annulus_eao;
  // Emissivity of outside tubin/insulation surface
  Real _Annulus_eai;
  // Density annulus fluid
  Real _Annulus_rho;
  // Dyn Viscosity annulus fluid
  Real _Annulus_nu;
  // Thermal conductivity annulus fluid
  Real _Annulus_lambda;
  // Heat capacity annulus fluid
  Real _Annulus_cp;
  //  Thermal volumetric expansion coefficient annulus fluid
  Real _Annulus_beta;
  // Annulus fluid thermal diffusivity = Temperaturleitfähigkeit
  Real _Annulus_alpha;
  // Thermal conductivoty formation (earth)
  Real _Rock_lambda;
  // Formation thermal diffusivity
  Real _Rock_alpha;
  // Thermal wellbore resistivity
  MaterialProperty<Real> & _Uto;
  //  Temperature at wellbore formation interface
  MaterialProperty<Real> & _Twb;
  // Diameter of the pipe
  const MaterialProperty<Real> & _diameter;
  // Conductivity vector
  std::vector<Real> _conductivity_vector;
  // The gravity acceleration as a vector
  const MaterialProperty<RealVectorValue> & _gravity;
  // Function for geothermal gradient, has to be defined in the Input
  const Function & _gradT;
  // mixing approach
    MooseEnum _hc;
  enum HC_case {Dropkin_Sommerscales, Raithby_Hollands, Churchill};
  // Ramey function parameter
  Real ft;
  // Definition of user time for steady state simulation or transient simulation time
  Real Timing;
  // Definition of calculation scheme for dimensionless time
  MooseEnum _Dim_time;
  enum Dim_time_case {Ramey_1962, Kutasov_2003_eq15, Kutasov_2003_eq16, Kutasov_2003_eq17, Kutasov_2003_eq18, Kutasov_1987_eq19, Kutun_2015_eq20 };
    // Annulus outer and inner Rasius, is determined within this material
  Real _rai = 0.0, _rao = 0.0;
  // Tolerance of finite difference derivation
  const Real _tol;
  // Independent gravity for calculation of Raleigh and Grashof numers in case gravity is set to "0"
  RealVectorValue _independ_gravity;
  // Combine all diameters in one vector
  std::vector<Real> _well_assembly;
  // Well direction for correction of gravity vector in terms of deviated well
  const MaterialProperty<RealVectorValue> & _well_dir;
  // Stefan Bolzmann Konstante in SI units (W/(m² * K⁴))
  Real Boltz = 0.00000005670367;
  Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;

  const std::array<Real, 14> column_vec{
    {0.01, 0.02, 0.05,
     0.1,  0.2,  0.5,
     1.0,  2.0,  5.0,
     10.0, 20.0, 50.0,
     100.0, 1000.0}};

  const std::array<Real, 10> row_vec{
    {0.1,  0.2,  0.5,
     1.0,  2.0,  5.0,
     10.0, 20.0, 50.0,
     100.0}};

  const std::array<std::array<Real, 14>, 10> Ramey{
   {{{0.313, 0.313, 0.314, 0.316, 0.318, 0.323, 0.330, 0.345, 0.373, 0.396, 0.417, 0.433, 0.438, 0.445}},
    {{0.423, 0.423, 0.424, 0.427, 0.430, 0.439, 0.452, 0.473, 0.511, 0.538, 0.568, 0.572, 0.578, 0.588}},
    {{0.616, 0.617, 0.619, 0.623, 0.629, 0.644, 0.666, 0.698, 0.745, 0.772, 0.790, 0.802, 0.806, 0.811}},
    {{0.802, 0.803, 0.806, 0.811, 0.820, 0.842, 0.872, 0.910, 0.958, 0.984,  1.00,  1.01,  1.01,  1.02}},
    {{ 1.02,  1.02,  1.03,  1.04,  1.05,  1.08,  1.11,  1.15,  1.20,  1.22,  1.24,  1.24,  1.25,  1.25}},
    {{ 1.36,  1.37,  1.37,  1.38,  1.40,  1.44,  1.48,  1.52,  1.56,  1.57,  1.58,  1.59,  1.59,  1.59}},
    {{ 1.65,  1.66,  1.66,  1.67,  1.69,  1.73,  1.77,  1.81,  1.84,  1.86,  1.86,  1.87,  1.87,  1.88}},
    {{ 1.96,  1.97,  1.97,  1.99,  2.00,  2.05,  2.09,  2.12,  2.15,  2.16,  2.16,  2.17,  2.17,  2.17}},
    {{ 2.39,  2.39,  2.40,  2.42,  2.44,  2.48,  2.51,  2.54,  2.56,  2.57,  2.57,  2.57,  2.58,  2.58}},
    {{ 2.73,  2.73,  2.74,  2.75,  2.77,  2.81,  2.84,  2.86,  2.88,  2.89,  2.89,  2.89,  2.89,  2.90}}}};
};
