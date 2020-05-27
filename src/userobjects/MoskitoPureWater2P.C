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

#include "MoskitoPureWater2P.h"

registerMooseObject("MoskitoApp", MoskitoPureWater2P);

template <>
InputParameters
validParams<MoskitoPureWater2P>()
{
  InputParameters params = validParams<MoskitoEOS2P>();
  return params;
}

MoskitoPureWater2P::MoskitoPureWater2P(const InputParameters & parameters)
  : MoskitoEOS2P(parameters)
{
  std::string eos_name = UserObjectName(name() + ":LiquidGas");
  {
    std::string class_name = "MoskitoWater97FluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, eos_name, params);
  }
  _eos_lg = &_fe_problem.getUserObject<MoskitoWater97FluidProperties>(eos_name);
}

void
MoskitoPureWater2P::VMFrac_T_from_p_h(
  const Real & pressure, const Real & enthalpy, Real & vmfrac, Real & temperature, Real & phase) const
{
  unsigned int region = _eos_lg->inRegionPH(pressure, enthalpy);
  switch (region)
  {
    case 1:
      vmfrac = 0.0;
      phase = 0;
      temperature = _eos_lg->temperature_from_ph(pressure, enthalpy);
      break;

    case 2:
      vmfrac = 1.0;
      phase = 1;
      temperature = _eos_lg->temperature_from_ph(pressure, enthalpy);
      break;

    case 3:
      vmfrac = 0.0;
      phase = 3;
      temperature = _eos_lg->temperature_from_ph(pressure, enthalpy);
      break;

    case 4:
      {
        phase = 2;
        Real hsat_l, hlat, hsat_g;
        h_lat(pressure, hlat, hsat_l, hsat_g);
        vmfrac  = enthalpy - hsat_l;
        vmfrac /= hlat;
        temperature = _eos_lg->vaporTemperature(pressure);
        break;
      }

    case 5:
      vmfrac = 1.0;
      phase = 1;
      temperature = _eos_lg->temperature_from_ph(pressure, enthalpy);
      break;

    default:
      mooseError(name(), ": inRegionPH() has given an incorrect region");
  }
}

void
MoskitoPureWater2P::rho_g_from_p_T(
  const Real &  pressure, const Real &  temperature, Real & rho, Real & drho_dp, Real & drho_dT, const Real & phase) const
{
  switch ((unsigned int)phase)
  {
    case 0:
      rho = 0.0;
      drho_dp = 0.0;
      drho_dT = 0.0;
      break;

    case 1:
      _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
      break;

    case 2:
      _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT,2);
      break;

    case 3:
      _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT,3);
      break;
  }
}

void
MoskitoPureWater2P::rho_l_from_p_T(
  const Real &  pressure, const Real &  temperature, Real & rho, Real & drho_dp, Real & drho_dT, const Real & phase) const
{
  switch ((unsigned int)phase)
  {
    case 0:
      _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT,1);
      break;

    case 1:
      rho = 0.0;
      drho_dp = 0.0;
      drho_dT = 0.0;
      break;

    case 2:
      _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT,1);
      break;

    case 3:
      _eos_lg->rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT,3);
      break;
  }
}

Real
MoskitoPureWater2P::rho_g_from_p_T(const Real & pressure, const Real & temperature, const Real & phase) const
{
  Real rho = 0.0;
  switch ((unsigned int)phase)
  {
    case 0:
      rho = 0.0;
      break;

    case 1:
      rho = _eos_lg->rho_from_p_T(pressure, temperature);
      break;

    case 2:
      rho = _eos_lg->rho_from_p_T(pressure, temperature,2);
      break;

    case 3:
      rho = _eos_lg->rho_from_p_T(pressure, temperature,3);
      break;
  }

  return rho;
}

Real
MoskitoPureWater2P::rho_l_from_p_T(const Real &  pressure, const Real &  temperature, const Real & phase) const
{
  Real rho = 0.0;
  switch ((unsigned int)phase)
  {
    case 0:
      rho = _eos_lg->rho_from_p_T(pressure, temperature,1);
      break;

    case 1:
      rho = 0.0;
      break;

    case 2:
      rho = _eos_lg->rho_from_p_T(pressure, temperature,1);
      break;

    case 3:
      rho = _eos_lg->rho_from_p_T(pressure, temperature,3);
      break;
  }

  return rho;
}

Real
MoskitoPureWater2P::cp_m_from_p_T(
      const Real & pressure, const Real & temperature, const Real & vmfrac, const Real & phase) const
{
  Real cp = 0.0;

  if (phase == 2.0)
  {
    cp  = _eos_lg->cp_from_p_T(pressure, temperature, 2) * vmfrac;
    cp += _eos_lg->cp_from_p_T(pressure, temperature, 1) * (1.0 - vmfrac);
  }
  else
    cp = _eos_lg->cp_from_p_T(pressure, temperature);

  return cp;
}

Real
MoskitoPureWater2P::rho_m_from_p_h(const Real & pressure, const Real & enthalpy) const
{
  Real rho,temperature,phase,vmfrac;

  VMFrac_T_from_p_h(pressure, enthalpy, vmfrac, temperature, phase);

  switch ((unsigned int)phase)
  {
    case 0:
      rho = rho_l_from_p_T(pressure, temperature, phase);
      break;

    case 1:
      rho = rho_g_from_p_T(pressure, temperature, phase);
      break;

    case 2:
      {
        Real rhol = rho_l_from_p_T(pressure, temperature, phase);
        Real rhog = rho_g_from_p_T(pressure, temperature, phase);

        rho  = rhol * rhog;
        rho /= vmfrac * (rhol - rhog) + rhog;
      }
      break;

    case 3:
      rho = rho_l_from_p_T(pressure, temperature, phase);
      break;
  }
  return rho;
}

void
MoskitoPureWater2P::h_lat(const Real & pressure, Real & hlat, Real & hsatl, Real & hsatg) const
{
  Real temperature = _eos_lg->vaporTemperature(pressure);
  hsatg = _eos_lg->h_from_p_T(pressure, temperature, 2);
  hsatl = _eos_lg->h_from_p_T(pressure, temperature, 1);
  hlat = hsatg - hsatl;
}
