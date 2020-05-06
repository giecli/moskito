# Over-all heat transfer coefficients in Steam and hot water injection wells
# Willhite, G. P. 19
# Appendix: Sample Calculation - Results of Uto differs slighly, because of rounding of the author

# Example has a lenght that is not of importance
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 100
  nx = 100
[]

# [Functions]
#   [./grad_func]
#     type = ParsedFunction
#     value = '360'
#   [../]
# []

[UserObjects]
  [./eos]
    type = MoskitoEOS1P_PureWater
  [../]
  [./viscosity]
    type = MoskitoViscosityWaterSmith
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_type = production
    well_direction = x
    well_diameter = 0.0890016
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0.0 0 0'
    outputs = exodus
    output_properties = 'temperature'
  [../]
  # [./Lateral]
  #   type = MoskitoLatHeatIterationXiong
  #    # Geometry of the well. As the example did not contain any tubing radius, which is required for teh material it was artificially set to a small radius
  #    well_diameter_vector = '0.0890016 0.089001602 0.216408 0.24384 0.3048'
  #    conductivities_vector = '80.42534006 0.0 80.42534006 0.346146923'
  #    # Rock parameters
  #    thermal_diffusivity_rock = 0.000000738063
  #    conductivity_rock = 1.7307346
  #    # Annulus parameters representing a stagnant gas at 14.7psia @ 296°C
  #    density_annulus = 0.6215162
  #    conductivity_annulus = 0.0441337
  #    dyn_viscosity_annulus = 0.0000285262
  #    capacity_annulus = 1025.76594
  #    thermal_expansion_annulus = 1.755e-3
  #    emissivity_annulus_outer = 0.9
  #    emissivity_annulus_inner = 0.9
  #    # Configuration of material
  #    geothermal_gradient = grad_func
  #    hc_calucation_model = Dropkin_Sommerscales
  #    DimTime_calculation_model = Ramey_1962
  #    user_defined_time = 1814400
  #    internal_solve_full_iteration_history = true
  #    outputs = exodus
  #    output_properties = 'thermal_resistivity_well'
  #  [../]
[]

[Variables]
  [./h]
    initial_condition = 350000
  [../]
  [./p]
    initial_condition = 5.0e6
  [../]
  [./q]
    initial_condition = 0.01
  [../]
[]

[BCs]
  [./hbc_top]
    type = DirichletBC
    variable = h
    boundary = left
    value = 300000
  [../]
[../]

# [BCs]
#   [./hbc]
#     type = MoskitoTemperatureToEnthalpy1P
#     variable = h
#     pressure = p
#     boundary = right
#     temperature = 333
#     eos_uo = eos
#   [../]
# [../]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy_1p1c
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./pkernel1]
    type = NullKernel
    variable = p
  [../]
  [./qkernel1]
    type = NullKernel
    variable = q
  [../]
  # [./heat]
  #   type = MoskitoLateralHeat
  #   variable = h
  # [../]
[]

[Preconditioning]
  active = pn1
  [./p2]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
  [./p1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                 '
  [../]
  [./pn1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -snes_type -snes_linesearch_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                   newtonls   basic               '
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 50
  l_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 50
  solve_type = NEWTON
  nl_abs_tol = 1e-7
[]

[Outputs]
  exodus = true
[]
