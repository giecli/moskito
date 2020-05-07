# Over-all heat transfer coefficients in Steam and hot water injection wells
# Willhite, G. P. 19
# Appendix: Sample Calculation - Results of Uto differs slighly, because of rounding of the author

# Example has a lenght that is not of importance
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1
  nx = 1
[]

[UserObjects]
  [./eos]
    type = MoskitoEOS1P_IdealFluid
  [../]
  [./viscosity]
    type = MoskitoViscosityWaterSmith
  [../]
[]

# No temperature gradient used in the example
[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '310.928'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell1P
    pressure = p
    temperature = T
    flowrate = q
    well_type = injection
    well_direction = x
    well_diameter = 0.0890016
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '9.8 0 0'
    # outputs = exodus
  [../]
  [./Lateral]
    type = MoskitoLatHeat_1p
     temperature = T
     # Geometry of the well. As the example did not contain any tubing radius, which is required for teh material it was artificially set to a small radius
     well_diameter_vector = '0.0890016 0.089001602 0.216408 0.24384 0.3048'
     conductivities_vector = '80.42534006 0.0 80.42534006 0.346146923'
     # Rock parameters
     thermal_diffusivity_rock = 0.000000738063
     conductivity_rock = 1.7307346
     # Annulus parameters representing a stagnant gas at 14.7psia @ 296°C
     density_annulus = 0.6215162
     conductivity_annulus = 0.0441337
     dyn_viscosity_annulus = 0.0000285262
     capacity_annulus = 1025.76594
     thermal_expansion_annulus = 1.755e-3
     emissivity_annulus_outer = 0.9
     emissivity_annulus_inner = 0.9
     # Configuration of material
     geothermal_gradient = grad_func
     hc_calucation_model = Dropkin_Sommerscales
     DimTime_calculation_model = Ramey_1962
     user_defined_time = 1814400
     internal_solve_full_iteration_history = true
     outputs = exodus
     output_properties = 'thermal_resistivity_well'
   [../]
[]

[Variables]
  [./T]
    initial_condition = 588.70555
  [../]
  [./p]
    initial_condition = 1.0e6
  [../]
  [./q]
    initial_condition = 0.01
  [../]
[]

[BCs]
  [./hbc]
    type = DirichletBC
    variable = T
    boundary = left
    value = 588.70555
  [../]
[../]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy_1p1c
    variable = T
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
  [./heat]
    type = MoskitoLateralHeat_1p
    variable = T
  [../]
[]

[Preconditioning]
  active = p1
  [./p1]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 2
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-4
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
