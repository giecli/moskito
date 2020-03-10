# Over-all heat transfer coefficients in Steam and hot water injection wells
# Willhite, G. P. 19
# Appendix: Sample Calculation - Results of Uto differs slighly, because of rounding of the author

# Example has a lenght that is not of importance
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 100
  nx = 50
[]

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
[]

[Variables]
  [./h]
    initial_condition = 3000000
  [../]
  [./p]
    initial_condition = 1.0e6
  [../]
  [./q]
    initial_condition = 0.01
  [../]
[]

[BCs]
  [./hbc_top]
    type = DirichletBC
    variable = h
    boundary = right
    value = 3000000
  [../]
[../]

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
