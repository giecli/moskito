[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 3000
  nx = 100
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 0.0001
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
  [./viscosity_2p]
    type = MoskitoViscosity2P
    ve_uo_gas = viscosity_gas
    ve_uo_liquid = viscosity_liqid
  [../]
  [./df]
    type = MoskitoDFHK
    surface_tension = 0.0288
  [../]
  [./eos]
    type = MoskitoPureWater2P
    derivative_tolerance = 1e-5
  [../]
[]
#
# [AuxVariables]
#   [./maz]
#   [../]
# []
#
# [AuxKernels]
#   [./maz1]
#     type =
#   [../]

[Materials]
  [./area]
    type = MoskitoFluidWell_2p1c
    well_diameter = 0.1
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_type = production
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
    output_properties = ''
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 1e5
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.01
  [../]
  [./hbc]
    type = DirichletBC
    variable = h
    boundary = left
    value = 7e5
  [../]
[]

[Variables]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      variable = p
      function = '1e5+10000*x'
    [../]
  [../]
  [./q]
    initial_condition = 0.01
    # scaling = 1e-5
  [../]
  [./h]
    initial_condition = 7e5
    # scaling = 1e-5
  [../]
[]

[Kernels]
  [./pkernel]
    type = NullKernel
    variable = p
  [../]
  [./qkernel]
    type = NullKernel
    variable = q
  [../]
  [./hkernel]
    type = NullKernel
    variable = h
  [../]
[]

[Preconditioning]
  active = pn1
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
  nl_max_its = 50
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-4
  solve_type = NEWTON
  automatic_scaling = true
[]

[Outputs]
  exodus = true
  [./test]
    type = VariableResidualNormsDebugOutput
  [../]
[]
