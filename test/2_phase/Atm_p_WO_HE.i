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
    output_properties = 'profile_mixture_density gas_velocity liquid_velocity void_fraction mass_fraction flow_pattern current_phase gas_density liquid_density density temperature well_velocity flow_type_c0 drift_velocity'
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 100000
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.01
  [../]
[]

[Variables]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      variable = p
      function = '90000+400*x'
    [../]
  [../]
  [./q]
    scaling = 1e-3
  [../]
[]

[AuxVariables]
  [./h]
    order = FIRST
    family = LAGRANGE
    initial_condition = 5e5
  [../]
[]

[Kernels]
  [./pkernel]
    type = MoskitoMass_2p1c
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./qkernel]
    type = MoskitoMomentum_2p1c
    variable = q
    pressure = p
    enthalpy = h
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
  nl_rel_tol = 1e-9
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
