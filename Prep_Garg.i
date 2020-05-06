[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 248.4
  nx = 300
[]

[MeshModifiers]
  [./openhole]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '159.6 -1 0'
    top_right =   '250 1 0'
  [../]
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
    derivative_tolerance = 1e-3
    execute_on = LINEAR
  [../]
[]

[GlobalParams]
  pressure = p
  enthalpy = h
  flowrate = q
  well_direction = x
  eos_uo = eos
  viscosity_uo = viscosity_2p
  drift_flux_uo = df
  roughness_type = smooth
  gravity = '9.8 0 0'
  outputs = exodus
  output_properties = 'gas_velocity liquid_velocity void_fraction mass_fraction flow_pattern current_phase gas_density liquid_density density temperature well_velocity'
[]


[Materials]
  [./area]
    type = MoskitoFluidWell2P
    well_diameter = 0.1
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = right
    value = 1.918e6
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.0008
  [../]
[]

[Variables]
  [./h]
    [./InitialCondition]
      type = FunctionIC
      variable = h
      function = 7e5
    [../]
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      variable = p
      # function = '1.918e6 - 6920 * (248.4 - x)'
      function = '5e5+x*6000'
    [../]
  [../]
  [./q]
    initial_condition = 0.0008
    scaling = 1e-7
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
  [../]
  [./pkernel]
    type = MoskitoMass
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./qkernel]
    type = MoskitoMomentum
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
    petsc_options_value = ' bjacobi  ilu          NONZERO                   newtontr   basic               '
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 50
  nl_max_its = 50
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  # [./test]
  #   type = VariableResidualNormsDebugOutput
  # [../]
[]
