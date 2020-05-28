[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 100
  nx = 50
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 0.0001
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityWaterSmith
  [../]
  [./viscosity_2p]
    type = MoskitoViscosity2P
    ve_uo_gas = viscosity_gas
    ve_uo_liquid = viscosity_liqid
  [../]
  [./df]
    type = MoskitoDFShi
    surface_tension = 0.0288
    Pan_param_cMax = 1.2
    Shi_param_Fv = 0.3
  [../]
  [./eos]
    type = MoskitoPureWater2P
  [../]
[]

[Materials]
  [./area]
    type = MoskitoFluidWell_2p1c
    pressure = p
    enthalpy = h
    flowrate = q
    well_type = production
    well_direction = x
    well_diameter = 0.0890016
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    roughness_type = smooth
    gravity = '10 0 0'
    outputs = exodus
    output_properties = 'void_fraction temperature'
  [../]
[]

[Variables]
  [./p]
    initial_condition = 2e6
    # scaling = 1e-3
  [../]
  [./h]
    initial_condition = 2e5
    # scaling = 1e-5
  [../]
  [./q]
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = '2e6'
  [../]
  [./qbc]
    type = FunctionDirichletBC
    variable = q
    boundary = left
    function = 'if(t>20,0.01,0.00001)'
  [../]
  [./hbc]
    type = DirichletBC
    variable = h
    boundary = right
    value = 5.3e5
  [../]
[]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy_2p1c
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./htimekernel]
    type = MoskitoTimeEnergy_2p1c
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./pkernel]
    type = MoskitoMass_2p1c
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./ptimekernel]
    type = MoskitoTimeMass_2p1c
    variable = p
    enthalpy = h
  [../]
  [./qkernel]
    type = MoskitoMomentum_2p1c
    variable = q
    pressure = p
    enthalpy = h
  [../]
  [./qtimekernel]
    type = MoskitoTimeMomentum_2p1c
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

[Dampers]
  # [./test]
  #   type = MaxIncrement
  #   variable = h
  #   max_increment = 0.1
  #   increment_type = fractional
  # [../]
  # [./test]
  #   type = BoundingValueNodalDamper
  #   variable = p
  #   min_value = 1e4
  # [../]
[]

[Executioner]
  type = Transient
  dt = 1
  end_time = 120
  l_max_its = 50
  nl_max_its = 50
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  solve_type = NEWTON
  automatic_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
