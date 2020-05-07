[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 200
  nx = 200
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 1e-3
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 1e-3
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
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_type = injection
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    well_diameter = 0.05
    gravity = '9.8 0 0'
    outputs = exodus
    output_properties = 'well_velocity void_fraction flow_pattern current_phase drift_velocity flow_type_c0 volFrac_HK superficial_gas_velocity superficial_liquid_velocity'
  [../]
[]

[Variables]
  [./h]
    initial_condition = 4.25e5
    # initial_condition = enthalpy?
  [../]
  [./p]
    initial_condition = 1e5
  [../]
  [./q]
    [./InitialCondition]
      type = FunctionIC
      variable = q
      function = 'if(x<100,1e-2*x,(x-99)*2)*3.85e-3'
    [../]
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
  [../]
  [./pkernel]
    type = NullKernel
    variable = p
  [../]
  [./qkernel]
    type = NullKernel
    variable = q
  [../]
[]

[Preconditioning]
  active = p3
  [./p1]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
  [./p2]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -sub_pc_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu newtonls basic NONZERO 51'
  [../]
  [./p3]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
  [../]
  [./p4]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu NONZERO 51'
  [../]
[]

[Executioner]
  type = Steady
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  [./test]
    type = VariableResidualNormsDebugOutput
  [../]
[]
