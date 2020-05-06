# Validation of Heat exchange using the Paper from Satman & Tureyen 2016
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 100
  nx = 20
[]

[MeshModifiers]
  [./openhole]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '50 -1 0'
    top_right = '100 1 0'
  [../]
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 1e-5
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 1e-4
  [../]
  [./viscosity_2p]
    type = MoskitoViscosity2P
    ve_uo_gas = viscosity_gas
    ve_uo_liquid = viscosity_liqid
  [../]
  [./eos21]
    type = MoskitoPureWater2P
  [../]
  # [./DF]
  #   type = MoskitoDFShi
  #   inclination_vector = '1.85 0.21 0.95'
  #   surface_tension = 0.0288
  # [../]
  [./DF]
    type = MoskitoDFHK
    surface_tension = 0.0288
  [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '373'
  [../]
  [./vm_func]
    type = ParsedFunction
    value = '0'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_diameter = 0.102
    eos_uo = eos21
    drift_flux_uo = DF
    viscosity_uo = viscosity_2p
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
    block = 0
    output_properties = 'temperature density well_velocity specific_heat well_reynolds_no well_moody_friction viscosity diameter flow_pattern gas_velocity liquid_velocity mass_fraction gas_density liquid_density void_fraction current_phase superficial_gas_velocity superficial_liquid_velocity flow_type_c0'
  [../]
  [./area1]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_diameter = 0.99
    eos_uo = eos21
    drift_flux_uo = DF
    viscosity_uo = viscosity_2p
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
    block = 1
    output_properties = 'temperature density well_velocity specific_heat well_reynolds_no well_moody_friction viscosity diameter flow_pattern gas_velocity liquid_velocity mass_fraction gas_density liquid_density void_fraction current_phase superficial_gas_velocity superficial_liquid_velocity flow_type_c0'
  [../]
  # [./Lateral]
  #   type = MoskitoLatHeatIterationXiong
  #    radius_wellbore = 0.0510000001
  #    radius_tubbing_outer = 0.0510000002
  #    conductivity_tubing = 40
  #    # Rock parameters
  #    Surface_temperature = 288.15
  #    thermal_diffusivity_rock = 1.102e-6
  #    conductivity_rock = 1.3
  #    # Configuration of material
  #    geothermal_gradient = grad_func
  #    hc_calucation_model = Dropkin_Sommerscales
  #    time_model = user_time
  #    # user_defined_time = 8640
  #    user_defined_time = 2592000
  #    DimTime_calculation_model = Kutun_2015_eq20
  #    internal_solve_full_iteration_history = true
  #    output_properties = 'Temperature_Rankine temperature_well_formation_interface heat_loss formation_temperature'
  #    outputs = exodus
  #  [../]
[]

[Variables]
  [./h]
    # scaling = 1e-2
    initial_condition = 1400000
  [../]
  [./p]
    initial_condition = 2e6
    scaling = 1e-3
  [../]
  [./q]
    initial_condition = 0.02
    scaling = 1e-4
  [../]
[]

[BCs]
  # [./qbc]
  #   type = MoskitoMassFlowRate
  #   enthalpy = h
  #   mass_flowrate = 7.08
  #   mixture_density = 1000
  #   boundary = left
  # [../]
  # [./qbc]
  #   type = MoskitoMassFlowRateCoupled
  #   variable = q
  #   pressure = p
  #   enthalpy = h
  #   eos_uo = eos21
  #   mass_flowrate = 7.1
  #   boundary = left
  # [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.02
  [../]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = right
    value = 2e6
  [../]
  # [./hbc]
  #   type = DirichletBC
  #   variable = h
  #   boundary = right
  #   value = 691431.822
  # [../]
  # [./hbc]
  #   type = MoskitoTemperatureToEnthalpy1P
  #   variable = h
  #   pressure = p
  #   boundary = right
  #   temperature = 500.45
  #   eos_uo = eos
  # [../]
[../]

# [ICs]
#   # [./hic]
#   #   type = ConstantIC
#   #   variable = h
#   #   value = 1100000
#   # [../]
#   [./hic2]
#     type = MoskitoTemp2Enthalpy2P
#     variable = h
#     pressure = p
#     geothermal_gradient = grad_func
#     vmfrac_gradient = vm_func
#   [../]
# []


[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./qkernel1]
    type = MoskitoMomentum
    variable = q
    pressure = p
    enthalpy = h
  [../]
  [./pkernel1]
    type = MoskitoMass
    variable = p
    enthalpy = h
    flowrate = q
  [../]
  # [./heat]
  #   type = MoskitoHeat
  #   variable = h
  # [../]
[]

[Preconditioning]
  active = p2
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
    petsc_options_value = 'asm ilu newtonls basic NONZERO 51'
  [../]
  [./p3]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
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
  print_linear_residuals = true
  [./test]
    type = VariableResidualNormsDebugOutput
  [../]
  console = true
[]
