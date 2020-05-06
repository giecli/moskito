# Validation of Heat exchange using the Paper from Satman & Tureyen 2016
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
    top_right = '250 1 0'
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
  [./eos]
    type = MoskitoEOS1P_PureWater
  [../]
  [./HK]
    type = MoskitoDFHK
    surface_tension = 0.0288
  [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '0.5695652 * x'
  [../]
[]

[GlobalParams]
  pressure = p
  enthalpy = h
  flowrate = q
  well_direction = x
  eos_uo = eos
  viscosity_uo = viscosity_2p
  # drift_flux_uo = HK
  roughness_type = smooth
  gravity = '9.8 0 0'
  outputs = exodus
  output_properties = 'temperature density well_velocity specific_heat well_reynolds_no well_moody_friction viscosity well_diameter flow_pattern gas_velocity liquid_velocity mass_fraction gas_density liquid_density void_fraction current_phase'
[]

[Materials]
  [./tubing]
    type = MoskitoFluidWell1P
    well_diameter = 0.102
    block = 0
  [../]
  [./openhole]
    type = MoskitoFluidWell1P
    well_diameter = 0.099
    block = 1
  [../]
  # [./Lateral_tubing]
  #   type = MoskitoLatHeatIterationXiong
  #    block = 0
  #    radius_wellbore = 0.0510000001
  #    radius_tubbing_outer = 0.0510000002
  #    conductivity_tubing = 40
  #    # Rock parameters
  #    Surface_temperature = 295.15
  #    thermal_diffusivity_rock = 1.102e-6
  #    conductivity_rock = 1.3
  #    # Configuration of material
  #    geothermal_gradient = grad_func
  #    hc_calucation_model = Dropkin_Sommerscales
  #    time_model = user_time
  #    # user_defined_time = 8640
  #    user_defined_time = 2592000
  #    DimTime_calculation_model = Kutun_2015_eq20
  #    internal_solve_full_iteration_history = false
  #    output_properties = 'Temperature_Rankine temperature_well_formation_interface formation_temperature'
  #    outputs = exodus
  #  [../]
   # [./Lateral_openhole]
   #   type = MoskitoLatHeatIterationXiong
   #    block = 1
   #    radius_wellbore = 0.0495000001
   #    radius_tubbing_outer = 0.0495000002
   #    conductivity_tubing = 40
   #    # Rock parameters
   #    Surface_temperature = 295.15
   #    thermal_diffusivity_rock = 1.102e-6
   #    conductivity_rock = 1.3
   #    # Configuration of material
   #    geothermal_gradient = grad_func
   #    hc_calucation_model = Dropkin_Sommerscales
   #    time_model = user_time
   #    # user_defined_time = 8640
   #    user_defined_time = 2592000
   #    DimTime_calculation_model = Kutun_2015_eq20
   #    internal_solve_full_iteration_history = false
   #    output_properties = 'Temperature_Rankine temperature_well_formation_interface formation_temperature'
   #    outputs = exodus
   #  [../]
[]

[Variables]
  [./h]
    # [./InitialCondition]
    #   type = FunctionIC
    #   variable = h
    #   function = 850000-(248.4-x)*500
    # [../]
    scaling = 1e-4
    # initial_condition = 1300000
    initial_condition = 691431.822
  [../]
  [./p]
    # [./InitialCondition]
    #   type = FunctionIC
    #   variable = p
    #   function = 1.918e6-(248.4-x)*6100
    # [../]
    initial_condition = 1.918e6
        # scaling = 1e4
  [../]
  [./q]
    initial_condition = 0.0088
    # scaling = 1e-2
  [../]
[]

[BCs]
  # [./qbc]
  #   # type = MoskitoMassFlowRateCoupled
  #   type = MoskitoMassFlowRate
  #   variable = q
  #   # pressure = p
  #   # enthalpy = h
  #   # eos_uo = eos
  #   mass_flowrate = 7.08
  #   mixture_density = 950
  #   boundary = left
  # [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = left
    value = 0.0088
  [../]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = right
    value = 1.918e6
  [../]
  # [./hbc1]
  #   type = DirichletBC
  #   variable = h
  #   boundary = right
  #   value = 691431.822
  # [../]
  # [./hbc2]
  #   type = DirichletBC
  #   variable = h
  #   boundary = left
  #   value = 800000
  # [../]
  # [./hbc]
  #   type = MoskitoTemperatureToEnthalpy1P
  #   variable = h
  #   pressure = p
  #   boundary = right
  #   temperature = 436.63
  #   eos_uo = eos
  # [../]
[../]

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
[]

[Executioner]
  type = Steady
  l_max_its = 50
  l_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 50
  solve_type = NEWTON
  nl_abs_tol = 1e-7
  # automatic_scaling = true
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  [./test]
    type = VariableResidualNormsDebugOutput
  [../]
  console = true
[]
