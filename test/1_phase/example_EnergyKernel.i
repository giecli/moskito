
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 100
  nx = 50
[]

[UserObjects]
  [./eos]
    type = MoskitoEOS1P_IdealFluid
  [../]
  [./viscosity]
    type = MoskitoViscosityWaterSmith
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell1P
    pressure = p
    temperature = T
    flowrate = q
    well_type = production
    well_direction = x
    well_diameter = 0.0890016
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '10.0 0 0'
  [../]
[]

[Variables]
  [./T]
    initial_condition = 300
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
    variable = T
    boundary = left
    value = 300
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
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
