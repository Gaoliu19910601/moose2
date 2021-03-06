[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 3
  nz = 0
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  active = 'u v'

  [./u]
    order = FIRST
    family = LAGRANGE
  [../]

  [./v]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff1 diff2'

  [./diff1]
    type = DiffMKernel
    variable = u
    mat_prop = diff1
  [../]

  [./diff2]
    type = DiffMKernel
    variable = v
    mat_prop = diff2
  [../]
[]

[BCs]
  active = 'left_u right_u left_v right_v'

  # Mesh Generation produces boundaries in counter-clockwise fashion
  [./left_u]
    type = DirichletBC
    variable = u
    boundary = 3
    value = 0
  [../]

  [./right_u]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 1
  [../]

  [./left_v]
    type = DirichletBC
    variable = v
    boundary = 3
    value = 1
  [../]

  [./right_v]
    type = DirichletBC
    variable = v
    boundary = 1
    value = 0
  [../]
[]

[Materials]
  active = 'dm1 dm2'

  [./dm1]
    type = GenericConstantMaterial
    block = 0
    prop_names  = 'diff1'
    prop_values = '2'
  [../]

  [./dm2]
    type = GenericConstantMaterial
    block = 0
    prop_names  = 'diff2'
    prop_values = '4'
  [../]
[]

[Executioner]
  type = Steady

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK'
[]

[Outputs]
  file_base = out
  output_initial = true
  exodus = true
  print_perf_log = true
[]
