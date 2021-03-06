# Use the exodus file for restarting the problem:
# - restart one variable
# - and have one extra variable
# - have PBP active to have more system in Equation system
#

[Mesh]
  file = out_xda_restart_part1.e
[]

[Functions]
  [./exact_fn]
    type = ParsedFunction
    value = t*((x*x)+(y*y))
  [../]

  [./forcing_fn]
    type = ParsedFunction
    value = -4+(x*x+y*y)
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]

  [./v]
    order = FIRST
    family = LAGRANGE
  [../]

  [./diffusivity]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./ie]
    type = TimeDerivative
    variable = u
  [../]

  [./diff]
    type = MatDiffusion
    variable = u
    prop_name = diffusivity
  [../]

  [./ffn]
    type = UserForcingFunction
    variable = u
    function = forcing_fn
  [../]

  [./diff_v]
    type = Diffusion
    variable = v
  [../]

  [./out_diffusivity]
    type = RealPropertyOutput
    variable = diffusivity
    prop_name = diffusivity
  [../]
[]

[Materials]
  [./mat]
    type = StatefulMaterial
    block = 0
  [../]
[]

[BCs]
  [./all]
    type = FunctionDirichletBC
    variable = u
    boundary = '0 1 2 3'
    function = exact_fn
  [../]

  [./left_v]
    type = DirichletBC
    variable = v
    boundary = '3'
    value = 0
  [../]

  [./right_v]
    type = DirichletBC
    variable = v
    boundary = '1'
    value = 1
  [../]
[]

[Preconditioning]
  [./PBP]
    type = PBP
    solve_order = 'u v diffusivity'
    preconditioner = 'AMG AMG AMG'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = JFNK

  restart_file_base = out_xda_restart_part1_cp/0005

  start_time = 1

  dt = 0.1
  reset_dt = true #NECESSARY to force a change in DT when using restart!

  num_steps = 3
[]

[Outputs]
  file_base = out_xda_restart_part2
  output_initial = true
  exodus = true
[]
