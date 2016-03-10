[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 3
  ny = 3
  nz = 1
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  zmin = 0.0
  zmax = 0.25
  elem_type = HEX8
[]

[XFEM]
  cut_type = 'square_cut_3d' # rectangular cut plane
  cut_data = ' -0.001 0.5 -0.001
                1.001 0.5 -0.001
                1.001 0.5  1.001
               -0.001 0.5  1.001'
  qrule = volfrac
  output_cut_plane = true
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[Constraints]
  [./xfem_constraint]
    type = XFEMEqualValueConstraint
    variable = u
    xfem_interface_id = 1
    jump = 0.2
    jump_flux = 0.1
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    boundary = top
    value = 2
  [../]

  [./right_u]
    type = DirichletBC
    variable = u
    boundary = bottom
    value = 0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
#  petsc_options_iname = '-pc_type -pc_hypre_type'
#  petsc_options_value = 'hypre boomeramg'
   petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'lu                      101'

  line_search = 'none'

  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 1.0
  end_time = 1.0
[]

[Outputs]
  file_base = diffusion_3d_out
  interval = 1
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]