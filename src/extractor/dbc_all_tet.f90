program monolis_dbc_all_tet
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_io
  use mod_monolis_io_arg
  use mod_monolis_dbc_all_util
  use mod_monolis_util_debug
  implicit none
  type(monolis_mesh) :: mesh
  integer(kint) :: n_block
  character :: output_dir*200, fnname*100, fename*100, foname*100
  logical :: is_format_id
  integer(kint), allocatable :: is_surf_node(:)
  real(kdouble), allocatable :: val(:)

  call monolis_global_initialize()

  call monolis_set_debug(.true.)
  call monolis_debug_header("monolis_dbc_all_tet")
  is_format_id = .false.

  call monolis_get_dbc_all_arg(n_block, val, fnname, fename, foname)

  call monolis_input_mesh(mesh, is_format_id)

  if(mesh%nbase_func /= 4)then
    call monolis_debug_header("ERROR: please input 1st order hex mesh")
    stop 1
  endif

  call monolis_get_surf_node(mesh, 4, 4, 3, is_surf_node)

  call output_dbc(foname, mesh, is_surf_node, n_block, val)

  call monolis_global_finalize()
end program monolis_dbc_all_tet

