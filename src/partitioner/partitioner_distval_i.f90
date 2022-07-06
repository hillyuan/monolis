program monolis_partitioner_distval
  use mod_monolis_mesh
  use mod_monolis_io
  use mod_monolis_io_arg
  implicit none
  integer(kint) :: n_domain, nnode, ndof
  type(monolis_mesh), allocatable :: mesh(:)
  integer(kint), allocatable :: val(:,:)
  character :: fname*100, label*100

  call monolis_get_part_bc_arg(n_domain, fname)

  call monolis_input_mesh_distval_i(fname, nnode, ndof, val, label)

  call monolis_par_input_node_id(n_domain, mesh)

  call monolis_par_output_distval_i(n_domain, mesh, fname, ndof, val, label)

end program monolis_partitioner_distval