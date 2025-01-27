program monolis_refiner
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_io
  use mod_monolis_io_arg
  implicit none
  type(monolis_mesh) :: mesh
  logical :: is_format_id

  call monolis_global_initialize()

  call monolis_set_debug(.true.)
  call monolis_debug_header("monolis_p_refiner for HEX element")
  call monolis_debug_header("*** function does not supported")

  !call monolis_get_refiner_arg(is_format_id)
  !is_format_id = .true.
  !call monolis_input_mesh(mesh, is_format_id)
  !call monolis_output_refined_mesh(mesh, is_format_id)
  call monolis_global_finalize()
end program monolis_refiner
