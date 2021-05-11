module mod_monolis_dbc_all_util
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_hash
  use mod_monolis_stdlib
  implicit none

  type(type_monolis_hash_tree) :: hash_tree

contains

  subroutine monolis_get_surf_node(mesh, nbase_func, nsurf, nsurf_node, is_surf_node)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: itable(3,4), nbase_func, conn(nbase_func), nsurf, nsurf_node
    integer(kint) :: i, j, in, eid, i1, i2, i3, tmp
    character :: ckey*27
    logical :: is_exist, is_pushed
    integer(kint), allocatable :: is_inner(:,:), is_surf_node(:)

!    itable(1,1) = 1; itable(2,1) = 2; itable(3,1) = 3
!    itable(1,2) = 1; itable(2,2) = 2; itable(3,2) = 4
!    itable(1,3) = 2; itable(2,3) = 3; itable(3,3) = 4
!    itable(1,4) = 3; itable(2,4) = 1; itable(3,4) = 4

    call monolis_hash_init(hash_tree)

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)

      do i = 1, nsurf
!        i1 = conn(itable(1,i))
!        i2 = conn(itable(2,i))
!        i3 = conn(itable(3,i))
        ckey = get_key(i1, i2, i3)

        is_exist = .false.
        call monolis_hash_get(hash_tree, ckey, in, is_exist)
        if(is_exist)then
          in = in + 1
        else
          in = 1
        endif
        call monolis_hash_push(hash_tree, ckey, in, is_pushed, is_exist)
      enddo
    enddo

    allocate(is_inner(nsurf,mesh%nelem), source = 0)

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)

      do i = 1, nsurf
!        i1 = conn(itable(1,i))
!        i2 = conn(itable(2,i))
!        i3 = conn(itable(3,i))
        ckey = get_key(i1, i2, i3)

        is_exist = .false.
        call monolis_hash_get(hash_tree, ckey, in, is_exist)
        if(.not. is_exist) stop "error: monolis_get_surf_node_tet"
        if(in == 2) is_inner(i,eid) = 1
      enddo
    enddo

    call monolis_hash_finalize(hash_tree)

    allocate(is_surf_node(mesh%nnode), source = 0)

    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, nbase_func, conn)
      do i = 1, nsurf
        if(is_inner(i,eid) == 1)then
          do j = 1, nsurf_node
            in = conn(itable(j,i))
            is_surf_node(in) = 1
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_get_surf_node
end module mod_monolis_dbc_all_util