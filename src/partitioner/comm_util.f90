module mod_monolis_comm_util
  use mod_monolis_util
  use mod_monolis_util_debug
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  implicit none

  type monolis_send_list
    integer(kint) :: nnode = 0
    integer(kint), allocatable :: domid(:)
    integer(kint), allocatable :: local_nid(:)
  end type monolis_send_list

contains

  subroutine get_overlap_domain(graph, n_domain)
    implicit none
    type(monolis_graph) :: graph
    integer(kint) :: nelem, n_domain
    integer(kint) :: i, j, in, id, maxid, minid, jS, jE, jn
    integer(kint), allocatable :: temp(:)

    call monolis_debug_header("get_overlap_domain")

    nelem = graph%nelem
    allocate(graph%elem_domid(nelem), source = 0)
    allocate(graph%elem_domid_uniq(nelem), source = 0)

    if(n_domain == 1)then
      graph%elem_domid = 1
      graph%elem_domid_uniq = 1
      return
    endif

    do i = 1, nelem
      jS = graph%ebase_func(i) + 1
      jE = graph%ebase_func(i+1)
      allocate(temp(jE-jS+1), source = 0)
      jn = 1
      do j = jS, jE
        in = graph%connectivity(j)
        temp(jn) = graph%node_domid_raw(in)
        jn = jn + 1
      enddo
      minid = minval(temp)
      maxid = maxval(temp)
      deallocate(temp)

      graph%elem_domid_uniq(i) = maxid
      graph%elem_domid(i) = maxid
      if(minid /= maxid) graph%elem_domid(i) = -1
    enddo

    do i = 1, nelem
      if(graph%elem_domid(i) == 0) stop "get_overlap_domain"
    enddo
  end subroutine get_overlap_domain

  subroutine get_nnode_and_nid_at_subdomain(graph, local_node, nid, avg)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node
    integer(kint) :: nnode, nid, in, j, jn, k, jS, jE
    integer(kint) :: count_in, count_out, avg
    logical, allocatable :: is_in(:)

    nnode = graph%nnode
    allocate(is_in(nnode), source = .false.)

    do j = 1, local_node%nelem
      in = local_node%eid(j)
      jS = graph%ebase_func(in) + 1
      jE = graph%ebase_func(in+1)
      do k = jS, jE
        jn = graph%connectivity(k)
        is_in(jn) = .true.
      enddo
    enddo

    count_in = 0
    count_out = 0
    do j = 1, nnode
      if(is_in(j))then
        if(graph%node_domid_raw(j) /= nid)then
          count_out = count_out + 1
        else
          count_in = count_in + 1
        endif
      endif
    enddo

    !> save to structure
    local_node%nnode     = count_in + count_out
    local_node%nnode_in  = count_in
    local_node%nnode_out = count_out

write(*,"(4i10)") nid, local_node%nnode, local_node%nnode_in, local_node%nnode_out

    avg = avg + local_node%nnode
    allocate(local_node%nid(local_node%nnode), source = 0)

    in = 1
    do j = 1, nnode
      if(is_in(j) .and. graph%node_domid_raw(j) == nid)then
        local_node%nid(in) = j
        in = in + 1
      endif
    enddo
    do j = 1, nnode
      if(is_in(j) .and. graph%node_domid_raw(j) /= nid)then
        local_node%nid(in) = j
        in = in + 1
      endif
    enddo

    call get_perm_node(local_node%nnode, local_node%nid, local_node%nid_perm)
  end subroutine get_nnode_and_nid_at_subdomain

  subroutine get_nnode_and_nid_at_subdomain_graph(graph, local_node, nid, avg)
    use mod_monolis_stdlib
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node
    integer(kint) :: nnode, nid, in, i, j, jn, k, jS, jE, id
    integer(kint) :: count_in, count_out, avg
    logical, allocatable :: is_in(:)

    nnode = graph%N
    allocate(is_in(nnode), source = .false.)

    do i = 1, nnode
      if(graph%node_domid_raw(i) == nid)then
        jS = graph%index(i) + 1
        jE = graph%index(i+1)
        do j = jS, jE
          in = graph%item(j)
          is_in(in) = .true.
        enddo
      endif
    enddo

    count_in = 0
    count_out = 0
    do j = 1, nnode
      if(is_in(j))then
        if(graph%node_domid_raw(j) /= nid)then
          count_out = count_out + 1
        else
          count_in = count_in + 1
        endif
      endif
    enddo

    !> save to structure
    local_node%nnode     = count_in + count_out
    local_node%nnode_in  = count_in
    local_node%nnode_out = count_out

write(*,"(4i10)") nid, local_node%nnode, local_node%nnode_in, local_node%nnode_out

    avg = avg + local_node%nnode
    allocate(local_node%nid(local_node%nnode), source = 0)

    in = 1
    do j = 1, nnode
      if(is_in(j) .and. graph%node_domid_raw(j) == nid)then
        local_node%nid(in) = j
        in = in + 1
      endif
    enddo
    do j = 1, nnode
      if(is_in(j) .and. graph%node_domid_raw(j) /= nid)then
        local_node%nid(in) = j
        in = in + 1
      endif
    enddo

    call get_perm_node(local_node%nnode, local_node%nid, local_node%nid_perm)

    !> reconstruct nodal graph
    allocate(local_node%index(local_node%nnode+1), source = 0)

    k = 0
    do i = 1, local_node%nnode
      in = local_node%nid(i)
      jS = graph%index(in) + 1
      jE = graph%index(in+1)
      jn = 0
      do j = jS, jE
        in = graph%item(j)
        if(is_in(in)) jn = jn + 1
      enddo
      k = k + 1
      local_node%index(k+1) = local_node%index(k) + jn
    enddo

    jn = local_node%index(local_node%nnode+1)
    allocate(local_node%item(jn), source = 0)

    k = 0
    do i = 1, local_node%nnode
      in = local_node%nid(i)
      jS = graph%index(in) + 1
      jE = graph%index(in+1)
      do j = jS, jE
        in = graph%item(j)
        if(is_in(in))then
          k = k + 1
          local_node%item(k) = in
        endif
      enddo
    enddo
  end subroutine get_nnode_and_nid_at_subdomain_graph

  subroutine get_perm_node(nnode, nid, perm)
    implicit none
    integer(kint) :: i, in, nenode
    integer(kint) :: imax, imin
    integer(kint) :: nnode, nid(:)
    integer(kint), allocatable :: perm(:)

    imax = maxval(nid)
    imin = minval(nid)
    allocate(perm(imin:imax))
    perm = -1

    in = 1
    do i = 1, nnode
      perm(nid(i)) = in
      in = in + 1
    enddo
  end subroutine get_perm_node

  subroutine get_nelem_and_eid_at_subdomain(graph, local_node, nid)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node
    integer(kint) :: nelem, nid, in, j, k, count_in, count_out, jS, jE, i
    logical, allocatable :: is_in(:), is_out(:)

    nelem = graph%nelem
    allocate(is_in(nelem), source = .false.)

    do i = 1, nelem
      jS = graph%ebase_func(i) + 1
      jE = graph%ebase_func(i+1)
      do j = jS, jE
        in = graph%connectivity(j)
        if(graph%node_domid_raw(in) == nid)then
          is_in(i) = .true.
        endif
      enddo
    enddo

    in = 0
    do j = 1, nelem
      if(is_in(j)) in = in + 1
    enddo

    local_node%nelem = in
    allocate(local_node%eid(in), source = 0)

    !> internal element (global unique element)
    in = 0
    count_in = 0
    do j = 1, nelem
      if(is_in(j) .and. graph%elem_domid_uniq(j) == nid)then
        in = in + 1
        count_in = count_in + 1
        local_node%eid(in) = j
      endif
    enddo

    !> external element
    count_out = 0
    do j = 1, nelem
      if(is_in(j) .and. graph%elem_domid_uniq(j) /= nid)then
        in = in + 1
        count_out = count_out + 1
        local_node%eid(in) = j
      endif
    enddo

    local_node%nelem_in = count_in
    local_node%nelem_out = count_out
  end subroutine get_nelem_and_eid_at_subdomain

  subroutine get_neib_PE(local_node, graph, comm, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node(:)
    type(monolis_com) :: comm(:)
    integer(kint) :: i, in, j, jn, n_neib, n_domain
    integer(kint), allocatable :: domain(:)

    allocate(domain(n_domain), source = 0)

    do i = 1, n_domain
      domain = 0
      do j = 1, local_node(i)%nnode
        in = local_node(i)%nid(j)
        jn = graph%node_domid_raw(in)
        if(jn /= i)then
          domain(jn) = domain(jn) + 1
        endif
      enddo

      n_neib = 0
      do j = 1, n_domain
        if(domain(j) /= 0) n_neib = n_neib + 1
      enddo
      comm(i)%recv_n_neib = n_neib
      if(n_neib == 0)then
        allocate(comm(i)%recv_neib_pe(1), source = 0)
        allocate(comm(i)%recv_index(0:1), source = 0)
      else
        allocate(comm(i)%recv_neib_pe(n_neib))
        allocate(comm(i)%recv_index(0:n_neib), source = 0)
      endif
      comm(i)%recv_index = 0

      in = 0
      do j = 1, n_domain
        if(domain(j) /= 0)then
          in = in + 1
          comm(i)%recv_neib_pe(in) = j
          comm(i)%recv_index(in) = comm(i)%recv_index(in-1) + domain(j)
        endif
      enddo
    enddo
  end subroutine get_neib_PE

  subroutine get_recv_table(local_node, graph, comm, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node(:)
    type(monolis_com) :: comm(:)
    integer(kint) :: i, in, j, jn, n_recv, n_domain, nnode
    integer(kint), allocatable :: node_master_localid(:)

    nnode = graph%nnode
    allocate(node_master_localid(nnode), source = 0)

    do i = 1, n_domain
      do j = 1, local_node(i)%nnode_in
        in = local_node(i)%nid(j)
        if(graph%node_domid_raw(in) == i)then
          node_master_localid(in) = local_node(i)%nid_perm(in)
        endif
      enddo
    enddo

    !> get recv table
    do i = 1, n_domain
      n_recv = comm(i)%recv_index(comm(i)%recv_n_neib)
!>>> if(n_recv == 0)
      allocate(local_node(i)%recv_item_perm(n_recv), source = 0)
      allocate(local_node(i)%recv_item_domid(n_recv), source = 0)

      !> set recv node
      n_recv = 0
      do j = local_node(i)%nnode_in + 1, local_node(i)%nnode
        in = local_node(i)%nid(j)
        jn = graph%node_domid_raw(in)
        if(node_master_localid(in) == 0) stop "get_recv_table"

        if(jn /= i)then
          n_recv = n_recv + 1
          local_node(i)%recv_item_perm(n_recv) = node_master_localid(in)
          local_node(i)%recv_item_domid(n_recv) = jn
        endif
      enddo
      !> get recv item
      call reorder_recv_node(comm(i), local_node(i), n_recv)
    enddo
  end subroutine get_recv_table

  subroutine reorder_recv_node(comm, local_node, n_recv)
    use mod_monolis_stdlib
    implicit none
    type(monolis_com) :: comm
    type(monolis_node_list) :: local_node
    integer(kint) :: n_recv, i, j, jS, jE, nnode_in, nnode_out
    integer(kint), allocatable :: temp(:)

    nnode_in = local_node%nnode_in
    nnode_out = local_node%nnode_out

    allocate(comm%recv_item(nnode_out), source = 0)
    do i = 1, nnode_out
      comm%recv_item(i) = nnode_in + i
    enddo

    allocate(temp(n_recv), source = 0)
    temp = local_node%recv_item_domid
    call monolis_qsort_int_with_perm(temp, 1, n_recv, local_node%recv_item_perm)
    call monolis_qsort_int_with_perm(local_node%recv_item_domid, 1, n_recv, comm%recv_item)

    do i = 1, comm%recv_n_neib
      jS = comm%recv_index(i-1) + 1
      jE = comm%recv_index(i)
      n_recv = jE - jS + 1
      call monolis_qsort_int_with_perm(local_node%recv_item_perm(jS:jE), 1, n_recv, comm%recv_item(jS:jE))
    enddo
  end subroutine reorder_recv_node

  subroutine get_send_table(local_node, graph, comm, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list) :: local_node(:)
    type(monolis_com) :: comm(:)
    type(monolis_send_list), allocatable :: send(:)
    integer(kint) :: i, in, j, k, jn, jS, jE, n_domain, nnode, nelem, n_bound_node
    logical, allocatable :: is_bound(:)

    nnode = graph%nnode
    !> get is_bound
    allocate(is_bound(nnode), source = .false.)
    do i = 1, graph%nelem
      if(graph%elem_domid(i) == -1)then
        jS = graph%ebase_func(i) + 1
        jE = graph%ebase_func(i+1)
        do j = jS, jE
          in = graph%connectivity(j)
          is_bound(in) = .true.
        enddo
      endif
    enddo

    !> get local id of bound node
    n_bound_node = 0
    do i = 1, nnode
      if(is_bound(i)) n_bound_node = n_bound_node + 1
    enddo

    !> get send table
    allocate(send(n_domain))
    do i = 1, n_domain
      send%nnode = 0
    enddo

    do i = 1, n_domain
      do j = 1, comm(i)%recv_n_neib
        in = comm(i)%recv_neib_pe(j)
        jS = comm(i)%recv_index(j-1) + 1
        jE = comm(i)%recv_index(j)
        call append_send_node(send(in), i, jS, jE, local_node(i)%recv_item_perm(jS:jE))
      enddo
    enddo

    !> rebuild send table
    call rebuild_send_table(comm, send, n_domain)
  end subroutine get_send_table

  subroutine append_send_node(send_list, domid, jS, jE, recv_item_perm)
    use mod_monolis_stdlib
    implicit none
    type(monolis_send_list) :: send_list
    integer(kint) :: i, j, k, in, domid, jS, jE, nnode
    integer(kint) :: recv_item_perm(:)
    integer(kint), allocatable :: temp(:)

    nnode = jE - jS + 1
    allocate(temp(nnode), source = 0)
    temp = domid

    call monolis_reallocate_integer(send_list%domid    , send_list%nnode, nnode, temp)
    call monolis_reallocate_integer(send_list%local_nid, send_list%nnode, nnode, recv_item_perm)
    send_list%nnode = send_list%nnode + nnode
  end subroutine append_send_node

  subroutine rebuild_send_table(comm, send, n_domain)
    implicit none
    type(monolis_com) :: comm(:)
    type(monolis_send_list) :: send(:)
    integer(kint) :: i, in, j, n_domain, n_neib
    integer(kint), allocatable :: domain(:)

    allocate(domain(n_domain), source = 0)

    do i = 1, n_domain
      domain = 0
      do j = 1, send(i)%nnode
        in = send(i)%domid(j)
        domain(in) = domain(in) + 1
      enddo

      n_neib = 0
      do j = 1, n_domain
        if(domain(j) /= 0) n_neib = n_neib + 1
      enddo
      comm(i)%send_n_neib = n_neib

      if(n_neib == 0)then
        allocate(comm(i)%send_neib_pe(1), source = 0)
        allocate(comm(i)%send_index(0:1), source = 0)
      else
        allocate(comm(i)%send_neib_pe(n_neib), source = 0)
        allocate(comm(i)%send_index(0:n_neib), source = 0)
      endif
      comm(i)%send_index = 0

      in = 0
      do j = 1, n_domain
        if(domain(j) /= 0)then
          in = in + 1
          comm(i)%send_neib_pe(in) = j
          comm(i)%send_index(in) = comm(i)%send_index(in-1) + domain(j)
        endif
      enddo

      call reorder_send_node(comm(i), send(i))
    enddo
  end subroutine rebuild_send_table

  subroutine reorder_send_node(comm, send)
    use mod_monolis_stdlib
    implicit none
    type(monolis_com) :: comm
    type(monolis_send_list) :: send
    integer(kint) :: i, in, jS, jE, n_send

    n_send = send%nnode
    if(n_send > 0)then
      allocate(comm%send_item(n_send), source = 0)
      comm%send_item = send%local_nid
    else
      return
    endif

    call monolis_qsort_int_with_perm(send%domid, 1, n_send, comm%send_item)

    do i = 1, comm%send_n_neib
      jS = comm%send_index(i-1) + 1
      jE = comm%send_index(i)
      n_send = jE - jS + 1
      call monolis_qsort_int(comm%send_item(jS:jE), 1, n_send)
    enddo
  end subroutine reorder_send_node
end module mod_monolis_comm_util
