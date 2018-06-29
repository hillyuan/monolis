module mod_monolis_hash
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  private
  public :: monolis_hash_init
  public :: monolis_hash_finalize
  public :: monolis_hash_get
  public :: monolis_hash_push

  integer(kind=kint), save :: monolis_current_hash_size = 1
  integer(kind=kint), parameter :: monolis_hash_size(22) = (/&
  &     1021,     2039,      4093,      8191,     16381, &
  &    32749,    65521,    131071,    262139,    524287, &
  &  1048573,  2097143,   4194301,   8388593,  16777213, &
  & 33554393, 67108859, 134217689, 268435399, 536870909, 1073741789, 2147483647/)

  type type_monolis_hash_list
    integer(kind=kint) :: hash = 0
    character :: key*27 = "                           "
    integer(kind=kint) :: val  = 0
  end type type_monolis_hash_list

  type type_monolis_hash_bin
    integer(kind=kint) :: n = 0
    type(type_monolis_hash_list), pointer :: list(:) => null()
  end type type_monolis_hash_bin

  type type_monolis_hash_tree
    integer(kind=kint) :: n_put = 0
    integer(kind=kint) :: tree_size = 0
    type(type_monolis_hash_bin), pointer :: bin(:) => null()
  end type type_monolis_hash_tree

  type(type_monolis_hash_tree), save :: monolis_hash_tree

contains

  subroutine monolis_hash_init()
    implicit none
    type(type_monolis_hash_bin), pointer :: bin(:)

    monolis_hash_tree%n_put = 0
    monolis_hash_tree%tree_size = monolis_hash_size(monolis_current_hash_size)
    allocate(bin(monolis_hash_tree%tree_size))
    monolis_hash_tree%bin => bin
    nullify(bin)
  end subroutine monolis_hash_init

  subroutine monolis_hash_finalize()
    implicit none
    integer(kind=kint) :: i
    type(type_monolis_hash_list), pointer :: list(:)
    do i = 1, monolis_hash_tree%tree_size
      if(0 < monolis_hash_tree%bin(i)%n)then
        list => monolis_hash_tree%bin(i)%list
        if(associated(list)) deallocate(list)
      endif
    enddo
    deallocate(monolis_hash_tree%bin)
    nullify(list)
  end subroutine monolis_hash_finalize

  subroutine monolis_hash_get(key, val, is_exist)
    implicit none
    integer(kind=kint) :: hash, val
    character :: key*27
    logical :: is_exist

    val = 0
    is_exist = .false.
    call monolis_hash_key(key, hash)
    call monolis_hash_list_get(key, hash, val, is_exist)
  end subroutine monolis_hash_get

  subroutine monolis_hash_push(key, val, is_pushed, is_exist)
    implicit none
    integer(kind=kint) :: hash, val
    character :: key*27
    logical :: is_exist, is_pushed

    is_exist  = .false.
    is_pushed = .false.
    if(0.75d0*dble(monolis_hash_size(monolis_current_hash_size)) < dble(monolis_hash_tree%n_put))then
      call monolis_hash_resize()
    endif

    call monolis_hash_key(key, hash)
    call monolis_hash_list_get(key, hash, val, is_exist)

    if(.not. is_exist)then
      call monolis_hash_list_push(key, hash, val)
      monolis_hash_tree%n_put = monolis_hash_tree%n_put + 1
      is_pushed = .true.
    endif
  end subroutine monolis_hash_push

  subroutine monolis_hash_resize()
    implicit none
    integer(kind=kint) :: i, j, n, hash, val
    integer(kind=kint) :: new_size_id, new_size, old_size
    type(type_monolis_hash_bin), pointer :: new_bin(:), old_bin(:), temp_bin
    type(type_monolis_hash_list), pointer :: list(:)
    character :: key*27

    old_size = monolis_hash_tree%tree_size
    if(monolis_current_hash_size < 22)then
      monolis_current_hash_size = monolis_current_hash_size + 1
    endif
    new_size_id = monolis_current_hash_size
    new_size = monolis_hash_size(new_size_id)
    monolis_hash_tree%tree_size = new_size

    allocate(new_bin(new_size))
    old_bin => monolis_hash_tree%bin
    monolis_hash_tree%bin => new_bin

    do i = 1, old_size
      temp_bin => monolis_hash_tree%bin(i)
      temp_bin%list => old_bin(i)%list
      do j = 1, old_bin(i)%n
        hash = temp_bin%list(j)%hash
        key  = temp_bin%list(j)%key
        val  = temp_bin%list(j)%val
        call monolis_hash_list_push(key, hash, val)
      enddo
    enddo

    do i = 1, old_size
      list => old_bin(i)%list
      if(associated(list)) deallocate(list)
    enddo
    deallocate(old_bin)
    nullify(temp_bin)
    nullify(old_bin)
    nullify(new_bin)
  end subroutine monolis_hash_resize

  subroutine monolis_hash_list_get(key, hash, val, is_exist)
    implicit none
    integer(kind=kint) :: n, i
    integer(kind=kint) :: index, hash, val
    character :: key*27
    logical :: is_exist

    is_exist = .false.
    call monolis_index_key(hash, index)
    n = monolis_hash_tree%bin(index)%n
    do i = 1, n
      if(monolis_hash_tree%bin(index)%list(i)%key == key)then
        val = monolis_hash_tree%bin(index)%list(i)%val
        is_exist = .true.
      endif
    enddo
  end subroutine monolis_hash_list_get

  subroutine monolis_hash_list_push(key, hash, val)
    implicit none
    integer(kind=kint) :: i, iold, inew
    integer(kind=kint) :: index, hash, val
    character :: key*27
    type(type_monolis_hash_list), pointer :: old_list(:), new_list(:)

    call monolis_index_key(hash, index)
    iold = monolis_hash_tree%bin(index)%n
    old_list => monolis_hash_tree%bin(index)%list

    inew = iold + 1
    allocate(new_list(inew))
    do i = 1, iold
      new_list(i)%hash = old_list(i)%hash
      new_list(i)%key  = old_list(i)%key
      new_list(i)%val  = old_list(i)%val
    enddo

    new_list(inew)%hash = hash
    new_list(inew)%key  = key
    new_list(inew)%val  = val

    monolis_hash_tree%bin(index)%n = inew
    monolis_hash_tree%bin(index)%list => new_list
    if(associated(old_list)) deallocate(old_list)
    nullify(old_list)
    nullify(new_list)
  end subroutine monolis_hash_list_push

  !> BJD2 hash function
  subroutine monolis_hash_key(key, hash)
    implicit none
    character :: key*27
    integer(kind=kint) :: hash, i, t
    hash = 5381
    do i = 1, 27
      t = mod(hash*33, 65536_4)
      hash = mod(t + iachar(key(i:i)), 65536_4)
    enddo
  end subroutine monolis_hash_key

  subroutine monolis_index_key(hash, index)
    implicit none
    integer(kind=kint) :: hash, index
    index = mod(hash, monolis_hash_tree%tree_size) + 1
  end subroutine monolis_index_key
end module mod_monolis_hash