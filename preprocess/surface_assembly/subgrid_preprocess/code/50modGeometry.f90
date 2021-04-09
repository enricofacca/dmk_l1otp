!>-------------------------------------------------------------
!> defines objects containing geometry and mesh information
!>
!> objects:
!> type star --> TODO Description
!> type mesh --> TODO Description
!> 
!>
!> methods:
!> 
!<------------------------------------------------------------- 
! TODO type submesh as extension of type mesh
module GeometrySurface
  use Globals
  implicit none
  private
  !> structure variable containing the
  !> set of elements (nodes, edges) sharing (around) a node (element) 
  !>
  type, public:: star
     !> Number of elements in star
     integer :: nel
     !> indices of elements in star
     !> Dimension(nel)
     integer, allocatable :: el(:)
   contains
     !> static constructor
     !> (procedure public for type star)
     procedure, public, pass :: init => init_star
     !> static destructor
     !> (procedure public for type star)
     procedure, public, pass :: kill => kill_star
     !> Info procedure.
     !> (public procedure for type star)
     procedure, public, pass :: info => write_star
  end type star
  !> structure variable containing mesh and geometrical information 
  type, public :: mesh
     !> Number of mesh nodes.
     !> **Input in the constructor**.
     integer :: nnode = 0  
     !> Number of mesh triangles.
     !> **Input in the constructor**.
     integer :: ncell = 0 
     !> Number of nodes in each cell.
     !> **Input in the constructor**.
     integer  :: nnodeincell = 3 
     !> Characteristic dimension of the grid
     !> **Input in the constructor**.
     integer :: ndim=3 
     !> Number of mesh edges (including boundary). 
     !> Calculated in the constructor.
     integer :: nedge = 0 
     !> Number of boundary edges.
     !> **Calculated in the constructor.**
     integer :: nedge_bc = 0 
     !> Number of boundary edges.
     !> **Calculated in the constructor.**
     integer :: nnode_bc = 0 
     !> level of mesh topological quantities optionally calculated
     !> (true/false)
     !> (default =0, =1 first level cnc, =2 second level cnc)
     !>
     !> Default level cnc vars:
     !> triang, coord, size_cell, bar_cell
     !> First level cnc vars:
     !> plist, side_cnc, iside, neigh, peri_cell,length_edge, bar_edge, normal
     !> Second level cnc vars:
     !> node_in_node, edge_in_node, tria_in_node, tria_in_tria
     integer :: cnc_built = 0
     !>-------------------------------------------------------------------------
     !> Topological info
     !> Dimension (4,ncell).
     !> **Input in the constructor.**
     integer, allocatable :: topol(:,:) 
     !> Dimension (3,ncell). 
     !> Edges (in the global enumeration) in each cell.
     !> **Calculated in the constructor.**
     integer, allocatable :: side_cnc(:,:)
     !> Dimension (3,ncell). 
     !> Neighboring triangles for each triangle.
     !> **Calculated in the constructor.**
     integer, allocatable :: neigh(:,:)
     !> Dimension (2,nedge). 
     !> Nodes for each edge. 
     !> **Calculated in the constructor.**
     integer, allocatable :: iside(:,:)
     !> Dimension (nnode_bc). 
     !> Nodes at boundary 
     !> **Calculated in the constructor.**
     integer, allocatable :: node_bc(:)
     !> Dimension (2,nedge). 
     !> East and West cell for each edge. 
     !> East and West are calculated with respect to the nodal ordering in CELL.
     !> **Calculated in the constructor.**
     integer, allocatable :: plist(:,:)
     !> Dimension (3,nnode). 
     !> Nodal coordinate. 
     !> **Input in the constructor.**
     real(kind=double), allocatable :: coord(:,:)
     !>-------------------------------------------------------------
     !> Normal to the triangles
     !> TODO: use extends of 2dgrid type ??
     !>-------------------------------------------------------------
     !> Label for the normal stored
     !> normal_storage='node' : normal_node is allocated and used
     !> normal_storage='cell' : normal_cell is allocated and used
     character(len=4) :: normal_storage
     !> Dimension (3,ncell). 
     !> Nodal coordinate. 
     !> **Input in the constructor.**
     real(kind=double), allocatable :: normal_cell(:,:)
     !> Dimension (3,nnode). 
     !> Nodal coordinate. 
     !> **Input in the constructor.**
     real(kind=double), allocatable :: normal_node(:,:)
     !>-------------------------------------------------------------
     !> Dimension (ncell).
     !> Area of the cell
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: size_cell(:)
     !> Dimension (ncell).
     !> Perimeter of the cell
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: peri_cell(:)
     !> Dimension (ncell).
     !> coordinates of center of mass of the cell
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: bar_cell(:,:)
     !> Dimension (nedge).
     !> length of the edge
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: length_edge(:)
     !> Dimension (3,nedge).
     !> normal vector to the edge (consistent with iside and plist)
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: normal(:,:)
     !> Dimension (3,nedge). 
     !> Coord of gravity center of edges
     !> **Computed in procedure Prop.**
     real(kind=double), allocatable :: bar_edge(:,:)
     !>-------------------------------------------------------------------------
     !> Optional array, not always allocated
     !>-------------------------------------------------------------------------
     !> Arrays allocated only if cnc2lev quantities needs to be calculated 
     !> Dimension (nnode)
     !> Star of the triangles for each node
     type(star), allocatable :: tria_in_node(:)
     !> Dimension (nnode)
     !> Star of the edges for each node
     type(star), allocatable :: edge_in_node(:)
     !> Dimension (nnode)
     !> Star of the nodes for each node
     type(star), allocatable :: node_in_node(:)
     !> Dimension (ncell)
     !> Star of the triangles touching  each triangles
     type(star), allocatable :: tria_in_tria(:)
     !> Dimension (nnode)
     !> permutation vector for nodes
     !> **Calculated in the constructor.**
     integer, allocatable :: perm_nodes(:)
     !>-------------------------------------------------------------------------
     !> Arrays allocated only for grid levels (using uniform refinements) 
     !> (grid level 0 is the coarsest)
     !>-------------------------------------------------------------------------
     !> Grid level 
     integer  :: grid_level=0
     !> Number of nodes of parent grid
     integer :: nnode_parent=0
     !> Number of triangles parent grid = ncell/4
     integer :: ncell_parent=0
     !> Number of edge of parent grid
     integer :: nedge_parent=0
     !> Parent nodes of each node (can be the same)
     !> Dimension (2,nnode)
     !> A node can be generated by a single parent
     !> (if it corresponds to a vertex node of the parent mesh)
     !> or by two parents (if it is an edge midpoint node of the
     !> parent mesh).
     !> E.G.: if (node_parent(1,inode) .eq. node_parent(2,inode) )
     !> inode has only one parent with index node_parent(1,inode)
     !> (is its direct son),
     !> otherwise inode is the mid point of the edge defined by
     !> node_parent(1,inode) and node_parent(2,inode)
     !> (this type of array and of reasoning is needed for nodes and
     !> edges) 
     integer, allocatable :: node_parent(:,:)
     !> Dimension (ncell)
     !> Parent triangles
     integer, allocatable :: cell_parent(:)
   contains
     !> static constructor
     !> (procedure public for type mesh)
     procedure, public, pass :: init => init_mesh
     !> static destructor
     !> (procedure public for type mesh)
     procedure, public, pass :: kill => kill_mesh
     !> Info procedure.
     !> (public procedure for type mesh)
     procedure, public, pass :: info => info_mesh
     !> Procedure to write mesh in read-ready format
     !> (public procedure for type mesh)
     procedure, public, pass :: write => write_mesh
     !> Reads mesh from files
     !> (private procedure for type mesh)
     procedure, private, pass :: read_mesh
     !> builds mesh refining a previous mesh 
     !> (private procedure for type mesh)
     procedure, private, pass :: refine => build_refined_triang_coord
     !> Build mesh connectivity data starting from triang
     !> (private procedure for type mesh, used in init)
     procedure, private, pass :: connection
     !> Calculate cell areas, edge lengths, normals, etc
     !> (private procedure for type mesh, used in init)
     procedure, private, pass :: cnc2lev
     !> Calculate cell areas, edge lengths, normals, etc
     !> (private procedure for type mesh, used in init)
     procedure, private, pass :: meshProp
          !> calculate Lebesgue p-norm of piecewise-constant array
     !> (public procedure for type mesh)
     procedure, public, pass :: normp
     !> calculate parameter of the spatial discretization
     !> (public procedure for type mesh)
     procedure, public, pass :: meshpar
     !> procedure to find the index of the node
     !> closer to a given point
     procedure, public, pass :: closer_node
     !> procedure to find the index of the node
     !> closer to a given point
     procedure, public, pass :: info_point
     !> Renumber node numbering for 
     !> (public procedure for type mesh)
     procedure, public, pass :: reorder
     !> Renumber node numbering for 
     !> (private procedure for type mesh)
     procedure, private, pass :: build_rcm_permutation
     !> calculate the projection along the edge outward unit normal
     !> (public procedure for type mesh)
     procedure, public, pass :: proj => normal_projection
     !>---------------------------------------------------------------
     !> Subgrid procedure
     !> calculate the average of real array defined over subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: avg
     !> calculate the average of real 2d-array defined over subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: avg_vec
     !> calculate the average of real array defined over subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: proj_subgrid
     !> Project an integrated values w.r.t to P1 basis
     !> functions to a subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: project_p1integrated
     !> calculate the average of real aray defined over subgrid
     !> (public procedure for type mesh)
     procedure, public, pass :: subnode2tria
     !> calculate ancestrors in series of refined meshes
     !> (public procedure for type mesh)
     procedure, public, pass :: heraldry
     !> Procedure to read subgrid cell_parent and node_parent
     !> (public procedure for type mesh)
     procedure, public, pass :: read_parent => read_grid_subgrid_relation
     !> Procedure to write subgrid cell_parent and node_parent
     !> (public procedure for type mesh)
     procedure, public, pass :: write_parent => write_grid_subgrid_relation
     !> Procedure to check grid subgrid dependence
     !> (public procedure for type mesh)
     procedure, public, nopass :: check_relations => check_grid_subgrid_relation
     !> Procedure to compute node2node_connection
     !> (public procedure for type mesh)
     procedure, public, pass :: nodenode_connection
  end type mesh
contains
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type star)
  !> Instantiate variable of type star
  !>
  !> usage:
  !>     call 'var'%init(lun_err,nel)
  !>
  !> where:
  !> \param[in] lun_err -> integer. IO unit for error msg
  !> \param[in] nel     -> integer. Number of element in star
  !<-------------------------------------------------------------
  subroutine init_star(this,lun_err,nel)
    use Globals
    implicit none
    class(star), intent(inout) :: this 
    integer,     intent(in   ) :: lun_err
    integer,     intent(in   ) :: nel
    !local 
    logical rc
    integer res

    this%nel= nel
    if (allocated(this%el)) then
       deallocate(this%el,stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_dealloc , 'init_star', &
            ' type star member el', res)
    end if
    allocate(this%el(nel),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_star', &
         '  type star member el',res)

  end subroutine init_star
  !>-------------------------------------------------------------
  !> Static destructor
  !> (procedure public for type star)
  !> Instantiate variable of type star
  !>
  !> usage:
  !>     call 'var'%init(lun_err,nel)
  !>
  !> where:
  !> \param[in] lun_err -> integer. IO unit for error msg
  !> \param[in] nel     -> integer. Number of element in star
  !<-------------------------------------------------------------
  subroutine kill_star(this,lun_err)
    use Globals
    implicit none
    class(star), intent(inout) :: this 
    integer,     intent(in   ) :: lun_err
    !local 
    logical rc
    integer res

    deallocate(this%el,stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill star', &
         '  type star member el',res)
  end subroutine kill_star
  !>-------------------------------------------------------------
  !> Info procedure
  !> (procedure public for type star)
  !> Instantiate variable of type star
  !>
  !> usage:
  !>     call 'var'%info(lun_out)
  !>
  !> where:
  !> \param[in] lun_out -> integer. Ouput unit msg
  !<-------------------------------------------------------------
  subroutine write_star(this,lun_out)
    use Globals
    implicit none
    class(star), intent(inout) :: this 
    integer,     intent(in   ) :: lun_out
    !local 
    integer i
    write(lun_out,*) 'Number of elements in star ',this%nel
    write(lun_out,'(10i8)') (this%el(i),i=1,this%nel)
  end subroutine write_star
  !>-------------------------------------------------------------
  !> Static constructor.
  !> (procedure public for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file if requested)
  !> variable of type mesh
  !>
  !> usage:
  !>     call 'var'%init(lun_err,&
  !>                flag_read_build, &
  !>                [file2read[=...], &
  !>                input_mesh[=...], &
  !>                flag_cnc[=...], &
  !>                flag_reorder[=...]
  !>                ])
  !>
  !> where:
  !> \param[in] lun_err                 -> integer. IO unit for error msg
  !> \param[in] flag_read_build         -> integer. flag=0 reads from file
  !>                                                flag=1 refine from input_mesh
  !>                                                flag=2 read refined mesh from file
  !> \param[in] (optional) grid2read    -> type(file). I/O file for mesh data
  !> \param[in] (optional) parent2read  -> type(file). I/O file for parent data
  !> \param[in] (optional) input_mesh   -> type(mesh). mesh to refine
  !>                                       (used only if flag_read_build=1)
  !> \param[in] (optional) flag_cnc     -> integer.
  !>                                       flag=0 build only triang, coord,
  !>                                              and related real properties
  !>                                       flag=1 in addition, build
  !>                                              level 1 connectivity vars, 
  !>                                              and related real properties
  !>                                       flag=2 in addition, build
  !>                                              level 2 connectivity vars, 
  !>                                              and related real properties
  !> \param[in] (optional) flag_reorder -> integer.
  !>                                       flag=0 build mesh only
  !>                                       flag=1 build mesh and do reorder
  !>                                              minimizing the bandwith
  !>                                              (Cathill-McKee) for P1 Galerkin FEM
  !<-------------------------------------------------------------
  subroutine init_mesh(this,lun_err,&
       flag_read_build, & ! next are optional arguments
       grid2read, parent2read, input_mesh,&
       flag_cnc, flag_reorder)
    use Globals
    implicit none
    !vars
    class(mesh),           intent(inout) :: this
    integer,               intent(in   ) :: lun_err
    integer,               intent(in   ) :: flag_read_build
    type(file),  optional, intent(in   ) :: grid2read, parent2read
    type(mesh),  optional, intent(in   ) :: input_mesh
    integer,     optional, intent(in   ) :: flag_cnc, flag_reorder

    !local 
    logical :: rc
    integer :: res
    integer, allocatable :: perm(:), inv_perm(:)
    
    ! Initialized triang and coord
    select case (flag_read_build)
    case (0)
       ! Read coord and triang from file
       call this%read_mesh(lun_err,grid2read)
    case (1)
       ! build uniformly refined mesh
       call this%refine(lun_err,input_mesh)
    case (2)
       ! read refined mesh+parent
       call this%read_mesh(lun_err,grid2read)
       call this%read_parent(lun_err,parent2read)
    case default
       write(lun_err,*) ' Wrong Flag in procedure init mesh'
       write(lun_err,*) ' Flag_read_build not equal to 0,1 or 2'
       stop
    end select

    if ( present(flag_reorder) .and. (flag_reorder .eq. 1)  ) then
       ! allocation of local working arrays
       allocate(perm(this%nnode),inv_perm(this%nnode),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err,err_alloc ,'init_mesh',&
            'work array perm, inv_perm',res)
       ! Compute node permutation 
       call this%build_rcm_permutation(lun_err,perm,inv_perm)
       ! Permute nodes
       call this%reorder(lun_err,perm,inv_perm)
       ! Create a copy 
       allocate(this%perm_nodes(this%nnode),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err,err_alloc,'init_mesh',&
            ' perm_nodes',res)
       this%perm_nodes = perm
       ! Free memory
       deallocate(perm,inv_perm,stat=res) 
       if(res .ne. 0) rc = IOerr(lun_err,err_dealloc,'init_mesh',&
            'work array perm',res) 
    end if

    ! Build connection 
    if ( present(flag_cnc) ) this%cnc_built=flag_cnc
    
    select case (this%cnc_built)
    case default
       call this%meshProp(lun_err)
    case (1)
       call this%connection(lun_err)
       call this%meshProp(lun_err)
    case (2)
       call this%connection(lun_err)
       call this%cnc2lev(lun_err)
       call this%meshProp(lun_err)
    end select

  end subroutine init_mesh
  !>-------------------------------------------------------------
  !> Part of the static constructor.
  !> (procedure private for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type mesh
  !>
  !> usage:
  !>     call this%read_mesh(lun_err, file2read)
  !>
  !> where:
  !> \param[in] lun_err   -> integer. IO unit for error msg
  !> \param[in] file2read -> type(file). I/O file information
  !<-------------------------------------------------------------
  subroutine read_mesh(this, lun_err, file2read)
    use Globals
    implicit none
    !vars
    class(mesh),    intent(out) :: this
    integer,        intent(in ) :: lun_err
    type(file),     intent(in ) :: file2read
    ! local vars
    logical :: rc
    integer :: res 
    integer :: u_number
    integer :: j, k,icell, ncoord,nnodeincell,nnormals
    integer :: n1,n2,n3
    character(len=256) :: rdwr,str,fname,first_line
    real(kind=double) :: vec12(3),vec13(3),cross_prod(3)
    real(kind=double) :: dnrm2



    u_number    = file2read%lun
    fname       = file2read%fn
    nnodeincell = 3

    read(u_number,*,iostat=res) this%nnode
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_mesh', &
         etb(fname) // ' vars nnode',res)

    read(u_number,*,iostat=res) this%ncell
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_mesh', &
         etb(fname) // ' vars ncell',res)

    allocate(this%coord(3,this%nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
         '  type mesh member coord (array)',res)
    allocate(this%topol(nnodeincell+1,this%ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
         '  type mesh member triangles (array)',res)

    ! read first line to check if the the third coordinate is present 
    read(u_number,'(a)',iostat=res) first_line
    ncoord=3
    read(first_line,*,iostat=res) (this%coord(k,1),k=1,ncoord)
    if (res .ne. 0) THEN
       this%coord(3,:) = zero
       ! read only first two column, the first is initialized to zero
       ncoord = 2
       read(first_line,*,iostat=res) (this%coord(k,1),k=1,ncoord)
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') 1
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_inp , 'read_mesh', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end if
    
    do j=2,this%nnode
       read(u_number,*,iostat=res) (this%coord(k,j),k=1,ncoord)
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') j
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_inp , 'read_mesh', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end do

    
    do j=1,this%ncell
       read(u_number,*,iostat=res) (this%topol(k,j), k=1,nnodeincell+1)
       if(res .ne. 0) then
          write(rdwr,'(i5)') j
          str=etb(rdwr)//'/'
          write(rdwr,'(i5)') k
          str=trim(str)//etb(rdwr)
          rc = IOerr(lun_err, err_inp , 'read_mesh', &
               trim(etb(fname)) // &
               ' type mesh member array triang at line/col '//trim(str),res)
       end if
    end do

    read(u_number,*,iostat=res) nnormals
    if ( res .eq. 0) then
       !
       ! read from file
       !
       if (nnormals .eq. this%ncell ) then
          !
          ! read normals to triangles
          !
          this%normal_storage = 'cell'
          allocate(this%normal_cell(3,this%ncell),stat=res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
               '  type mesh member normal_cell',res)
          
          do j=1,this%ncell
             read(u_number,*,iostat=res) (this%normal_cell(k,j),k=1,3)
             if(res .ne. 0) THEN
                write(rdwr,'(i5)') j
                str=etb(rdwr)//'/'
                rc = IOerr(lun_err, err_inp , 'read_mesh', &
                     trim(etb(fname)) //&
                     ' type mesh member array normal_cell at line '//&
                     trim(str),res)
                
             end if
          end do
       end if

       if (nnormals .eq. this%nnode ) then
          !
          ! read normals to nodes
          !
          this%normal_storage = 'node'
          allocate(this%normal_node(3,this%nnode),stat=res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
               '  type mesh member normal_node ',res)
          do j=1,this%nnode
             read(u_number,*,iostat=res) (this%normal_node(k,j),k=1,3)
             if(res .ne. 0) THEN
                write(rdwr,'(i5)') j
                str=etb(rdwr)//'/'
                rc = IOerr(lun_err, err_inp , 'read_mesh', &
                     trim(etb(fname)) //&
                     ' type mesh member array normal_node at line '//&
                     trim(str),res)
                
             end if
          end do
       end if
    else
       !
       ! normal to triangles (the difference w.r.t. reading from file
       ! is that the normal can be not oriented, if an orientation exists)
       !
       this%normal_storage = 'cell'
       allocate(this%normal_cell(3,this%ncell),stat=res)
          if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_mesh', &
               '  type mesh member normal_node ',res)
       do icell=1,this%ncell
          n1 = this%topol(1,icell)
          n2 = this%topol(2,icell)
          n3 = this%topol(3,icell)
          vec12 = this%coord(:,n2)-this%coord(:,n1)
          vec13 = this%coord(:,n3)-this%coord(:,n1)
          cross_prod = cross(vec12,vec13)
          this%normal_cell(:,icell) = cross_prod/dnrm2(3,cross_prod,1)
       end do
    end if

  end subroutine read_mesh
  !>-------------------------------------------------------------
  !> Procedure to read parent grid information from file
  !> (procedure public for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type mesh
  !>
  !> usage:
  !>     call 'var'%read_parent(lun_err, file2read)
  !>
  !> where:
  !> \param[in] lun_err    -> integer. IO unit for error msg
  !> \param[in] file2read  -> type(file). I/O file information
  !>
  !> Remark: the first number read in is this%grid_level. This refers
  !> to the grid level of the parent mesh contained in file2read.
  !<-------------------------------------------------------------
  subroutine read_grid_subgrid_relation(this, lun_err, file2read)
    use Globals
    implicit none
    !vars
    class(mesh),    intent(inout) :: this
    integer,        intent(in   ) :: lun_err
    type(file),     intent(in   ) :: file2read

    !local 
    logical :: rc
    integer :: res
    integer :: k, inode, icell, parent_level
    integer :: u_number
    character(len=256) :: rdwr,str,fname
    integer :: nnode_parent_read,ncell_parent_read
    integer :: subnnode_read,subncell_read

    u_number = file2read%lun
    fname = file2read%fn

    read(u_number,*,iostat=res) parent_level
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars nnode',res)
    this%grid_level = parent_level + 1
    ! read  array dimension and check that the mesh is consistent
    read(u_number,*,iostat=res) nnode_parent_read,  subnnode_read 
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars nnode_parent_read,  subnnode_read ',res)
    if ( subnnode_read .ne. this%nnode) then
       write(lun_err,*) ' Mismatch node number'
       write(lun_err,*) ' subgrid nnode = ', this%nnode
       write(lun_err,*) ' subnnode read  = ', subnnode_read 
       stop
    end if

    read(u_number,*,iostat=res) ncell_parent_read,  subncell_read 
    if(res .ne. 0) rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
         etb(fname) // ' vars ncell_parent_read, subncell_read',res)
    if ( subncell_read .ne. this%ncell) then
       write(lun_err,*) ' Mismatch node number'
       write(lun_err,*) ' subgrid ncell = ', this%ncell
       write(lun_err,*) ' subncell read = ', subncell_read 
       stop
    end if


    ! set dimension arrays
    this%nnode_parent = nnode_parent_read
    this%ncell_parent = ncell_parent_read

    ! set allocation
    allocate(this%cell_parent(this%ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_grid_subgrid_relation', &
         '  type mesh member cell_parent',res)

    allocate(this%node_parent(2,this%nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'read_grid_subgrid_relation', &
         '  type mesh member node_parent',res)

    do inode=1,this%nnode
       read(u_number,*,iostat=res) (this%node_parent(k,inode),k=1,2)
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') inode
          str=etb(rdwr)//'/'
          rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end do

    do icell=1,this%ncell
       read(u_number,*,iostat=res) this%cell_parent(icell)
       if(res .ne. 0) then
          write(rdwr,'(i5)') icell
          str=etb(rdwr)
          rc = IOerr(lun_err, err_inp , 'read_grid_subgrid_relation', &
               trim(etb(fname)) // &
               ' type mesh member array cell_parent at line/col '//trim(str),res)
       end if
    end do
  end subroutine read_grid_subgrid_relation
  !>-------------------------------------------------------------
  !> Procedure to write grid parent information to file
  !> (procedure public for type mesh)
  !> Instantiate (allocate if necessary)
  !> and initilize (by also reading from input file)
  !> variable of type mesh
  !>
  !> usage:
  !>     call 'var'%write_parent(lun_err, file2write)
  !>
  !> where:
  !> \param[in] lun_err    -> integer. IO unit for error msg
  !> \param[in] file2write -> type(file). I/O file information
  !>
  !> Remark: the first number written out is this%grid_level. This refers
  !> to the grid level of the parent mesh contained in file2write.
  !<-------------------------------------------------------------
  subroutine write_grid_subgrid_relation(this, lun_err, file2write)
    use Globals
    implicit none
    !vars
    class(mesh),    intent(in) :: this
    integer,        intent(in) :: lun_err
    type(file),     intent(in) :: file2write
    !local 

    integer :: k, inode, icell
    integer :: u_number
    character(len=256) :: rdwr,str,fname
    logical :: rc
    integer :: res


    u_number = file2write%lun
    fname = file2write%fn

    write(u_number,*,iostat=res) this%grid_level
    if(res .ne. 0) rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
         etb(fname) // ' vars grid_level',res)

    if (this%grid_level.eq.0) then
       rc=IOerr(lun_err, wrn_out , &
            'write_grid_subgrid_relation', &
            etb(fname) // ' grid_level 0: no parent info',res)
       return
    end if

    ! Write array dimensions
    write(u_number,*,iostat=res) this%nnode_parent, this%nnode
    if(res .ne. 0) rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
         etb(fname) // ' vars nnode',res)

    write(u_number,*,iostat=res) this%ncell_parent, this%ncell
    if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write_grid_subgrid_relation', &
         etb(fname) // ' vars ncell',res)

    ! Write node parent
    do inode=1,this%nnode
       write(u_number,*,iostat=res) (this%node_parent(k,inode),k=1,2)
       if(res .ne. 0) THEN
          write(rdwr,'(i5)') inode
          str=trim(adjustl(rdwr))//'/'
          rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
               trim(etb(fname)) //&
               ' type mesh member array coord at line '//&
               trim(str),res)
       end if
    end do

    ! Write tria parent
    do icell=1,this%ncell
       write(u_number,*,iostat=res) this%cell_parent(icell)
       if(res .ne. 0) then
          write(rdwr,'(i5)') icell
          str=trim(adjustl(rdwr))
          str=trim(str)//trim(adjustl(rdwr))
          rc = IOerr(lun_err, err_out , 'write_grid_subgrid_relation', &
               trim(etb(fname)) // &
               ' type mesh member array this%cell_parent at line/col '//trim(str),res)
       end if
    end do

  end subroutine write_grid_subgrid_relation
  !>-------------------------------------------------------------
  !> Procedure to check grid/subgrid relations
  !> (procedure public for type mesh)
  !> Return ierr not equal to zero if something is not right
  !> 
  !> usage:
  !>     ierr=check_relations(grid,subgrid,lun_err,flag)
  !>
  !> where:
  !> \param[in] grid     -> type(mesh). 
  !> \param[in] subgrid  -> type(mesh). 
  !> \param[in] flag     -> integer. flag=1 check grid_level
  !>                                 flag=2 check coordinates
  !> \param[out] ierr    -> integer. ierr=1 wrong level of subgrid
  !>                                 ierr=2 wrong cordinate
  !<-------------------------------------------------------------
  function check_grid_subgrid_relation(grid,subgrid,flag) result (ierr)
    use Globals
    implicit none
    ! External variables
    type(mesh), intent(in)  :: grid
    type(mesh), intent(in)  :: subgrid
    integer,     intent(in) :: flag
    integer :: ierr
    !local 
    integer :: inode,n1,n2
    real(kind=double) :: temp(3),dnrm2
      
    ierr= 0
    if ( flag .ge. 1 )  then
       if ( (grid%grid_level .ne. subgrid%grid_level -1) .or. &
            ( subgrid%ncell/4 .ne. grid%ncell) ) then
          ierr=1
       end if
       if ( ( ierr .eq. 0) .and. (flag .eq.  2) ) then
          do inode = 1, subgrid%nnode
             n1 = subgrid%node_parent(1,inode)
             n2 = subgrid%node_parent(2,inode)
             
             temp(:) = subgrid%coord(:,inode) - &
                  onehalf*(grid%coord(:,n1) +grid%coord(:,n2))
             if ( dnrm2(2,temp,1) .gt. 1.0d-13 ) then
                ierr=2
                exit
             end if
          end do
       end if
    end if
  end function check_grid_subgrid_relation
  !>-------------------------------------------------------------
  !> Part of the build_mesh procedure
  !> (procedure private for type mesh)
  !> Build the triang and coord array of a conformally refined
  !> grid from a given coarser grid.
  !> Builds also the parent relationships.
  !>
  !> usage:
  !>     call 'var'%refine(lun_err, input_mesh)
  !>
  !> where:
  !> \param[in ] lun_err      -> integer. IO unit for err. msg. 
  !> \param[in ] input_mesh -> type(mesh). mesh to be refined
  !>
  !<-------------------------------------------------------------
  subroutine build_refined_triang_coord(this, lun_err , input_mesh)
    use Globals
    implicit none
    ! External variables
    class(mesh), intent(out) :: this
    integer,     intent(in)  :: lun_err
    type(mesh),  intent(in)  :: input_mesh
    ! Local vars
    type(mesh) :: local_mesh 
    integer :: vert(3),edg(3),midnod(2),nod(3,2)
    integer :: iedge,inode,icell,ivert,izone
    integer :: res
    integer :: iedgeloc,jedgeloc,inodloc,i,ntri1
    logical :: rc
    
    local_mesh=input_mesh
    if (local_mesh%cnc_built .lt. 1) then
       call local_mesh%connection(lun_err)
       local_mesh%cnc_built = 1 
       call local_mesh%MeshProp(lun_err)
    end if

    ! assign subgrid dimension
    this%nnode        = local_mesh%nnode + local_mesh%nedge
    this%ncell        = 4 * local_mesh%ncell
    this%grid_level   = local_mesh%grid_level + 1
    this%nnode_parent = local_mesh%nnode
    this%ncell_parent = local_mesh%ncell

    ! allocate additional array for subgrid
    allocate (&
         this%topol(4,this%ncell),&
         this%coord(3,this%nnode),&
         this%cell_parent(this%ncell),&
         this%node_parent(2,this%nnode),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc,&
         'build_refined_triang_coord', &
         ' type mesh member triang, coord, cell_parent'//&
         'node_parent',res)

    ! assign coordinates
    this%coord(:,1:local_mesh%nnode)=local_mesh%coord
    do iedge=1,local_mesh%nedge
       this%coord(1:3,local_mesh%nnode+iedge)=local_mesh%bar_edge(1:3,iedge)
    end do

    ! subnode parent
    do inode = 1, local_mesh%nnode
       this%node_parent(1,inode) = inode  
       this%node_parent(2,inode) = inode
    end do
    do iedge = 1, local_mesh%nedge 
       inode = local_mesh%nnode + iedge
       this%node_parent(1,inode) = local_mesh%iside(1,iedge)
       this%node_parent(2,inode) = local_mesh%iside(2,iedge)
    end do

    ! build new triang
    ntri1 = 0
    do icell=1,local_mesh%ncell
       ! store vertices
       do ivert = 1,3
          vert(ivert) = local_mesh%topol(ivert,icell)
       end do
       izone = local_mesh%topol(4,icell)
       ! store edges in triangle
       do iedgeloc=1,3
          edg(iedgeloc) = local_mesh%side_cnc(iedgeloc,icell)
          ! store nodes in edges
          do inodloc = 1,2
             nod(iedgeloc,inodloc) = local_mesh%iside(inodloc,edg(iedgeloc))
          end do
       end do
       ! build triangles
       do ivert = 1,3
          inode = vert(ivert)
          !find edges including inode
          iedgeloc = 0
          do jedgeloc = 1,3
             if ((inode .eq. nod(jedgeloc,1)) .or. &
                  (inode .eq. nod(jedgeloc,2))) then 
                ! edg(jedgeloc) includes inod
                iedgeloc = iedgeloc + 1
                midnod(iedgeloc) = edg(jedgeloc)+local_mesh%nnode
             endif
          end do              
          ! jedgeloc
          ntri1 = ntri1 + 1
          this%topol(1,ntri1) = inode
          do i = 1,2
             this%topol(i+1,ntri1) = midnod(i) 
          end do
          this%topol(4,ntri1) = izone
          this%cell_parent(ntri1) =icell
       end do                 ! iedgeloc
       ! build inner triangle
       ntri1 = ntri1 + 1
       do inodloc = 1,3
          iedgeloc = inodloc
          this%topol(inodloc,ntri1) = edg(iedgeloc)+local_mesh%nnode
       end do
       this%topol(4,ntri1) = izone
       this%cell_parent(ntri1) = icell
    end do

    call local_mesh%kill(lun_err)
    
  end subroutine build_refined_triang_coord

  !>----------------------------------------------------------------------------
  !> Procedure to project to subgrid a quantity integrated w.r.t to p1
  !> basin function
  !> TODO EF move it to P1Galerkin
  !>----------------------------------------------------------------------------
  subroutine project_p1integrated(subgrid,&
       p1integrated_grid,p1integrated_subgrid)  
    use Globals
    implicit none
    class(mesh),       intent(in ) :: subgrid
    real(kind=double), intent(in ) :: p1integrated_grid(subgrid%nnode_parent)
    real(kind=double), intent(out) :: p1integrated_subgrid(subgrid%nnode)
    ! local
    integer :: inode_sub, n1, n2

    do inode_sub = 1,subgrid%nnode
       n1 = subgrid%node_parent(1,inode_sub)
       n2 = subgrid%node_parent(2,inode_sub)
       p1integrated_subgrid(inode_sub) = (p1integrated_grid(n1) + p1integrated_grid(n2)) / 8.0d0
    end do

  end subroutine project_p1integrated

  !>-------------------------------------------------------------
  !> Static destructor.
  !> (procedure public for type mesh)
  !> deallocate all arrays for a var of type mesh
  !>
  !> usage:
  !>     call 'var'%kill(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-----------------------------------------------------------
  subroutine kill_mesh(this, lun)
    implicit none
    ! vars
    class(mesh) :: this
    integer, intent(in) :: lun
    ! local vars
    integer :: res
    logical :: rc
    integer :: inode, icell

    select case (this%cnc_built)
    case default
       deallocate(this%topol,this%coord,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members triang, coord')
       
       deallocate(this%size_cell,this%bar_cell,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members size_cell, bar_cell')
    case (1)
       deallocate(this%topol,this%coord,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members triang, coord')
       
       deallocate(this%size_cell,this%bar_cell,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members size_cell, bar_cell')

       if (allocated(this%normal)) then
          deallocate(this%peri_cell,this%length_edge,this%normal,this%bar_edge,stat=res)
          if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
               'mesh members peri_cell, length_edge, normal, bar_edge')
       end if

       deallocate(this%plist,this%side_cnc,this%iside,this%neigh,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members plist,side_cnc,iside,neigh')
       
       deallocate(this%node_bc,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members node_bc')
    case (2)
       deallocate(this%topol,this%coord,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members triang, coord')
       
       deallocate(this%size_cell,this%bar_cell,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members size_cell, bar_cell')
       
       if (allocated(this%normal)) then
          deallocate(this%peri_cell,this%length_edge,this%normal,this%bar_edge,stat=res)
          if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
               'mesh members peri_cell, length_edge, normal, bar_edge')
       end if
       
       deallocate(this%plist,this%side_cnc,this%iside,this%neigh,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members plist,side_cnc,iside,neigh')

       deallocate(this%node_bc,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members node_bc')

       do inode = 1, this%nnode
          call this%node_in_node(inode)%kill(lun)
          call this%tria_in_node(inode)%kill(lun)
          call this%edge_in_node(inode)%kill(lun)
       end do
       deallocate(this%node_in_node,this%tria_in_node,this%edge_in_node,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members nodein_node, tria_in_node,edge_in_node')

       do icell = 1, this%ncell
          call this%tria_in_tria(inode)%kill(lun)
       end do
       deallocate(this%tria_in_tria,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members tria_in_tria')
    end select

    this%nnode = 0
    this%ncell = 0

    if (this%grid_level.gt.0) then
       deallocate(this%cell_parent,this%node_parent,stat=res)
       if (res.ne.0) rc=IOerr(lun, err_dealloc, 'kill_mesh', &
            'mesh members cell_parent,node_parent')
    end if
    this%grid_level=0
    
  end subroutine kill_mesh
  !>-------------------------------------------------------------
  !> Info procedure.
  !> (public procedure for type mesh)
  !> Prints content of a variable of type mesh
  !> 
  !> usage:
  !>     call 'var'%info(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !<-------------------------------------------------------------
  subroutine info_mesh(this, lun)
    use Globals
    implicit none
    ! vars
    class(mesh), intent(in) :: this
    integer :: lun
    ! local vars
    integer :: j, k
    
    write(lun,*) ' '
    write(lun,*) ' Info: grid structure definition:'

    write(lun,*) 'number of nodes       ',this%nnode
    write(lun,*) 'number of triangles   ',this%ncell
    write(lun,*) 'number of total edges ',this%nedge
    write(lun,*) 'number of bnd edges   ',this%nedge_bc
    write(lun,*) 'number of int edges   ',this%nedge-this%nedge_bc
    write(lun,*) 'number of bnd nodes   ',this%nnode_bc


    write(lun,*) ' '
    write(lun,*) 'node coordinates: (type mesh member coord)'
    write(lun,*) '    node   x           y           z'
    do j=1,this%nnode
       write(lun,'(i8,2x,3(1pe12.4))') j,(this%coord(k,j),k=1,2)
    end do

    write(lun,*) ' '
    write(lun,*) 'triangles: (type mesh memeber Triangle)'
    write(lun,*) '    tri    node 1    node 2    node 3   material'
    do j=1,this%ncell
       write(lun,'(i8,2x,i8,2x,i8,2x,i8,2x,i5)') j,(this%topol(k,j),k=1,4)
    end do

    write(lun,*) ' '    
    write(lun,*) 'element Area: (type mesh member area)'
    do j=1,this%ncell
       write(lun,'(i8,2x,(1pe12.4))') j,this%size_cell(j)
    end do

    if (this%cnc_built.ge.1) then
       write(lun,*) ' '    
       write(lun,*) 'Boundary edges: (type mesh member iside, bar_edge)'
       write(lun,*) '   edge    node 1    node 2  s1_G       s2_G'
       do j=1,this%nedge_bc
          write(lun,'(i8,2x,i8,2x,i8,2(1pe12.4))') j,(this%iside(k,j),k=1,2), &
               (this%bar_edge(k,j),k=1,2)
       end do

       write(lun,*) ' '    
       write(lun,*) 'Internal edges: (type mesh member iside, bar_edge)'
       write(lun,*) '   edge    node 1    node 2  s1_G       s2_G'
       do j=this%nedge_bc+1,this%nedge
          write(lun,'(i8,2x,i8,2x,i8,2(1pe12.4))') j,(this%iside(k,j),k=1,2), &
               (this%bar_edge(k,j),k=1,2)
       end do
       
       write(lun,*) ' '    
       write(lun,*) 'neighbors: (type mesh member neigh)'
       write(lun,*) '    tri     tri 1     tri 2     tri 3 '
       do j=1,this%ncell
          write(lun,'(i8,2x,i8,2x,i8,2x,i8)') j,(this%neigh(k,j),k=1,3)
       end do
       
       write(lun,*) ' '    
       write(lun,*) 'edge numbering: (type mesh member side_cnc)'
       write(lun,*) '    tri    edge 1    edge 2    edge 3 '
       do j=1,this%ncell
          write(lun,'(i8,2x,i8,2x,i8,2x,i8)') j,(this%side_cnc(k,j),k=1,3)
       end do
       
       write(lun,*) ' '    
       write(lun,*) 'west-east elements: (type mesh member plist)'
       write(lun,*) '   edge     tri 1     tri 2 '
       do j=1,this%nedge
          write(lun,'(i8,2x,i8,2x,i8)') j,(this%plist(k,j),k=1,2)
       end do

       write(lun,*) ' '    
       write(lun,*) 'element Perimeter: (type mesh member perim)'
       do j=1,this%ncell
          write(lun,'(i8,2x,(1pe12.4))') j,this%peri_cell(j)
       end do

       write(lun,*) ' '    
       write(lun,*) 'coordinates of cell baricenters: (type mesh members bar_cell)'
       write(lun,*) '   elem    s1c         s2c'
       do j=1,this%ncell
          write(lun,'(i8,2x,2(1pe12.4))') j,(this%bar_cell(k,j),k=1,2)
       end do
       
       write(lun,*) ' '    
       write(lun,*) 'coordinates length of edges: (type mesh member length_edge)'
       write(lun,*) '   edge    length'
       do j=1,this%nedge
          write(lun,'(i8,2x,(1pe12.4))') j,this%length_edge(j)
       end do
       write(lun,*) ' '    
       write(lun,*) 'components of edge normals: (type mesh member normal)'
       write(lun,*) '   edge    n(1)        n(2)'
       do j=1,this%nedge
          write(lun,'(i8,2x,2(1pe12.4))') j,this%normal(1,j),this%normal(2,j)
       end do
    end if
    
    if (this%grid_level.gt.0) then
       write(lun,*) ' '    
       write(lun,*) 'parent information'
       write(lun,*) ' '    
       write(lun,*) ' node_parent'    
       write(lun,*) '   node     par 1    par 2'
       do j=1,this%nnode
          write(lun,'(i8,2x,i8,2x,i8)') j,(this%node_parent(k,j),k=1,2)
       end do
       write(lun,*) ' '    
       write(lun,*) ' cell_parent'    
       write(lun,*) '   tria     parent'
       do j=1,this%ncell
          write(lun,'(i8,2x,i8)') j,this%cell_parent(j)
       end do
    end if
  end subroutine info_mesh

  !>--------------------------------------------------------------
  !> Procedure write the mesh in the format
  !> nnode
  !> nncell
  !> coord1,coord2
  !> triang1,triang2,triang3,triang4
  !>-----------------------------------------------------------
  subroutine write_mesh(this, lun_err, file2write)
    use Globals
    implicit none
    ! vars
    class(mesh), intent(in) :: this
    integer,     intent(in) :: lun_err
    type(file),  intent(in) :: file2write
    ! local vars
    logical :: rc
    integer :: res,lun_write
    integer :: inode, icell,j
    
    lun_write=file2write%lun
    write(lun_write,*,iostat=res) this%nnode
    if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
            'type mesh member nnode',res)

    write(lun_write,*,iostat=res) this%ncell
    if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
         'type mesh member ncell',res)

    do inode = 1,this%nnode
       write(lun_write,*,iostat=res) this%coord(:,inode)
       if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
            'type mesh meber coord',res)
    end do
    do icell = 1, this%ncell
       write(lun_write,*,iostat=res) ( this%topol(j,icell), j=1,4 )
       if(res .ne. 0) rc = IOerr(lun_err, err_out, 'write mesh', &
            'type mesh meber triang',res)
    end do
  end subroutine write_mesh

  !>-------------------------------------------------------------
  !> calculate the projection along the edge outward unit normal
  !> (public procedure for type mesh)
  !>
  !>
  !> usage:
  !>     call 'var'%proj(vec, proj)
  !>
  !> where:
  !> \param[in] vec -> real(double), dimension (2,nedge)
  !>                   input vector to be projected along the edge normal
  !> \param[out] proj-> real(double), dimension (nedge)
  !>                    array of projected velocity along the normal to the edge
  !>                    It is a scalar for each edge as it is assumed to
  !>                    be directed along the normal.
  !>
  !<-------------------------------------------------------------
  subroutine normal_projection(this, vec, proj)
    use Globals
    implicit none

    class(mesh), intent(in) :: this
    real(kind=double), intent(in) :: vec(:,:)
    real(kind=double), intent(out) :: proj(:)

    ! local vars
    integer :: iedge
    real(kind=double) :: locvec(3),locnorm(3)
    ! external functions (BLAS lvl 1)
    real(kind=double) :: ddot

    do iedge=1,this%nedge
       locvec = vec(:,iedge)
       locnorm = this%normal(:,iedge)
       proj(iedge)=ddot(3,locvec,1,locnorm,1)
    end do
  end subroutine normal_projection
  !>-------------------------------------------------------------
  !> Build level 1 mesh connectivity data starting from 
  !> previous mesh%topol mesh%coord
  !> (private procedure for type mesh, used in init)
  !> 
  !> usage:
  !>     call 'var'%connection(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !> part of the static constructor
  !> at the end, we permutate all the edge structure and we
  !> renumber the edges so that the first nedge_bc edges are at
  !> the boundary
  !<-------------------------------------------------------------
  !! we work with local copies, and do proper assignments
  !! at the beginning and at the end
  subroutine connection(this, lun)
    use Globals
    implicit none
    !  External variables 
    class(mesh), intent(inout) :: this
    integer, intent(in) :: lun
    ! Local variables
    !  work arrays
    integer, allocatable :: triang(:,:)
    integer, allocatable :: neigh(:,:),side_cnc(:,:)
    integer, allocatable :: iside(:,:),plist(:,:)
    integer, allocatable :: perm(:),perminv(:)
    integer, allocatable :: mark_node(:)
    !  scalars and static arrays
    integer :: nnode,nedgemax,nedge,nedge_bc,ncell
    integer :: res
    integer :: icell,jtria,m
    integer :: kedg, iedg, iloc,inode
    integer :: il1,il2,in(2),n_ind,k3,k2,node,k4,k1
    logical :: rc,bound_edge

    ! Assign local scalar variabls
    ncell = this%ncell 
    nedgemax=3*ncell+1
    nnode = this%nnode
    ! Allocate local and work arrays 
    allocate(triang(4,ncell), &
         side_cnc(3,ncell), &
         plist(2,nedgemax), &
         iside(2,nedgemax), &
         neigh(3,ncell),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection', &
         '  local arrays triang,side_cnc,plist,iside,neigh',res)

    ! Assign local copies of arrays
    triang=this%topol

    nedge = 0
    do k1 = 1,3
       do  icell = 1,ncell
          side_cnc(k1,icell) = 0
       end do
    end do
    !
    ! * global edge enumeration
    !
    do icell = 1,ncell
       do k1 = 1,3
          !
          ! side_cnc(k1,icell)=0 means that this edge has
          ! already been
          ! considered, in the cell opposite to icell
          !
          if(side_cnc(k1,icell).eq.0) then
             k2 = mod(k1,3)+1
             !
             ! 'local' egde k1 shares 'local' nodes k1 and k2
             !
             il1 = triang(k1,icell)
             il2 = triang(k2,icell)
             jtria = icell 
             bound_edge = .true.

             do while (jtria.lt.ncell.and.bound_edge)  
                jtria = jtria + 1
                n_ind = 0
                k3 = 0
                do while (k3.lt.3.and.k3-n_ind.lt.2.and.bound_edge)
                   k3 = k3 + 1
                   node = triang(k3,jtria)
                   if(node.eq.il1 .or. node.eq.il2) then
                      n_ind = n_ind+1
                      in(n_ind) = k3
                      !
                      ! n_ind can be equal to either 1 or 2.
                      ! In this latter case 
                      ! the edge is an internal one.
                      !
                      bound_edge = (n_ind.eq.1)
                   endif
                end do
             end do

             if (bound_edge) then
                !
                ! * the edge connecting node il1 with node il2
                !   is a boundary one. 
                ! * by convention, the normal is outward
                !
                nedge = nedge + 1

                side_cnc(k1,icell)    = nedge
                neigh(k1,icell)       = 0

                iside(1,nedge) = triang(k1,icell)
                iside(2,nedge) = triang(k2,icell)

                plist(1,nedge) = icell
                plist(2,nedge) = 0
             else       
                !
                ! * the edge connecting node il1 with node il2
                !   is shared by trinlges
                ! * icell and jtria. By assumption the normal
                !   is oriented from icell
                ! * to jtria.
                ! * The local index of the new edge is given
                !   by the local index
                ! * of the 'first' node  that is
                ! * nodes(1,2) -> edge 1,
                !   nodes(2,3) -> edge 2,
                !   nodes(3,1) -> edge 3
                !
                if(in(1).eq.1.and.in(2).eq.2) k4 = 1
                if(in(1).eq.2.and.in(2).eq.3) k4 = 2
                if(in(1).eq.1.and.in(2).eq.3) k4 = 3

                nedge = nedge + 1

                side_cnc(k1,icell)    = nedge
                neigh(k1,icell)       = jtria
                side_cnc(k4,jtria)    = nedge
                neigh(k4,jtria)       = icell

                iside(1,nedge) = triang(k1,icell)
                iside(2,nedge) = triang(k2,icell)

                plist(1,nedge) = icell
                plist(2,nedge) = jtria

             endif
          endif
       end do
    end do

    allocate(perm(nedge),perminv(nedge),stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection', &
         ' temporary array perm, perminv ',res)
    kedg=0
    do iedg=1,nedge
       if(plist(1,iedg)*plist(2,iedg).eq.0)then
          ! for the edges at the boundary
          kedg=kedg+1
          perm(kedg)=iedg
       end if
    end do
    nedge_bc=kedg
    do iedg=1,nedge
       if(plist(1,iedg)*plist(2,iedg).ne.0)then
          ! for the edges at the boundary
          kedg=kedg+1
          perm(kedg)=iedg
       end if
    end do
    do iedg=1,nedge
       !perm(iedg) = iedg
       perminv(perm(iedg))=iedg
    end do

    ! Allocate global arrays of the structure "mesh"
    ! if not yet allocated 
    !
    if ( .not. allocated(this%neigh) ) then
       allocate(this%neigh(3,ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection', &
            ' type mesh member array neigh',res)
    end if
    if ( .not. allocated(this%side_cnc) ) then
       allocate(this%side_cnc(3,ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection', &
            ' type mesh member array side_cnc',res)
    end if
    if ( .not. allocated(this%iside) ) then
       allocate(this%iside(2,nedge),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection', &
            ' type mesh member array iside',res)
    end if
    if ( .not. allocated(this%plist) ) then
       allocate(this%plist(2,nedge),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'connection', &
            ' type mesh member array plist',res)
    end if
    ! copy local scalars and arrays into output
    ! global scalars and arrays
    this%nedge = nedge
    this%nedge_bc = nedge_bc
    do iedg=1,nedge
       this%iside(:,iedg)=iside(:,perm(iedg))
       this%plist(:,iedg)=plist(:,perm(iedg))       
    end do
    do icell=1,ncell
       do iloc=1,3         
          this%side_cnc(iloc,icell)=perminv(side_cnc(iloc,icell))
       end do
    end do
    this%neigh=neigh

    
    ! boundary_nodes
    this%nnode_bc = this%nedge_bc
    allocate(&
         this%node_bc(this%nnode_bc),&
         mark_node(nnode),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection', &
         ' local arrays mark_node',res)
    mark_node = 0
    do iedg=1,this%nedge_bc
       mark_node(this%iside(1:2,iedg)) = 1
    end do
    m=0
    do inode=1,nnode
       if( mark_node(inode) .eq. 1 ) then
          m=m+1
          this%node_bc(m) = inode
       end if          
    end do

    ! dellocate local arrays
    deallocate(triang,side_cnc,plist,iside,neigh,perm,perminv,mark_node,stat=res)
    if(res .ne. 0) rc = IOerr(lun, err_dealloc, 'connection', &
         ' local arrays triang, side_cnc, iside, plist, neigh, perm, perminv',res)

  end subroutine connection
  !>-------------------------------------------------------------
  !> Build level 2 mesh connectivity data starting from 
  !> previous mesh%topol mesh%coord and data built in connection
  !> (private procedure for type mesh, used in init)
  !> 
  !> usage:
  !>     call 'var'%cnc2lev(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !> part of the static constructor
  !<-------------------------------------------------------------
  !! we work with local copies, and do proper assignments
  !! at the beginning and at the end
  subroutine cnc2lev(this, lun_err)
    use Globals
    implicit none
    !vars
    class(mesh),   intent(inout) :: this
    integer, intent(in) :: lun_err
    ! local vars
    logical :: rc, add
    integer :: res    
    integer :: max_conn,nel
    integer :: i,j,k,k1,m, icell, inode, nneigh,tria2add
    integer, allocatable :: ncell_in_node(:)
    integer, allocatable :: icell_in_node(:,:)
    integer, allocatable :: nedge_in_node(:)
    integer, allocatable :: iedge_in_node(:,:)
    integer, allocatable :: nnode_in_node(:)
    integer, allocatable :: inode_in_node(:,:)    

    max_conn=30
    allocate(&
         ncell_in_node(this%nnode),&
         icell_in_node(max_conn,this%nnode),&
         nedge_in_node(this%nnode),&
         iedge_in_node(max_conn,this%nnode),&
         nnode_in_node(this%nnode),&
         inode_in_node(max_conn,this%nnode),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc , 'cnc2lev', &
         ' local arrays ncell_in_tria etc.', res)

    call make_stars(this%nnode,this%ncell,this%nedge,&
         max_conn,&
         this%topol,&
         this%iside,&
         ncell_in_node,&
         icell_in_node,&
         nedge_in_node,&
         iedge_in_node,&
         nnode_in_node,&
         inode_in_node)

    if (.not. allocated(this%tria_in_node)) then
       allocate(this%tria_in_node(this%nnode),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc , 'cnc2lev', &
         'mesh member tria_in_node',res)
    end if
    if (.not. allocated(this%tria_in_node)) then
       allocate(this%edge_in_node(this%nnode),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc , 'cnc2lev', &
            'mesh member edge_in_node',res)
    end if
    if (.not. allocated(this%node_in_node)) then
       allocate(this%node_in_node(this%nnode),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc , 'cnc2lev', &
            'mesh member node_in_node',res)
    end if

    do inode = 1 ,this%nnode
       nel=ncell_in_node(inode)
       call this%tria_in_node(inode)%init(lun_err,nel)
       this%tria_in_node(inode)%el(1:nel) = icell_in_node(1:nel,inode)

       nel = nedge_in_node(inode)
       call this%edge_in_node(inode)%init(lun_err,nel)
       this%edge_in_node(inode)%el(1:nel) = iedge_in_node(1:nel,inode)

       nel = nnode_in_node(inode)
       call this%node_in_node(inode)%init(lun_err,nel)
       this%node_in_node(inode)%el(1:nel) = inode_in_node(1:nel,inode)
    end do

    if (.not. allocated(this%tria_in_tria)) then
       allocate(this%tria_in_tria(this%ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun_err, err_alloc , 'cnc2lev', &
            'mesh member tria_in_tria',res)
    end if

    do icell = 1 , this%ncell
       nneigh = 0
       do j = 1,3
          inode = this%topol(j,icell) 
          nneigh = nneigh + this%tria_in_node(inode)%nel - 1
       end do
       ! remove triangles in side_cnc counted 2 times
       ! and count it
       m = 0
       do j = 1,3
          if ( this%neigh(j, icell ) .ne. 0 ) then
             nneigh = nneigh - 1
             m = m + 1
          end if
       end do

       call this%tria_in_tria(icell)%init(lun_err, nneigh)
       ! add element neigh
       i=0
       do j = 1, 3
          if ( this%neigh(j, icell ) .ne. 0 ) then
             i = i +1
             this%tria_in_tria(icell)%el(i) = this%neigh(j,icell)
          end if
       end do

       ! add the element connected by a node
       k1 = 0
       do j = 1, 3
          inode = this%topol(j,icell) 
          do k = 1, this%tria_in_node(inode)%nel
             tria2add = this%tria_in_node(inode)%el(k)
             ! check if is the central or one of its neighs
             add = ( &
                  tria2add  .ne. icell               .and. &
                  tria2add  .ne. this%neigh(1,icell) .and. & 
                  tria2add  .ne. this%neigh(2,icell) .and. &
                  tria2add  .ne. this%neigh(3,icell) )
             if ( add ) then
                k1=k1+1
                this%tria_in_tria(icell)%el(m+k1) = tria2add
             end if
          end do
       end do
       call isort(this%tria_in_tria(icell)%nel,this%tria_in_tria(icell)%el)
    end do

    deallocate(&
         ncell_in_node,&
         icell_in_node,&
         nedge_in_node,&
         iedge_in_node,&
         nnode_in_node,&
         inode_in_node,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc , 'cnc2lev', &
         ' local var stars',res)

  contains
    subroutine make_stars(nnode,ncell,nedge,&
         max_connectivity,&
         triang,&
         iside,&
         nmb_tria_in_node,&
         ids_tria_in_node,&
         nmb_edge_in_node,&
         ids_edge_in_node,&
         nmb_node_in_node,&
         ids_node_in_node)
      implicit none
      !     input      
      integer :: nnode,ncell,nedge
      integer :: max_connectivity
      integer :: triang(4,ncell)
      integer :: iside(2,nedge)
      !     ouput
      integer :: nmb_tria_in_node(nnode)
      integer :: ids_tria_in_node(max_connectivity,nnode)
      integer :: nmb_edge_in_node(nnode)
      integer :: ids_edge_in_node(max_connectivity,nnode)
      integer :: nmb_node_in_node(nnode)
      integer :: ids_node_in_node(max_connectivity,nnode)
      !     local
      integer :: j,k,icell,iedge,inode,inode_loc,pos,other
      logical :: in_list
      !
      pos=0

      nmb_tria_in_node = 0
      nmb_edge_in_node = 0
      nmb_node_in_node = 0

      ids_tria_in_node = 0
      ids_edge_in_node = 0
      ids_node_in_node = 0

      ! triangles in node
      do icell=1,ncell
         do inode_loc=1,3
            inode=triang(inode_loc,icell)
            nmb_tria_in_node(inode)=nmb_tria_in_node(inode)+1
            pos=nmb_tria_in_node(inode)
            ids_tria_in_node(pos,inode)=icell
         end do
      end do
      ! edges in node
      do iedge = 1,nedge
         do j = 1,2
            inode = iside(j,iedge)
            if ( nmb_edge_in_node(inode) == 0) then 
               ids_edge_in_node(1,inode) = iedge
               nmb_edge_in_node(inode)   = nmb_edge_in_node(inode) + 1
            else
               in_list = .false.
               do k=1,nmb_edge_in_node(inode)
                  if(ids_edge_in_node(k,inode) == iedge ) then
                     in_list = .true.
                  end if
               end do
               if (.not. in_list ) then
                  nmb_edge_in_node(inode)  = nmb_edge_in_node(inode) + 1
                  pos= nmb_edge_in_node(inode)
                  ids_edge_in_node(pos,inode)   = iedge
               end if
            end if
         end do
      end do
      ! nodes in node
      do iedge=1,nedge
         do j=1,2
            inode = iside(j,iedge)
            other = iside(mod(j,2)+1,iedge)
            if (nmb_node_in_node(inode) == 0) then 
               ids_node_in_node(1,inode) = other
               nmb_node_in_node(inode)   = nmb_node_in_node(inode) + 1
            else
               in_list = .false.
               do k=1,nmb_node_in_node(inode)
                  if(ids_node_in_node(k,inode) == other ) then
                     in_list = .true.
                  end if
               end do
               if (.not. in_list ) then
                  nmb_node_in_node(inode)  = nmb_node_in_node(inode) + 1
                  pos= nmb_node_in_node(inode)
                  ids_node_in_node(pos,inode)   = other
               end if
            end if
         end do
      end do
      do inode=1,nnode
         call isort(nmb_tria_in_node(inode),&
              ids_tria_in_node(1:nmb_tria_in_node(inode),inode))
         call isort(nmb_edge_in_node(inode),&
              ids_edge_in_node(1:nmb_edge_in_node(inode),inode))
         call isort(nmb_node_in_node(inode),&
              ids_node_in_node(1:nmb_node_in_node(inode),inode))
      end do

    end subroutine make_stars
  end subroutine cnc2lev

  !>-------------------------------------------------------------
  !> Calculate cell areas, edge lengths, etc
  !> completing the instantiation of a var of type mesh
  !> (private procedure for type mesh, used in init)
  !> 
  !> usage:
  !>     call 'var'%meshProp(lun)
  !>
  !> where:
  !> \param[in] lun -> integer. I/O unit for error message output
  !>
  !> part of the static constructor
  !<-------------------------------------------------------------
  subroutine meshProp(this, lun)
    use Globals
    implicit none
    ! vars
    class(mesh), intent(inout) :: this
    integer, intent(in) :: lun    
    ! local variables
    logical :: rc
    integer :: res
    integer :: icell,iedg,iloc
    integer :: n1,n2,n3
    integer :: celL,celR
    real(kind=double) :: s1n1,s2n1,s3n1,s1n2,s2n2,s3n2,s1n3,s2n3,s3n3
    real(kind=double) :: areaL,areaR
    real(kind=double) :: vec12(3),vec13(3),vec1L(3),vec1R(3),cross_prod(3)
    real(kind=double) :: coordL(3), coordR(3)
    real(kind=double) :: dnrm2
    
    if (.not. allocated(this%size_cell)) then
       allocate(this%size_cell(this%ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
            '  type mesh member size_cell (array)',res)
    end if
    if(.not. allocated(this%bar_cell)) then
       allocate(this%bar_cell(3,this%ncell),stat=res)
       if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
            '  type mesh member bar_cell (array)',res)
    end if
    if (this%cnc_built .ge. 1) then
       if (.not. allocated(this%length_edge)) then
          allocate(this%length_edge(this%nedge),stat=res)
          if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
               '  type mesh member length_edge (array)',res)
       end if
       if (.not. allocated(this%peri_cell)) then
          allocate(this%peri_cell(this%ncell),stat=res)
          if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
               '  type mesh member peri_cell (array)',res)    
       end if
       if (.not. allocated(this%bar_edge)) then
          allocate(this%bar_edge(3,this%nedge),stat=res)
          if(res .ne. 0) rc = IOerr(lun, err_alloc, 'meshProp', &
               '  type mesh member bar_edge (array)',res)
       end if
    endif
    if (this%cnc_built .ge. 2) then
       ! nothing to be one here at this moment
    end if

    do icell=1,this%ncell
       n1 = this%topol(1,icell)
       n2 = this%topol(2,icell)
       n3 = this%topol(3,icell)
       s1n1 = this%coord(1,n1)
       s2n1 = this%coord(2,n1)
       s3n1 = this%coord(3,n1)
       s1n2 = this%coord(1,n2)
       s2n2 = this%coord(2,n2)
       s3n2 = this%coord(3,n2)
       s1n3 = this%coord(1,n3)
       s2n3 = this%coord(2,n3)
       s3n3 = this%coord(3,n3)
       ! bar_cell
       this%bar_cell(1,icell) = (s1n1+s1n2+s1n3)*onethird
       this%bar_cell(2,icell) = (s2n1+s2n2+s2n3)*onethird
       this%bar_cell(3,icell) = (s3n1+s3n2+s3n3)*onethird
       ! size_cell
       vec12 = (/s1n2-s1n1,s2n2-s2n1,s3n2-s3n1/)
       vec13 = (/s1n3-s1n1,s2n3-s2n1,s3n3-s3n1/)    
       cross_prod = cross(vec12,vec13)
       this%size_cell(icell) = onehalf*dnrm2(3,cross_prod,1)
    end do

    if (this%cnc_built .ge. 1) then

       do iedg=1,this%nedge
          n1 = this%iside(1,iedg)
          n2 = this%iside(2,iedg)
          s1n1 = this%coord(1,n1)
          s2n1 = this%coord(2,n1)
          s3n1 = this%coord(3,n1)
          s1n2 = this%coord(1,n2)
          s2n2 = this%coord(2,n2)
          s3n2 = this%coord(3,n2)
          ! length_edge
          this%length_edge(iedg) = sqrt((s1n2-s1n1)*(s1n2-s1n1)+&
               (s2n2-s2n1)*(s2n2-s2n1)+(s3n2-s3n1)*(s3n2-s3n1))
          ! bar_edge
          this%bar_edge(1,iedg) = (s1n1+s1n2)*onehalf
          this%bar_edge(2,iedg) = (s2n1+s2n2)*onehalf
          this%bar_edge(3,iedg) = (s3n1+s3n2)*onehalf
       end do
        
       do icell=1,this%ncell
          this%peri_cell(icell) = zero
          do iloc=1,3
             iedg = this%side_cnc(iloc,icell)
             this%peri_cell(icell) = this%peri_cell(icell) + this%length_edge(iedg)
          end do
       end do

    end if
  
  end subroutine meshProp

  !>------------------------------------------------------------------
  !> procedure for averaging in  coarser mesh a quantity defined on a 
  !> refined mesh.
  !>      $(data_avg(i)=sum_{i=1,4}data(subtriangles(i))/4$
  !> It works  only if triang_sub is the topology 
  !> of a conformally refined mesh
  !>
  !> usage:    call  avg(subgrid, data, data_avg)
  !>
  !> where:
  !> \param[in ] data       -> real(ncell_sub) data to be averaged
  !> \param[out] data_avg   -> real(ncell_sub) data averaged
  !>------------------------------------------------------------------
  subroutine avg(subgrid, grid,data, data_avg)
    use Globals
    implicit none
    class(mesh),       intent(inout) :: subgrid
    class(mesh),       intent(inout) :: grid
    real(kind=double), intent(in   ) :: data(subgrid%ncell)
    real(kind=double), intent(inout) :: data_avg(subgrid%ncell_parent)
    !local 
    integer :: icell,icell_sub
    
    data_avg=zero
    do icell_sub = 1, subgrid%ncell
       icell = subgrid%cell_parent(icell_sub)
       data_avg(icell) = data_avg(icell) + &
            data(icell_sub)* onefourth
    end do
    

  end subroutine avg

  !>------------------------------------------------------------------
  !> procedure for averaging in  coarser mesh a quantity defined on a 
  !> refined mesh.
  !>      $(data_avg(i)=sum_{i=1,4}data(subtriangles(i))/4$
  !> It works  only if triang_sub is the topology 
  !> of a conformally refined mesh
  !>
  !> usage:    call  avg_vec(subgrid, data, data_avg)
  !>
  !> where:
  !> \param[in ] ndata      -> integer. Number of data to be averaged
  !> \param[in ] data       -> real Dimension (ndata,ncell_sub) 
  !>                             data to be averaged
  !> \param[out] data_avg   -> real Dimension(ndata,ncell_sub) 
  !>                             data averaged
  !>------------------------------------------------------------------
  subroutine avg_vec(subgrid, grid,ndata,data, data_avg)
    use Globals
    implicit none
    class(mesh),       intent(inout) :: subgrid
    class(mesh),       intent(inout) :: grid
    integer,           intent(in   ) :: ndata
    real(kind=double), intent(in   ) :: data(ndata,subgrid%ncell)
    real(kind=double), intent(inout) :: data_avg(ndata,subgrid%ncell_parent)
    !local 
    integer :: icell,icell_sub
    
    data_avg=zero
    do icell_sub = 1, subgrid%ncell
       icell = subgrid%cell_parent(icell_sub)
       data_avg(:,icell) = data_avg(:,icell) + &
            data(:,icell_sub)* onefourth
    end do
    

  end subroutine avg_vec

  !>------------------------------------------------------------------
  !> procedure for projection over a rifened mesh vreal var. data
  !> (public procedure for type mesh)
  !>
  !> usage:    call  avg(subgrid,data,data_avg)
  !>
  !> where:
  !> \param[in ] data       -> real(ncell_sub) data to be averaged
  !> \param[out] data_avg   -> real(ncell_sub) data averaged
  !>------------------------------------------------------------------
  subroutine proj_subgrid(subgrid, data, data_prj)
    use Globals
    implicit none
    class(mesh),       intent(in   ) :: subgrid
    real(kind=double), intent(in   ) :: data(subgrid%ncell_parent)
    real(kind=double), intent(inout) :: data_prj(subgrid%ncell)
    !local 
    integer :: icell_sub

    do icell_sub = 1, subgrid%ncell 
       data_prj(icell_sub) = data(subgrid%cell_parent(icell_sub))
    end do

  end subroutine proj_subgrid

  !>------------------------------------------------------------------
  !> procedure for projection data defined of subgrid nodes
  !> to triangles of not refined grid
  !> (public procedure for type mesh
  !>------------------------------------------------------------------
  subroutine subnode2tria(subgrid,data,data_prj)
    implicit none
    class(mesh),        intent(in ) :: subgrid
    real(kind=double) , intent(in ) :: data(subgrid%nnode)
    real(kind=double) , intent(out) :: data_prj(subgrid%ncell/4)
    !local
    integer :: i,j, icell, inode, sub_ncell

    sub_ncell = subgrid%ncell

    data_prj = zero
    do i = 1, sub_ncell
       icell = subgrid%topol(4,i)
       do j=1,3
          inode = subgrid%topol(j,i)
          data_prj(icell) = data_prj(icell) + data(inode)/1.2d1 
       end do
    end do
  end subroutine subnode2tria


  !>------------------------------------------------------------------
  !> procedure 
  !> (public procedure for type mesh)
  !>
  !> usage:    call  
  !>
  !> where:
  !> \param[in ]   -> real(ncell_sub) data to be averaged
  !> \param[out]   -> real(ncell_sub) data averaged
  !>------------------------------------------------------------------
  subroutine heraldry(this, nend, ngrids, grids, tria_ancestors)
    use Globals
    implicit none
    class(mesh),  intent(in   ) :: this
    integer   ,   intent(in   ) :: nend, ngrids
    type(mesh),   intent(in   ) :: grids(ngrids)
    integer,      intent(out  ) :: tria_ancestors(grids(nend)%ncell)
    !local 
    integer :: icell, icell_father,grid_level, ncell_max
    
    ! max level of refined
    ncell_max = grids(nend)%ncell    
    do icell = 1, ncell_max
       tria_ancestors(icell) = icell
    end do

    grid_level = grids(nend)%grid_level
    do while ( grid_level .ne. this%grid_level )
       do icell = 1, ncell_max
          icell_father = & 
               grids(grid_level+1)%cell_parent(tria_ancestors(icell))
          tria_ancestors(icell) = icell_father 
       end do
       grid_level = grid_level - 1
    end do
  end subroutine heraldry

  !>------------------------------------------------------------------
  !> Evaluates the p-norm of functions piecewise constant on triangles
  !>------------------------------------------------------------------
  function normp(this,power,data) result(total)
    use Globals
    implicit none
    class(mesh),       intent(in) :: this
    real(kind=double), intent(in) :: power
    real(kind=double), intent(in) :: data(this%ncell)
    real(kind=double) :: total
    !local 
    integer :: icell
    real(kind=double) :: ddot

    if ( power .eq.  0 ) then
       total = maxval( abs( data ) )
    else if ( power .eq. 1 ) then
       total = ddot(this%ncell, abs(data) ,1,this%size_cell,1)
    else
       total = zero
       do icell = 1, this%ncell
          total = total + abs(data(icell)) ** power * this%size_cell(icell)
       end do
       total = total ** ( one / power )
    end if
  end function normp

  !>-----------------------------------------------------
  !> Find the index of the closer node to given point
  !>-----------------------------------------------------
  function closer_node(this,point) result(id_closernode)
    use Globals
    implicit none
    class(mesh), intent(in) :: this
    real(kind=double) ,intent(in) :: point(2)
    integer :: id_closernode
    !local 
    integer i
    real(kind=double) :: dist,dist_min
    real(kind=double) :: temp(2)
    dist_min=1.0d30
    do i = 1, this%nnode
       temp = point - this%coord(1:2,i)
       dist = sqrt(temp(1)**2+temp(2)**2)
       if ( dist < dist_min ) then
          id_closernode=i
          dist_min = dist
       end if
    end do
  end function closer_node


  !>-----------------------------------------------------
  !> Minimal mesh parameter
  !>-----------------------------------------------------
  function meshpar(this,flag) result(h)
    use Globals
    implicit none
    class(mesh), intent(in) :: this
    integer,     intent(in) :: flag 
    real(kind=double) :: h
    !local 
    if (flag .eq. 0) h=sqrt(2*sum(this%size_cell)/this%ncell)
    if (flag .eq. 1) h=minval(sqrt(this%size_cell)/2.0d0)
    if (flag .eq. 2) h=maxval(sqrt(this%size_cell)/2.0d0)
  end function meshpar
  !>-----------------------------------------------------
  !> Given a point gives the index of the closer
  !> point, the minimal_distance, the index of the supporting triangle
  !> (index equal = 0 if dist_min < small ) 
  !>-----------------------------------------------------
  subroutine info_point(this,point,&
       index_closer_node, index_closer_tria, dist_min)
    use Globals
    implicit none
    class(mesh),        intent(in ) :: this
    real(kind=double) , intent(in ) :: point(2)
    integer,            intent(out) :: index_closer_node
    integer,            intent(out) :: index_closer_tria
    real(kind=double) , intent(out) :: dist_min
    !local 
    integer j,inode,icell
    real(kind=double) :: dist
    real(kind=double) :: temp(2)
    dist_min=1.0d30
    do icell = 1, this%ncell
       do j=1,3
          inode = this%topol(j,icell)
          temp = point - this%coord(:,inode)
          dist = sqrt(temp(1)**2+temp(2)**2)
          if ( dist < dist_min ) then
             index_closer_node=inode
             index_closer_tria=icell
             dist_min = dist
             if (dist_min < small) index_closer_tria=0
          end if
       end do
    end do
  end subroutine info_point
  !>-------------------------------------------------------------------
  !> Procedure that reorders mesh nodes
  !> given a permutation vector.
  !> It works on members triang and coord
  !> and rebuilds all other member if necessary. 
  !>
  !> usage:
  !>     call 'var'%reorder(lun_err,perm_nodes,inv_perm_nodes)
  !>
  !> where:
  !> \param[in] lun_err        -> integer. I/O unit for error message output
  !> \param[in] perm_nodes     -> integer. dimension(nnode)
  !>                                       node permutation array
  !> \param[in] inv_perm_nodes -> integer. dimension(nnode)
  !>                                       inverse of perm_nodes
  !>------------------------------------------------------------------
  subroutine reorder(this,lun_err,perm_nodes,inv_perm_nodes)
    use Globals
    implicit none
    class(mesh), intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    integer,     intent(in   ) :: perm_nodes(this%nnode) 
    integer,     intent(in   ) :: inv_perm_nodes(this%nnode) 
    ! local 
    integer :: inode,icell,i
    
    !  Permute the nodes according to the permutation vector.
    call double_col_permute(3, this%nnode, perm_nodes, this%coord)
    !
    !  Permute the node indices in the triangle array.
    do icell = 1, this%ncell
       do i = 1, 3
          inode = this%topol(i,icell)
          this%topol(i,icell) = inv_perm_nodes ( inode )
       end do
    end do

    ! By UNANIMOUS DECISION we check in connection
    ! if arrays are allocated
    
    ! Create local working mesh
    select case (this%cnc_built)
    case (0)
       ! Rebuild real Properties
       call this%meshProp(lun_err)
    case (1)
       ! Rebuild 1st level connections
       call this%connection(lun_err)
       ! Rebuild real Properties
       call this%meshProp(lun_err)
    case (2)
       call this%connection(lun_err)
       ! Rebuild 2nd level connections
       call this%cnc2lev(lun_err) 
       ! Rebuild real Properties
       call this%meshProp(lun_err)
    end select

    ! Renumber also node_parent if this is a subgrid
    if ( this%grid_level .gt. 0 ) then
       call integer_col_permute(2,this%nnode,perm_nodes, &
            this%node_parent)
    end if

  end subroutine reorder
  !>-------------------------------------------------------------------
  !> Procedure to build the permutation of the nodes to
  !> minimize the matrix bandwidth for Galerkin P1 FEM
  !>
  !> usage:
  !>     call 'var'%build_rcm_permutation(lun_err,
  !>                                      perm_nodes,
  !>                                      inv_perm_nodes)
  !>
  !> where:
  !> \param[in] lun_err        -> integer. I/O unit for error message output
  !> \param[in] perm_nodes     -> integer. dimension(nnode)
  !>                                       node permutation array
  !> \param[in] inv_perm_nodes -> integer. dimension(nnode)
  !>                                       inverse of perm_nodes
  !>-------------------------------------------------------------------
  subroutine build_rcm_permutation(this,lun_err,perm_nodes,iperm_nodes)
    use Globals
    use SparseMatrix
    implicit none
    class(mesh),        intent(inout) :: this
    integer,            intent(in   ) :: lun_err
    integer, optional,  intent(inout) :: perm_nodes(this%nnode) 
    integer, optional,  intent(inout) :: iperm_nodes(this%nnode) 
    !local 
    logical :: rc
    integer :: res,i,j,inode
    type(spmat) :: adj_matrix
    
    perm_nodes=0
    iperm_nodes=0

    call this%nodenode_connection(lun_err,.False.,adj_matrix) 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ADJ bandwidth = ', adj_matrix%bandwidth()
    !
    !  Compute the RCM permutation.
    !
    call adj_matrix%genrcm ( lun_err, perm_nodes,iperm_nodes )
    !
    !write(*,*) iperm_nodes
    
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Permuted ADJ bandwidth = ', &
         adj_matrix%bandwidth(perm_nodes,iperm_nodes)
    write ( *, '(a)' ) ' '
    
   
    !
    !  Free memory.
    !
    call adj_matrix%kill(lun_err)   
  end subroutine build_rcm_permutation

  subroutine nodenode_connection(this, lun_err,include_diagonal,adj_matrix)
    use Globals
    use SparseMatrix
    implicit none
    class(mesh),  intent(inout) :: this
    integer,      intent(in   ) :: lun_err
    logical,      intent(in   ) :: include_diagonal
    type(spmat),  intent(inout) :: adj_matrix
    !local
    logical :: rc
    logical :: clean_memory
    integer :: res
    integer :: iedge,inode,i1,i2,nnode,nterm,j
    integer, allocatable :: iaw(:)
    integer, allocatable :: jaw(:)
    integer, allocatable :: count(:)
    
    if ( .not. allocated(this%iside))  then
       call this%connection(lun_err)
       clean_memory = .true.
    end if

    nnode = this%nnode
    if ( include_diagonal) then
       nterm = 2*this%nedge+this%nnode
    else
       nterm = 2*this%nedge
    end if

    allocate(iaw(nnode+1),jaw(nterm),count(nnode),stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_star', &
         'work array ia ja count',res)
    
    ! build ia pointer
    if ( include_diagonal) then
       count=1
    else
       count=0
    end if
    do iedge = 1,this%nedge
       i1=this%iside(1,iedge)
       i2=this%iside(2,iedge)
       count(i1) = count(i1)+1
       count(i2) = count(i2)+1
       iaw(this%iside(1,iedge)) = iaw(this%iside(1,iedge))+1 
       iaw(this%iside(2,iedge)) = iaw(this%iside(2,iedge))+1
    end do
    iaw(1) = 1 
    do inode = 1, nnode
       iaw(inode+1) =  iaw(inode) + count(inode) 
    end do
       
    if ( include_diagonal ) then
       ! add self connection
       do inode = 1, this%nnode
          jaw(iaw(inode)) = inode
       end do
       count = 1
    else
       count = 0
    end if
    
    do iedge = 1, this%nedge
       i1 = this%iside(1,iedge)
       i2 = this%iside(2,iedge)
       jaw( iaw(i1)+count(i1) ) = i2 
       count( i1 ) = count( i1 ) + 1 
       jaw( iaw(i2)+count(i2) ) = i1
       count( i2 ) = count( i2 ) + 1 
    end do

   
    call adj_matrix%init(lun_err,&
         nnode,nnode, nterm,&
         storage_system='csr',&
         is_symmetric=.true.)

    adj_matrix%is_sorted = .False.    
    adj_matrix%ia    = iaw
    adj_matrix%ja    = jaw
    adj_matrix%coeff = zero
   
    call adj_matrix%sort()


    ! add extra diagonal terms
    deallocate(iaw,jaw,count,stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_star', &
         'work array ia ja count',res)
    
    if (clean_memory) then
       ! clean array built in connection
       deallocate(this%plist,this%side_cnc,this%iside,this%neigh,stat=res)
       if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_mesh', &
            'mesh members plist,side_cnc,iside,neigh')

       deallocate(this%node_bc,stat=res)
       if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_mesh', &
            'mesh members node_bc')
    end if

  end subroutine nodenode_connection

!!$  subroutine nodenode_connection2(this, lun_err,include_diagonal,&
!!$       nterm,ia,ja)
!!$    use Globals
!!$    use SparseMatrix
!!$    implicit none
!!$    class(mesh),  intent(inout) :: this
!!$    integer,      intent(in   ) :: lun_err
!!$    integer,      intent(in   ) :: nterm_max
!!$    logical,      intent(in   ) :: include_diagonal
!!$    integer,      intent(inout) :: nterm
!!$    integer,      intent(inout) :: ia(this%nnode)
!!$    type(spmat),  intent(inout) :: adj_matrix(nterm_max)
!!$    !local
!!$    logical :: rc
!!$    logical :: clean_memory
!!$    integer :: res
!!$    integer :: iedge,inode,i1,i2,nnode,nterm,j
!!$    integer, allocatable :: iaw(:)
!!$    integer, allocatable :: jaw(:)
!!$    integer, allocatable :: count(:)
!!$    
!!$    if ( .not. allocated(this%iside))  then
!!$       call this%connection(lun_err)
!!$       clean_memory = .true.
!!$    end if
!!$
!!$    nnode = this%nnode
!!$    if ( include_diagonal) then
!!$       nterm = 2*this%nedge+this%nnode
!!$    else
!!$       nterm = 2*this%nedge
!!$    end if
!!$
!!$    allocate(count(nnode),stat=res)
!!$    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_star', &
!!$         'work array ia ja count',res)
!!$    
!!$    ! build ia pointer
!!$    if ( include_diagonal) then
!!$       count=1
!!$    else
!!$       count=0
!!$    end if
!!$    do iedge = 1,this%nedge
!!$       i1=this%iside(1,iedge)
!!$       i2=this%iside(2,iedge)
!!$       count(i1) = count(i1)+1
!!$       count(i2) = count(i2)+1
!!$       ia(this%iside(1,iedge)) = ia(this%iside(1,iedge))+1 
!!$       ia(this%iside(2,iedge)) = ia(this%iside(2,iedge))+1
!!$    end do
!!$    ia(1) = 1 
!!$    do inode = 1, nnode
!!$       ia(inode+1) =  ia(inode) + count(inode) 
!!$    end do
!!$       
!!$    if ( include_diagonal ) then
!!$       ! add self connection
!!$       do inode = 1, this%nnode
!!$          jaw(ia(inode)) = inode
!!$       end do
!!$       count = 1
!!$    else
!!$       count = 0
!!$    end if
!!$    
!!$    do iedge = 1, this%nedge
!!$       i1 = this%iside(1,iedge)
!!$       i2 = this%iside(2,iedge)
!!$       jaw( ia(i1)+count(i1) ) = i2 
!!$       count( i1 ) = count( i1 ) + 1 
!!$       jaw( ia(i2)+count(i2) ) = i1
!!$       count( i2 ) = count( i2 ) + 1 
!!$    end do
!!$
!!$    m=0
!!$    do i=1,nterm_max
!!$       if ( ja(i) .ne. 0 )  then
!!$          m=m+1
!!$          ja(m) = ja(i)
!!$       end if
!!$    end do
!!$
!!$    do i=1,nrow
!!$       start = ia%(i)
!!$       end   = ia%(i+1)-1
!!$       nel   = end - start +1 
!!$       call sort(nel, ja(start:end) 
!!$    
!!$
!!$    ! add extra diagonal terms
!!$    deallocate(iaw,jaw,count,stat=res)
!!$    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'init_star', &
!!$         'work array ia ja count',res)
!!$    
!!$    if (clean_memory) then
!!$       ! clean array built in connection
!!$       deallocate(this%plist,this%side_cnc,this%iside,this%neigh,stat=res)
!!$       if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_mesh', &
!!$            'mesh members plist,side_cnc,iside,neigh')
!!$
!!$       deallocate(this%node_bc,stat=res)
!!$       if (res.ne.0) rc=IOerr(lun_err, err_dealloc, 'kill_mesh', &
!!$            'mesh members node_bc')
!!$    end if
!!$
!!$  end subroutine nodenode_connection2


end module GeometrySurface

   
