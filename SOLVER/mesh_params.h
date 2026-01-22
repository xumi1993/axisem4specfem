! Proc   4: Header for mesh information to run static solver
! created by the mesher on 01/13/2026, at 15h 52min
 
!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::
!   Background model     :            external
!   Dominant period [s]  :   50.0000
!   Elements/wavelength  :    1.5000
!   Courant number       :    0.6000
!   Coarsening levels    :         3
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
 integer, parameter ::         npol =         4  !            polynomial order
 integer, parameter ::        nelem =      1272  !                   proc. els
 integer, parameter ::       npoint =     31800  !               proc. all pts
 integer, parameter ::    nel_solid =      1068  !             proc. solid els
 integer, parameter ::    nel_fluid =       204  !             proc. fluid els
 integer, parameter :: npoint_solid =     26700  !             proc. solid pts
 integer, parameter :: npoint_fluid =      5100  !             proc. fluid pts
 integer, parameter ::  nglob_fluid =      3554  !            proc. flocal pts
 integer, parameter ::     nel_bdry =        72  ! proc. solid-fluid bndry els
 integer, parameter ::        ndisc =         3  !   # disconts in bkgrd model
 integer, parameter ::   nproc_mesh =         4  !        number of processors
 integer, parameter :: lfbkgrdmodel =         8  !   length of bkgrdmodel name
 
!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::
!   Time step [s]        :    0.3337
!   Min(h/vp),dt/courant :    3.4259    2.2249
!   max(h/vs),T0/wvlngth :   32.5166   33.3333
!   Inner core r_min [km]:  928.0628
!   Max(h) r/ns(icb) [km]:   79.6525
!   Max(h) precalc.  [km]:   89.2368
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
