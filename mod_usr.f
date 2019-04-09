module mod_usr
  use mod_mhd
  implicit none


double precision, dimension(0:1300) :: data_t, qn, qvx, qvy, qvz, qBx1, qBy, qBz, qTemp,qe, qp, qrho
double precision :: qdt, qlength_unit, qmass_unit
integer :: qtc

contains

  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr

    call set_coordinate_system('Cartesian_1.75D')
    call mhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    integer :: i, stat, l
    
    save_physical_boundary = .true.

    mhd_gamma=1.66d0 
    mhd_eta = 0
    qtc = 0
	
    qlength_unit =1.d3 !m
    qmass_unit = 1.d-12 !kg 1.d-12
    !Read in data
    open(unit = 10, file = 'data.csv')
    i = 0
    do i=0,1300
      read(10,*, iostat = stat) data_t(i),  qn(i), qvx(i), qvy(i), qvz(i), qBx1(i), qBy(i), qBz(i), qTemp(i)
    end do
    
    close(unit = 10)
    !n = n/10d0
    !print*, n(0:30)

    qrho = qn *1.d6 * (1.6726d-27/qmass_unit) * qlength_unit**3  !in /cc to start
    qvx = qvx/1.d3 * qlength_unit !in km to start
    qvy = qvy/1.d3 * qlength_unit
    qvz = qvz /1.d3 * qlength_unit
    qBx1 = qBx1 * 1.d-9  /dsqrt(1.257d-6 * qmass_unit/qlength_unit)!in nT to start
    qBy = qBy * 1.d-9 /dsqrt(1.257d-6 * qmass_unit/qlength_unit)
    qBz = qBz * 1.d-9  /dsqrt(1.257d-6 * qmass_unit/qlength_unit)
    qp = qTemp * qn *1.d6 * qlength_unit**3 * 1.38065d-23 / (qmass_unit * qlength_unit**2) !in Kelvin
    !p = 30000.0
    !e = p/(mhd_gamma-1) + half*n*(vx**2+vy**2+vz**2) + half*(Bx1**2+By**2+Bz**2) 

    print*, 1.d6 * (1.6726d-27/qmass_unit) * qlength_unit**3
    print*,  1.d-9  /dsqrt(1.257d-6 * qmass_unit/qlength_unit)
    print*, 1.d6 * qlength_unit**3 * 1.38065d-23 / (qmass_unit * qlength_unit**2)/(1.d6*(1.6726d-27/qmass_unit)*qlength_unit**3)
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    ! initialize one grid 
    integer, intent(in) :: ixGmin1,ixGmax1,ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)

    double precision:: rholeft,rhoright,slocx1
    double precision:: vleft(3),pleft,vright(3),pright
    double precision:: bx,byleft,bzleft,byright,bzright
    logical,save :: first=.true.



    w(ixGmin1:ixGmax1,rho_)     = qrho(0)
    w(ixGmin1:ixGmax1,mom(1))   = qvx(0)
    w(ixGmin1:ixGmax1,mom(2))   = qvy(0)
    w(ixGmin1:ixGmax1,mom(3))   = qvz(0)
    w(ixGmin1:ixGmax1,p_ )      = qp(0)
    w(ixGmin1:ixGmax1,mag(1) )  = qBx1(0)
    w(ixGmin1:ixGmax1,mag(2) )  = qBy(0)
    w(ixGmin1:ixGmax1,mag(3) )  = qBz(0)

    !print*, qrho(0)
    !print*, qvx(0)
    !print*, qvy(0)
    !print*, qvz(0)
    !print*, qBx1(0)
    !print*, qBy(0)
    !print*, qBz(0)
    !print*, qp(0)
	
    !print*, e(0)
    !print*, w(ixGmin1,p_)
    !w(ixGmin1:ixGmax1,rho_) = exp(-10.*(x(ixGmin1:ixGmax1,1)-0.5)**2)+0.2


    call mhd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
  


  end subroutine initonegrid_usr


  subroutine specialbound_usr(qt,ixGmin1,ixGmax1,ixBmin1,ixBmax1,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixGmin1,ixGmax1, ixBmin1,ixBmax1, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
    double precision :: rho(ixGmin1:ixGmax1)
    integer :: ix
    double precision :: progress    
    
    if (qt > data_t(qtc+1)) then
            qtc = qtc+1
            !print*, p(qtc)
            !print*, w(17, p_)
            !print*, w(18, p_)
            !print*, w(19, p_)
            !print*, ' '


    end if
    progress = ((qt - data_t(qtc))/(data_t(qtc+1)-data_t(qtc)))
    do ix=ixBmin1,ixBmax1
      w(ix,rho_)   = (qrho(qtc)+(qrho(qtc+1)-qrho(qtc))*progress)
      w(ix,mom(1)) = (qvx(qtc)+(qvx(qtc+1)-qvx(qtc))*progress)
      w(ix,mom(2)) = (qvy(qtc)+(qvy(qtc+1)-qvy(qtc))*progress)
      w(ix,mom(3)) = (qvz(qtc)+(qvz(qtc+1)-qvz(qtc))*progress)
      w(ix,mag(1)) = qBx1(0)!(qBx1(qtc)+(qBx1(qtc+1)- qBx1(qtc))*progress)
      w(ix,mag(2)) = (qBy(qtc)+(qBy(qtc+1)-qBy(qtc))*progress)
      w(ix,mag(3)) = (qBz(qtc)+(qBz(qtc+1)-qBz(qtc))*progress)
      w(ix, p_)    = (qp(qtc)+(qp(qtc+1)-qp(qtc))*progress)
    end do
    !print*, qt
    !print*, 0.5*sin(4*3.14159*qt)

    !print*, (qp(qtc)+(qp(qtc+1)-qp(qtc))*progress)

    call mhd_to_conserved(ixGmin1,ixGmax1,ixBmin1,ixBmax1,w,x)

    !print*, "HERE"
    !print*, unit_magneticfield
	
  end subroutine specialbound_usr



end module mod_usr
