c===================================================================
c
c  VUMAT subroutine for dislocation density crystal plasticity
c
c-------------------------------------------------------------------
c This code has being modified to incorporate the effects of FCC crystal
c plasticity
c-------------------------------------------------------------------

      SUBROUTINE VUMAT(
     & NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME,
     & TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY,
     & STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRAdoLD, FIELdoLD,
     & STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW,
     & STRETCHNEW, DEFGRADNEW, FIELDNEW,
     & STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)

      INCLUDE 'VABA_PARAM.INC'

c      implicit double precision (a-h,o-z)

      parameter(zero = 0., max_loops = 10, num_slip_sys = 12)

      character*(*) CMNAME
      logical Converged, Improved

c-------------------------------------------------------------------
c  Dimension arrays passed into the VUMAT subroutine in ABAQUS
c-------------------------------------------------------------------

      dimension
     &  props(nprops), ! Material properties passed from ABAQUS
     &  DENSITY(NBLOCK), ! Density of the declared block
     &  COORDMP(NBLOCK,*),  ! Coordinates of Gauss pt. being evaluated
     &  STRAININC(NBLOCK,NDIR+NSHR),  ! Strain increment (used for energy calculations)
     &  RELSPININC(*), 
     &  TEMPOLD(*),   ! Old temperature
     &  FIELdoLD(NBLOCK,NFIELDV),
     &  FIELDNEW(NBLOCK,NFIELDV),
     &  STRESSOLD(NBLOCK,NDIR+NSHR),  ! Previous time step: Cauchy stress tensor stored in vector form
     &  STATEOLD(NBLOCK,NSTATEV), ! Previous time step: State Variables in vector form
     &  ENERINTERNOLD(NBLOCK), ! Old internal energy
     &  ENERINELASOLD(NBLOCK), ! Old inelastic energy
     &  TEMPNEW(*),   ! New temperature
     &  STRETCHOLD(NBLOCK,NDIR+NSHR), ! Old stretch tensor stored in vector form
     &  DEFGRADOLD(NBLOCK,NDIR+2*NSHR), ! Deformation gradient at the beginning of step
     &  STRETCHNEW(NBLOCK,NDIR+NSHR), 
     &  DEFGRADNEW(NBLOCK,NDIR+2*NSHR), ! Deformation gradient at the end of step
     &  STRESSNEW(NBLOCK,NDIR+NSHR), ! Current time step: Cauchy stress tensor stored in vector form
     &  STATENEW(NBLOCK,NSTATEV), ! Current time step: State Variables
     &  ENERINTERNNEW(NBLOCK), ! Current time step: Internal energy
     &  ENERINELASNEW(NBLOCK) ! Current time step: inelastic energy
 
      ! MOOSE
      DIMENSION stressM(3,3),dstranM(3,3),ddsddeM(3,3,3,3),Cel(3,3,3,3)

      DIMENSION AUX33(3,3),AUX3333(3,3,3,3),AUX6(6),AUX66(6,6)

c-------------------------------------------------------------------
c  Dimension other arrays used in this VUMAT
c-------------------------------------------------------------------

      dimension
     &  array1(3,3), ! Dummy array
     &  array2(3,3), ! Dummy array
     &  array3(num_slip_sys,num_slip_sys), ! A_matrix of Newton Raphson (df_a/dgamma_dot_b)
     &  array4(6,6), ! Dummy array used in Voigt notation
     &  array5(6,6), ! Inverse of array4()
     &  array6(3,3,3,3), ! 4th rank dummy array
     &  array7(3,3,3,3), ! Another 4th rank dummy array
     &  cross(3), ! Dummy variable
     &  s_t (num_slip_sys), ! Thermal slip resistance
     &  s_a (num_slip_sys), ! Athermal slip resistance
     &  C0(3,3,3,3), ! Local  4th rank elastic stiffness tensor
     &  C(3,3,3,3), ! Global 4th rank elastic stiffness tensor
     &  C_avg(3,3,3,3), ! Average over all grains
     &  del(3,3), ! Kronecker delta tensor
     &  dir_cos(3,3),! Updated direction cosines
     &  dir_cos0(3,3), ! original direction cosines
     &  E_el(3,3), ! Elastic Green stran tensor
     &  E_tot(3,3), ! Green stran tensor E_el + E_p
     &  DEFGRADNEW_M(3,3),  ! New deformation gradient
     &  DEFGRAdoLD_M(3,3),  ! Old deformation gradient
     &  dev_stress(3,3),  !Deviatoric stress tensor
     &  F0(3,3), ! F at beginning of sub increment
     &  F1(3,3), ! F at end of sub increment
     &  F_el(3,3), ! Elastic part of F
     &  F_el_inv(3,3), ! Inverse of elastic part of F
     &  F_dot(3,3), ! F-dot
     &  F_inv(3,3), ! Inverse of deformation gradient
     &  F_p_inv_0(3,3), ! Inverse of plastic part of F at beginning
     &  F_p_inv(3,3), ! Inverse of plastic part of F at end of step
     &  func(num_slip_sys), ! Function that is solved to get gamma_dot(k)
     &  rho_m0(num_slip_sys), ! Mobile density at beginning of step
     &  rho_m(num_slip_sys), ! Mobile density at end of step
     &  rho_i0(num_slip_sys), ! Immobile density at beginning of step
     &  rho_i(num_slip_sys), ! Immobile density at end of step
     &  gamma_dot(num_slip_sys),! Shear strain rate on system
     &  gamma_dot_g(num_slip_sys),! Dislocation glide rate
     &  gamma_try(num_slip_sys),! Candidate gamma_dots for line search
     &  gamma_try_g(num_slip_sys),! Candidate gamma_dot_gs for line search
     &  gradient(num_slip_sys),! Gradient of sum-sqares-error
     &  psi(3), ! Updated euler angles
     &  psi0(3), ! Original euler angles
     &  sig0(3,3),! Stress tensor at beginning of step
     &  sig(3,3),! Stress tensor
     &  sig_avg(3,3),! Rate of stress change for a grain
     &  Spk2(3,3),! 2nd Piola Kirkhoff stress
     &  tau(num_slip_sys), ! Resolved shear stress for slip system
     &  xL(3,3),! Eulerian velocity gradient
     &  xL_p(3,3),! Plastic vel grad in current configuration
     &  xs0(3,num_slip_sys),! Inter config slip directions in global coords
     &  xs(3,num_slip_sys),! Current config slip directions in global coords
     &  xm0(3,num_slip_sys),! Inter config plane normals in global coords
     &  xm(3,num_slip_sys),! Current config plane normals in global coords
     &  y(3,num_slip_sys), ! Miller indices of slip plane normals
     &  z(3,num_slip_sys), ! Miller indices of slip directions
     &  delta_E_p(3,3),         !increment of plastic stran tensor
     &  delta_gamma(num_slip_sys), !increment of shear strain on each slip system
     &  E_p(3,3),
     &  A (num_slip_sys,num_slip_sys), ! Interaction coefficients between slip systems for slip resistance
     &  H (num_slip_sys,num_slip_sys), ! Geometric interaction coefficients between slip systems for mean free path
     &  eta(num_slip_sys),! Non-Schmid tensor
     &  ddpdsig(3,3,3,3),! deriv of D_p wrt sig * dt
     &  ddsdde_4th(3,3,3,3),! 4th rank tangent stiffness tensor
     &  dsadgb (num_slip_sys,num_slip_sys), ! deriv of s wrt del_gamma_beta
     &  dstdgb (num_slip_sys,num_slip_sys), ! deriv of g wrt del_gamma_beta
     &  drhomdgb(num_slip_sys,num_slip_sys), ! deriv of rho_m wrt del_gamma_beta
     &  drhoidgb(num_slip_sys,num_slip_sys), ! deriv of rho_i wrt del_gamma_beta
     &  d_gamma_dot(num_slip_sys),! Gamma_dot step from Newton-Raphson
     &  dlpdsd(3,3,3,3),! deriv of xL_p wrt sigma_dot
     &  dTaudGd(num_slip_sys,num_slip_sys),! deriv of Tau wrt gamma_dot_beta
     &  temp1(num_slip_sys,num_slip_sys),
     &  temp2(num_slip_sys,num_slip_sys),
     &  temp1_inv(num_slip_sys,num_slip_sys),
     &  shr_g(num_slip_sys),
     &  d_disl(num_slip_sys),
     &  dldrhom(num_slip_sys,num_slip_sys),
     &  orient(3,2500),
     &  time(2) ! Step Time and Total Time
!     &  ddsdde(6,6) ! 6x6 stiffness tensor (to be used only if we are using implicit)

      Pi = 4.*atan(1.)
      !print *, 'Value of Pi: ', Pi

      ! These values come from UMAT
      ndi = NDIR
      nshr = NSHR

      isingular = 0
c      OPEN(11,FILE='GAMMA_DOT.OUT', STATUS='UNKNOWN',
c     &    POSITION='APPEND')
      ! Convert from ABAQUS to MOOSE type of variable
      time(1) = STEPTIME
      time(2) = TOTALTIME

      ! since ABAQUS does not pass KINC, we have to manually add this
      if((steptime.eq.0.0).or.(totaltime.eq.DT)) then
        KINC = 1
        !write(11, *) 'Increment: ', KINC
      else
        KINC = 2
        !write(11, *) 'Increment: ', KINC+1
      endif
c-------------------------------------------------------------------
c  Assign props() array to logical variable names
c-------------------------------------------------------------------

c     We are directly assigning C11, C12 and C44 in VUMAT
c     So this is not needed
      !E_1 = props(1)  ! Young's modulus
      !E_2 = props(2)
      !G_1 = props(3)  ! Shear modulus
      !G_2 = props(4)
      !xnu = props(5)  ! Poisson's ratio

      gammadot0g = props(6)
      p0 = props(7)   ! Dislcoation barrier strength
      rho_m_zero = props(8)  ! Initial mobile dislocation density
      rho_i_zero = props(9)  ! Initial immobile dislocation density
      d_disl_zero = props(10) ! Initial dislocation line length
      R_c = props(11)  ! Critical capture radius for dislocation annihilation
      q_ann = props(12)  ! Dislocation evolution annihilation constant
      q_dyn = props(13)   ! Immobile dislocation evolution dynamic recovery constant
      b = props(14)    ! Burgers vector
      q_gen = props(15)  ! Dislocation line generation constant
      drag_stress = props(16)  ! Peierls barrier resistance
      q_t = props(17) ! Taylor hardening parameter
      x_d = props(18) ! Mean free path constant (dislocation network)
      p = props(19) ! Shape parameter (dislocation glide)
      q = props(20) ! Shape parameter (dislocation glide)
      Alatent = props(21) ! Latent hardening coefficient
      enthalpy_const = props(22) ! Multiplication constant for enthalpy term
      tau0 = props(23) ! Athermal slip resistance constant
      hp_coeff = props(24) ! Hall-Petch coefficient
      grain_size = props(25) ! Grain size (for Hall-Petch term)

      euler_1 = props(26) ! phi_1
      euler_2 = props(27) ! phi
      euler_3 = props(28) ! phi_2

      ! Some additional constants
      B_k = 1.3806503e-23  ! Boltzmann constant
      freq = 1e13  ! Debye frequency
      omega = b**3  ! Atomic volume
      dtime = DT !timestep increment

      ! gammadot0g = gammadot0g*rho_m_zero*d_disl_zero
      ! refgdg = gammadot0g*rho_m_zero*d_disl_zero

      ! Hard coded for now
      xtemp = 300.0e0  ! TEMPOLD
      TOLERANCE = 5.0e-7  ! this is manual (~1e-3 x strain rate)

c     ANISOTROPIC elastic constants for ferrite
c     From Patra and Tome 2024
      C11 = 231.4e3
      C12 = 134.7e3
      C44 = 116.7e3
      G = 75.02e3

c     adjusted to change units of product with stress to SI units
      act_vol = (b**3)*1.0e-3

c     delF0 scaled to SI units (Ngan Zhang 1999)
      delF0 = enthalpy_const*(G*1.0e6)*(b**3)*1.0e-9

c-------------------------------------------------------------------
c  Initialize Kronecker delta tensor
c-------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          del(i,j) = 0.0
        end do
        del(i,i) = 1.0
      end do

c-------------------------------------------------------------------
c  Assign slip system normals and slip directions for an FCC.
c  These slip systems are paired according to the cross-slip direction
c-------------------------------------------------------------------

c     !   plane      dir
c     /  1, 1, 1, -1, 1, 0 /
c     / -1,-1, 1, -1, 1, 0 /

c     /  1, 1, 1,  0,-1, 1 /
c     / -1, 1, 1,  0,-1, 1 /

c     /  1, 1, 1,  1, 0,-1 /
c     /  1,-1, 1,  1, 0,-1 /

c     / -1, 1, 1, -1,-1, 0 /
c     /  1,-1, 1, -1,-1, 0 /

c     / -1, 1, 1,  1, 0, 1 /
c     / -1,-1, 1,  1, 0, 1 /

c     / -1,-1, 1,  0, 1, 1 /
c     /  1,-1, 1,  0, 1, 1 /

C.... system 1 (111)[-110]
      y(1,1)=1
      y(2,1)=1
      y(3,1)=1
      z(1,1)=-1
      z(2,1)=1
      z(3,1)=0
C.... system 2 (-1-11)[-110]
      y(1,2)=-1
      y(2,2)=-1
      y(3,2)=1
      z(1,2)=-1
      z(2,2)=1
      z(3,2)=0
C.... system 3 (111)[0-11]
      y(1,3)=1
      y(2,3)=1
      y(3,3)=1
      z(1,3)=0
      z(2,3)=-1
      z(3,3)=1
C.... system 4 (-111)[0-11]
      y(1,4)=-1
      y(2,4)=1
      y(3,4)=1
      z(1,4)=0
      z(2,4)=-1
      z(3,4)=1
C.... system 5 (111)[10-1]
      y(1,5)=1
      y(2,5)=1
      y(3,5)=1
      z(1,5)=1
      z(2,5)=0
      z(3,5)=-1
C.... system 6 (1-11)[10-1]
      y(1,6)=1
      y(2,6)=-1
      y(3,6)=1
      z(1,6)=1
      z(2,6)=0
      z(3,6)=-1
C.... system 7 (-111)[-1-10]
      y(1,7)=-1
      y(2,7)=1
      y(3,7)=1
      z(1,7)=-1
      z(2,7)=-1
      z(3,7)=0
C.... system 8 (1-11)[-1-10]
      y(1,8)=1
      y(2,8)=-1
      y(3,8)=1
      z(1,8)=-1
      z(2,8)=-1
      z(3,8)=0
C.... system 9 (-111)[101]
      y(1,9)=-1
      y(2,9)=1
      y(3,9)=1
      z(1,9)=1
      z(2,9)=0
      z(3,9)=1
C.... system 10 (-1-11)[101]
      y(1,10)=-1
      y(2,10)=-1
      y(3,10)=1
      z(1,10)=1
      z(2,10)=0
      z(3,10)=1
C.... system 11 (-1-11)[011]
      y(1,11)=-1
      y(2,11)=-1
      y(3,11)=1
      z(1,11)=0
      z(2,11)=1
      z(3,11)=1
C.... system 12 (1-11)[011]
      y(1,12)=1
      y(2,12)=-1
      y(3,12)=1
      z(1,12)=0
      z(2,12)=1
      z(3,12)=1
c... end of slip system assignment
c-------------------------------------------------------------------
c  Normalize the Miller indices to length one.
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
        call normalize_vector(y(1,i),y(2,i),y(3,i))
        call normalize_vector(z(1,i),z(2,i),z(3,i))
      end do

c-------------------------------------------------------------------
c  Initialize interaction coefficient matrix
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
        do j = 1,num_slip_sys
          A(i,j) = Alatent
          H(i,j) = 1.0
        end do
        A(i,i) = 1.0
        H(i,i) = 1.0
      end do

c-------------------------------------------------------------------
c   Start looping over all material points in the geometry
c-------------------------------------------------------------------

      do ia_block = 1, NBLOCK

c-------------------------------------------------------------------
c  Initialize internal variables for initial time step
c-------------------------------------------------------------------
!       ! ABAQUS
!       if (time(1) .eq. 0.0) then
!       ! MOOSE
      if (kinc .eq. 1) then
          !print*, kinc
!         checkpoint("Initializing ISVs, initial time step")

! ! MPS comment
! c-------------------------------------------------------------------
! c  Store element position in global array
! c-------------------------------------------------------------------
!         xpos(noel) = coords(1)
!         ypos(noel) = coords(2)

c-------------------------------------------------------------------
c  Check for normality of Miller indices
c-------------------------------------------------------------------
        do k = 1,num_slip_sys
          sum1 = 0.0
          do i = 1,3
            sum1 = sum1 + y(i,k)*z(i,k)
          end do
          if (abs(sum1) .gt. tolerance) then
            write(6,*) 'The Miller indices are WRONG!!!'
            write(6,*) 'on slip system # ',k
            STOP
          end if
        end do

c-------------------------------------------------------------------
c  Generate Euler angles.  1-3
c-------------------------------------------------------------------

        psi0(1) = euler_1 * Pi/180
        psi0(2) = euler_2 * Pi/180
        psi0(3) = euler_3 * Pi/180

! c-------------------------------------------------------------------
! c  Initialize F_el.    4-12 ! Not needed.
! c-------------------------------------------------------------------
!         do i = 1,3
!           do j = 1,3
!             F_el(i,j) = 0.0
!           end do
!           F_el(i,i) = 1.0
!         end do

c-------------------------------------------------------------------
c  Initialize F_p_inv_0.    13-21
c-------------------------------------------------------------------
        do i = 1,3
          do j = 1,3
            F_p_inv_0(i,j) = 0.0
          end do
          F_p_inv_0(i,i) = 1.0
        end do

c-------------------------------------------------------------------
c  Initialize E_p.  22-30
c-------------------------------------------------------------------
        do i = 1,3
          do j = 1,3
            E_p(i,j) = 0.
          end do
        end do

c-------------------------------------------------------------------
c  Initialize E_eff.  31
c-------------------------------------------------------------------
        E_eff = 0.

c-------------------------------------------------------------------
c  Initialize E_p_eff.   32
c-------------------------------------------------------------------
        E_p_eff = 0.

c-------------------------------------------------------------------
c  Initialize E_p_eff_cum.   33
c-------------------------------------------------------------------
        E_p_eff_cum = 0.

c-------------------------------------------------------------------
c  Initialize average values of defect densities and the
c  corresponding loop sizes.   34-35
c-------------------------------------------------------------------
        rho_m_avg0 = rho_m_zero
        rho_i_avg0 = rho_i_zero

c-------------------------------------------------------------------
c  Initialize average dislocation glide and climb rates.   36
c-------------------------------------------------------------------
        gamma_dot_g_avg0 = gammadot0g

c-------------------------------------------------------------------
c  Initialize mobile density for each slip system of
c  each grain.  37-48
c-------------------------------------------------------------------
        do n = 1,num_slip_sys
          rho_m0(n) = rho_m_zero
        end do

c-------------------------------------------------------------------
c  Initialize immobile density for each slip system of
c  each grain.  49-60
c-------------------------------------------------------------------
        do n = 1,num_slip_sys
          rho_i0(n) = rho_i_zero
        end do
      end if
c-------------------------------------------------------------------
c  End of initializations.  Read in internal variables.
c-------------------------------------------------------------------

!       ! ABAQUS
!       if(time(1) .ne. 0.0) then
!       ! MOOSE

      if (kinc .ne. 1) then
!         checkpoint("Initializing ISVs, nonzero time step")
        n = 0

c-------------------------------------------------------------------
c  Read in Euler Angles 1-3
c-------------------------------------------------------------------
        do i = 1,3
          n = n + 1
          psi0(i) = STATEOLD(ia_block, n)
        end do

c-------------------------------------------------------------------
c  Read the elastic part of F   4-12
c-------------------------------------------------------------------
        do i = 1,3
          do j = 1,3
            n = n + 1
            ! F_el(i,j) = STATEOLD(ia_block, n) ! Not needed. This is only for post-processing
          end do
        end do

c-------------------------------------------------------------------
c  Read inverse of the plastic part of F    13-21
c-------------------------------------------------------------------
        do i = 1,3
          do j = 1,3
            n = n + 1
            F_p_inv_0(i,j) = STATEOLD(ia_block, n)
          end do
        end do

c-------------------------------------------------------------------
c  Read E_p   22-30
c-------------------------------------------------------------------
        do i = 1,3
          do j = 1,3
            n = n + 1
            E_p(i,j) = STATEOLD(ia_block, n)
          end do
        end do

c-------------------------------------------------------------------
c  Read E_eff      31
c-------------------------------------------------------------------
        n=n+1
        E_eff = STATEOLD(ia_block, n)
c-------------------------------------------------------------------
c  Read E_p_eff      32
c-------------------------------------------------------------------
        n=n+1
        E_p_eff = STATEOLD(ia_block, n)
c-------------------------------------------------------------------
c  Read E_p_eff_cum   33
c-------------------------------------------------------------------
        n = n+1
        E_p_eff_cum = STATEOLD(ia_block, n)

c-------------------------------------------------------------------
c  Read average values of defect densities and loop sizes  34-35
c-------------------------------------------------------------------
        n=n+1
        rho_m_avg0 = STATEOLD(ia_block, n)

        n=n+1
        rho_i_avg0 = STATEOLD(ia_block, n)

c-------------------------------------------------------------------
c  Read average dislocation glide and climb rates  36
c-------------------------------------------------------------------

        n=n+1
        gamma_dot_g_avg0 = STATEOLD(ia_block, n)

c-------------------------------------------------------------------
c  Read mobile dislocation density values       37-48
c-------------------------------------------------------------------
        do i = 1,num_slip_sys
          n = n + 1
          rho_m0(i) = STATEOLD(ia_block, n)
          !print*, rho_m0(i), n, kinc
        end do

c-------------------------------------------------------------------
c  Read immobile dislocation density values       49-60
c-------------------------------------------------------------------
        do i = 1,num_slip_sys
          n = n + 1
          rho_i0(i) = STATEOLD(ia_block, n)
          !print*, STATEOLD(ia_block, n), n
        end do

c-------------------------------------------------------------------
c  End of initializations
c-------------------------------------------------------------------
      end if

c     n = 61, 62, 63 are the euler angles

c-------------------------------------------------------------------
c  Calculate direction cosines based on Euler angles
c  Only being calculated for the initial timestep
c  Imported from state variables for the remaining steps
c-------------------------------------------------------------------

      n = 63
      if (kinc .eq. 1) then
            s1 = sin(psi0(1))
            c1 = cos(psi0(1))
            s2 = sin(psi0(2))
            c2 = cos(psi0(2))
            s3 = sin(psi0(3))
            c3 = cos(psi0(3))

            ! Bunge notation
            ! From VPSC
            ! Also verified from: http://engineering.snu.ac.kr/lecture/texture&anisotropy/Texture%20&%20Anisotropy%2010(Representation).pdf
            ! 88-96
            dir_cos0(1,1) = c1*c3-s1*s3*c2
            dir_cos0(1,2) = s1*c3+c1*s3*c2
            dir_cos0(1,3) = s3*s2
            dir_cos0(2,1) = -c1*s3-s1*c3*c2
            dir_cos0(2,2) = -s1*s3+c1*c3*c2
            dir_cos0(2,3) = c3*s2
            dir_cos0(3,1) = s1*s2
            dir_cos0(3,2) = -c1*s2
            dir_cos0(3,3) = c2

            ! since its the initial timestep
            do ia = 1,3
                  do ja = 1,3
                        dir_cos(ia,ja) = dir_cos0(ia,ja)
                  end do !ja
            end do !ia
      end if

      if (kinc .ne. 1) then

            ! Original orientation matrix (64-72)
            do ia = 1,3
                  do ja=1,3
                        n = n + 1
                        dir_cos0(ia,ja) = STATEOLD(ia_block, n)
                  end do !ja
            end do !ia

            ! Updated orientation matrix (73-81)
            do ia = 1,3
                  do ja=1,3
                        n = n + 1
                        dir_cos(ia,ja) = STATEOLD(ia_block, n)
                  end do !ja
            end do !ia

      end if

c-------------------------------------------------------------------
c  Initialize ANISOTROPIC 4th rank elastic stiffness tensor
c-------------------------------------------------------------------
      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          C0(i,j,k,l) = C12*del(i,j)*del(k,l) +
     &     C44*(del(i,k)*del(j,l)+del(i,l)*del(k,j))
         end do
        end do
       end do
      end do
      C0(1,1,1,1) = C11
      C0(2,2,2,2) = C11
      C0(3,3,3,3) = C11

c--------------------------------------------------------------------
c  Initialize arrays for averaging over grains
c--------------------------------------------------------------------
      do i = 1,3
       do j = 1,3
         sig_avg(i,j) = 0.0
       end do
      end do

c====================================================================
c  Begin calculations over the element
c====================================================================

c--------------------------------------------------------------------
c  Rotate local anisotropic elasticity tensor to
c  global coordinates
c--------------------------------------------------------------------
      call rotate_4th(dir_cos,C0,C)

c--------------------------------------------------------------------
c  Convert Miller Indices to global coordinates
c--------------------------------------------------------------------
      do n = 1,num_slip_sys
        do i = 1,3
          xs0(i,n) = 0.0
          xm0(i,n) = 0.0
          do j = 1,3
            xs0(i,n) = xs0(i,n) + dir_cos0(i,j)*z(j,n)
            xm0(i,n) = xm0(i,n) + dir_cos0(i,j)*y(j,n)
          end do
        end do
      end do

!       ! MOOSE
!       if (time(2) .eq. 0.0) goto 999

c--------------------------------------------------------------------
c  Initialize number of subincrements.  Note that the
c  subincrement initially equals the total increment.  This
c  remains the case unless the process starts to diverge.
c--------------------------------------------------------------------
      N_incr       = 1
      N_incr_total = 1

c====================================================================
c  Top of Subincrement Time Step Integration Loop.
c  Calculate subincrement time step size.
c====================================================================

  100 dt_incr = dtime / N_incr_total
c      print*, dt_incr

c-------------------------------------------------------------------
c  Initialize mobile dislocation density
c-------------------------------------------------------------------
      do n = 1,num_slip_sys
        rho_m(n) = rho_m0(n)
      end do

c-------------------------------------------------------------------
c  Initialize immobile dislocation density
c-------------------------------------------------------------------
      do n = 1,num_slip_sys
        rho_i(n) = rho_i0(n)
      end do

c-------------------------------------------------------------------
c  Initialize average defect densities and loop sizes
c-------------------------------------------------------------------
      rho_m_avg = rho_m_avg0
      rho_i_avg = rho_i_avg0

c-------------------------------------------------------------------
c  Convert deformation gradient from vector form to matrix form
c-------------------------------------------------------------------

      ! non-symmetric tensor
      ! I have taken it from https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/sub/default.htm?startat=ch01s02asb04.html
      DEFGRAdoLD_M(1,1) = DEFGRAdoLD(ia_block,1)
      DEFGRAdoLD_M(2,2) = DEFGRAdoLD(ia_block,2)
      DEFGRAdoLD_M(3,3) = DEFGRAdoLD(ia_block,3)

      DEFGRAdoLD_M(1,2) = DEFGRAdoLD(ia_block,4)
      DEFGRAdoLD_M(2,3) = DEFGRAdoLD(ia_block,5)
      DEFGRAdoLD_M(3,1) = DEFGRAdoLD(ia_block,6)

      DEFGRAdoLD_M(2,1) = DEFGRAdoLD(ia_block,7)
      DEFGRAdoLD_M(3,2) = DEFGRAdoLD(ia_block,8)
      DEFGRAdoLD_M(1,3) = DEFGRAdoLD(ia_block,9)
      
      ! update deformation gradient
      DEFGRADNEW_M(1,1) = DEFGRADNEW(ia_block,1)
      DEFGRADNEW_M(2,2) = DEFGRADNEW(ia_block,2)
      DEFGRADNEW_M(3,3) = DEFGRADNEW(ia_block,3)

      DEFGRADNEW_M(1,2) = DEFGRADNEW(ia_block,4)
      DEFGRADNEW_M(2,3) = DEFGRADNEW(ia_block,5)
      DEFGRADNEW_M(3,1) = DEFGRADNEW(ia_block,6)

      DEFGRADNEW_M(2,1) = DEFGRADNEW(ia_block,7)
      DEFGRADNEW_M(3,2) = DEFGRADNEW(ia_block,8)
      DEFGRADNEW_M(1,3) = DEFGRADNEW(ia_block,9)      

c--------------------------------------------------------------------
c  Initialize deformation gradients for beginning and
c  end of subincrement.
c--------------------------------------------------------------------
      do i = 1,3
        do j = 1,3
        
          F0(i,j) = DEFGRAdoLD_M(i,j) + (DEFGRADNEW_M(i,j)
     &     - DEFGRAdoLD_M(i,j)) * (N_incr - 1)/N_incr_total

          F1(i,j) = DEFGRAdoLD_M(i,j) + (DEFGRADNEW_M(i,j)
     &     - DEFGRAdoLD_M(i,j))*N_incr/N_incr_total

        end do
      end do

c--------------------------------------------------------------------
c  Multiply F() by F_p_inv() to get F_el()
c--------------------------------------------------------------------
      call aa_dot_bb(3,F0,F_p_inv_0,F_el)
      if (determinant(F_el) .eq. 0.0e0) then
        F_el(1,1) = 1.0e0
        F_el(2,2) = 1.0e0
        F_el(3,3) = 1.0e0
      end if
      call inverse_3x3(F_el,F_el_inv,isingular)
      if (isingular .eq. 1) go to 500

c--------------------------------------------------------------------
c  Rotate xs0 and xm0 to current coordinates, called xs and xm
c--------------------------------------------------------------------
      do n = 1,num_slip_sys
        do i = 1,3
          xs(i,n) = 0.0
          xm(i,n) = 0.0
          do j = 1,3
            xs(i,n) = xs(i,n) + F_el(i,j)*xs0(j,n)
            xm(i,n) = xm(i,n) + xm0(j,n)*F_el_inv(j,i)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate elastic Green strain
c--------------------------------------------------------------------
      array1(:,:) = 0.0
      array2(:,:) = 0.0
      call transpose1(3,F_el,array1)
      call aa_dot_bb(3,array1,F_el,E_el)
      do i = 1,3
        E_el(i,i) = E_el(i,i) - 1
        do j = 1,3
          E_el(i,j) = E_el(i,j)/2
        end do
      end do

c--------------------------------------------------------------------
c  Multiply the anisotropic stiffness tensor by the Green strain
c  to get the 2nd Piola Kirkhhoff stress
c--------------------------------------------------------------------
      call aaaa_dot_dot_bb(3,C,E_el,Spk2)

c--------------------------------------------------------------------
c  Convert from PK2 stress to Cauchy stress
c--------------------------------------------------------------------
      det = 0.0
      det = determinant(F_el)
      call transpose1(3,F_el,array2)
      call aa_dot_bb(3,F_el,Spk2,array1)
      call aa_dot_bb(3,array1,array2,sig)
      do i = 1,3
        do j = 1,3
          sig(i,j) = sig(i,j)/det
        end do
      end do

c--------------------------------------------------------------------
c  Calculate resolved shear stress for each slip system
c--------------------------------------------------------------------
      do k = 1,num_slip_sys
        tau(k) = 0.0
        do j = 1,3
          do i = 1,3
            tau(k) = tau(k) + xs(i,k)*xm(j,k)*sig(i,j)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate athermal slip resistance for each slip system.
c--------------------------------------------------------------------
      do ia = 1,num_slip_sys
        sum1 = 0.0
        sum2 = 0.0
        do ib = 1,num_slip_sys
          sum1 = sum1 + p0 * A(ia,ib) * (rho_m(ib) + rho_i(ib))
        end do
        s_a(ia) = tau0 + hp_coeff/sqrt(grain_size) + q_t*G*b*sqrt(sum1)
      end do

c--------------------------------------------------------------------
c  Calculate reference shear stress for each slip system.
c--------------------------------------------------------------------
      do ia = 1,num_slip_sys
        s_t(ia) = drag_stress
      end do

c--------------------------------------------------------------------
c  Calculate 1st estimate of gamma_dot for each slip system.
c--------------------------------------------------------------------
!       checkpoint("Calculating 1st estimate of gamma_dot")

      do ia = 1,num_slip_sys
        sum1 = 0.0
        sum2 = 0.0
        do ib = 1,num_slip_sys
          sum1 = sum1 + H(ia,ib)*(rho_m(ib) + rho_i(ib))
        end do
        d_disl(ia) = 1.e0/sqrt(sum1)
      end do

      do k = 1,num_slip_sys
!         refgdg = gammadot0g*(rho_m(k)/rho_m_zero)*
!      &   (d_disl(k)/d_disl_zero)
        refgdg = 1.0e8 * rho_m(k) * b * d_disl(k)

        if (abs(tau(k)) .lt. s_a(k)) then

          gamma_dot_g(k) = refgdg * exp((-delF0/B_k/xtemp) *
     &    power(1.0 - power(abs(tau(k))/s_a(k), p),q)) *
     &    sgn(tau(k))
!          print *, 'Slip System: ', k, gamma_dot_g(k)

        else if (abs(tau(k)) .ge. s_a(k)) then

          gamma_dot_g(k) = refgdg * sgn(tau(k))

        end if

        gamma_dot(k) = gamma_dot_g(k)

      end do

c--------------------------------------------------------------------
c  Calculate d(Tau)/d(Gamma_dot)
c--------------------------------------------------------------------
!       checkpoint("Calculating d(Tau) / d(gamma_dot)")
      do ia = 1,num_slip_sys
       do ib = 1,num_slip_sys
        dTaudGd(ia,ib) = 0.0

        do i = 1,3
         do j = 1,3
          array1(i,j) = xs0(i,ib)*xm0(j,ib)
         end do
        end do
        call aaaa_dot_dot_bb(3,C0,array1,array2)

        do i = 1,3
         do j = 1,3
          array1(i,j) = xs0(i,ia)*xm0(j,ia)
         end do
        end do
        call aa_dot_dot_bb(3,array1,array2,dTaudGd(ia,ib))

        dTaudGd(ia,ib) = (-1)*dTaudGd(ia,ib)*dt_incr
       end do !ib
      end do !ia

c====================================================================
c  Begin Newton-Raphson iterative loops.
c====================================================================
      converged = .false.

      do while (.not.converged)

      converged = .true.

c--------------------------------------------------------------------
c  Calculate g, F_p_inv, F_el, sig, tau, func, sse.
c--------------------------------------------------------------------
!       checkpoint("Calling eval_func() based on initial estimate")
      call eval_func(xs0, xm0, dt_incr,
     &  gamma_dot, gamma_dot_g,
     &  F1, F_p_inv, F_p_inv_0, F_el, num_slip_sys, C,
     &  rho_m0, rho_m, rho_i0, rho_i,
     &  Alatent, p0, B_k, R_c, q_dyn,x_d, drag_stress, q_gen, q_ann,
     &  delF0, G, b, xtemp, sse, sig, tau,
     &  s_t, s_a, A, H, p, q, q_t, ! shr_g, shr_c,
     &  func, xL_p, gammadot0g, refgdg, rho_m_zero, d_disl_zero,
     &  tau0, hp_coeff, grain_size, isingular)
      if (isingular .eq. 1) go to 500

!       checkpoint("1st run of eval_func() complete")
      sse_ref = sse
c--------------------------------------------------------------------
c  Begin calculation of the partial derivatives needed for
c  the Newton-Raphson step!!!
c--------------------------------------------------------------------

!       checkpoint("Calculating derivatives")

c--------------------------------------------------------------------
c  Calculate derivative of the mobile dislocation density, rho_m and
c  immobile dislocation density, rho_i w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------
c     mobile dislocations

      do ia = 1,num_slip_sys
        do ib = 1,num_slip_sys
          temp1(ia,ib) = 0.0
          temp2(ia,ib) = 0.0
          drhomdgb(ia,ib) = 0.0
          drhoidgb(ia,ib) = 0.0
          dldrhom(ia,ib) = 0.0
        end do
      end do

      do ia = 1,num_slip_sys
        sum1 = 0.0
        do ib = 1,num_slip_sys
          sum1 = sum1 + H(ia,ib)*(rho_m(ib) + rho_i(ib))
        end do
        dldrhom(ia,ia) = -0.5*d_disl(ia)*d_disl(ia)*H(ia,ia)/
     &   sqrt(sum1)
!       end do
!
!       do ia = 1,num_slip_sys
        temp2(ia,ia) = temp2(ia,ia) + (q_gen/b/d_disl(ia)
     &   - 1.0/b/d_disl(ia) - (q_ann*R_c*rho_m(ia)/b))*
     &   sgn(gamma_dot(ia))*dt_incr

        temp1(ia,ia) = temp1(ia,ia) + 1.0 +
     &   (q_gen/b/(d_disl(ia)**2.0))*
     &   dldrhom(ia,ia)*abs(gamma_dot(ia))*dt_incr +
     &   q_ann*R_c*abs(gamma_dot(ia))*dt_incr/b +
     &   0.5*x_d*H(ia,ia)*abs(gamma_dot(ia))*dt_incr/
     &   sqrt(sum1)/b
      end do

      call inverse(num_slip_sys,temp1,temp1_inv,isingular)
      if (isingular .eq. 1) go to 500

      do ia = 1,num_slip_sys
        do ib = 1,num_slip_sys
          do ic = 1,num_slip_sys
            drhomdgb(ia,ib) = drhomdgb(ia,ib) + temp1_inv(ia,ic)*
     &       temp2(ic,ib)
          end do
        end do
      end do

c     immobile dislocations
      do ia = 1,num_slip_sys
        sum1 = 0.0
        do ib = 1,num_slip_sys
          sum1 = sum1 + H(ia,ib)*(rho_m(ib) + rho_i(ib))
        end do
        t3 = x_d*sqrt(sum1)/b - q_dyn*rho_i(ia)

        t4 = 1.0e0 - dt_incr*abs(gamma_dot(ia))*0.5*x_d*H(ia,ia)/b/
     &   sqrt(sum1) + q_dyn*dt_incr*abs(gamma_dot(ia))

        drhoidgb(ia,ia) = t3*sgn(gamma_dot(ia))*dt_incr/t4
      end do

c--------------------------------------------------------------------
c  Calculate derivative of the thermal slip resistance, s_t
c  w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------
      do ia = 1,num_slip_sys
        do ib = 1,num_slip_sys
          dstdgb(ia,ib) = 0.0
        end do
      end do

c--------------------------------------------------------------------
c  Calculate derivative of the athermal slip resistance, s_a
c  w.r.t. gamma-dot-beta.
c--------------------------------------------------------------------
      do ia = 1,num_slip_sys
        sum1 = 0.0
        do ib = 1,num_slip_sys
          sum1 = sum1 + p0*A(ia,ib)*(rho_m(ib) + rho_i(ib))

          dsadgb(ia,ib) = 0.5*(q_t*G*b)*(p0*A(ia,ib)*(drhomdgb(ia,ib)
     &     + drhoidgb(ia,ib))/sqrt(sum1))
        end do
      end do

c--------------------------------------------------------------------
c  Form "A-matrix" of derivatives wrt d_gamma_beta.
c--------------------------------------------------------------------
!       checkpoint("Calculating A-matrix")

      do ia = 1,num_slip_sys
        do ib = 1,num_slip_sys

          if(abs(tau(ia)) .lt. s_a(ia)) then
            array3(ia,ib) = - gamma_dot_g(ia)*(p*q*(delF0/B_k/xtemp)*
     &      power(1.0 - power(abs(tau(ia))/s_a(ia), p), q - 1)*
     &      power(abs(tau(ia))/s_a(ia), p - 1))*
     &      (dTaudGd(ia,ib)*sgn(tau(ia))/s_a(ia)
     &      - abs(tau(ia))*(dsadgb(ia,ib) + dstdgb(ia,ib))/
     &      (s_a(ia)*s_a(ia))) - drhomdgb(ia,ib) *
     &      gamma_dot_g(ia)/rho_m(ia)
          else if(abs(tau(ia)) .ge. s_a(ia)) then
            array3(ia,ib) = 0.0;
         end if
        end do
        array3(ia,ia) = array3(ia,ia) + 1
      end do

!       checkpoint("All derivatives calculated.")
c--------------------------------------------------------------------
c  Calculate the gradient of sse wrt gamma_dot(). Will be used
c  later to ensure that line search is in the correct direction.
c--------------------------------------------------------------------
      do j = 1,num_slip_sys
        gradient(j) = 0.0
        do i = 1,num_slip_sys
          gradient(j) = gradient(j) + func(i)*array3(i,j)
        end do
        gradient(j) = 2 * gradient(j)
      end do

c--------------------------------------------------------------------
c  Solve for increment of gamma_dot. Solution is returned
c  in the func() array.
c--------------------------------------------------------------------
      call simeq(num_slip_sys,array3,func,isingular)
      if (isingular .eq. 1) go to 500

c--------------------------------------------------------------------
c  Store offset in d_gamma_dot(k)
c--------------------------------------------------------------------
      do k = 1,num_slip_sys
        d_gamma_dot(k) = func(k)
      end do

c--------------------------------------------------------------------
c  Check to make sure that N-R step leads 'down hill' the
c  sse surface.
c--------------------------------------------------------------------
      sum1 = 0.0
      do i = 1,num_slip_sys
        sum1 = sum1 - gradient(i)*d_gamma_dot(i)
      end do

      if (sum1 .gt. 0.0) then
        do i = 1,num_slip_sys
          d_gamma_dot(i) = -d_gamma_dot(i)
        end do
      end if

c--------------------------------------------------------------------
c  Multiply step size by two 'cause next loop will divide it by 2.
c--------------------------------------------------------------------
      do k = 1,num_slip_sys
        d_gamma_dot(k) = d_gamma_dot(k)*2
      end do

c====================================================================
c  Begin line search.
c====================================================================

      improved = .false.

      do N_ctr = 1,max_loops

      sse_old = sse

c--------------------------------------------------------------------
c  Divide step size by two.
c--------------------------------------------------------------------
c     Weighted estimate of new values of gamma_try_g and gamma_try_c
!       checkpoint("Estimating new gamma_dots")
      do k = 1,num_slip_sys
        d_gamma_dot(k) = d_gamma_dot(k)/2.0e0

        gamma_try_g(k) = gamma_dot_g(k) + d_gamma_dot(k)
        gamma_try(k) = gamma_try_g(k)
      end do

c--------------------------------------------------------------------
c  Calculate g, F_p_inv, F_el, sig, tau, func, and sse based
c  on gamma_try(k)
c--------------------------------------------------------------------
!       checkpoint("Call eval_func() based on new estimates")
      call eval_func(xs0, xm0, dt_incr,
     &  gamma_try, gamma_try_g,
     &  F1, F_p_inv, F_p_inv_0, F_el, num_slip_sys, C,
     &  rho_m0, rho_m, rho_i0, rho_i,
     &  Alatent, p0, B_k, R_c, q_dyn,x_d, drag_stress, q_gen, q_ann,
     &  delF0, G, b, xtemp, sse, sig, tau,
     &  s_t, s_a, A, H, p, q, q_t, ! shr_g, shr_c,
     &  func, xL_p, gammadot0g, refgdg, rho_m_zero, d_disl_zero,
     &  tau0, hp_coeff, grain_size,isingular)
      if (isingular .eq. 1) go to 500

!       checkpoint("2nd run of eval_func() complete.")

c--------------------------------------------------------------------
c  Check for 'convergence' of the line search.  Note that the line
c  search doesn't waste time converging closely.  This is 'cause
c  the N-R step does so much better.
c--------------------------------------------------------------------

      if ((sse_old.le.sse_ref).and.(sse.ge.sse_old).and.
     & (N_ctr.gt.1)) improved=.true.

      if (improved) go to 200

      end do ! Linear Search

c--------------------------------------------------------------------
c  Add "d_gamma_dot" to gamma_dot to get new values for
c  this iteration.
c--------------------------------------------------------------------

  200 do k = 1,num_slip_sys
        gamma_dot(k) = gamma_dot(k) + d_gamma_dot(k)*2.0
        gamma_dot_g(k) = gamma_dot_g(k) + d_gamma_dot(k)*2.0
      end do

c--------------------------------------------------------------------
c  If (sse_old > tolerance) then this step has not converged.
c--------------------------------------------------------------------
      if (sse_old .gt. tolerance) converged = .false.

c--------------------------------------------------------------------
c  If (sse_old > sse_ref/2) then convergence is too slow and
c  increment is divided into two subincrements.
c--------------------------------------------------------------------
      if ((sse_old.gt.sse_ref/2.0).and.(.not.converged)) then
        N_incr = 2 * N_incr - 1
        N_incr_total = 2 * N_incr_total
        if (N_incr_total.gt.100000)  go to 500
c        print*, N_incr_total
        go to 100
      end if

c--------------------------------------------------------------------
c  End iterative loop.
c--------------------------------------------------------------------

      end do ! 'Newton Raphson Iterative Loop'
!       checkpoint("End of Newton-Rhapson step")
c--------------------------------------------------------------------
c  If another subincrement remains to be done, then reinitialize
c  F_p_inv_0 and g0.  F0 and F1 gets reinitialized back at the
c  top of this loop.
c--------------------------------------------------------------------
      if (N_incr .lt. N_incr_total) then
        if (N_incr .eq. (N_incr/2)*2) then ! N_incr is 'even'
          N_incr = N_incr/2 + 1
          N_incr_total = N_incr_total/2
        else ! N_incr is 'odd'
          N_incr = N_incr + 1
        end if

        do i = 1,3
          do j = 1,3
            F_p_inv_0(i,j) = F_p_inv(i,j)
          end do
        end do

        do i = 1,num_slip_sys
          rho_m0(i) = rho_m(i)
          rho_i0(i) = rho_i(i)
        end do

        go to 100
      end if

! c--------------------------------------------------------------------
! c  Store shearing rates.
! c--------------------------------------------------------------------
!
!       do ia = 1,num_slip_sys
!         gamma_dot_g(ia) = shr_g(ia)
!         gamma_dot(ia) = shr_g(ia)
!       end do

c--------------------------------------------------------------------
c  Calculate the average Cauchy stress.
c--------------------------------------------------------------------
      do i = 1,3
        do j = 1,3
          sig_avg(i,j) = sig_avg(i,j) + sig(i,j)
        end do
      end do

c--------------------------------------------------------------------
c  Write out Euler Angles in Roe/Bunge Convention.
c--------------------------------------------------------------------

      ! Note: Update dir_cos0 to dir_cos using F_el
      ! This dir_cos is directly used in next timestep
      call aa_dot_bb(3,F_el,dir_cos0,dir_cos)
      call aa_dot_bb(3,F_el,dir_cos0,array1)

      ! updated euler angles is psi
      psi(1) = psi0(1)
      psi(2) = psi0(2)
      psi(3) = psi0(3)
      call bunge_angles(time(2),array1,psi) ! This is only a post-processing step

c--------------------------------------------------------------------
c  Calculate the average elasticity tensor.
c--------------------------------------------------------------------
      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          C_avg(i,j,k,l) = C(i,j,k,l)
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c--------------------------------------------------------------------
c  Rotate xs0 and xm0 to current coordinates, called xs and xm.
c--------------------------------------------------------------------
      call inverse_3x3(F_el,F_el_inv,isingular)
      if (isingular .eq. 1) go to 500
      do n = 1,num_slip_sys
        do i = 1,3
          xs(i,n) = 0.0
          xm(i,n) = 0.0
          do j = 1,3
            xs(i,n) = xs(i,n) + F_el(i,j)*xs0(j,n)
            xm(i,n) = xm(i,n) + xm0(j,n)*F_el_inv(j,i)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate the derivative of the plastic part of the rate of
c  deformation tensor in the current configuration wrt sigma.
c--------------------------------------------------------------------
      ddpdsig(:,:,:,:) = 0.0

c      do i = 1,3
c       do j = 1,3
c        do k = 1,3
c         do l = 1,3
c          ! if (m .eq. 1) ddpdsig(i,j,k,l) = 0.0

c          do n = 1,num_slip_sys

c           refgdg = gammadot0g*rho_m(n)/rho_m_zero

c           ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l) + (xs(i,n)*xm(j,n) +
c     &      xm(i,n)*xs(j,n))*(((xs(k,n)*xm(l,n) + xm(k,n)*xs(l,n))*
c     &      (gamma_dot_g(n)*(p*q*(delF0/B_k/xtemp)*
c     &      power(1 - power((abs(tau(n)) - s_a(n))/s_t(n),p)
c     &      ,q-1)*power((abs(tau(n)) - s_a(n))/s_t(n),p-1))*
c     &      sgn(tau(n))/s_t(n))))

c          end do ! num_slip_sys

c         end do ! l
c        end do ! k
c       end do ! j
c      end do ! i

c      write(7,'(12f7.3)')(gamma_dot(i)/0.0004,i=1,12)

c--------------------------------------------------------------------
c  End calculations over the integration point
c--------------------------------------------------------------------

c--------------------------------------------------------------------
c  Calculate Green strain
c--------------------------------------------------------------------
      call transpose1(3,F1,array1)
      call aa_dot_bb(3,array1,F1,E_tot)

      do i = 1,3
        E_tot(i,i) = E_tot(i,i) - 1
        do j = 1,3
          E_tot(i,j) = E_tot(i,j)/2
        end do
      end do

c====================================================================
c  Begin calculation of the Jacobian (the tangent stiffness matrix)
c====================================================================

c--------------------------------------------------------------------
c  Store sig_avg() into sig()
c--------------------------------------------------------------------
      do i = 1,3
        do j = 1,3
          sig(i,j) = sig_avg(i,j)
        end do
      end do

! c--------------------------------------------------------------------
! c  Divide C_avg by num_grains to get correct average.
! c--------------------------------------------------------------------
!
!       do i = 1,3
!        do j = 1,3
!         do k = 1,3
!          do l = 1,3
!           C_avg(i,j,k,l) = C_avg(i,j,k,l)
!          end do ! l
!         end do ! k
!        end do ! j
!       end do ! i

c--------------------------------------------------------------------
c  Output stress, strain, etc. for post processing.
c--------------------------------------------------------------------
      sum1 = 0.0
      do i = 1,3
       do j = 1,3
        sum1 = sum1 + sig(i,j) * sig(i,j)
       end do
      end do
      trace = sig(1,1) + sig(2,2) + sig(3,3)
      sig_eff=sqrt(1.5 * sum1 - 0.5 * trace**2)

c--------------------------------------------------------------------
c  Calculate the inverse of F_el
c--------------------------------------------------------------------
      call inverse_3x3(F_el,F_el_inv,isingular)
      if (isingular .eq. 1) go to 500

c--------------------------------------------------------------------
c  Scale by appropriate constants and divide by num_grains to
c  get the average value.  ALSO multiply by 'dtime' which is
c  d(sig)/d(sig_dot).
c--------------------------------------------------------------------
c      do i = 1,3
c       do j = 1,3
c        do k = 1,3
c         do l = 1,3
c          ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l)*dtime
c     &     /4.0
c         end do ! l
c        end do ! k
c       end do ! j
c      end do ! i

c--------------------------------------------------------------------
c  Multiply the 4th rank elastic stiffness tensor by the derivative
c  of the plastic part of the rate of deformation tensor wrt sig_dot.
c--------------------------------------------------------------------
c      call aaaa_dot_dot_bbbb(3,C_avg,ddpdsig,array6)

c--------------------------------------------------------------------
c  Add 4th rank identity tensor to array6()
c--------------------------------------------------------------------
c      do i = 1,3
c       do j = 1,3
c        do k = 1,3
c         do l = 1,3
c          array6(i,j,k,l) = array6(i,j,k,l) + 0.5*
c     &     (del(i,k)*del(j,l) + del(i,l)*del(j,k))
c         end do
c        end do
c       end do
c      end do

c--------------------------------------------------------------------
c  Need to take the inverse of Array4.  Since it relates two 2nd
c  rank tensors that are both symmetric, Array4 can be xformed to
c  Voigt notation style in order to do the inverse, then xformed back.
c--------------------------------------------------------------------
c      call forth_to_Voigt(Array6,Array4)
c      call inverse(6,Array4,Array5,isingular)
c      if (isingular .eq. 1) go to 500
c      call Voigt_to_forth(Array5,Array6)

c--------------------------------------------------------------------
c  Multiply Array6 by C, the elastic stiffness matrix to
c  finally get the Jacobian, but as a 4th rank tensor.
c--------------------------------------------------------------------
c      call aaaa_dot_dot_bbbb (3,Array6,C_avg,ddsdde_4th)

c--------------------------------------------------------------------
c  Store the stress tensor in the ABAQUS stress 'vector'
c--------------------------------------------------------------------
c      do i = 1,ndi
c        STRESSNEW(ia_block,i) = sig(i,i)
c      end do
c      if (nshr .eq. 1) STRESSNEW(ia_block,ndi+1) = sig(1,2)
c      if (nshr .eq. 3) then
c        STRESSNEW(ia_block,4) = sig(1,2)
c        STRESSNEW(ia_block,5) = sig(1,3)
c        STRESSNEW(ia_block,6) = sig(2,3)
c      end if

      STRESSNEW(ia_block,1) = sig(1,1)
      STRESSNEW(ia_block,2) = sig(2,2)
      STRESSNEW(ia_block,3) = sig(3,3)
      STRESSNEW(ia_block,4) = sig(1,2)
      STRESSNEW(ia_block,5) = sig(2,3)
      STRESSNEW(ia_block,6) = sig(3,1)

c--------------------------------------------------------------------
c  Store the Jacobian in Voigt notation form.
c--------------------------------------------------------------------
c      do i = 1,3
c       do j = i,3   ! not 1 !!!
c        ia = i
c        if (i.ne.j) ia=i+j+1
c        do k = 1,3
c         do l = k,3 ! not 1 !!!
c          ib = k
c          if (k.ne.l) ib=k+l+1
c            array4(ia,ib) = ddsdde_4th(i,j,k,l)
c            if(ib.ge.4) array4(ia,ib) = 2 * array4(ia,ib)
c         end do
c        end do
c       end do
c      end do

c      do i =1,6
c        do j = 1,6
c          ddsdde(i,j) = 0.0
c        end do
c      end do

c      if (ndi .eq. 1) then ! 1-D
c         ddsdde(1,1) = array4(1,1)
c      else if (ndi .eq. 2) then ! 2-D plane stress & axi
c         do i = 1,2
c            do j = 1,2
c               ddsdde(i,j) = array4(i,j)
c            end do
c         end do
c         ddsdde(1,3) = array4(1,4)
c         ddsdde(2,3) = array4(2,4)
c         ddsdde(3,1) = array4(4,1)
c         ddsdde(3,2) = array4(4,2)
c         ddsdde(3,3) = array4(4,4)
c      else if (ndi .eq. 3 .and. nshr .eq. 1) then ! plane strain
c         do i = 1,4
c            do j = 1,4
c               ddsdde(i,j) = array4(i,j)
c            end do
c         end do
c      else ! Full 3-D
c         do i = 1,6
c            do j = 1,6
c               ddsdde(i,j) = array4(i,j)
c            end do
c         end do
c      end if

c===================================================================
c  Store the internal variables in the STATENEW() array
c===================================================================

      n = 0
c-------------------------------------------------------------------
c  Store the Original Euler Angles  1-3
c-------------------------------------------------------------------
      do i = 1,3
            n = n + 1
            STATENEW(ia_block,n) = psi0(i)
            !print*, STATENEW(ia_block,n)
      end do

c-------------------------------------------------------------------
c  Store the elastic part of F  4-12
c-------------------------------------------------------------------
      do i = 1,3
        do j = 1,3
          n = n + 1
          STATENEW(ia_block,n) = F_el(i,j)
c          print*,STATENEW(ia_block,n)
        end do
      end do

c-------------------------------------------------------------------
c  Store the plastic part of F   13-21
c-------------------------------------------------------------------
      do i = 1,3
        do j = 1,3
          n = n + 1
          STATENEW(ia_block,n) = F_p_inv(i,j)
c          print*,STATENEW(ia_block,n)
        end do
      end do

c-------------------------------------------------------------------
c  Plastic strain calculations
c-------------------------------------------------------------------
c  Increment of plastic shear strain accumulated on each slip system
c  over this time step.
      do i = 1,num_slip_sys
        delta_gamma(i) = gamma_dot(i)*dtime
      end do

      do j = 1,3
        do k = 1,3
          delta_E_p(j,k) =  0.0
        end do
      end do

c  Increment of the plastic strain tensor
      do j = 1,3
        do k = 1,3
          do l = 1,num_slip_sys
            delta_E_p(j,k) = delta_E_p(j,k) + 0.5*delta_gamma(l)*
     &       (xs0(j,l)*xm0(k,l)+xs0(k,l)*xm0(j,l))
          end do
        end do
      end do

c  Plastic strain tensor
      do i = 1,3
        do j = 1,3
          E_p(i,j) = E_p(i,j)+delta_E_p(i,j)
        end do
      end do

      ! SDV 22-30
      do i = 1,3
        do j = 1,3
          n = n + 1
          STATENEW(ia_block,n) = E_p(i,j)
c          print*,STATENEW(ia_block,n)
        end do
      end do

c  Effective total strain
      call aa_dot_dot_bb(3,E_tot,E_tot,sum1)
      E_eff = sqrt(2./3.*sum1)

      ! SDV 31
      n = n + 1
      STATENEW(ia_block,n)= E_eff
c      print*,STATENEW(ia_block,n)

c  Effective plastic strain
      call aa_dot_dot_bb(3,E_p,E_p,sum1)
      E_p_eff = sqrt(2./3.*sum1)

      if (E_p_eff.lt.0) E_p_eff = 0.
      ! SDV 32
      n = n + 1
      STATENEW(ia_block,n)= E_p_eff
c      print*,STATENEW(ia_block,n)

c  Effective plastic strain increment
      sum2 = 0.0
      call aa_dot_dot_bb(3,delta_E_p,delta_E_p,sum2)
      sum2=sqrt(2./3.*sum2)
      E_p_eff_cum = sum2 + E_p_eff_cum

      ! SDV 33
      n = n + 1
      STATENEW(ia_block,n) = E_p_eff_cum
c      print*,STATENEW(ia_block,n)

c---------------------------------------------------------------------
c  Store average values of defect densities and loop sizes   34-35
c---------------------------------------------------------------------
      rho_m_avg = 0.0
      rho_i_avg = 0.0
      gamma_dot_g_avg = 0.0

      do i = 1,num_slip_sys
        rho_m_avg = rho_m_avg + rho_m(i)
        rho_i_avg = rho_i_avg + rho_i(i)
        gamma_dot_g_avg = gamma_dot_g_avg + gamma_dot_g(i)
      end do
      rho_tot = rho_m_avg + rho_i_avg

      n=n+1
      STATENEW(ia_block,n) = rho_m_avg/num_slip_sys
      !print*, STATENEW(ia_block,n)
      n=n+1
      STATENEW(ia_block,n) = rho_i_avg/num_slip_sys

c---------------------------------------------------------------------
c  Store average shearing rates   36
c---------------------------------------------------------------------
      n=n+1
      STATENEW(ia_block,n) = gamma_dot_g_avg/num_slip_sys

c-------------------------------------------------------------------
c  Store the mobile dislocation density     37-48
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
        n = n + 1
        STATENEW(ia_block,n) = rho_m(i)
        !print*,STATENEW(ia_block,n), n, KINC
      end do

c-------------------------------------------------------------------
c  Store the immobile dislocation density     49-60
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
        n = n + 1
        STATENEW(ia_block,n) = rho_i(i)
        !print*, STATENEW(ia_block,n)
      end do

c-------------------------------------------------------------------
c  Store the updated Euler Angles  61, 62, 63
c-------------------------------------------------------------------
      do i = 1,3
            n = n + 1
            STATENEW(ia_block,n) = psi(i)
            !print*, STATENEW(ia_block,n)
      end do

c-------------------------------------------------------------------
c  Store the dir_cos0 matrix   64-72
c-------------------------------------------------------------------
      do i = 1,3
        do j = 1,3
          n = n + 1
          STATENEW(ia_block,n) = dir_cos0(i,j)
c          print*,STATENEW(ia_block,n)
        end do
      end do

c-------------------------------------------------------------------
c  Store the dir_cos matrix   73-81
c-------------------------------------------------------------------
      do i = 1,3
        do j = 1,3
          n = n + 1
          STATENEW(ia_block,n) = dir_cos(i,j)
c          print*,STATENEW(ia_block,n)
        end do
      end do

c-------------------------------------------------------------------
c  Store the slip rate  82-93
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
        n = n + 1
        STATENEW(ia_block,n) = tau(i)
        !print*,STATENEW(ia_block,n), n, KINC
      end do

c-------------------------------------------------------------------
c  Store the gamma_dot_g  94-105
c-------------------------------------------------------------------
      do i = 1,num_slip_sys
        n = n + 1
        STATENEW(ia_block,n) = gamma_dot_g(i)
        !print*,STATENEW(ia_block,n), n, KINC
      end do

c      write(11, *) gamma_dot_g(1), gamma_dot_g(2), gamma_dot_g(3),
c     & gamma_dot_g(4), gamma_dot_g(5), gamma_dot_g(6),
c     & gamma_dot_g(7), gamma_dot_g(8), gamma_dot_g(9),
c     & gamma_dot_g(10), gamma_dot_g(11), gamma_dot_g(12)
c-------------------------------------------------------------------
c  SDV controlling element deletion     106
c-------------------------------------------------------------------
      n = n + 1
      STATENEW(ia_block,n) = 1

c-------------------------------------------------------------------
c  Update internal energy
c-------------------------------------------------------------------

      enerInternNew(ia_block) = enerInternOld(ia_block)
     &     +((stressOld(ia_block,1)+stressNew(ia_block,1))
     &     *strainInc(ia_block,1)+(stressOld(ia_block,2)+
     &     stressNew(ia_block,2))*strainInc(ia_block,2)
     &     +(stressOld(ia_block,3)+stressNew(ia_block,3))
     &     *strainInc(ia_block,3)+2.0*(stressOld(ia_block,4)+
     &     stressNew(ia_block,4))*strainInc(ia_block,4)
     &     +2.0*(stressOld(ia_block,5)+stressNew(ia_block,5))
     &     *strainInc(ia_block,5)+2.0*(stressOld(ia_block,6)+
     &     stressNew(ia_block,6))*strainInc(ia_block,6))
     &     *0.5/DENSITY(ia_block)

c-------------------------------------------------------------------
c  Update inelastic energy
c-------------------------------------------------------------------
      dev_stress(:,:) = 0.0
      sigma_mean = 0.0
      sum1 = 0.0
      sigma_mean = (stressNew(ia_block,1) + stressNew(ia_block,2)
     &     + stressNew(ia_block,3))/3
      do i=1,3
            do j=1,3
                  dev_stress(i,j) = sig(i,j)
            end do !j
            dev_stress(i,i) = dev_stress(i,i) - sigma_mean
      end do !i

      call aa_dot_dot_bb(3,dev_stress,dev_stress,sum1)
      eff_dev_stress = sqrt(3/2 * sum1)

      ENERINELASNEW(ia_block) = ENERINELASOLD(ia_block)
     &      + eff_dev_stress * E_p_eff_cum / DENSITY(ia_block)

c-------------------------------------------------------------------
!       checkpoint("UMAT complete")

      if (N_incr_total .gt. 100000) then
        write(*,*) N_incr_total
  500   DTIME = 0.5*DTIME
      end if

c      ! elasticity tensor
c      do i = 1,3
c        do j = 1,3
c          do k = 1,3
c            do l = 1,3
c              Cel(i,j,k,l) = C_avg(i,j,k,l)
c            end do
c          end do
c        end do
c      end do

      end do !  end loop over all material points (NBLOCK)

      return
      end


c====================================================================
c====================================================================
c====================== S U B R O U T I N E S =======================
c====================================================================
c====================================================================
c
c  Calculate a vector cross product.
c
c  c = a X b
c
c--------------------------------------------------------------------

      subroutine cross_product(a1,a2,a3,b1,b2,b3,c1,c2,c3)

      implicit double precision (a-h,o-z)

      c1 = a2*b3 - a3*b2
      c2 = a3*b1 - a1*b3
      c3 = a1*b2 - a2*b1

      return
      end

c====================================================================
c====================================================================
c
c  Calculate sin(theta) between two vectors.
c
c--------------------------------------------------------------------

      subroutine calc_sin(x1,x2,x3,y1,y2,y3,z)

      implicit double precision (a-h,o-z)

      z = sqrt(1 - (x1*y1 + x2*y2 + x3*y3)**2)

      return
      end

c====================================================================
c====================================================================
c
c  Normalize the length of a vector to one.
c
c--------------------------------------------------------------------

      subroutine normalize_vector(x,y,z)

      implicit double precision (a-h,o-z)

      xlength = sqrt(x*x+y*y+z*z)
      x = x/xlength
      y = y/xlength
      z = z/xlength

      return
      end

c====================================================================
c====================================================================
c
c  Transpose a (n x n) tensor.
c
c--------------------------------------------------------------------

      subroutine transpose1(n,a,b)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n)

      do i = 1,n
        do j = 1,n
          b(i,j) = a(j,i)
        end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Calculate the dot product of two 2nd rank tensors.
c  Result is stored in cc(i,j)
c
c--------------------------------------------------------------------

      subroutine aa_dot_bb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n), c(n,n)

      do i = 1,n
        do j = 1,n
          c(i,j) = 0.0e0
          do k = 1,n
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          end do
        end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of two 2nd rank tensors.
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bb(n,a,b,sum1)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n)

      sum1 = 0.0
      do i = 1,n
        do j = 1,n
          sum1 = sum1 + a(i,j)*b(i,j)
        end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of two 4th rank tensors.
c  Result is stored in c(i,j,k,l)
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bbbb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n,n,n), b(n,n,n,n), c(n,n,n,n)

      do i = 1,n
       do j = 1,n
        do k = 1,n
         do l = 1,n
          c(i,j,k,l) = 0
          do m1 = 1,n
           do m2 = 1,n
            c(i,j,k,l) = c(i,j,k,l) + a(i,j,m1,m2)*b(m1,m2,k,l)
           end do !m2
          end do !m1
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of a 4th rank tensor and
c  a 2nd rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aaaa_dot_dot_bb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n,n,n), b(n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(i,j,k,l)*b(k,l)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of a 2nd rank tensor and
c  a 4th rank tensor.  Result is stored in c(i,j).
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bbbb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n,n,n), c(n,n)

      do i = 1,n
       do j = 1,n
        c(i,j) = 0
        do k = 1,n
         do l = 1,n
          c(i,j) = c(i,j) + a(k,l)*b(k,l,i,j)
         end do !l
        end do !k
       end do !j
      end do !i

      return
      end

c====================================================================
c====================================================================
c
c  Rotates any 3x3x3x3 tensor by a rotation matrix.
c
c  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
c
c--------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3,3,3), c(3,3,3,3), d(3,3,3,3)

      do m = 1,3
       do n = 1,3
        do k = 1,3
         do l = 1,3
          d(m,n,k,l) = a(k,1)*(a(l,1)*b(m,n,1,1) +
     &     a(l,2)*b(m,n,1,2) + a(l,3)*b(m,n,1,3)) +
     &     a(k,2)*(a(l,1)*b(m,n,2,1) +
     &     a(l,2)*b(m,n,2,2) + a(l,3)*b(m,n,2,3)) +
     &     a(k,3)*(a(l,1)*b(m,n,3,1) +
     &     a(l,2)*b(m,n,3,2) + a(l,3)*b(m,n,3,3))
         end do
        end do
       end do
      end do

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          c(i,j,k,l) = a(i,1)*(a(j,1)*d(1,1,k,l) +
     &     a(j,2)*d(1,2,k,l) + a(j,3)*d(1,3,k,l)) +
     &     a(i,2)*(a(j,1)*d(2,1,k,l) +
     &     a(j,2)*d(2,2,k,l) + a(j,3)*d(2,3,k,l)) +
     &     a(i,3)*(a(j,1)*d(3,1,k,l) +
     &     a(j,2)*d(3,2,k,l) + a(j,3)*d(3,3,k,l))
         end do
        end do
       end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      subroutine inverse_3x3(a,b,isingular)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3)

      isingular = 0

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      b(2,1) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      b(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      b(2,3) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(3,2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)

      det = a(1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1)

      if (det .eq. 0.) then
        isingular = 1
        return
      end if

      do i = 1,3
         do j = 1,3
            b(i,j) = b(i,j)/det
         end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Solve simultaneous equations using LU decomposition (Crout's method)
c  Result is stored in b(i)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine simeq(n,a,b,isingular)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n), index(n)

      isingular = 0

      call LU_Decomp(n,a,index,isingular)
      if (isingular .eq. 1) return
      call LU_BackSub(n,a,index,b)

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a matrix using
c  LU decomposition (Crout's method)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine inverse(n,a,b,isingular)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n), c(n,n), index(n)

      isingular = 0
      do i = 1,n
         do j = 1,n
            c(i,j) = a(i,j)
         end do
      end do

      do i = 1,n
         do j = 1,n
            b(i,j) = 0.0
         end do
         b(i,i) = 1.0
      end do

      call LU_Decomp(n,c,index,isingular)
      if (isingular .eq. 1) return
      do j = 1,n
         call LU_BackSub(n,c,index,b(1,j))
      end do

      return
      end

c====================================================================
c====================================================================
c
c  This sub performs an LU Decomposition (Crout's method) on the
c  matrix "a". It uses partial pivoting for stability. The index()
c  vector is used for the partial pivoting.  The v() vector is
c  a dummy work area.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine LU_Decomp(n,a,index,isingular)

      implicit double precision (a-h,o-z)

      dimension a(n,n), index(n), v(n)

      tiny = 1.0e-20

c--------------------------------------------------------------------
c  Loop over the rows to get the implicit scaling info.
c--------------------------------------------------------------------

      do i = 1,n
         a_max = 0.0
         do j = 1,n
            a_max = max(a_max,abs(a(i,j)))
         end do !j
         if(a_max.eq.0.) then
          isingular=1
          return
         endif
         v(i) = 1.0/a_max
      end do !i

c--------------------------------------------------------------------
c  Begin big loop over all the columns.
c--------------------------------------------------------------------

      do j = 1,n

         do i = 1,j-1
            sum1 = a(i,j)
            do k = 1,i-1
               sum1 = sum1 - a(i,k)*a(k,j)
            end do
            a(i,j) = sum1
         end do

         a_max = 0.0
         do i = j,n
            sum1 = a(i,j)
            do k = 1,j-1
               sum1 = sum1 - a(i,k)*a(k,j)
            end do
            a(i,j) = sum1
            dummy = v(i) * abs(sum1)
            if (dummy .gt. a_max) then
               imax = i
               a_max = dummy
            end if
         end do

c--------------------------------------------------------------------
c  Pivot rows if necessary.
c--------------------------------------------------------------------

         if (j .ne. imax) then
            do k = 1,n
               dummy = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dummy
            end do
            v(imax) = v(j)
         end if
         index(j) = imax

c--------------------------------------------------------------------
c  Divide by the pivot element.
c--------------------------------------------------------------------

         if (a(j,j) .eq. 0.0) a(j,j) = tiny
         if (j .ne. n) then
            dummy = 1.0/a(j,j)
            do i = j+1,n
               a(i,j) = a(i,j)*dummy
            end do
         end if

      end do !j
      isingular=0

      return
      end

c====================================================================
c====================================================================
c
c  Solves a set of simultaneous equations by doing back substitution.
c  The answer in returned in the b() vector.  The a(,) matrix
c  must have already been "LU Decomposed" by the above subroutine.
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine LU_BackSub(n,a,index,b)

      implicit double precision (a-h,o-z)

      dimension a(n,n), index(n), b(n)

      ii = 0

c--------------------------------------------------------------------
c  Do the forward substitution.
c--------------------------------------------------------------------

      do i = 1,n
         m = index(i)
         sum1 = b(m)
         b(m) = b(i)
         if (ii .ne. 0) then
            do j = ii,i-1
               sum1 = sum1 - a(i,j)*b(j)
            end do
         else if (sum1 .ne. 0.0) then
            ii = i
         end if
         b(i) = sum1
      end do

c--------------------------------------------------------------------
c  Do the back substitution.
c--------------------------------------------------------------------

      do i = n,1,-1
         sum1 = b(i)
         if (i .lt. n) then
            do j = i+1,n
               sum1 = sum1 - a(i,j)*b(j)
            end do
         end if
         b(i) = sum1/a(i,i)
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Restore a symmetric 4th rank tensor stored in Voigt notation
c  back to its 4th rank form.
c
c--------------------------------------------------------------------

      subroutine Voigt_to_forth(b,a)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = 1,3
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = 1,3
          ib = k
          if (k.ne.l) ib=9-k-l
          a(i,j,k,l) = b(ia,ib)
          if (ia.gt.3) a(i,j,k,l) = a(i,j,k,l)/2
          if (ib.gt.3) a(i,j,k,l) = a(i,j,k,l)/2
         end do
        end do
       end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Store a SYMMETRIC 4th rank tensor in Voigt notation.
c
c--------------------------------------------------------------------

      subroutine forth_to_Voigt(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=9-k-l
          b(ia,ib) = a(i,j,k,l)
         end do
        end do
       end do
      end do

      return
      end



c====================================================================
c====================================================================
c
c  Perform x**y but while retaining the sign of x.
c
c--------------------------------------------------------------------

      function power(x,y)

      implicit double precision (a-h,o-z)

      if (x.eq.0.0) then
        if (y.gt.0.0) then
          power = 0.0
        else if (y .lt. 0.0) then
          power = 1.0d+300
        else
          power = 1.0
        end if
      else
         power = y*log10(abs(x))
         if (power .gt. 300.) then
           power = 1.d+300
         else
           power = 10.d0 ** power
         end if
         if (x .lt. 0.0) power = -power
      end if

      return
      end

c====================================================================
c====================================================================
c
c  Return the sign of a number.
c
c--------------------------------------------------------------------

      function sgn(a)

      implicit double precision (a-h,o-z)

      sgn = 1.0
      if (a .lt. 0.0) sgn = -1.0

      return
      end

! c===================================================================
! c===================================================================
! c
! c  Print out an array.
! c
! c-------------------------------------------------------------------
!
!       subroutine print_array(n,a)
!
!       implicit double precision (a-h,o-z)
!       dimension a(n,n)
!
!       do i = 1,n
!          write(6,'(10f12.3)')(a(i,j),j=1,n)
!       end do
!       print*,' '
!
!       return
!       end

c===================================================================
c===================================================================
c
c  Print out Euler angles in Roe notation.
c  This is modified from the Kocks subroutine
c  In the Roe notation: the third angle, phi = 180 - phi(Kocks)
c  Ref.: Slide 31, http://engineering.snu.ac.kr/lecture/texture&anisotropy/Texture%20&%20Anisotropy%2010(Representation).pdf
c
c-------------------------------------------------------------------

      subroutine roe_angles(time,array1,psi)

      implicit double precision (a-h,o-z)

      dimension array1(3,3),psi(3)

      pi = 4.0 * atan(1.0)

      if (abs(array1(3,3)) .gt. 0.99999) then
        psix   = atan2(array1(2,1),array1(1,1))
        thetay = 0.0
        phiz   = 0.0
      else
        psix   = atan2(array1(2,3),array1(1,3))
        thetay = acos(array1(3,3))
        phiz   = atan2(array1(3,2),-array1(3,1))
      end if

      psi(1) = 180*psix/pi
      psi(2) = 180*thetay/pi
      psi(3) = 180 - 180*phiz/pi


c      write(7,'(a3,2i5,4f10.3)')'BM2',ninc,m,time,psi,theta,phi

      return
      end

c===================================================================
c===================================================================
c
c  Print out Euler angles in Bunge notation.
c
c-------------------------------------------------------------------

      subroutine bunge_angles(time,array1,psi)

      implicit double precision (a-h,o-z)

      dimension array1(3,3),psi(3)

      PI = 4.0*atan(1.0)

      ! From VPSC
      if (array1(3,3) .ge. 1.0e0) array1(3,3) = 1.0e0
      if (array1(3,3) .le. -1.0e0) array1(3,3) = -1.0e0
      APhi = acos(array1(3,3))

      if(abs(array1(3,3)).ge.0.9999) then
        Aphi_2 = 0.
        Aphi_1 = atan2(array1(1,2),array1(1,1))
      else
        sth=sin(APhi)
        Aphi_2 = atan2(array1(1,3)/sth,array1(2,3)/sth)
        Aphi_1 = atan2(array1(3,1)/sth,-array1(3,2)/sth)
      endif

      psi(1) = 180.0d0*Aphi_1/PI
      psi(2) = 180.0d0*APhi/PI
      psi(3) = 180.0d0*Aphi_2/PI


c      write(7,'(a3,2i5,4f10.3)')'BM2',ninc,m,time,psi,theta,phi

      return
      end

c=======================================================================
c=======================================================================
c
c  Evaluate function to be minimized to zero.  gamma_dot()'s are
c  input and several things are output.
c
c-----------------------------------------------------------------------

      subroutine eval_func(xs0, xm0, dtime,
     &  gamma_dot, gamma_dot_g,
     &  F1, F_p_inv, F_p_inv_0, F_el, num_slip_sys, C,
     &  rho_m0, rho_m, rho_i0, rho_i,
     &  Alatent, p0, B_k, R_c, q_dyn, x_d, drag_stress, q_gen, q_ann,
     &  delF0, G, b, xtemp, sse, sig, tau,
     &  s_t, s_a, A, H, p, q, q_t, ! shr_g, shr_c,
     &  func, xL_p, gammadot0g, refgdg, rho_m_zero, d_disl_zero,
     &  tau0, hp_coeff, grain_size, isingular)

      implicit double precision (a-h,o-z)

      dimension
     &  xm(3,num_slip_sys), xm0(3,num_slip_sys),
     &  xs(3,num_slip_sys), xs0(3,num_slip_sys),
     &  F1(3,3), F_p_inv(3,3), C(3,3,3,3), F_el_inv(3,3), xL_p(3,3),
     &  xL_p_inter(3,3), F_el(3,3), F_p_inv_0(3,3), E_el(3,3),
     &  Spk2(3,3), sig(3,3), tau(num_slip_sys),
     &  array1(3,3), array2(3,3), func(num_slip_sys),
     &  gamma_dot(num_slip_sys),
     &  gamma_dot_g(num_slip_sys),
     &  shr_g(num_slip_sys),
     &  A(num_slip_sys,num_slip_sys),
     &  H(num_slip_sys,num_slip_sys),
     &  s_t(num_slip_sys), s_a(num_slip_sys),
     &  rho_m0(num_slip_sys), rho_m(num_slip_sys),
     &  rho_i0(num_slip_sys), rho_i(num_slip_sys),
     &  drhomdt(num_slip_sys), drhoidt(num_slip_sys),
     &  d_disl(num_slip_sys)


      Pi = 4.*atan(1.)

c*** Note that xs0 and xm0 are in INTERMEDIATE configuration!!!
!       checkpoint("Begin eval_func()")
c--------------------------------------------------------------------
c  Calculate the plastic part of the
c  velocity gradient in the intermediate configuration.
c--------------------------------------------------------------------
!       checkpoint("Calculate plastic velocity gradient")
      do i = 1,3
        do j = 1,3
          xL_p_inter(i,j) = 0.0
          do k = 1,num_slip_sys
            xL_p_inter(i,j) = xL_p_inter(i,j) +
     &       xs0(i,k)*xm0(j,k)*gamma_dot(k)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Begin calculation process of F_p_n+1 = exp(xL_p_inter*dt).F_p_n
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          array1(i,j) = xL_p_inter(i,j)*dtime
        end do
      end do

c--------------------------------------------------------------------
c  Calculate omega.
c--------------------------------------------------------------------

      sum1 = 0
      do i = 1,3
        do j = 1,3
          sum1 = sum1 + array1(i,j)*array1(i,j)
        end do
      end do
      omega = sqrt(0.5 * sum1)

c--------------------------------------------------------------------
c  Continue calculating intermediate stuff needed for F_p_n+1
c--------------------------------------------------------------------

      call aa_dot_bb(3,array1,array1,array2)

      do i = 1,3
        if (omega .ne. 0.0) then   ! if omega=0 then no need.
          do j = 1,3
            array1(i,j) = array1(i,j)*sin(omega)/omega +
     &       array2(i,j)*(1.0-cos(omega))/omega**2
          end do
        end if
        array1(i,i) = 1.0 + array1(i,i)
      end do

c--------------------------------------------------------------------
c   Finally multiply arrays to get F_p_inv at end of time step.
c--------------------------------------------------------------------

      call inverse_3x3(array1,array2,isingular)
      if (isingular .eq. 1) return
      call aa_dot_bb(3,F_p_inv_0,array2,F_p_inv)

c--------------------------------------------------------------------
c  Multiply F() by F_p_inv() to get F_el()
c--------------------------------------------------------------------

      call aa_dot_bb(3,F1,F_p_inv,F_el)
      if (determinant(F_el) .eq. 0.0e0) then
        F_el(1,1) = 1.0e0
        F_el(2,2) = 1.0e0
        F_el(3,3) = 1.0e0
      end if
      call inverse_3x3(F_el,F_el_inv,isingular)
      if (isingular .eq. 1) return

c--------------------------------------------------------------------
c  Rotate director vectors from intermediate configuration to
c  the current configuration.
c--------------------------------------------------------------------

      do n = 1,num_slip_sys
        do i = 1,3
          xs(i,n) = 0.0
          xm(i,n) = 0.0
          do j = 1,3
            xs(i,n) = xs(i,n) + F_el(i,j) * xs0(j,n)
            xm(i,n) = xm(i,n) + xm0(j,n)  * F_el_inv(j,i)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate elastic Green Strain
c--------------------------------------------------------------------

      call transpose1(3,F_el,array1)
      call aa_dot_bb(3,array1,F_el,E_el)
      do i = 1,3
        E_el(i,i) = E_el(i,i) - 1
        do j = 1,3
          E_el(i,j) = E_el(i,j)/2
        end do
      end do

c--------------------------------------------------------------------
c  Multiply the stiffness tensor by the Green strain to get
c  the 2nd Piola Kirkhhoff stress
c--------------------------------------------------------------------

      call aaaa_dot_dot_bb(3,C,E_el,Spk2)

c--------------------------------------------------------------------
c  Convert from PK2 stress to Cauchy stress
c--------------------------------------------------------------------

      det = determinant(F_el)
      call transpose1(3,F_el,array2)
      call aa_dot_bb(3,F_el,Spk2,array1)
      call aa_dot_bb(3,array1,array2,sig)
      do i = 1,3
        do j = 1,3
          sig(i,j) = sig(i,j)/det
        end do
      end do

c--------------------------------------------------------------------
c  Calculate resolved shear stress for each slip system.
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        tau(k) = 0.0
        do j = 1,3
          do i = 1,3
            tau(k) = tau(k) + xs(i,k)*xm(j,k)*sig(i,j)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate dislocation densities, rho_m, rho_i
c--------------------------------------------------------------------
!       checkpoint("Update dislocation densities")

      do ia = 1,num_slip_sys
        sum1 = 0.0e0
        do ib = 1,num_slip_sys
          sum1 = sum1 + H(ia,ib)*(rho_m0(ib) + rho_i0(ib))
        end do
        d_disl(ia) = 1.0e0/sqrt(sum1)
      end do

      do ia = 1,num_slip_sys
        sum1 = 0.0e0
        do ib = 1,num_slip_sys
          sum1 = sum1 + H(ia,ib)*(rho_m0(ib) + rho_i0(ib))
        end do
        drhomdt(ia) = (q_gen/b/d_disl(ia) - q_ann*R_c*rho_m0(ia)/b
     &   - (x_d*sqrt(sum1))/b)*abs(gamma_dot(ia))

        drhoidt(ia) = ((x_d*sqrt(sum1))/b - q_dyn*
     &   rho_i0(ia))*abs(gamma_dot(ia))
      end do

      do ia = 1,num_slip_sys
        rho_m(ia) = rho_m0(ia) + drhomdt(ia)*dtime
        !print*, drhomdt(ia)
        if (rho_m(ia) .le. 1.0e4) then
          rho_m(ia) = 1.0e4
c          print*, "mobile", ia
        end if
        rho_i(ia) = rho_i0(ia) + drhoidt(ia)*dtime
        if (rho_i(ia) .le. 1.0e5) then
          rho_i(ia) = 1.0e5
c          print*, "immobile", ia
        end if
      end do

c--------------------------------------------------------------------
c  Calculate the athermal slip resistance for each slip system.
c--------------------------------------------------------------------
!       checkpoint("Update slip resistances")
      do ia = 1,num_slip_sys
        sum1 = 0.0
        do ib = 1,num_slip_sys
          sum1 = sum1 + p0*A(ia,ib)*(rho_m(ib) + rho_i(ib))
        end do
        s_a(ia) = tau0 + hp_coeff/sqrt(grain_size) + q_t*G*b*sqrt(sum1)
      end do

c--------------------------------------------------------------------
c  Calculate thermal slip resistance for each slip system.
c--------------------------------------------------------------------

      do ia = 1,num_slip_sys
        s_t(ia) = drag_stress
      end do

c--------------------------------------------------------------------
c  Calculate new slip rate for each slip system.
c--------------------------------------------------------------------
!       checkpoint("Calculate new shearing rates")
      do k = 1,num_slip_sys
!         refgdg = gammadot0g*(rho_m(k)/rho_m_zero)*
!      &   (d_disl(k)/d_disl_zero)
        refgdg = 1.0e8 * rho_m(k) * b * d_disl(k)

        if (abs(tau(k)) .lt. s_a(k)) then
          shr_g(k) = refgdg * exp((-delF0/B_k/xtemp)*
     &     power(1.0 - power(abs(tau(k))/s_a(k),p),q))*
     &     sgn(tau(k))

        else if (abs(tau(k)) .ge. s_a(k)) then
          shr_g(k) = refgdg * sgn(tau(k))
        end if

      end do

c--------------------------------------------------------------------
c  Calculate function values
c--------------------------------------------------------------------
!       checkpoint("Calculating function values")
      do k = 1,num_slip_sys
        func(k) = gamma_dot(k) - shr_g(k)
c        print*, func(k)
c        gamma_dot(k) = shr_g(k)
      end do

c--------------------------------------------------------------------
c  Calculate root mean square error that is to be minimized to zero
c--------------------------------------------------------------------
!       checkpoint("Calculating rms error")
      sse = 0.0
      gamma_dot_max = abs(gamma_dot(1))
      do k = 1,num_slip_sys
        sse = sse + abs(gamma_dot(k))*(func(k)**2)
c        print*, k
        if (abs(gamma_dot(k)) .gt. gamma_dot_max) then
          gamma_dot_max = abs(gamma_dot(k))
        end if
c        print*, k, "end"
      end do
!       checkpoint("Updating sse")
      if (gamma_dot_max .gt. 0.0) then
        sse = sqrt(sse/gamma_dot_max)/num_slip_sys
      end if
!       checkpoint("End eval_func()")
      return
      end

c====================================================================
c====================================================================
c
c  Calculate the determinant of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      function determinant(a)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

      determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3

      return
      end

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     subroutine voigt_abq   ---->   version of
c
c     transforms 6x1 matrix t1 into second order tensor t2 if iopt=1
c     and viceversa if iopt=2.
c     transforms 6x6 matrix c2 into fourth order tensor c4 if iopt=3
c     and viceversa if iopt=4.
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine voigt_abq(t1,t2,c2,c4,iopt)

      implicit double precision (a-h,o-z)

      dimension t1(6),t2(3,3),c2(6,6),c4(3,3,3,3)
      dimension ijv(6,2)
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/

      if(iopt.eq.1) then
      do 30 i=1,6
      i1=ijv(i,1)
      i2=ijv(i,2)
      t2(i1,i2)=t1(i)
   30 t2(i2,i1)=t1(i)
      endif
c
      if(iopt.eq.2) then
      do 40 i=1,6
      i1=ijv(i,1)
      i2=ijv(i,2)
   40 t1(i)=t2(i1,i2)
      endif
c
      if (iopt.eq.3) then
      do 10 i=1,6
      i1=ijv(i,1)
      i2=ijv(i,2)
      do 10 j=1,6
      j1=ijv(j,1)
      j2=ijv(j,2)
      c4(i1,i2,j1,j2)=c2(i,j)
      c4(i2,i1,j1,j2)=c2(i,j)
      c4(i1,i2,j2,j1)=c2(i,j)
   10 c4(i2,i1,j2,j1)=c2(i,j)
      endif
c
      if(iopt.eq.4) then
      do 20 i=1,6
      i1=ijv(i,1)
      i2=ijv(i,2)
      do 20 j=1,6
      j1=ijv(j,1)
      j2=ijv(j,2)
   20 c2(i,j)=c4(i1,i2,j1,j2)
      endif
c
      return
      end  ! end voigt_abq

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     subroutine voigt_abq   ---->   version of
c
c     transforms 6x1 matrix t1 into second order tensor t2 if iopt=1
c     and viceversa if iopt=2.
c     transforms 6x6 matrix c2 into fourth order tensor c4 if iopt=3
c     and viceversa if iopt=4.
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine voigt_abq_strain(t1,t2,c2,c4,iopt)

      implicit double precision (a-h,o-z)

      dimension t1(6),t2(3,3),c2(6,6),c4(3,3,3,3)
      dimension ijv(6,2),fact(6,6)
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,1,2,1,3,2,3/
      data ((fact(n,m),n=1,6),m=1,6)/1,1,1,2,2,2,1,1,1,2,2,2,
     1  1,1,1,2,2,2,2,2,2,4,4,4,2,2,2,4,4,4,2,2,2,4,4,4/
      if(iopt.eq.1) then
         do 30 i=1,3
            i1=ijv(i,1)
            i2=ijv(i,2)
            t2(i1,i2)=t1(i)
            t2(i2,i1)=t1(i)
 30      continue
         do 31 i=4,6
            i1=ijv(i,1)
            i2=ijv(i,2)
            t2(i1,i2)=t1(i)/2.
            t2(i2,i1)=t1(i)/2.
 31      continue
      endif
c
      if(iopt.eq.2) then
         do 40 i=1,3
            i1=ijv(i,1)
            i2=ijv(i,2)
            t1(i)=t2(i1,i2)
 40      continue
         do 41 i=4,6
            i1=ijv(i,1)
            i2=ijv(i,2)
            t1(i)=t2(i1,i2)*2.
 41      continue
      endif
c
      if (iopt.eq.3) then
         do 10 i=1,6
            i1=ijv(i,1)
            i2=ijv(i,2)
      do 10 j=1,6
         j1=ijv(j,1)
         j2=ijv(j,2)
            c4(i1,i2,j1,j2)=c2(i,j)*(1/fact(i,j))
c            write(*,*)'cijkl=cij',i1,i2,j1,j2
c            write(*,*)'cij',i,j
c            write(*,*)'fac',(1/fact(i,j))
            c4(i2,i1,j1,j2)=c2(i,j)*(1/fact(i,j))
            c4(i1,i2,j2,j1)=c2(i,j)*(1/fact(i,j))
            c4(i2,i1,j2,j1)=c2(i,j)*(1/fact(i,j))
 10   continue
      endif
c
      if(iopt.eq.4) then
         do 20 i=1,6
            i1=ijv(i,1)
            i2=ijv(i,2)
         do 20 j=1,6
            j1=ijv(j,1)
            j2=ijv(j,2)
            c2(i,j)=fact(i,j)*c4(i1,i2,j1,j2)
c            write(*,*)'cij',i,j
c            write(*,*)'cijkl=cij',i1,i2,j1,j2
c            write(*,*)'fac',fact(i,j)

 20      continue
      endif

      return
      end  ! voigt_abq_strain
