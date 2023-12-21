MODULE K_EPSILON_MODELS

    implicit none

    contains

        !!!*************************************************
        !!!*						         	             *
        !!!*                    ddeta                        *
        !!!*								                 *
        !!!*************************************************
            
        subroutine  ddeta(nr,A,DA,deta)
            implicit none
            integer, intent(in) :: nr
            real*8, intent(in) :: A(1:nr),deta
            real*8, intent(out) :: DA(1:nr)
            integer j

            DA(1) =  -(3.d0*A(1) -4.d0*A(2) +A(3))/(2.d0*deta)

            do j=2,nr-1
                DA(j) = (A(j+1) -A(j-1))/(2.d0*deta)
            enddo

            DA(nr) = (3.d0*A(nr) -4.d0*A(nr-1) +A(nr-2))/(2.d0*deta)

            end

        !!!*************************************************
        !!!*						         	             *
        !!!*                    d2deta2                        *
        !!!*								                 *
        !!!*************************************************

        subroutine  d2deta2(nr,A,D2A,deta)
            implicit none
            integer, intent(in) :: nr
            real*8, intent(in) :: A(1:nr),deta
            real*8, intent(out) :: D2A(1:nr)
            real*8 deta2
            integer j

            deta2 = deta*deta
            
            D2A(1) =  (12.d0*a(1) -30.d0*a(2) +24.d0*a(3) -6.d0*a(4))/(6.d0*deta2)
            
            do j=2,nr-1
                D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/deta2
            enddo
            
            D2A(nr) = (12.d0*a(nr) -30.d0*a(nr-1) +24.d0*a(nr-2) -6.d0*a(nr-3))/(6.d0*deta2)
            
            end

        !!!*************************************************
        !!!*						         	             *
        !!!*            initialize K - epsilon 1D                       *
        !!!*								                 *
        !!!*************************************************

        subroutine initialize_k_epsilon_1D(y_plus,U,Kt,eps,Re_tau,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Re_tau
            real*8, intent(out) :: y_plus(1:ny),U(1:ny),kt(1:ny),eps(1:ny)
            integer j
            real*8 Kappa, Cmu,y_plus_mid

            Kappa = 4.d-1
            Cmu = 9.d-2
            y_plus_mid = 11.635

            do j=1,ny
                if(y_plus(j)<=y_plus_mid) then
                    U(j) = y_plus(j)
                    eps(j) = (y_plus(j)*Kappa -y_plus(j)*y_plus(j)*Kappa/Re_tau)/(Kappa*y_plus_mid)**2.d0 ! 0.1335d0 / ( 1.d0 + ( ( y_plus(j) - 15.515d0 )**2.d0 ) / 166.7634d0 ) !
                    Kt(j) = dsqrt( (y_plus(j)*Kappa -y_plus(j)*y_plus(j)*Kappa/Re_tau)*eps(j)/Cmu ) ! 0.019678d0 * y_plus(j) * y_plus(j) / ( 1.d0 + ( ( y_plus(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
                else
                    U(j) = (1.d0/Kappa)*dlog(y_plus(j)) +5.5d0 
                    eps(j) = (y_plus(j)*Kappa -y_plus(j)*y_plus(j)*Kappa/Re_tau)/(Kappa*y_plus(j))**2.d0 !
                    Kt(j) = dsqrt( (y_plus(j)*Kappa -y_plus(j)*y_plus(j)*Kappa/Re_tau)*eps(j)/Cmu ) ! 0.019678d0 * y_plus(j) * y_plus(j) / ( 1.d0 + ( ( y_plus(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
                endif
            enddo

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Turbulent Reynolds Number K epsilon	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),eps(1:ny)
            real*8, intent(out) :: Ret(1:ny)
            real*8 eps_min
            integer j

            eps_min = 1.d-12

            do j=1,ny
                if (eps(j) <= eps_min) then
                    Ret(j)= dabs(Kt(j)*Kt(j)/eps_min)
                else
                    Ret(j)= dabs(Kt(j)*Kt(j)/eps(j))
                endif
            enddo

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Turbulent Reynolds Number K y_plus	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  turbulent_reynolds_number_k_y_plus(Rey,kt,y_plus,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),y_plus(1:ny)
            real*8, intent(out) :: Rey(1:ny)
            integer j

            do j=1,ny
                Rey(j)= dsqrt(dabs(Kt(j)))*y_plus(j)
            enddo

            end


        !!!***************************************************
        !!!*						         	               *
        !!!*       Nagano and Takawa K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.3d0
            Ce1=1.45d0
            Ce2=1.9d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Abid K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  abid_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.4d0
            Ce1=1.45d0
            Ce2=1.83d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Lam and Bremhorst K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  lam_bremhorst_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.3d0
            Ce1=1.44d0
            Ce2=1.92d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Jones and Launder K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  jones_launder_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.3d0
            Ce1=1.45d0
            Ce2=2.d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Launder and Sharma K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  launder_sharma_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.3d0
            Ce1=1.44d0
            Ce2=1.92d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Yang and Shih K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  yang_shih_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.3d0
            Ce1=1.44d0
            Ce2=1.92d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Abe, Kondoh and Nagano K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  abe_kondoh_nagano_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.4d0
            sigmae= 1.4d0
            Ce1=1.5d0
            Ce2=1.9d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*    Chang, Hsieh and Chen K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  chang_hsieh_chen_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.3d0
            Ce1=1.44d0
            Ce2=1.92d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*    Chien K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  chien_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.0d0
            sigmae= 1.3d0
            Ce1=1.35d0
            Ce2=1.8d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*    Myong and Kasagi K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  myong_kasagi_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.4d0
            sigmae= 1.3d0
            Ce1=1.4d0
            Ce2=1.8d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*    Rodi and Mansour K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  rodi_mansour_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            sigmaK= 1.3d0
            sigmae= 1.3d0
            Ce1=1.44d0
            Ce2=1.92d0
            Cmu=0.09d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  k_epsilon_constants(model,sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            character(len = 17), intent(in) :: model
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            select case (model)
                case ("AB")
                    call abid_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("AKN")
                    call abe_kondoh_nagano_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("CH")
                    call chien_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("CHC")
                    call chang_hsieh_chen_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("JL")
                    call jones_launder_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("LB")
                    call lam_bremhorst_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("LS")
                    call launder_sharma_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("MK")
                    call myong_kasagi_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("NT")
                    call nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("RMM")
                    call rodi_mansour_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case ("YS")
                    call yang_shih_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
                case default
                    print*, ' Model Not Recognised. Defaulting to NT. '
                    call nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            end select

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Nagano Takawa K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  nagano_takawa_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)

            Ret_min = 1.d-12
            do j=1,ny
                if (Ret(j) <= Ret_min) then
                    fmu(j)= (1.d0 +4.1d0/Ret_min**0.75d0)*(1.d0 -dexp(-y_plus(j)/26.d0))**2.d0
                else
                    fmu(j)= (1.d0 +4.1d0/Ret(j)**0.75d0)*(1.d0 -dexp(-y_plus(j)/26.d0))**2.d0
                endif
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -0.3d0*dexp(-(Ret(j)/6.5d0)**2.d0))*(1.d0 -dexp(-y_plus(j)/6.d0))**2.d0
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Abid K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  abid_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), Rey(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)
            call turbulent_reynolds_number_k_y_plus(Rey,kt,y_plus,ny)

            Ret_min = 1.d-12
            do j=1,ny
                if (Ret(j) <= Ret_min) then
                    fmu(j)= dtanh(0.008d0*Rey(j)) * (1.d0 +4.d0/Ret(j)**0.75d0)
                else
                    fmu(j)= dtanh(0.008d0*Rey(j)) * (1.d0 +4.d0/Ret_min**0.75d0)
                endif
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -(2.d0/9.d0)*dexp(-(Ret(j)/6.d0)**2.d0))*(1.d0 -dexp(-y_plus(j)/12.d0))
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Lam and Bremhorst K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  lam_bremhorst_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), Rey(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)
            call turbulent_reynolds_number_k_y_plus(Rey,kt,y_plus,ny)

            Ret_min = 1.d-12
            do j=1,ny
                fmu(j)= ( ( 1.d0 - dexp(-0.0165d0*Rey(j)) )**2.d0 ) * (1.d0 +20.5d0/Ret(j))
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -dexp(-Ret(j)**2.d0))
            enddo

            f1 = 1.d0 + (0.05d0 / fmu)**2.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Jones and Launder K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  jones_launder_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)

            Ret_min = 1.d-12
            do j=1,ny
                fmu(j)= dexp(-2.5d0/(1.d0 + Ret(j)/50.d0))
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -0.3d0*dexp(-Ret(j)**2.d0))
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Launder and Sharma K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  launder_sharma_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)

            Ret_min = 1.d-12
            do j=1,ny
                fmu(j)= dexp(-3.4d0/(1.d0 + Ret(j)/50.d0))
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -0.3d0*dexp(-Ret(j)**2.d0))
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Yang and Shih K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  yang_shih_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), Rey(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)
            call turbulent_reynolds_number_k_y_plus(Rey,kt,y_plus,ny)

            Ret_min = 1.d-12
            do j=1,ny
                fmu(j)= dsqrt( 1.d0 - dexp( -1.5d-04 * Rey(j) -5.d-07 * Rey(j)**3.d0 &
                        -1.5-10 * Rey(j)**5.d0 ) ) * ( 1.d0 + 1.d0 / dsqrt(Ret(j)) )
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= dsqrt( Ret(j) )/( 1.d0 + dsqrt( Ret(j) ) )
            enddo

            f1 = f2

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*   Abe, Kondoh and Nagano K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  abe_kondoh_nagano_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)

            Ret_min = 1.d-12
            do j=1,ny
                fmu(j)= (1.d0 + 5.d0 * dexp(-(Ret(j)/2.d+02)**2.d0) / Ret(j)**0.75d0)*(1.d0 -dexp(-y_plus(j)/14.d0))**2.d0
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -0.3d0*dexp(-(Ret(j)/6.5d0)**2.d0))*(1.d0 -dexp(-y_plus(j)/3.1d0))**2.d0
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*   Chang, Hsieh and Chen K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  chang_hsieh_chen_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), Rey(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)
            call turbulent_reynolds_number_k_y_plus(Rey,kt,y_plus,ny)

            Ret_min = 1.d-12
            do j=1,ny
                if(Ret(j)<=Ret_min) then
                    fmu(j)= ( (1.d0 - dexp( -0.0215d0*Rey(j)))**2.d0 ) * (1.d0 + 31.66d0 / Ret_min**1.25d0)
                else
                    fmu(j)= ( (1.d0 - dexp( -0.0215d0*Rey(j) ))**2.d0 ) * (1.d0 + 31.66d0 / Ret(j)**1.25d0)
                endif
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -0.01d0*dexp(-Ret(j)**2.d0))*(1.d0 -dexp(-0.0631d0*Rey(j)))
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*        Chien K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  chien_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)

            Ret_min = 1.d-12
            do j=1,ny
                fmu(j)= (1.d0 - dexp( -0.0115d0*y_plus(j) ))
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -(2.d0/9.d0)*dexp(-(Ret(j)/6.d0)**2.d0))
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*   Myong and Kasagi K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  myong_kasagi_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)

            Ret_min = 1.d-12
            do j=1,ny
                if(Ret(j)<=Ret_min) then
                    fmu(j)= (1.d0 - dexp( -y_plus(j)/70.d0 )) * (1.d0 + 3.45d0/dsqrt(Ret_min) )
                else
                    fmu(j)= (1.d0 - dexp( -y_plus(j)/70.d0 )) * (1.d0 + 3.45d0/dsqrt(Ret(j)) )
                endif
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            do j=1,ny
                f2(j)= (1.d0 -(2.d0/9.d0)*dexp(-(Ret(j)/6.d0)**2.d0)) * (1.d0 - dexp( -y_plus(j)/5.d0) )**2.d0
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*   Rodi and Mansour K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  rodi_mansour_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,U,detady,Cmu,deta)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),U(1:ny),detady(1:ny),Cmu,deta
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), Rey(1:ny), fmu(1:ny), Ret_min, eps_min, dUdeta(1:ny), dUdy(1:ny), Rep(1:ny)
            integer j

            call turbulent_reynolds_number_k_epsilon(Ret,kt,eps,ny)
            call turbulent_reynolds_number_k_y_plus(Rey,kt,y_plus,ny)

            Ret_min = 1.d-12
            do j=1,ny
                if(y_plus(j)<=100.d0) then
                    fmu(j)= (1.d0 - dexp( -2.d-04*y_plus(j) -6.d-04*y_plus(j)**2.d0 &
                        +2.5d-07*y_plus(j)**3.d0 )) / (1.d0 - dexp(-0.095d0*Rey(j)) )
                else
                    fmu(j)= 1.d0 / (1.d0 - dexp(-0.095d0*Rey(j)) )
                endif
            enddo

            do j=1,ny
                nuT(j)= Cmu*fmu(j)*Ret(j)
            enddo

            call ddeta(ny,U,dUdeta,deta)
            dUdy = dUdeta*detady

            do j=1,ny
                Rep(j) = ( nut(j)*(dUdy(j))**2.d0 ) / (Kt(j) *dsqrt(Cmu*eps(j)))
                f2(j)= (1.d0 -0.22d0*dexp( -0.03357d0*dsqrt(Ret(j)) )) * (1.d0 - dexp( -0.095d0*Rey(j) ) ) + &
                    dexp(1.8d0*Rep(j)**3.d0) - 1.d0
            enddo

            f1 = 1.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  k_epsilon_functions(model,nut,f1,f2,ny,y_plus,kt,eps,U,detady,Cmu,deta)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),U(1:ny),detady(1:ny),Cmu,deta
            character(len = 17), intent(in) :: model
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            select case (model)
                case ("AB")
                    call abid_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("AKN")
                    call abe_kondoh_nagano_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("CH")
                    call chien_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("CHC")
                    call chang_hsieh_chen_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("JL")
                    call jones_launder_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("LB")
                    call lam_bremhorst_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("LS")
                    call launder_sharma_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("MK")
                    call myong_kasagi_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("NT")
                    call nagano_takawa_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case ("RMM")
                    call rodi_mansour_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,U,detady,Cmu,deta)
                case ("YS")
                    call yang_shih_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
                case default
                    print*, ' Model Not Recognised. Defaulting to NT. '
                    call nagano_takawa_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            end select

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Nagano Takawa K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  nagano_takawa_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            D = 0.d0
            E = 0.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Abid K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  abid_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            D = 0.d0
            E = 0.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Lam and Bremhorst K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  lam_bremhorst_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            D = 0.d0
            E = 0.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Jones and Launder K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  jones_launder_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,nut,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: U(1:ny),Kt(1:ny),nut(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 dUdeta(1:ny), d2Udeta2(1:ny), uk(1:ny), dukdeta(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 0.d0
            eps_wall_ny = 0.d0

            uk = dsqrt(Kt)
            call ddeta(ny,uk,dukdeta,deta)
            
            D = 2.d0 * (dukdeta * detady)**2.d0

            call ddeta(ny,U,dUdeta,deta)
            call d2deta2(ny,U,d2Udeta2,deta)

            E = 2.d0 * nut * ( d2Udeta2*detady**2.d0 + dUdeta*d2etady2 )**2.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Launder and Sharma K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  launder_sharma_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,nut,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: U(1:ny),Kt(1:ny),nut(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 dUdeta(1:ny), d2Udeta2(1:ny), uk(1:ny), dukdeta(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 0.d0
            eps_wall_ny = 0.d0

            uk = dsqrt(Kt)
            call ddeta(ny,uk,dukdeta,deta)
            
            D = 2.d0 * (dukdeta * detady)**2.d0

            call ddeta(ny,U,dUdeta,deta)
            call d2deta2(ny,U,d2Udeta2,deta)

            E = 2.d0 * nut * ( d2Udeta2*detady**2.d0 + dUdeta*d2etady2 )**2.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Yang and Shih K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  yang_shih_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,nut,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: U(1:ny),Kt(1:ny),nut(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 dUdeta(1:ny), d2Udeta2(1:ny), uk(1:ny), dukdeta(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            D = 0.d0

            call ddeta(ny,U,dUdeta,deta)
            call d2deta2(ny,U,d2Udeta2,deta)

            E = nut * ( d2Udeta2*detady**2.d0 + dUdeta*d2etady2 )**2.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*  Abe, Kondoh and Nagano K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  abe_kondoh_nagano_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            D = 0.d0
            E = 0.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Chang, Hsieh and Chen K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  chang_hsieh_chen_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            D = 0.d0
            E = 0.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Chien K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  chien_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,eps,y_plus,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),eps(1:ny),y_plus(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 0.d0
            eps_wall_ny = 0.d0

            D = 2.d0 * Kt / y_plus**2.d0
            E = -( 2.d0 * eps / y_plus**2.d0 ) * dexp(-0.5d0 * y_plus)

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       Myong and Kasagi K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  myong_kasagi_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            D = 0.d0
            E = 0.d0

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*   Rodi and Mansour K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  rodi_mansour_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,eps,nut,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: U(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),detady(1:ny),d2etady2(1:ny),deta
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 dUdeta(1:ny), d2Udeta2(1:ny), dKtdeta(1:ny), Ret_min, eps_min
            integer j

            eps_wall_1 = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2))) &
                        -dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
            eps_wall_ny = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1))) &
                        -dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) )*detady(ny) )**2.d0

            call ddeta(ny,Kt,dKtdeta,deta)
            call ddeta(ny,U,dUdeta,deta)
            
            D = 0.d0

            call ddeta(ny,U,dUdeta,deta)
            call d2deta2(ny,U,d2Udeta2,deta)

            E = 1.2d0 * nut * ( d2Udeta2*detady**2.d0 + dUdeta*d2etady2 )**2.d0 + &
                7.5d-03 * (Kt / eps) *(dKtdeta * detady) * (dUdeta * detady) * &
                ( d2Udeta2*detady**2.d0 + dUdeta*d2etady2 )

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*       K - Epsilon D,E,eps_wall (1D)	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  D_E_epsilon_wall_1D(model,D,E,eps_wall_1,eps_wall_ny,Kt,eps,f2,y_plus,nut,U,detady,d2etady2,deta,Ce2,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny),y_plus(1:ny),detady(1:ny),d2etady2(1:ny),f2(1:ny),Ce2,deta
            character(len = 17), intent(in) :: model
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            select case (model)
                case ("AB")
                    call abid_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
                case ("AKN")
                    call abe_kondoh_nagano_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
                case ("CH")
                    call chien_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,eps,y_plus,detady,d2etady2,deta,ny)
                case ("CHC")
                    call chang_hsieh_chen_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
                case ("JL")
                    call jones_launder_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,nut,detady,d2etady2,deta,ny)
                case ("LB")
                    call lam_bremhorst_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
                case ("LS")
                    call launder_sharma_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,nut,detady,d2etady2,deta,ny)
                case ("MK")
                    call myong_kasagi_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
                case ("NT")
                    call nagano_takawa_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
                case ("RMM")
                    call rodi_mansour_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,eps,nut,detady,d2etady2,deta,ny)
                case ("YS")
                    call yang_shih_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,U,Kt,nut,detady,d2etady2,deta,ny)
                case default
                    print*, ' Model Not Recognised. Defaulting to NT. '
                    call nagano_takawa_D_E_epsilon_wall_1D(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            end select

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*                U coefficients	       	   *
        !!!*								                   *
        !!!***************************************************
        subroutine  u_coefficients(aU_w,aU_e,sU,nut,dnutdy,deta,Re_tau,d2etady2,detady,r,is_cartesian)
            implicit none
            real*8, intent(in) :: nut,dnutdy,deta,Re_tau,d2etady2,detady,r
            LOGICAL, INTENT(IN) :: is_cartesian
            real*8, intent(out) :: aU_w,aU_e,sU
            real*8 dev

            if(is_cartesian) then
                dev = deta*( d2etady2 + dnutdy*detady/(1.d0+nut) )/(4.d0*detady**2.d0)
                sU = (deta*deta)/(2.d0*Re_tau*(1.d0+nut)*(detady)**2.d0)
            else 

                if(r>0.d0) then
                    dev = deta*( d2etady2 + (1.d0/r + dnutdy/(1.d0+nut))*detady )/(4.d0*detady**2.d0)
                    sU = (deta*deta)/(Re_tau*(1.d0+nut)*(detady)**2.d0)
                else 
                    dev = (deta / 4.d0) * d2etady2 / (detady)**2.d0
                    sU = (deta*deta)/(2.d0*Re_tau*(1.d0+nut)*detady**2.d0)
                endif

            endif

            aU_w = 5.d-1 - dev
            aU_e = 5.d-1 + dev
            
            end

        !!!***************************************************
        !!!*						         	               *
        !!!*                K coefficients	       	   *
        !!!*								                   *
        !!!***************************************************
        subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdy,dUdy,D,deta,sigmak,d2etady2,detady,r,is_cartesian)
            implicit none
            real*8, intent(in) :: eps,nut,dnutdy,dUdy,D,deta,sigmak,d2etady2,detady,r
            LOGICAL, INTENT(IN) :: is_cartesian
            real*8, intent(out) :: aK_w,aK_e,sK
            real*8 dev

            if(is_cartesian) then
                dev = deta*( d2etady2 + dnutdy/(sigmak+nut)*detady )/(4.d0*detady**2.d0)
                sK = (nut*dUdy*dUdy - D - eps)*(deta*deta)/(2.d0*(1.d0+nut/sigmak)*detady**2.d0)
            else
                if(r>0.d0) then
                    dev = deta*( d2etady2 + (1.d0/r + dnutdy/(sigmak+nut))*detady )/(4.d0*detady**2.d0)
                    sK = (nut*dUdy*dUdy - eps)*(deta*deta)/(2.d0*(1.d0+nut/sigmak)*detady**2.d0)
                else
                    dev = (deta / 4.d0) * d2etady2 / (detady)**2.d0
                    sK = - (D+eps)*(deta*deta)/(4.d0*(1.d0+nut/sigmak)*detady**2.d0)
                endif
            endif

            aK_w = 5.d-1 - dev
            aK_e = 5.d-1 + dev
            

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*                E coefficients	       	   *
        !!!*								                   *
        !!!***************************************************
        subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdy,dUdy,E,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2, &
                    detady,r,is_cartesian,in_source)
            implicit none
            real*8, intent(in) :: eps,Kt,nut,dnutdy,dUdy,E,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady,r
            real*8, intent(out) :: aE_w,aE_e,sE
            logical, INTENT(IN) :: is_cartesian,in_source
            real*8 K_min, Kb, dev
            
            if(is_cartesian) then
                dev = deta*( (sigmae+nut)*d2etady2 + dnutdy*detady )/(4.d0*(sigmae+nut)*detady**2.d0)
                Kb = (Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps)*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detady**2.d0))
            else 
                if(r>0.d0) then
                    dev = deta*( d2etady2 + (1.d0/r + dnutdy/(sigmae+nut))*detady )/(4.d0*detady**2.d0)
                    Kb = (Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps)*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detady**2.d0))
                else
                    dev = (deta / 4.d0) * d2etady2 / (detady)**2.d0
                    Kb = -Ce2*f2*eps*(deta*deta/(4.d0*(1.d0+nut/sigmae)*detady**2.d0))
                endif

            endif

            K_min = 1.d-60

            if (in_source) then
                aE_w = 5.d-1 - dev
                aE_e = 5.d-1 + dev
                if (Kt<=K_min) then
                    sE = Kb*eps/K_min 
                else
                    sE = Kb*eps/Kt
                endif
            else
                aE_w = (5.d-1 - dev)/(1.d0 - Kb/Kt)
                aE_e = (5.d-1 + dev)/(1.d0 - Kb/Kt)
                sE = 0.d0 
            endif

            sE = sE + E*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detady**2.d0))

            end
        !!!*************************************************
        !!!*						         	           *
        !!!*       Calculate U and nut Gradient 1D                 *
        !!!*								               *
        !!!*************************************************

        subroutine  u_nut_gradient_1D(dUdy,dnutdy,U,nut,deta,detady,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: deta,nut(1:ny),U(1:ny),detady(1:ny)
            real*8, INTENT(OUT) :: dUdy(1:ny),dnutdy(1:ny)
            real*8 dUdeta(1:ny),dnutdeta(1:ny)

            call ddeta(ny,U,dUdeta,deta)
            dUdy = dUdeta*detady
            call ddeta(ny,nut,dnutdeta,deta)
            dnutdy = dnutdeta*detady

            end

        !!!*************************************************
        !!!*						         	           *
        !!!*          Budget k - epsilon                 *
        !!!*								               *
        !!!*************************************************

        subroutine  budget_k_epsilon(ny,U,Kt,eps,nut,r,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps, &
                    epseps, detady, d2etady2, is_cartesian)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: deta,sigmaK,sigmaE,Ce1,Ce2,f1(1:ny),f2(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny)
            real*8, intent(in) :: detady(1:ny),d2etady2(1:ny),r(1:ny)
            LOGICAL, INTENT(IN) :: is_cartesian
            real*8, INTENT(OUT) :: tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny)
            real*8, INTENT(OUT) :: Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny)
            real*8 dUdy(1:ny),d2Ktdeta2(1:ny),d2epsdeta2(1:ny),dKtdeta(1:ny),depsdeta(1:ny),dnutdy(1:ny)
            real*8 dUdeta(1:ny),dnutdeta(1:ny)

            call u_nut_gradient_1D(dUdy,dnutdy,U,nut,deta,detady,ny)

            tau_mu = dUdy
            tau_R = nut*dUdy
            Pk = tau_R*dUdy

            call d2deta2(ny,Kt,d2Ktdeta2,deta)
            call ddeta(ny,Kt,dKtdeta,deta)

            Dk = d2Ktdeta2*detady**2.d0 + dKtdeta*d2etady2
            
            if(.NOT. is_cartesian) then
                if(r(1)==0.d0) then
                    Dk(1) = Dk(1) + d2Ktdeta2(1)*detady(1)**2.d0 + dKtdeta(1)*d2etady2(1)
                    Dk(2:ny) = Dk(2:ny) + dKtdeta(2:ny) * detady(2:ny) / r(2:ny)
                    Deps(1) = Deps(1) + d2epsdeta2(1)*detady(1)**2.d0 + depsdeta(1)*d2etady2(1)
                    Deps(2:ny) = Deps(2:ny) + depsdeta(2:ny) * detady(2:ny) / r(2:ny)
                else if(r(ny)==0.d0) then
                    Dk(ny) = Dk(ny) + d2Ktdeta2(ny)*detady(ny)**2.d0 + dKtdeta(ny)*d2etady2(ny)
                    Dk(1:ny-1) = Dk(1:ny-1) + dKtdeta(1:ny-1) * detady(1:ny-1) / r(1:ny-1)
                    Deps(ny) = Deps(ny) + d2epsdeta2(ny)*detady(ny)**2.d0 + depsdeta(ny)*d2etady2(ny)
                    Deps(1:ny-1) = Deps(1:ny-1) + depsdeta(1:ny-1) * detady(1:ny-1) / r(1:ny-1)
                else
                    Dk = Dk + dKtdeta * detady / r
                    Deps = Deps + depsdeta * detady / r
                endif
            endif

            Tk = (nut/sigmaK)*Dk + (dKtdeta*detady/sigmaK)*dnutdy
            Teps = (nut/sigmaE)*Deps + (depsdeta*detady/sigmaE)*dnutdy

            if(Kt(1)==0.d0 .AND. Kt(ny)/=0.d0) then
                Peps(2:ny) = f1(2:ny)*Ce1*(eps(2:ny)/Kt(2:ny))*Pk(2:ny)
                Peps(1) = Peps(2)
                epsEps(2:ny) = -f2(2:ny)*Ce2*(eps(2:ny)/Kt(2:ny))*eps(2:ny)
                epsEps(1) = epsEps(2)
            else if(Kt(ny)==0.d0 .AND. Kt(1)/=0.d0) then
                Peps(1:ny-1) = f1(1:ny-1)*Ce1*(eps(1:ny-1)/Kt(1:ny-1))*Pk(1:ny-1)
                Peps(ny) = Peps(ny-1)
                epsEps(1:ny-1) = -f2(1:ny-1)*Ce2*(eps(1:ny-1)/Kt(1:ny-1))*eps(1:ny-1)
                epsEps(ny) = epsEps(ny-1)
            else if(Kt(ny)==0.d0 .AND. Kt(1)==0.d0) then
                Peps(2:ny-1) = f1(2:ny-1)*Ce1*(eps(2:ny-1)/Kt(2:ny-1))*Pk(2:ny-1)
                Peps(1) = Peps(2)
                Peps(ny) = Peps(ny-1)
                epsEps(2:ny-1) = -f2(2:ny-1)*Ce2*(eps(2:ny-1)/Kt(2:ny-1))*eps(2:ny-1)
                epsEps(1) = epsEps(2)
                epsEps(ny) = epsEps(ny-1)
            else
                Peps = f1*Ce1*(eps/Kt)*Pk
                epsEps = -f2*Ce2*(eps/Kt)*eps
            endif

            call d2deta2(ny,eps,D2epsdeta2,deta)
            call ddeta(ny,eps,depsdeta,deta)

            end

        !!!*************************************************
        !!!*						         	           *
        !!!*                 Solve U  1D                      *
        !!!*								               *
        !!!*************************************************
            
        subroutine  solve_u_1D(U,nut,dnutdy,detady,d2etady2,deta,Re_tau,r,ny,is_cartesian,is_thomas,is_wall_1,is_wall_ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),r(1:ny),deta,Re_tau
            LOGICAL, INTENT(IN) :: is_cartesian,is_thomas,is_wall_1, is_wall_ny
            real*8, intent(inout) :: U(1:ny)
            real*8 aU_w,aU_e,sU
            real*8 A(1:ny),C_apex(1:ny),denominator
            integer j

            if(is_wall_1) then
                U(1) = 0.d0
            else
                j=1
                call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j),r(j),is_cartesian)
                U(1) =  sU + (aU_e + aU_w)*U(2)
            endif

            if(is_wall_ny) then
                U(ny) = 0.d0
            else
                j=ny
                call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j),r(j),is_cartesian)
                U(ny) =  sU + (aU_e + aU_w)*U(ny-1)
            endif

            if(is_thomas) then

                if(is_wall_1) then
                    call u_coefficients(aU_w,aU_e,sU,nut(2),dnutdy(2),deta,Re_tau,d2etady2(2),detady(2),r(2),is_cartesian)
                    A(2) = aU_e
                    C_apex(2) = sU + aU_w * U(1)
                else
                    call u_coefficients(aU_w,aU_e,sU,nut(1),dnutdy(1),deta,Re_tau,d2etady2(1),detady(1),r(2),is_cartesian)                   
                    A(1) = aU_e + aU_w
                    C_apex(1) = sU
                    call u_coefficients(aU_w,aU_e,sU,nut(2),dnutdy(2),deta,Re_tau,d2etady2(2),detady(2),r(2),is_cartesian)
                    denominator = ( 1.d0 - aU_w * A(1) )
                    A(2) = aU_e / denominator
                    C_apex(2) = ( aU_w * C_apex(1) + sU ) / denominator
                endif

                do j =3,ny
                    call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j),r(j),is_cartesian)
                    denominator = ( 1.d0 - aU_w * A(j-1) )
                    A(j) = aU_e / denominator
                    C_apex(j) = ( aU_w * C_apex(j-1) + sU ) / denominator
                enddo

                do j =ny-1,2,-1
                    U(j) = A(j) * U(j+1) + C_apex(j)
                enddo

            else

                do j =2,ny-1
                    call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j),r(j),is_cartesian)
                    U(j) =  sU + aU_e*U(j+1) + aU_w*U(j-1)
                enddo

            endif


            end

        !!!*************************************************
        !!!*						         	           *
        !!!*                 Solve Kt 1D                       *
        !!!*								               *
        !!!*************************************************
            
        subroutine  solve_Kt_1D(model,Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,D,r,deta,sigmak,ny,is_cartesian,is_thomas, &
            is_wall_1,is_wall_ny)
            implicit none
            integer, intent(in) :: ny
            character(len = 17), intent(in) :: model
            real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),r(1:ny),D(1:ny)
            real*8, intent(in) :: deta,sigmak
            LOGICAL, INTENT(IN) :: is_cartesian,is_thomas,is_wall_1, is_wall_ny
            real*8, intent(inout) :: Kt(1:ny)
            real*8 aK_w,aK_e,sK
            real*8 A(1:ny),C_apex(1:ny),denominator
            integer j

            if(is_wall_1) then
                Kt(1) = 0.d0
            else
                j=1
                call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),D(j),deta,sigmak,d2etady2(j),detady(j),r(j), &
                is_cartesian)
                Kt(1) = sK + (aK_e + aK_w)*Kt(2)
            endif

            if(is_wall_ny) then
                Kt(ny) = 0.d0
            else
                j=ny
                call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),D(j),deta,sigmak,d2etady2(j),detady(j),r(j), &
                is_cartesian)
                Kt(ny) = sK + (aK_e + aK_w)*Kt(ny-1)
            endif

            if(is_thomas) then

                if(is_wall_1) then
                    call K_coefficients(aK_w,aK_e,sK,eps(2),nut(2),dnutdy(2),dUdy(2),D(2),deta,sigmak,d2etady2(2),detady(2),r(2),&
                    is_cartesian)
                    A(2) = aK_e
                    C_apex(2) = sK + aK_w * Kt(1)
                else 
                    call K_coefficients(aK_w,aK_e,sK,eps(1),nut(1),dnutdy(1),dUdy(1),D(1),deta,sigmak,d2etady2(1),detady(1),r(1), &
                    is_cartesian)                   
                    A(1) = aK_e + aK_w
                    C_apex(1) = sK
                    call K_coefficients(aK_w,aK_e,sK,eps(2),nut(2),dnutdy(2),dUdy(2),D(2),deta,sigmak,d2etady2(2),detady(2),r(2), &
                    is_cartesian)
                    denominator = ( 1.d0 - aK_w * A(1) )
                    A(2) = aK_e / denominator
                    C_apex(2) = ( aK_w * C_apex(1) + sK ) / denominator
                endif
                do j =3,ny
                    call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),D(j),deta,sigmak,d2etady2(j),detady(j),r(j), &
                    is_cartesian)
                    denominator = ( 1.d0 - aK_w * A(j-1) )
                    A(j) = aK_e / denominator
                    C_apex(j) = ( aK_w * C_apex(j-1) + sK ) / denominator
                    !print*, ' A(j) = ', A(j),' C_apex(j) = ', C_apex(j)
                enddo

                do j =ny-1,2,-1
                    Kt(j) = A(j) * Kt(j+1) + C_apex(j)
                enddo

            else
                
                do j =2,ny-1
                    call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),D(j),deta,sigmak,d2etady2(j),detady(j),r(j), &
                    is_cartesian)
                    Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
                enddo
                
            endif

            end

        !!!*************************************************
        !!!*						         	           *
        !!!*                 Solve Eps  1D                      *
        !!!*								               *
        !!!*************************************************
            
        subroutine  solve_eps_1D(model,Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,E,r,deta,sigmae,ce1,ce2,f1,f2,eps_wall_1, &
                    eps_wall_ny, ny, is_cartesian,is_thomas,is_wall_1,is_wall_ny,in_source)
            implicit none
            integer, intent(in) :: ny
            character(len = 17), intent(in) :: model
            real*8, intent(in) :: Kt(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmae
            real*8, intent(in) :: E(1:ny),r(1:ny),ce1,ce2,f1(1:ny),f2(1:ny),eps_wall_1,eps_wall_ny
            LOGICAL, INTENT(IN) :: is_cartesian,is_thomas,is_wall_1, is_wall_ny, in_source
            real*8, intent(inout) :: eps(1:ny)
            real*8 aE_w,aE_e,sE
            real*8 A(1:ny),C_apex(1:ny),denominator
            integer j

            if(is_wall_1) then
                eps(1) = eps_wall_1
            else 
                j=1
                call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),E(j),deta,sigmae,Ce1,f1(j),Ce2,f2(j), &
                d2etady2(j), detady(j),r(j),is_cartesian, in_source)
                eps(1) = sE + (aE_e + aE_w)*eps(2)
            endif

            if(is_wall_ny) then
                eps(ny) = eps_wall_ny 
            else 
                j=ny
                call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),E(j),deta,sigmae,Ce1,f1(j),Ce2,f2(j), &
                d2etady2(j), detady(j),r(j),is_cartesian, in_source)
                eps(ny) = sE + (aE_e + aE_w)*eps(ny-1)
            endif

            if(is_thomas) then

                if(is_wall_1) then
                    call E_coefficients(aE_w,aE_e,sE,eps(2),Kt(2),nut(2),dnutdy(2),dUdy(2),E(2),deta,sigmae,Ce1,f1(2),Ce2,f2(2), &
                        d2etady2(2),detady(2),r(2),is_cartesian, in_source)
                    A(2) = aE_e
                    C_apex(2) = sE + aE_w * eps(1)
                else 
                    call E_coefficients(aE_w,aE_e,sE,eps(1),Kt(1),nut(1),dnutdy(1),dUdy(1),E(1),deta,sigmae,Ce1,f1(1),Ce2,f2(1), &
                        d2etady2(1),detady(1),r(1),is_cartesian, in_source)                   
                    A(1) = aE_e + aE_w
                    C_apex(1) = sE
                    call E_coefficients(aE_w,aE_e,sE,eps(2),Kt(2),nut(2),dnutdy(2),dUdy(2),E(2),deta,sigmae,Ce1,f1(2),Ce2,f2(2), &
                        d2etady2(2),detady(2),r(2),is_cartesian, in_source)
                    denominator = ( 1.d0 - aE_w * A(1) )
                    A(2) = aE_e / denominator
                    C_apex(2) = ( aE_w * C_apex(1) + sE ) / denominator
                endif
                do j =3,ny
                    call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),E(j),deta,sigmae,Ce1,f1(j),Ce2,f2(j), &
                    d2etady2(j),detady(j),r(j),is_cartesian, in_source)
                    denominator = ( 1.d0 - aE_w * A(j-1) )
                    A(j) = aE_e / denominator
                    C_apex(j) = ( aE_w * C_apex(j-1) + sE ) / denominator
                enddo
            
                do j =ny-1,2,-1
                    eps(j) = A(j) * eps(j+1) + C_apex(j)
                enddo  

            else

                do j =2,ny-1
                    call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),E(j),deta,sigmae,Ce1,f1(j),Ce2,f2(j), &
                    d2etady2(j), detady(j),r(j),is_cartesian, in_source)
                    eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
                enddo 

            endif

            end


    end module 