MODULE K_EPSILON_MODELS

    implicit none

    contains

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
        !!!*       Turbulent Reynolds Number 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  turbulent_reynolds_number(Ret,kt,eps,ny)
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
        !!!*       Nagano Takawa K - Epsilon Constants 	       	   *
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
        !!!*       K - Epsilon Constants 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  k_epsilon_constants(model,sigmak,sigmae,Ce1,Ce2,Cmu)
            implicit none
            character(len = 17), intent(in) :: model
            real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu

            select case (model)
            case ("NT")
                call nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
            case ("AB")
                call abid_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu)
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

            call turbulent_reynolds_number(Ret,kt,eps,ny)

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
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            call turbulent_reynolds_number(Ret,kt,eps,ny)

            Ret_min = 1.d-12
            do j=1,ny
                fmu(j)= dtanh(0.008d0*Ret(j)) * (1.d0 +4.d0/Ret(j)**0.75d0)
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
        !!!*       K - Epsilon Functions 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  k_epsilon_functions(model,nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: y_plus(1:ny),Kt(1:ny),eps(1:ny),Cmu
            character(len = 17), intent(in) :: model
            real*8, intent(out) :: nut(1:ny),f1(1:ny),f2(1:ny)
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            select case (model)
            case ("NT")
                call nagano_takawa_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
            case ("AB")
                call abid_k_epsilon_functions(nut,f1,f2,ny,y_plus,kt,eps,Cmu)
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

        subroutine  nagano_takawa_D_E_epsilon_wall(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
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
        !!!*       Nagano Takawa K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  abid_D_E_epsilon_wall(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
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
        !!!*       K - Epsilon D,E,eps_wall 	       	   *
        !!!*								                   *
        !!!***************************************************

        subroutine  D_E_epsilon_wall(model,D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),detady(1:ny),d2etady2(1:ny),deta
            character(len = 17), intent(in) :: model
            real*8, intent(out) :: D(1:ny),E(1:ny),eps_wall_1,eps_wall_ny
            real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
            integer j

            select case (model)
            case ("NT")
                call nagano_takawa_D_E_epsilon_wall(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            case ("AB")
                call abid_D_E_epsilon_wall(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
            case default
                print*, ' Model Not Recognised. Defaulting to NT. '
                call nagano_takawa_D_E_epsilon_wall(D,E,eps_wall_1,eps_wall_ny,Kt,detady,d2etady2,deta,ny)
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

            call ddeta(ny,U,dUdeta,deta)
            dUdy = dUdeta*detady
            call ddeta(ny,nut,dnutdeta,deta)
            dnutdy = dnutdeta*detady

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
            
        subroutine  solve_Kt_1D(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,D,r,deta,sigmak,ny,is_cartesian,is_thomas, &
            is_wall_1,is_wall_ny)
            implicit none
            integer, intent(in) :: ny
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
            
        subroutine  solve_eps_1D(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,E,r,deta,sigmae,ce1,ce2,f1,f2,eps_wall_1,eps_wall_ny, &
                    ny, is_cartesian,is_thomas,is_wall_1,is_wall_ny,in_source)
            implicit none
            integer, intent(in) :: ny
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