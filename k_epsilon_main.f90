MODULE K_EPSILON_MODELS

    implicit none

    contains

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
        !!!*                U coefficients	       	   *
        !!!*								                   *
        !!!***************************************************
        subroutine  u_coefficients(aU_w,aU_e,sU,nut,dnutdy_plus,deta,Re_tau,d2etady2,detady)
            implicit none
            real*8, intent(in) :: nut,dnutdy_plus,deta,Re_tau,d2etady2,detady
            real*8, intent(out) :: aU_w,aU_e,sU
            real*8 dev

            dev = deta*( (1.d0+nut)*d2etady2 + dnutdy_plus*detady )/(4.d0*(1.d0+nut)*(detady)**2.d0)

            aU_w = 5.d-1 - dev
            aU_e = 5.d-1 + dev
            sU = (deta*deta)/(2.d0*Re_tau*(1.d0+nut)*(detady)**2.d0)

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*                K coefficients	       	   *
        !!!*								                   *
        !!!***************************************************
        subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdy_plus,dUdy_plus,D,deta,sigmak,d2etady2,detady)
            implicit none
            real*8, intent(in) :: eps,nut,dnutdy_plus,dUdy_plus,D,deta,sigmak,d2etady2,detady
            real*8, intent(out) :: aK_w,aK_e,sK
            real*8 dev

            dev = deta*( (sigmak+nut)*d2etady2 + dnutdy_plus*detady )/(4.d0*(sigmak+nut)*detady**2.d0)

            aK_w = 5.d-1 - dev
            aK_e = 5.d-1 + dev
            sK = (nut*dUdy_plus*dUdy_plus - D - eps)*(deta*deta)/(2.d0*(1.d0+nut/sigmak)*detady**2.d0)

            end

        !!!***************************************************
        !!!*						         	               *
        !!!*                E coefficients	       	   *
        !!!*								                   *
        !!!***************************************************
        subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdy_plus,dUdy_plus,E,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2, &
                    detady,in_source)
            implicit none
            real*8, intent(in) :: eps,Kt,nut,dnutdy_plus,dUdy_plus,E,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady
            real*8, intent(out) :: aE_w,aE_e,sE
            logical, INTENT(IN) :: in_source
            real*8 K_min, Kb, dev
            
            dev = deta*( (sigmae+nut)*d2etady2 + dnutdy_plus*detady )/(4.d0*(sigmae+nut)*detady**2.d0)

            K_min = 1.d-12

            Kb = (Ce1*f1*nut*dUdy_plus*dUdy_plus + E -Ce2*f2*eps)*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detady**2.d0))

            if (in_source) then
                aE_w = 5.d-1 - dev
                aE_e = 5.d-1 + dev
                if (Kt<=K_min) then
                    sE = 0.d0*Kb*eps/K_min
                else
                    sE = Kb*eps/Kt
                endif
            else
                aE_w = (5.d-1 - dev)/(1.d0 - Kb/Kt)
                aE_e = (5.d-1 + dev)/(1.d0 - Kb/Kt)
                sE = 0.d0 
            endif

            end

        !!!*************************************************
        !!!*						         	           *
        !!!*          Budget k - epsilon                 *
        !!!*								               *
        !!!*************************************************

        subroutine  budget_k_epsilon(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps, &
                    epseps, detady, d2etady2)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: deta,sigmaK,sigmaE,Ce1,Ce2,f1,f2(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny)
            real*8, intent(in) :: detady(1:ny),d2etady2(1:ny)
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
            Tk = (nut/sigmaK)*Dk + (dKtdeta*detady/sigmaK)*dnutdy

            Peps(2:ny) = f1*Ce1*(eps(2:ny)/Kt(2:ny))*Pk(2:ny)
            Peps(1) = Peps(2)

            call d2deta2(ny,eps,D2epsdeta2,deta)
            call ddeta(ny,eps,depsdeta,deta)

            Deps = d2epsdeta2*detady**2.d0 + depsdeta*d2etady2
            Teps = (nut/sigmaE)*Deps + (depsdeta*detady/sigmaE)*dnutdy
            epsEps(2:ny-1) = -f2(2:ny-1)*Ce2*(eps(2:ny-1)/Kt(2:ny-1))*eps(2:ny-1)
            epsEps(1) = epsEps(2)
            epsEps(ny) = epsEps(ny-1)

            end

        !!!*************************************************
        !!!*						         	           *
        !!!*                 Solve U  1D                      *
        !!!*								               *
        !!!*************************************************
            
        subroutine  solve_u_1D(U,nut,dnutdy,detady,d2etady2,deta,Re_tau,ny,is_thomas,is_wall_1,is_wall_ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,Re_tau
            LOGICAL, INTENT(IN) :: is_thomas,is_wall_1, is_wall_ny
            real*8, intent(inout) :: U(1:ny)
            real*8 aU_w,aU_e,sU
            real*8 A(1:ny),C_apex(1:ny),denominator
            integer j

            if(is_thomas) then

                U(1) = 0.d0
                call u_coefficients(aU_w,aU_e,sU,nut(2),dnutdy(2),deta,Re_tau,d2etady2(2),detady(2))
                A(2) = aU_e
                C_apex(2) = sU + aU_w * U(1)
                do j =3,ny-1
                    call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j))
                    denominator = ( 1.d0 - aU_w * A(j-1) )
                    A(j) = aU_e / denominator
                    C_apex(j) = ( aU_w * C_apex(j-1) + sU ) / denominator
                enddo
                !call u_coefficients(aU_w,aU_e,sU,nut(ny),dnutdy(ny),deta,Re_tau,d2etady2(ny),detady(ny))
                !denominator = ( 1.d0 - aU_w * A(ny-1) )
                !A(ny) = aU_e / denominator
                !C_apex(ny) = ( aU_w * C_apex(ny-1) + sU ) / denominator
                A(ny-1) = 0.d0
                U(ny) = 0.d0

                do j =ny-1,2,-1
                    U(j) = A(j) * U(j+1) + C_apex(j)
                enddo

            else

                U(1) = 0.d0
                do j =2,ny-1
                    call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j))
                    U(j) =  sU + aU_e*U(j+1) + aU_w*U(j-1)
                    !!!U(j) = A(j) * U(j+1) + C_apex(j)
                enddo
                U(ny) = 0.d0

            endif


            end

        !!!*************************************************
        !!!*						         	           *
        !!!*                 Solve Kt 1D                       *
        !!!*								               *
        !!!*************************************************
            
        subroutine  solve_Kt_1D(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,D,deta,sigmak,ny,is_thomas,is_wall_1,is_wall_ny)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: eps(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),D(1:ny),deta,sigmak
            LOGICAL, INTENT(IN) :: is_thomas,is_wall_1, is_wall_ny
            real*8, intent(inout) :: Kt(1:ny)
            real*8 aK_w,aK_e,sK
            real*8 A(1:ny),C_apex(1:ny),denominator
            integer j

            if(is_wall_1) then

                Kt(1) = 0.d0
                call K_coefficients(aK_w,aK_e,sK,eps(2),nut(2),dnutdy(2),dUdy(2),D(2),deta,sigmak,d2etady2(2),detady(2))
                A(2) = aK_e
                C_apex(2) = sK + aK_w * Kt(1)
                do j =3,ny-1
                    call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),D(j),deta,sigmak,d2etady2(j),detady(j))
                    denominator = ( 1.d0 - aK_w * A(j-1) )
                    A(j) = aK_e / denominator
                    C_apex(j) = ( aK_w * C_apex(j-1) + sK ) / denominator
                    !print*, ' A(j) = ', A(j),' C_apex(j) = ', C_apex(j)
                enddo
                !call K_coefficients(aK_w,aK_e,sK,eps(ny),nut(ny),dnutdy(ny),dUdy(ny),D(ny),deta,sigmak,d2etady2(ny),detady(ny))
                !denominator = ( 1.d0 - aK_w * A(ny-1) )
                !A(ny) = aK_e / denominator
                !C_apex(ny) = ( aK_w * C_apex(ny-1) + sK ) / denominator
                A(ny-1) = 0.d0
                Kt(ny) = 0.d0

                do j =ny-1,2,-1
                    Kt(j) = A(j) * Kt(j+1) + C_apex(j)
                enddo

            else
                
                Kt(1) = 0.d0
                do j =2,ny-1
                    call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),D(j),deta,sigmak,d2etady2(j),detady(j))
                    Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
                enddo
                Kt(ny) = 0.d0
                
            endif

            end

        !!!*************************************************
        !!!*						         	           *
        !!!*                 Solve Eps  1D                      *
        !!!*								               *
        !!!*************************************************
            
        subroutine  solve_eps_1D(Kt,eps,dUdy,nut,dnutdy,detady,d2etady2,E,deta,sigmae,ce1,ce2,f1,f2,eps_wall_1,eps_wall_ny,ny, &
                    is_thomas,is_wall_1,is_wall_ny,in_source)
            implicit none
            integer, intent(in) :: ny
            real*8, intent(in) :: Kt(1:ny),dUdy(1:ny),nut(1:ny),dnutdy(1:ny),detady(1:ny),d2etady2(1:ny),deta,sigmae
            real*8, intent(in) :: E(1:ny),ce1,ce2,f1(1:ny),f2(1:ny),eps_wall_1,eps_wall_ny
            LOGICAL, INTENT(IN) :: is_thomas,is_wall_1, is_wall_ny, in_source
            real*8, intent(inout) :: eps(1:ny)
            real*8 aE_w,aE_e,sE
            real*8 A(1:ny),C_apex(1:ny),denominator
            integer j

            if(is_thomas) then

                eps(1) = eps_wall_1
                call E_coefficients(aE_w,aE_e,sE,eps(2),Kt(2),nut(2),dnutdy(2),dUdy(2),E(2),deta,sigmae,Ce1,f1(2),Ce2,f2(2), &
                    d2etady2(2),detady(2), in_source)
                A(2) = aE_e
                C_apex(2) = sE + aE_w * eps(1)
                do j =3,ny-1
                    call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),E(j),deta,sigmae,Ce1,f1(j),Ce2,f2(j), &
                    d2etady2(j),detady(j), in_source)
                    denominator = ( 1.d0 - aE_w * A(j-1) )
                    A(j) = aE_e / denominator
                    C_apex(j) = ( aE_w * C_apex(j-1) + sE ) / denominator
                enddo
                !call E_coefficients(aE_w,aE_e,sE,eps(ny),Kt(ny),nut(ny),dnutdy(ny),dUdy(ny),E(ny),deta,sigmae,Ce1,f1(j),Ce2,f2(ny), & 
                !d2etady2(ny),detady(ny), in_source)
                !denominator = ( 1.d0 - aE_w * A(ny-1) )
                A(ny-1) = 0.d0
                eps(ny) = eps_wall_ny
            
                do j =ny-1,2,-1
                    eps(j) = A(j) * eps(j+1) + C_apex(j)
                enddo  

            else

                eps(1) = eps_wall_1
                do j =2,ny-1
                    call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),E(j),deta,sigmae,Ce1,f1(j),Ce2,f2(j), &
                    d2etady2(j), detady(j), in_source)
                    eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
                enddo 
                eps(ny) = eps_wall_ny 

            endif

            end


    end module 