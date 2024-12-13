!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE poissonSolver
        USE global
        IMPLICIT NONE
        INTEGER(KIND=8) :: i, j,k, n
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        REAL (KIND = 8)    :: dalt, derr1, derr2, div, dab, dudt, dvdt, dwdt, derrStdSt
                  
         ALLOCATE ( b(nx+2, ny+2, nz+2) )
               
         derrStdSt = 0._rk

         !initialize variables
         nIterPcor = 0  
         divmax = 0._rk
         dalt   = 0._rk         
         derr1  = 0._rk

        !$acc parallel loop gang vector collapse (3) present (b, pc, pco)
        DO k = 1, nz+2
        DO j = 1, ny+2
        DO i = 1, nx+2
           b(i,j,k)  = 0.
           pc(i,j,k) = 0.
           pco(i,j,k)= 0.
        END DO
        END DO
        END DO
        !$acc end parallel
              
         CALL computeDiv    !divergence vector  

         CALL cpu_time(dStart)
       
         CALL REDBLACKSOR(eps1, nIterPcor, derr1, derr2) 

         CALL cpu_time(dfinish)

         CALL correctPressure !pressure correction 
           
         CALL correctVelocity !velocity correction   

         CALL velocityBC      !correct velocity at boundaries

        !$acc parallel loop gang vector firstprivate (deltat)   &
        !$acc private (i, j, k, dudt, dvdt, dwdt)               &
        !$acc present (fluidIndexPtr, ut, u, vt, v, wt, w)      &
        !$acc reduction (max: derrStdSt)
         DO n = 1, fluidCellCount
           i = fluidIndexPtr(n, 1)
           j = fluidIndexPtr(n, 2)  
           k = fluidIndexPtr(n, 3)  
           dudt = dabs((ut(i,j,k) - u(i,j,k)))/deltat
           dvdt = dabs((vt(i,j,k) - v(i,j,k)))/deltat
	    dwdt = dabs((wt(i,j,k) - w(i,j,k)))/deltat
           derrStdSt = dmax1(derrStdSt, dudt, dvdt, dwdt)    
         ENDDO
	 !$acc end parallel

         OPEN(111,FILE='iter.dat',ACCESS='Append',STATUS='unknown')
         WRITE(111,16) ita, nIterPcor, divmax, derrStdSt, dfinish-dstart
 16      FORMAT(' ',I8, I10, 4E15.6)
         CLOSE(111)	 
         
         WRITE(*,16) ita, nIterPcor, divmax, derrStdSt, dfinish-dstart

        !$acc parallel loop gang vector present (u, ut, v, vt, w, wt) collapse (3)
         DO k = 1, nz+2
         DO j = 1, ny+2
         DO i = 1, nx+2
           u(i,j,k) = ut(i,j,k)
           v(i,j,k) = vt(i,j,k)
           w(i,j,k) = wt(i,j,k)
         END DO
         END DO
         END DO
        !$acc end parallel

         DEALLOCATE (b)
      END SUBROUTINE poissonSolver
!***********************************************************************
      SUBROUTINE computeDiv
         USE global
         IMPLICIT NONE
         INTEGER :: n, i, j, k

         divmax=0.

        !$acc parallel loop gang vector private (i, j, k)   &
        !$acc present (fluidIndexPtr, b, ut, vt, wt, deltax, deltay, deltaz) reduction (max:divmax)
         DO 70 n = 1, fluidCellCount
           i = fluidIndexPtr(n, 1)
           j = fluidIndexPtr(n, 2)  
           k = fluidIndexPtr(n, 3)  		   
           b(i,j,k) =  (ut(i,j,k) - ut(i-1,j,k))/deltax(i) +  &
                       (vt(i,j,k) - vt(i,j-1,k))/deltay(j) +  &
                       (wt(i,j,k) - wt(i,j,k-1))/deltaz(k)
            divmax=dmax1(abs(b(i,j,k)), divmax)
 70      CONTINUE
        !$acc end parallel
      END SUBROUTINE computeDiv

!***********************************************************************

      SUBROUTINE correctPressure
         USE global
         IMPLICIT NONE
         INTEGER(KIND=8) :: n, i, j, k
         REAL :: r1p, r2p

        !$acc parallel loop gang vector private (i, j, k)   &
        !$acc present (fluidIndexPtr, p, pc)
         DO 30 n = 1, fluidCellCount 
            i = fluidIndexPtr(n, 1)
            j = fluidIndexPtr(n, 2) 
            k = fluidIndexPtr(n, 3) 			
            p(i,j,k) = p(i,j,k) + pc(i,j,k) 
 30      CONTINUE
        !$acc end parallel   
      END SUBROUTINE correctPressure


!***********************************************************************

      SUBROUTINE correctVelocity
         USE global
         IMPLICIT NONE
         INTEGER(KIND=8) :: n, i, j, k
   
        !$acc parallel loop gang vector private (i, j, k) firstprivate (deltat) &
        !$acc present (fluidIndexPtr, ut, vt, wt, deltax, deltay, deltaz, pc)
         DO 30 n = 1, fluidCellCount 
            i = fluidIndexPtr(n, 1)
            j = fluidIndexPtr(n, 2) 
            k = fluidIndexPtr(n, 3) 			
			
            ut(i,j,k) = ut(i,j,k) - deltat/(0.5d0*(deltax(i+1)+deltax(i)))*(pc(i+1,j,k)-pc(i,j,k))
            vt(i,j,k) = vt(i,j,k) - deltat/(0.5d0*(deltay(j+1)+deltay(j)))*(pc(i,j+1,k)-pc(i,j,k))  
            wt(i,j,k) = wt(i,j,k) - deltat/(0.5d0*(deltaz(k+1)+deltaz(k)))*(pc(i,j,k+1)-pc(i,j,k))  			
 30      CONTINUE
         !$acc end parallel
      END SUBROUTINE correctVelocity
      
!***********************************************************************

!***********************************************************************
          
      SUBROUTINE REDBLACKSOR(epsi, isum, derr, derr2)      
         USE global
         IMPLICIT NONE
         INTEGER, PARAMETER :: rk = selected_real_kind(8)   
         INTEGER(KIND=8) :: n, i, j, k
         REAL (KIND = 8) :: omega, derr3, errSum
         REAL (KIND = 8), INTENT(IN)     :: epsi
         REAL (KIND = 8), INTENT(OUT)    :: derr, derr2
         INTEGER (KIND = 8), INTENT(OUT) :: isum
         isum = 0   
         omega = 1.99 ! 1.667_rk      
 3       isum = isum + 1  
         derr = 0._rk 
         derr2 = 0._rk
         errSum = 0._rk
		 
        !$acc parallel loop gang vector present(Acx, Acy, Acz, pc, b, pco, redCellIndexPtr) firstprivate(deltat, omega) private (i, j, k) 
         DO 10 n = 1, redCellCount 
             i = redCellIndexPtr(n, 1)
             j = redCellIndexPtr(n, 2)  
             k = redCellIndexPtr(n, 3)  

            pc(i,j,k) = (b(i,j,k)/deltat &
     &        -Acx(i-1,1)*pco(i-1,j,k) -Acx(i-1,3)*pco(i+1,j,k) &
     &        -Acy(j-1,1)*pco(i,j-1,k) -Acy(j-1,3)*pco(i,j+1,k) &
     &        -Acz(k-1,1)*pco(i,j,k-1) -Acz(k-1,3)*pco(i,j,k+1))/(Acx(i-1,2)+Acy(j-1,2)+Acz(k-1,2)) 
            pc(i,j,k) = (1._rk-omega)*pco(i,j,k) + omega*pc(i,j,k)
 10      CONTINUE 
        !$acc end parallel  
		 
        !$acc parallel loop gang vector present(Acx, Acy, Acz, pc, b, pco, blackCellIndexPtr) firstprivate(deltat, omega) private (i, j, k)
         DO 20 n = 1, blackCellCount   
             i = blackCellIndexPtr(n, 1)
             j = blackCellIndexPtr(n, 2)  
             k = blackCellIndexPtr(n, 3) 

            pc(i,j,k) = (b(i,j,k)/deltat &
     &        -Acx(i-1,1)*pc(i-1,j,k) -Acx(i-1,3)*pc(i+1,j,k) &
     &        -Acy(j-1,1)*pc(i,j-1,k) -Acy(j-1,3)*pc(i,j+1,k) &
     &        -Acz(k-1,1)*pc(i,j,k-1) -Acz(k-1,3)*pc(i,j,k+1))/(Acx(i-1,2)+Acy(j-1,2)+Acz(k-1,2)) 
            pc(i,j,k) = (1._rk-omega)*pco(i,j,k) + omega*pc(i,j,k)

 20      CONTINUE          
        !$acc end parallel  
              
        !$acc parallel loop gang vector reduction(max:derr2) present(pc, pco, fluidIndexPtr) private (i, j, k) 
         DO 30 n = 1, fluidCellCount   
             i = fluidIndexPtr(n, 1)
             j = fluidIndexPtr(n, 2) 
             k = fluidIndexPtr(n, 3) 			 
            derr2 = dmax1(derr2,abs(pc(i,j,k)-pco(i,j,k)))     
            !errSum = errSum + (pc(i,j,k)-pco(i,j,k))**2       
            pco(i,j,k) = pc(i,j,k)  
 30      CONTINUE 
        !$acc end parallel
		
         !derr = dsqrt(errSum/fluidCellCount) 
         !derr3 = dmin1(derr, derr2)
         !WRITE(*,*) isum, derr2
         !IF (mod(isum,2000).EQ.0)WRITE(*,*) isum, derr, derr2
         !IF (ita.LE.2.AND.isum.LT.50000) GOTO 3
         !IF (derr.gt.epsi) GOTO 3
         IF (derr2.GE.epsi) GOTO 3
         !IF (ita.lt.15.AND.isum.lt.100) GOTO 3         
      END SUBROUTINE REDBLACKSOR
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss    
