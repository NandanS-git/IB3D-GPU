      SUBROUTINE coefficientMatrix      
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) :: n, i, j, k, iType
        REAL (KIND = 8)    ::  rx1, rx2, rxsum, ry1, ry2, rysum, rz1, rz2, rzsum
         
         DO j = 2, ny+1
            ry1   = yp(j)   - yp(j-1)       
            ry2   = yp(j+1) - yp(j)     
            rysum = ry1 + ry2       
            Acy(j-1, 1) =   2._rk/(ry1*rysum)
            Acy(j-1, 2) =  -2._rk/(ry1*ry2)
            Acy(j-1, 3) =   2._rk/(ry2*rysum) 
         END DO 

         DO i = 2, nx+1   
            rx1   = xp(i)   - xp(i-1)       
            rx2   = xp(i+1) - xp(i)     
            rxsum = rx1 + rx2                  
            Acx(i-1, 1) =   2._rk/(rx1*rxsum)
            Acx(i-1, 2) =  -2._rk/(rx1*rx2)
            Acx(i-1, 3) =   2._rk/(rx2*rxsum)  
         END DO  

         DO i = 2, nz+1   
            rz1   = zp(i)   - zp(i-1)       
            rz2   = zp(i+1) - zp(i)     
            rzsum = rz1 + rz2                  
            Acz(i-1, 1) =   2._rk/(rz1*rzsum)
            Acz(i-1, 2) =  -2._rk/(rz1*rz2)
            Acz(i-1, 3) =   2._rk/(rz2*rzsum)    
         END DO

		 
         !inlet, i = 1
         Acx(1, 2)   =  Acx(1, 2)  + Acx(1, 1)
         Acx(1, 1)   =  0._rk
         !outlet, i = nx
         Acx(nx, 2)  =  Acx(nx, 2) - Acx(nx, 3) 
         Acx(nx, 3)  =  0._rk 
         !bottom, j = 1
         Acy(1, 2)   =  Acy(1, 2)  + Acy(1, 1)
         Acy(1, 1)   =  0._rk
         !top, i = ny
         Acy(ny, 2)  =  Acy(ny, 2) + Acy(ny, 3) 
         Acy(ny, 3)  =  0._rk 
		 !front, k = 1
         Acz(1, 2)   =  Acz(1, 2) + Acz(1, 1)
         Acz(1, 1)   =  0._rk
         !back, k = nz
         Acz(nz, 2)  =  Acz(nz, 2) + Acz(nz, 3) 
         Acz(nz, 3)  =  0._rk 

        !$acc update device (Acx, Acy, Acz)         
        
        print*, "Coefficient Matrix generated"

      END SUBROUTINE
