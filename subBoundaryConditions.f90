!***********************************************************************            
SUBROUTINE velocityBC
      USE global
      IMPLICIT NONE
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER (KIND = 8):: i, j, k
      REAL (KIND = 8) :: a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4, &
                         a5, b5, c5, a6, b6, c6, a7, b7, c7, a8, b8, c8, &
                         x3, uinlet
      REAL (KIND = 8) :: a00, qinlet, inletArea
      
      GO TO 11
      !new pulsatile fourier 
      
       a00 =       6.133  
       a1 =       1.212  
       b1 =       3.089  
       a2 =      -2.471  
       b2 =       3.284  
       a3 =      -2.453  
       b3 =      -1.821 
       a4 =      -1.082  
       b4 =      -1.669  
       a5 =         1.9 
       b5 =      -1.705  
       a6 =      0.9329  
       b6 =      0.8529  
       a7 =    -0.09047 
       b7 =      0.4561 
       a8 =     -0.3148  
       b8 =      0.1686 
       c1 =       7.181  
       
       x3 = totime-int(totime)

       IF(x3.GT.0.87) x3 = 0.87

       qinlet =a00 + a1*cos(x3*c1) + b1*sin(x3*c1) +                                      &
               a2*cos(2*x3*c1) + b2*sin(2*x3*c1) + a3*cos(3*x3*c1) + b3*sin(3*x3*c1) +   &
               a4*cos(4*x3*c1) + b4*sin(4*x3*c1) + a5*cos(5*x3*c1) + b5*sin(5*x3*c1) +   &
               a6*cos(6*x3*c1) + b6*sin(6*x3*c1) + a7*cos(7*x3*c1) + b7*sin(7*x3*c1) +   &
               a8*cos(8*x3*c1) + b8*sin(8*x3*c1)

       inletArea = 50.26548245743669181540229413247
    
       uinlet    = qinlet*1000/inletArea

     !$acc parallel loop gang vector collapse (2) present (ut, vt, wt, u, v, w, xp, yp, zp, deltax, deltay, deltaz)  &
     !$acc firstprivate (uinlet, uc11, uc22, deltat, nx, y11, y12, z11, z12, y21, y22, z21, z22)
      DO 30 k = 2, nz+1 
      DO 30 j = 2, ny+1 
        !uniform inlet
        ut(1,j,k) = uinlet     
        !ut(1,j,k) = 1._rk         
        vt(1,j,k) = -vt(2,j,k)  
        wt(1,j,k) = -wt(2,j,k)
        !Neumann - low Re convective flow (outlet)
        !ut(nx+1,j,k)  =  ut(nx,j,k)                                              
        !vt(nx+2,j,k)  =  vt(nx+1,j,k)
        !wt(nx+2,j,k)  =  wt(nx+1,j,k)	
        !Orlanski - vortex shedding Re Convective flow (outlet)		
        !ut(nx+1,j,k)  =  u(nx+1,j,k)-(deltat/deltax(nx+1))*uc*(u(nx+1,j,k)-u(nx,j,k))       
        !vt(nx+2,j,k)  = -vt(nx+1,j,k) + v(nx+2,j,k) + v(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc*(v(nx+2,j,k)-v(nx+1,j,k))
        !wt(nx+2,j,k)  = -wt(nx+1,j,k) + w(nx+2,j,k) + w(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc*(w(nx+2,j,k)-w(nx+1,j,k))

        !outlet1
        IF (yp(j).gt.y11.and.yp(j).lt.y12)THEN
        IF (zp(k).gt.z11.and.zp(k).lt.z12)THEN
        ut(nx+1,j,k)  =  u(nx+1,j,k)-(deltat/deltax(nx+1))*uc11*(u(nx+1,j,k)-u(nx,j,k))       
        vt(nx+2,j,k)  = -vt(nx+1,j,k) + v(nx+2,j,k) + v(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc11*(v(nx+2,j,k)-v(nx+1,j,k))
        wt(nx+2,j,k)  = -wt(nx+1,j,k) + w(nx+2,j,k) + w(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc11*(w(nx+2,j,k)-w(nx+1,j,k))
        END IF
        END IF

        !outlet2
        IF (yp(j).gt.y21.and.yp(j).lt.y22)THEN
        IF (zp(k).gt.z21.and.zp(k).lt.z22)THEN
        ut(nx+1,j,k)  =  u(nx+1,j,k)-(deltat/deltax(nx+1))*uc22*(u(nx+1,j,k)-u(nx,j,k))       
        vt(nx+2,j,k)  = -vt(nx+1,j,k) + v(nx+2,j,k) + v(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc22*(v(nx+2,j,k)-v(nx+1,j,k))
        wt(nx+2,j,k)  = -wt(nx+1,j,k) + w(nx+2,j,k) + w(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc22*(w(nx+2,j,k)-w(nx+1,j,k))
        END IF
        END IF
 30   CONTINUE   
     !$acc end parallel

 

       a1 = 0.5585
	b1 = 3.116
	c1 = 0.1628
	a2 = 0.06637
	b2 = 15.48
	c2 = -1.376
	a3 = 0.06423
	b3 = 27.03
	c3 = -4.713
	a4 = 0.2312
	b4 = 8.127
	c4 = 0.2692
	a5 = 0.1823
	b5 = 31.97
	c5 = 1.055
	a6 = 0.112
	b6 = 33.8
	c6 = 2.843
	a7 = 0.04077
	b7 = 48.69
	c7 = -3.027
	a8 = 0.03597	
	b8 = 51.21
	c8 = 4.455
       x3 = totime-int(totime)   
       uinlet = (a1*sin(b1*x3+c1) + a2*sin(b2*x3+c2) + a3*sin(b3*x3+c3) +  &
              &  a4*sin(b4*x3+c4) + a5*sin(b5*x3+c5) + a6*sin(b6*x3+c6) +  &
              &  a7*sin(b7*x3+c7) + a8*sin(b8*x3+c8))*1000._rk 
     
     !$acc parallel loop gang vector collapse (2) present (ut, vt, wt, u, v, w, xp, yp, zp, deltax, deltay, deltaz)  &
     !$acc firstprivate (uinlet, uc, deltat, nx)
      DO 40 k = 2, nz+1 
      DO 40 j = 2, ny+1 
        !uniform inlet
        ut(1,j,k) = uinlet     
        !ut(1,j,k) = 1._rk         
        vt(1,j,k) = -vt(2,j,k)  
        wt(1,j,k) = -wt(2,j,k)
        !Neumann - low Re convective flow (outlet)
        !ut(nx+1,j,k)  =  ut(nx,j,k)                                              
        !vt(nx+2,j,k)  =  vt(nx+1,j,k)
        !wt(nx+2,j,k)  =  wt(nx+1,j,k)	
        !Orlanski - vortex shedding Re Convective flow (outlet)		
        ut(nx+1,j,k)  =  u(nx+1,j,k)-(deltat/deltax(nx+1))*uc*(u(nx+1,j,k)-u(nx,j,k))       
        vt(nx+2,j,k)  = -vt(nx+1,j,k) + v(nx+2,j,k) + v(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc*(v(nx+2,j,k)-v(nx+1,j,k))
        wt(nx+2,j,k)  = -wt(nx+1,j,k) + w(nx+2,j,k) + w(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc*(w(nx+2,j,k)-w(nx+1,j,k))
 40   CONTINUE   
     !$acc end parallel

 11  CONTINUE
 
     !$acc parallel loop gang vector collapse (2) present (ut, vt, wt, u, v, w, xp, yp, zp, deltax, deltay, deltaz)  &
     !$acc firstprivate (uc, deltat, nx)
      DO 50 k = 2, nz+1 
      DO 50 j = 2, ny+1 
        !uniform inlet
        ut(1,j,k) = 1.    
        !ut(1,j,k) = 1._rk         
        vt(1,j,k) = -vt(2,j,k)  
        wt(1,j,k) = -wt(2,j,k)
        !Neumann - low Re convective flow (outlet)
        ut(nx+1,j,k)  =  ut(nx,j,k)                                              
        vt(nx+2,j,k)  =  vt(nx+1,j,k)
        wt(nx+2,j,k)  =  wt(nx+1,j,k)	
        !Orlanski - vortex shedding Re Convective flow (outlet)		
        !ut(nx+1,j,k)  =  u(nx+1,j,k)-(deltat/deltax(nx+1))*uc*(u(nx+1,j,k)-u(nx,j,k))       
        !vt(nx+2,j,k)  = -vt(nx+1,j,k) + v(nx+2,j,k) + v(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc*(v(nx+2,j,k)-v(nx+1,j,k))
        !wt(nx+2,j,k)  = -wt(nx+1,j,k) + w(nx+2,j,k) + w(nx+1,j,k) - (2._rk*deltat/deltax(nx+2))*uc*(w(nx+2,j,k)-w(nx+1,j,k))
 50   CONTINUE   
     !$acc end parallel

     !$acc parallel loop gang vector collapse (2) present (ut, vt, wt) &
     !$acc firstprivate (ny)
      DO 10 k = 2, nz+1                                                              
      DO 10 i = 2, nx+1
         !ut(i,1,k) = -ut(i,2,k) 
         !wt(i,1,k) = -wt(i,2,k)		
         ut(i,1,k) =  ut(i,2,k)     
         wt(i,1,k) =  wt(i,2,k) 
         vt(i,1,k) =  0._rk  
		 
         !ut(i,ny+2,k) = -ut(i,ny+1,k)                                      !wall no slip - closed channel
         !wt(i,ny+2,k) = -wt(i,ny+1,k)        
         ut(i,ny+2,k) =  ut(i,ny+1,k)                                      !symmetric bc - open channel        
         wt(i,ny+2,k) =  wt(i,ny+1,k)
         vt(i,ny+1,k) =  0._rk    
 10   CONTINUE
     !$acc end parallel

     !$acc parallel loop gang vector collapse (2) present (ut, vt, wt) &
     !$acc firstprivate (nz)      
      DO 20 j = 2, ny+1
      DO 20 i = 2, nx+1	  
        !ut(i,j,1) = -ut(i,j,2) 
        !vt(i,j,1) = -vt(i,j,2)		
         ut(i,j,1) =  ut(i,j,2)     
         vt(i,j,1) =  vt(i,j,2) 
         wt(i,j,1) =  0._rk  
		 
        !ut(i,j,nz+2) = -ut(i,j,nz+1)                                      !wall no slip - closed channel
        !vt(i,j,nz+2) = -vt(i,j,nz+1)        
         ut(i,j,nz+2) =  ut(i,j,nz+1)                                      !symmetric bc - open channel        
         vt(i,j,nz+2) =  vt(i,j,nz+1)
         wt(i,j,nz+1) =  0._rk    
 20   CONTINUE
      !$acc end parallel

      END SUBROUTINE velocityBC

      
!***********************************************************************
      SUBROUTINE solidCellBC
         USE global
         IMPLICIT NONE
         INTEGER, PARAMETER :: rk = selected_real_kind(8)
         INTEGER (KIND = 8):: i, j, k, n    
        !$acc parallel loop gang vector &
        !$acc present(solidIndexPtr, ut, vt, wt, p) &
        !$acc private (i, j, k)
         DO n = 1, solidCellCount
            i = solidIndexPtr(n,1)
            j = solidIndexPtr(n,2)
	     k = solidIndexPtr(n,3)
            ut(i,j,k) = 0._rk
            vt(i,j,k) = 0._rk
            wt(i,j,k) = 0._rk
            p(i,j,k)  = 0._rk
         END DO
        !$acc end parallel
      END SUBROUTINE solidCellBC
!***********************************************************************






