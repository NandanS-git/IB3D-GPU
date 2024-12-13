!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE readInput
       USE global 
       IMPLICIT NONE
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER (KIND=8) :: i, j , k  
        OPEN(60, FILE = 'inputdata', FORM = 'formatted')
        READ(60,*) nx, ny, nz,  &
                   re,      &
                   itamax, eps1, pcItaMax, &                           
                   u0, v0,  w0, &
                   surGeoPoints, a0, freq, aoa, piv_pt,&
                   xShift, yShift, zShift,           &
                   istart, dt_order, inor, &
                   i_startSearch, i_endSearch, &
                   j_startSearch, j_endSearch, &
                   k_startSearch, k_endSearch
        CLOSE(60)  

        !rev = 3.3018868943711498
        !rev = 3.285714  
        rev = 1/re 
        pi = 4.D0*ATAN(1.D0)
        deltat = dt_order*1e-3
        !uc11 = 116.236233
        !uc22 = 116.236233
        !uc = 396.33054782262406
        uc = 1.
        ita = 0
        ita1 = 0
        totime = 0._rk 
        alpha  = 1._rk 
 
        !outlet1 boundaries
        y11 = 0 
        y12 = 8.82
        z11 = 0
        z12 = 10

        !outlet2 boundaries
        y21 = 16.5
        y22 = 23.53
        z21 = 0
        z22 = 10  
  
        print*, 'dt =',  deltat, 'ita = ', ita, 'totime = ', totime
        print*, 'total cells =', nx*ny*nz        
        ALLOCATE( x1(nx+3), y1(ny+3), z1(nz+3), &
                  deltax(nx+2), deltay(ny+2), deltaz(nz+2))
        ALLOCATE( xu(nx+3), yu(ny+2), zu(nz+2), &
                  xv(nx+2), yv(ny+3), zv(nz+2), &
                  xw(nx+2), yw(ny+2), zw(nz+3), &
                  xp(nx+2), yp(ny+2), zp(nz+2))
        
        OPEN(61, FILE = 'gridx.dat', FORM = 'formatted')
        DO i = 2, nx+2
           READ(61, *) x1(i)
        END DO
        CLOSE(61)  

        DO i = 2, nx+1
           deltax(i) = x1(i+1) - x1(i)
        END DO  
        
        deltax(1)    = deltax(2)
        deltax(nx+2) = deltax(nx+1)
        x1(1)        = x1(2) - deltax(1)
        x1(nx+3)     = x1(nx+2) + deltax(nx+2)

        OPEN(62, FILE = 'gridy.dat', FORM = 'formatted')
        DO i = 2, ny+2
           READ(62, *) y1(i)
        END DO
        CLOSE(62)  
        
        DO i = 2, ny+1
           deltay(i) = y1(i+1) - y1(i)
        END DO  
        
        deltay(1)    = deltay(2)
        deltay(ny+2) = deltay(ny+1)
        y1(1)        = y1(2) - deltay(1)
        y1(ny+3)     = y1(ny+2) + deltay(ny+2)
		
		        
        OPEN(63, FILE = 'gridz.dat', FORM = 'formatted')
        DO i = 2, nz+2
           READ(63, *) z1(i)
        END DO
        CLOSE(63)  

        DO i = 2, nz+1
           deltaz(i) = z1(i+1) - z1(i)
        END DO  
        
        deltaz(1)    = deltaz(2)
        deltaz(nz+2) = deltaz(nz+1)
        z1(1)        = z1(2) - deltaz(1)
        z1(nz+3)     = z1(nz+2) + deltaz(nz+2)
        
        DO i = 1, nx+3
           xu(i) = x1(i)
        ENDDO
		
        DO i = 1, ny+3
           yv(i) = y1(i)
        ENDDO
		
	 DO i = 1, nz+3
           zw(i) = z1(i)
        ENDDO
		
        DO i = 1, ny+2
           yu(i) = 0.5_rk*(y1(i)+y1(i+1))
		   yw(i) = yu(i)
           yp(i) = yu(i)
        END DO  
		
        DO i = 1, nx+2
           xv(i) = 0.5_rk*(x1(i)+x1(i+1))
           xw(i) = xv(i)
           xp(i) = xv(i)
        END DO	
		
	 DO i = 1, nz+2
           zu(i) = 0.5_rk*(z1(i)+z1(i+1))
	    zv(i) = zu(i)
           zp(i) = zu(i)
        END DO
        
       !$acc update device (deltax, deltay, deltaz)
       !$acc update device (x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp)
      END SUBROUTINE readInput

!************************************************************************************************
      
      SUBROUTINE readSurfaceMeshGmsh
       USE global 
       IMPLICIT NONE           
       INTEGER (KIND = 8) :: n, i1, i2, i3, i4, i5, i6, i7, gPoints
       CHARACTER (LEN = 72) :: cLine          
             
       OPEN(121, FILE = 'geometries/sphere.msh', form = 'formatted')               !READ SURFACE MESH FILE  
        DO n = 1, 4
           READ (121,*) cLine
        END DO
        READ (121,*) ibNodes  !nsurf=total no. of points in file  
        ALLOCATE ( xnode(ibNodes), ynode(ibNodes), znode(ibNodes) )
        DO n = 1, ibNodes
           READ (121,*) i1, xnode(n), ynode(n), znode(n)
        END DO
        DO n = 1, 2
          READ (121,*) line
        END DO        
        READ (121,*) ibElems   !no. of elements 
        DO n = 1, surGeoPoints
          READ (121,*) cLine
        END DO
        ibElems = ibElems-surGeoPoints
        ALLOCATE ( ibElP1(ibElems), ibElP2(ibElems), ibElP3(ibElems) )
        ibElP1 = 0
        ibElP2 = 0
        ibElP3 = 0
        DO n = 1, ibElems
          READ (121,*) i1, i2, i3, i4, i5, ibElP1(n), ibElP2(n), ibElP3(n)
          !Print*, n, ibElP1(n), ibElP2(n)
        END DO
       !do n = 1, ibNodes
         !print*, n, xnode(n), ynode(n)
       !end do
       CLOSE(121)     
       PRINT *, 'SURFACE MESH READING COMPLETE'
       PRINT *, 'ibNodes =', ibNodes, 'ibElems =', ibElems

       !$acc update device (xnode, ynode, znode, ibElP1, ibElP2, ibElP3)
      END SUBROUTINE readSurfaceMeshGmsh
!***********************************************************************

!***********************************************************************      
      SUBROUTINE readSurfaceMeshGambit
       USE global 
       IMPLICIT NONE           
       INTEGER (KIND = 8) :: n, i1, i2, i3, i4
       CHARACTER (LEN = 72) :: cLine          
             
       OPEN(121, FILE = 'geometries/bend_ellipse_ha.neu', form = 'formatted')               !READ SURFACE MESH FILE  
        DO n = 1, 6
           READ (121,*) cLine
        END DO
        READ (121,*) ibNodes, ibElems, i1, i2, i3, i4
        ALLOCATE ( xnode(ibNodes), ynode(ibNodes), znode(ibNodes) )
        DO n = 1, 2
           READ (121,*) cLine
        END DO
        DO n = 1, ibNodes
           READ (121,*) i1, xnode(n), ynode(n), znode(n)
        END DO
        DO n = 1, 2
          READ (121,*) cLine
        END DO        
        ALLOCATE ( ibElP1(ibElems), ibElP2(ibElems), ibElP3(ibElems) )
        ibElP1 = 0
        ibElP2 = 0
        ibElP3 = 0
        DO n = 1, ibElems
          READ (121,*) i1, i2, i3, ibElP1(n), ibElP2(n), ibElP3(n)
        END DO
       CLOSE(121)     
       PRINT *, 'SURFACE MESH READING COMPLETE'
       PRINT *, 'ibNodes =', ibNodes, 'ibElems =', ibElems

       !$acc update device (xnode, ynode, znode, ibElP1, ibElP2, ibElP3)       
      END SUBROUTINE readSurfaceMeshGambit
!*******************************************************************
