!***********************************************************************
     SUBROUTINE shiftSurfaceNodesInitial
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  i
        REAL (KIND=8)      ::  xr, yr, zr
        ALLOCATE ( xnode1(ibNodes), ynode1(ibNodes), znode1(ibNodes) )
	 
	  !xnode1 = xnode
         !ynode1 = ynode
         !znode1 = znode 
         
        !IF(ita.eq.0) THEN !!For ita equal to 0
	  u_init = 0.
	  u_final = 0.
	  v_init = 0.
	  v_final = 0.
	  w_init = 0.
	  w_final = 0.
	  xmove = 0.
	  ymove = 0.
	  zmove = 0.
	 !END IF

       !$acc parallel loop gang vector present (xnode1, ynode1, znode1, xnode, ynode, znode) &
       !$acc private (xr, yr, zr) firstprivate (xshift, yshift, zshift)
        DO i = 1, ibnodes
           xr = xnode(i)
           yr = ynode(i)
           zr = znode(i)
           xnode1(i) = xr + xShift
           ynode1(i) = yr + yShift 
           znode1(i) = zr + zshift  
        ENDDO
       !$acc end parallel 

       !$acc update host (xnode1, ynode1, znode1)
      END SUBROUTINE shiftSurfaceNodesInitial
      
     SUBROUTINE shiftSurfaceNodesMovingbody
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  i
        REAL (KIND=8)      ::  xr, yr, zr
        REAL (KIND=8)      :: vol
             
	 !IF(ita.gt.0) THEN !!For ita greter than 0
	  vol = (pi*(1.0)**3)/6. 
	
	  u_init = u_final
	  v_init = v_final
	  w_init = w_final
	
	  u_final = (Total_Force_X/(vol*1.))*deltat + u_init
	  v_final = (Total_Force_Y/(vol*1.))*deltat + v_init
	  w_final = (Total_Force_Z/(vol*1.))*deltat + w_init
	
         xmove = u_final*deltat
         ymove = v_final*deltat
         zmove = w_final*deltat
	 !END IF
       !$acc parallel loop gang vector present (xnode1, ynode1, znode1) firstprivate(xmove, ymove, zmove) private(xr, yr, zr)
        DO i = 1, ibnodes
           xr = xnode1(i)
           yr = ynode1(i)
           zr = znode1(i)
           xnode1(i) = xr + xmove 
           ynode1(i) = yr + ymove  
           znode1(i) = zr + zmove   
        ENDDO
       !$acc end parallel

       !$acc update host (xnode1, ynode1, znode1)
      END SUBROUTINE shiftSurfaceNodesMovingbody
!*********************************************************************** 
     
!***********************************************************************  
      SUBROUTINE computeSurfaceNorm
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n  !c1, c2, c3, c4
        REAL (KIND=8)      :: p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL
                             
        ALLOCATE (xcent(ibElems), ycent(ibElems), zcent(ibElems), &
                  cosAlpha(ibElems), cosBeta(ibElems), cosGamma(ibElems), &
                  alpha3(ibElems), beta3(ibElems), gamma3(ibElems))
        !inor = -1._rk
        !compute centroid and direction cosines 
       !$acc parallel loop gang vector       &                                                                                                   
       !$acc present (xnode1, ynode1, znode1, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, alpha3, beta3, gamma3, ibElP1, ibElP2, ibElP3)   &
       !$acc private (p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, lenEL)  firstprivate (inor)           
        DO n = 1, ibElems
           p1x = xnode1(ibElP1(n))                       !x coordinate element node 1
           p1y = ynode1(ibElP1(n))                       !y coordinate element node 1
           p1z = znode1(ibElP1(n))                       !z coordinate element node 1
        
           p2x = xnode1(ibElP2(n))                       !x coordinate element node 2
           p2y = ynode1(ibElP2(n))                       !y coordinate element node 2
           p2z = znode1(ibElP2(n))                       !z coordinate element node 2	   
        
           p3x = xnode1(ibElP3(n))                       !x coordinate element node 3
           p3y = ynode1(ibElP3(n))                       !y coordinate element node 3
           p3z = znode1(ibElP3(n))                       !z coordinate element node 3		   
           
           xcent(n) =  (p2x+p1x+p3x)/3._rk                  !centroid x coordinate element
           ycent(n) =  (p2y+p1y+p3y)/3._rk                  !centroid y coordinate element    
           zcent(n) =  (p2z+p1z+p3z)/3._rk                  !centroid z coordinate element   	

           cosAlpha(n) = (p2y-p1y)*(p3z-p1z)-(p3y-p1y)*(p2z-p1z)  
           cosBeta(n)  = (p2z-p1z)*(p3x-p1x)-(p3z-p1z)*(p2x-p1x)  
           cosGamma(n) = (p2x-p1x)*(p3y-p1y)-(p3x-p1x)*(p2y-p1y)
           
           alpha3(n) = cosAlpha(n)
           beta3(n)  = cosBeta(n)
           gamma3(n) = cosGamma(n)

           lenEL       = dsqrt(cosAlpha(n)**2 + cosBeta(n)**2 + cosGamma(n)**2)   !length of element
  
           cosAlpha(n) = cosAlpha(n)/lenEl*inor               !direction cosine unit normal along x
           cosBeta(n)  = cosBeta(n)/lenEl*inor                !direction cosine unit normal along y
           cosGamma(n) = cosGamma(n)/lenEl*inor               !direction cosine unit normal along z
        ENDDO
       !$acc end parallel

       !$acc update host (xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, alpha3, beta3, gamma3)
        print*, 'SurfaceNorm done, inor =', inor     
       
     END SUBROUTINE computeSurfaceNorm
!******************************************************************

     SUBROUTINE tagging
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, m, i, j, k, n1, n2, n3, n4, nel2n
                                
        INTEGER            :: iPt, iPt1
        REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, &
                              n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
                              cent_x, cent_y, cent_z
                              
        ALLOCATE(minElemcell(nx+2,ny+2,nz+2))

        ibCellCount = 0  
        fluidCellCount = 0
        solidCellCount = 0
        cell = 0
        minElemcell = 0
        n1dotn = 0
        n2dotn = 0
        n3dotn = 0
        n4dotn = 0
        n5dotn = 0
        n6dotn = 0
        n7dotn = 0
        n8dotn = 0
        n9dotn = 0
        
       !$acc parallel loop gang vector collapse(3)        &                                                               
       !$acc present (cosAlpha, cosBeta, cosGamma, xcent, ycent, zcent, cell, minElemcell, x1, y1, z1, xp, yp, zp)       &
       !$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, nel2n, m, cent_x, cent_y, cent_z,         &
       !$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn)                                           
        DO 10 k = k_startSearch, k_endSearch
        DO 10 j = j_startSearch, j_endSearch
        DO 10 i = i_startSearch, i_endSearch
          
           m = 0
           minDis = 1e14
           
           n1x = xp(i)
           n1y = yp(j)
           n1z = zp(k)
		   
           n2x = x1(i)
           n3x = x1(i+1)		   
		   
           n2y = y1(j)
           n3y = y1(j+1)
		   
           n2z = z1(k)
           n3z = z1(k+1)
           
           !$acc loop seq
           DO m = 1, ibElems
              cent_x = xcent(m)
              cent_y = ycent(m)
              cent_z = zcent(m)
              dis  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2) 
              IF (dis.LT.minDis) THEN
                 minDis   = dis                 
                 nel2n    = m
              ENDIF                            
           ENDDO 
           
           minElemcell(i,j,k) = nel2n
           
	     !n1dotn  = (n1x - xcent(nel2n))*cosAlpha(nel2n) + &
                     !(n1y - ycent(nel2n))*cosBeta(nel2n)  + &
                     !(n1z - zcent(nel2n))*cosGamma(nel2n)

	     n2dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n2z - zcent(nel2n))*cosGamma(nel2n)

	     n3dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n2z - zcent(nel2n))*cosGamma(nel2n)

	     n4dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n2z - zcent(nel2n))*cosGamma(nel2n)	
               
	     n5dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n2z - zcent(nel2n))*cosGamma(nel2n)	

            n6dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n3z - zcent(nel2n))*cosGamma(nel2n)

	     n7dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n3z - zcent(nel2n))*cosGamma(nel2n)

	     n8dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n3z - zcent(nel2n))*cosGamma(nel2n)	  
              
	     n9dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                      (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                      (n3z - zcent(nel2n))*cosGamma(nel2n)		
 
           !n       = i-1  + nx*(j-2)  + nx*ny*(k-2)
	         
           IF ( n2dotn.LE.-1e-16 .AND. n3dotn.LE.-1e-16 .AND. n4dotn.LE.-1e-16 .AND. n5dotn.LE.-1e-16  &
               .AND. n6dotn.LE.-1e-16 .AND. n7dotn.LE.-1e-16 .AND. n8dotn.LE.-1e-16 .AND. n9dotn.LE.-1e-16) THEN
		    cell(i,j,k) = 1
           ELSEIF (n2dotn.GT.-1e-16 .AND. n3dotn.GT.-1e-16 .AND. n4dotn.GT.-1e-16 .AND. n5dotn.GT.-1e-16 &
		    .AND. n6dotn.GT.-1e-16 .AND. n7dotn.GT.-1e-16 .AND. n8dotn.GT.-1e-16 .AND. n9dotn.GT.-1e-16) THEN
                  cell(i,j,k) = 0   
           ELSE 
                  cell(i,j,k) = 2               
                             
           ENDIF 
          
 10      CONTINUE  
         !$acc end parallel
	
         !$acc update host(cell)
  	 
         ibCellCount = 0
         solidCellCount = 0
         fluidCellCount = 0
         
        !$acc parallel loop gang vector collapse(3) present(cell) reduction(+: solidCellCount, fluidCellCount, ibCellCount)
         DO 20 k = 2, nz+1
         DO 20 j = 2, ny+1
         DO 20 i = 2, nx+1
            !n       = i-1  + nx*(j-2)  + nx*ny*(k-2)
            IF (cell(i,j,k).eq.1) THEN
               solidCellCount = solidCellCount + 1
            ELSEIF (cell(i,j,k).eq.0) THEN
               fluidCellCount  = fluidCellCount + 1
            ELSEIF (cell(i,j,k).eq.2) THEN
               ibCellCount = ibCellCount + 1
            ENDIF
 20      CONTINUE
        !$acc end parallel
         
!     GOTO 1000
      open(82,file='inter.dat',status='unknown')
      write(82,*)'variables = "x", "y","z", "var"'
	do k = 2, nz+1
       do j = 2, ny+1
       do i = 2, nx+1
	!n = i-1  + nx*(j-2)  + nx*ny*(k-2)
       if(cell(i,j,k).eq.2)then
       write(82,*) xp(i),yp(j), zp(k), 2
       endif
       end do
       end do
	enddo
      close(82)
      
      open(83,file='fluid.dat',status='unknown')
      write(83,*)'variables = "x", "y","z","var"'
       do k = 2, nz+1
       do j = 2, ny+1
       do i = 2, nx+1
	!n = i-1  + nx*(j-2)  + nx*ny*(k-2)
       if(cell(i,j,k).eq.0)then
       write(83,*)xp(i),yp(j), zp(k), 0
       endif
       end do
       end do
	enddo
      close(83)
      
      open(84,file='solid.dat',status='unknown')
      write(84,*)'variables = "x", "y","z","var"'
       do k = 2, nz+1
       do j = 2, ny+1
       do i = 2, nx+1
	!n = i-1  + nx*(j-2)  + nx*ny*(k-2)
       if(cell(i,j,k).eq.1)then
       write(84,*)xp(i),yp(j), zp(k), 1
       endif
       end do
       end do
	enddo
      close(84)

!1000 CONTINUE       
        print*, 'search done'
        Print*, 'imms. cells=', ibCellCount
        Print*, 'fluid cells=', fluidCellCount
        Print*, 'solid cells=', solidCellCount            
     END SUBROUTINE tagging
     
!***********************************************************************

!***********************************************************************

     SUBROUTINE selectiveRetagging
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, m, i, j, k, n1, n2, n3, n4, i1, j1, k1, nn, &
                               nel2n, sumId
        INTEGER            :: iPt, iPt1
        REAL (KIND=8)      :: n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis , minDis, &
                              n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn,  &
                              cent_x, cent_y, cent_z

	 ALLOCATE(minElemcell(nx+2,ny+2,nz+2))
	 
	 minElemcell(:,:,:) = 0
	 n1dotn = 0
        n2dotn = 0
        n3dotn = 0
        n4dotn = 0
        n5dotn = 0
        n6dotn = 0
        n7dotn = 0
        n8dotn = 0
        n9dotn = 0

        !$acc parallel loop gang vector       &                                                                                         
        !$acc present(interceptedIndexPtr, xp, yp, zp, x1, y1, z1, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, cell, minElemcell)  &
        !$acc private(n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, dis, minDis, m, nel2n, cent_x, cent_y, cent_z,   &
        !$acc         n1dotn, n2dotn, n3dotn, n4dotn, n5dotn, n6dotn, n7dotn, n8dotn, n9dotn, i1, j1, k1, i, j, k)
        DO nn = 1, ibCellCount  
            i1 = interceptedIndexPtr(nn, 1)
            j1 = interceptedIndexPtr(nn, 2)
            k1 = interceptedIndexPtr(nn, 3)
            
           !$acc loop collapse(3) seq
            DO 10 k = k1-2, k1+3
            DO 10 j = j1-2, j1+3
            DO 10 i = i1-2, i1+3
               m = 0
	        minDis = 1e14
		 n1x = xp(i)
		 n1y = yp(j)
		 n1z = zp(k)
		  
		 n2x = x1(i)
		 n3x = x1(i+1)		   
			   
		 n2y = y1(j)
		 n3y = y1(j+1)
			   
		 n2z = z1(k)
		 n3z = z1(k+1)
		
		 !$acc loop seq	   
		 DO m = 1, ibElems
		    cent_x = xcent(m)
		    cent_y = ycent(m)
		    cent_z = zcent(m)
		    dis  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2  + (n1z-cent_z)**2) 
		    IF (dis.LT.minDis) THEN
			 minDis   = dis                 
			 nel2n    = m
	           ENDIF                            
	        ENDDO 
	        
	        minElemcell(i,j,k) = nel2n
	        
		 !n1dotn  = (n1x - xcent(nel2n))*cosAlpha(nel2n) + &
               !          (n1y - ycent(nel2n))*cosBeta(nel2n)  + &
               !          (n1z - zcent(nel2n))*cosGamma(nel2n)

		 n2dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n2z - zcent(nel2n))*cosGamma(nel2n)

		 n3dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n2z - zcent(nel2n))*cosGamma(nel2n)

		 n4dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n2z - zcent(nel2n))*cosGamma(nel2n)	  
      	        
		 n5dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n2z - zcent(nel2n))*cosGamma(nel2n)	

		 n6dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n3z - zcent(nel2n))*cosGamma(nel2n)

		 n7dotn  = (n2x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n3z - zcent(nel2n))*cosGamma(nel2n)

		 n8dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n2y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n3z - zcent(nel2n))*cosGamma(nel2n)
	        	        
		 n9dotn  = (n3x - xcent(nel2n))*cosAlpha(nel2n) + &
                         (n3y - ycent(nel2n))*cosBeta(nel2n)  + &
                         (n3z - zcent(nel2n))*cosGamma(nel2n)				
				
		 !n       = i-1  + nx*(j-2)  + nx*ny*(k-2)
				 
	        IF (n2dotn.LE.-1e-16 .AND. n3dotn.LE.-1e-16 .AND. n4dotn.LE.-1e-16 .AND. n5dotn.LE.-1e-16  &
		    .AND. n6dotn.LE.-1e-16 .AND. n7dotn.LE.-1e-16 .AND. n8dotn.LE.-1e-16 .AND. n9dotn.LE.-1e-16) THEN
			   cell(i,j,k) = 1
	        ELSEIF (n2dotn.GT.-1e-16 .AND. n3dotn.GT.-1e-16 .AND. n4dotn.GT.-1e-16 .AND. n5dotn.GT.-1e-16 &
		    .AND. n6dotn.GT.-1e-16 .AND. n7dotn.GT.-1e-16 .AND. n8dotn.GT.-1e-16 .AND. n9dotn.GT.-1e-16) THEN
			   cell(i,j,k) = 0  
	        ELSE 
			   cell(i,j,k) = 2               
		 ENDIF 
 10         CONTINUE
	  ENDDO
        !$acc end parallel

        !$acc update host(cell)

	  ibCellCount = 0  
	  solidCellCount = 0
	  fluidCellCount = 0
	  
        !$acc parallel loop collapse(3) present(cell) reduction(+: solidCellCount, fluidCellCount, ibCellCount)  
         DO 20 k = 2, nz+1		 
         DO 20 j = 2, ny+1
         DO 20 i = 2, nx+1
	     !n       = i-1  + nx*(j-2)  + nx*ny*(k-2)
            IF (cell(i,j,k).eq.1) THEN
		 solidCellCount = solidCellCount + 1
            ELSEIF (cell(i,j,k).eq.0) THEN
               fluidCellCount  = fluidCellCount + 1 
            ELSEIF (cell(i,j,k).eq.2) THEN         
               ibCellCount = ibCellCount + 1 
            ENDIF 
 20      CONTINUE
        !$acc end parallel
 
       GOTO 1001
       
       IF (ita.eq.500) THEN  
	 open(82,file='inter.dat',status='unknown')
        write(82,*)'variables = "x", "y","z", "var"'
	  do k = 2, nz+1
         do j = 2, ny+1
         do i = 2, nx+1
	  !n = i-1  + nx*(j-2)  + nx*ny*(k-2)
         if(cell(i,j,k).eq.2)then
         write(82,*) xp(i),yp(j), zp(k), 2
         endif
         end do
         end do
	  enddo
        close(82)
      
        open(83,file='fluid.dat',status='unknown')
        write(83,*)'variables = "x", "y","z","var"'
         do k = 2, nz+1
         do j = 2, ny+1
         do i = 2, nx+1
	  !n = i-1  + nx*(j-2)  + nx*ny*(k-2)
         if(cell(i,j,k).eq.0)then
         write(83,*)xp(i),yp(j), zp(k), 0
         endif
         end do
         end do
	  enddo
        close(83)
      
        open(84,file='solid.dat',status='unknown')
        write(84,*)'variables = "x", "y","z","var"'
         do k = 2, nz+1
         do j = 2, ny+1
         do i = 2, nx+1
	  !n = i-1  + nx*(j-2)  + nx*ny*(k-2)
         if(cell(i,j,k).eq.1)then
         write(84,*)xp(i),yp(j), zp(k), 1
         endif
         end do
         end do
	  enddo
        close(84)
       ENDIF   
 1001  CONTINUE     
       print*, 'selective retagging', ibCellCount, fluidCellCount, solidCellCount    
     END SUBROUTINE selectiveRetagging
!***********************************************************************

     SUBROUTINE cellCount
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j, k
        
        print*, "cellCount started"
        
         iPt  = 0
         iPt1 = 0
         iPt2 = 0

         ALLOCATE(interceptedIndexPtr(ibCellCount,3),fluidIndexPtr(fluidCellCount, 3), &
         solidIndexPtr(solidCellCount, 3))

        !$acc update host(cell)
         
         DO 30 k = 2, nz+1
         DO 30 j = 2, ny+1
         DO 30 i = 2, nx+1
            !n  = i-1  + nx*(j-2)  + nx*ny*(k-2)
            
            IF (cell(i,j,k).eq.0) THEN
               iPt1 = iPt1 + 1
               fluidIndexPtr(iPt1, 1) = i
               fluidIndexPtr(iPt1, 2) = j
               fluidIndexPtr(iPt1, 3) = k
            ELSEIF (cell(i,j,k).eq.1) THEN
               iPt2 = iPt2 + 1
               solidIndexPtr(iPt2, 1) = i
               solidIndexPtr(iPt2, 2) = j  
               solidIndexPtr(iPt2, 3) = k  			   
            ELSEIF (cell(i,j,k).eq.2) THEN
               iPt = iPt + 1
               interceptedIndexPtr(iPt, 1) = i
               interceptedIndexPtr(iPt, 2) = j
               interceptedIndexPtr(iPt, 3) = k
            ENDIF 
 30      CONTINUE  
         redCellCount = 0
         blackCellCount  = 0

         DO n = 1, fluidCellCount
            i = fluidIndexPtr(n, 1)
            j = fluidIndexPtr(n, 2) 
	     k = fluidIndexPtr(n, 3) 

            IF (mod(i+j+k,2).EQ.1) THEN
               redCellCount = redCellCount + 1
            ELSE
               blackCellCount = blackCellCount + 1
            ENDIF
         ENDDO
         ALLOCATE (redCellIndexPtr(redCellCount,3) ,blackCellIndexPtr(blackCellCount,3))
         ipt1 = 0
         iPt = 0 
         DO n = 1, fluidCellCount
            i = fluidIndexPtr(n, 1)
            j = fluidIndexPtr(n, 2) 
            k = fluidIndexPtr(n, 3) 

            IF (mod(i+j+k,2).EQ.1) THEN
               iPt = iPt + 1
               redCellIndexPtr(iPt, 1) = i
               redCellIndexPtr(iPt, 2) = j   
               redCellIndexPtr(iPt, 3) = k               
			   
            ELSE
               iPt1 = iPt1 + 1
               blackCellIndexPtr(iPt1, 1) = i
               blackCellIndexPtr(iPt1, 2) = j 
               blackCellIndexPtr(iPt1, 3) = k  
			   
            ENDIF 
         ENDDO  

        !$acc update device(interceptedIndexPtr, fluidIndexPtr, solidIndexPtr) 
        !$acc update device(redCellIndexPtr, blackCellIndexPtr)  

            print*, fluidCellCount, redCellCount, blackCellCount
            print*, "cellCount done"
     END SUBROUTINE cellCount
!**************************************************************************

     SUBROUTINE computeNormDistance
        USE global
        IMPLICIT NONE        
        INTEGER, PARAMETER :: rk = selected_real_kind(8)                                
        INTEGER            :: i, j, k, n, nel2p
        REAL (KIND=8)      :: n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z
                                                                         
        CHARACTER*70  filename1
        
    	ALLOCATE(pNormDis(ibCellCount), nelp(ibCellCount), &
        nelu1(ibCellCount), nelu2(ibCellCount), nelv1(ibCellCount), &
        nelv2(ibCellCount), nelw1(ibCellCount), nelw2(ibCellCount), &
        u1NormDis(ibCellCount), u2NormDis(ibCellCount) , v1NormDis(ibCellCount), &
        v2NormDis(ibCellCount) , w1NormDis(ibCellCount), w2NormDis(ibCellCount))
         
        print*, 'computeNormDistance started'

       !$acc parallel loop gang vector  &
       !$acc private (i, j, k, n1x, n2x, n3x, n1y, n2y, n3y, n1z, n2z, n3z, nel2p)          &
       !$acc present (xp, yp, zp, x1, y1, z1, interceptedIndexPtr, minElemcell, nelp, nelu1,          &
       !$acc          nelu2, nelv1, nelv2, nelw1, nelw2, pNormDis, u1NormDis, u2NormDis, v1NormDis,   &
       !$acc          v2NormDis, w1NormDis, w2NormDis, xcent, ycent, zcent, cosAlpha, cosBeta, cosGamma, cell)                                                 
        DO 10 n = 1, ibCellCount
           
           n1x = xp(interceptedIndexPtr(n, 1))
           n2x = x1(interceptedIndexPtr(n, 1))
           n3x = x1(interceptedIndexPtr(n, 1)+1)
           n1y = yp(interceptedIndexPtr(n, 2))           
           n2y = y1(interceptedIndexPtr(n, 2))          
           n3y = y1(interceptedIndexPtr(n, 2)+1)   
           n1z = zp(interceptedIndexPtr(n, 3))           
           n2z = z1(interceptedIndexPtr(n, 3))          
           n3z = z1(interceptedIndexPtr(n, 3)+1)  
                    
           i = interceptedIndexPtr(n, 1) 
           j = interceptedIndexPtr(n, 2) 
           k = interceptedIndexPtr(n, 3) 
           nel2p = minElemcell(i,j,k) 
           nelp(n)  = nel2p
           nelu1(n)  = nel2p
           nelu2(n)  = nel2p
           nelv1(n)  = nel2p
           nelv2(n)  = nel2p
           nelw1(n)  = nel2p
           nelw2(n)  = nel2p
          
           pNormDis(n)  = (n1x - xcent(nel2p))*cosAlpha(nel2p) + &
                          (n1y - ycent(nel2p))*cosBeta(nel2p)  + &
                          (n1z - zcent(nel2p))*cosGamma(nel2p)

           u1NormDis(n) = (n2x - xcent(nel2p))*cosAlpha(nel2p) + &
                          (n1y - ycent(nel2p))*cosBeta(nel2p)  + &
                          (n1z - zcent(nel2p))*cosGamma(nel2p)    
      
           u2NormDis(n) = (n3x - xcent(nel2p))*cosAlpha(nel2p) + &
                          (n1y - ycent(nel2p))*cosBeta(nel2p)  + &
                          (n1z - zcent(nel2p))*cosGamma(nel2p)

           v1NormDis(n) = (n1x - xcent(nel2p))*cosAlpha(nel2p) + &
                          (n2y - ycent(nel2p))*cosBeta(nel2p)  + &
                          (n1z - zcent(nel2p))*cosGamma(nel2p)
           
           v2NormDis(n) = (n1x - xcent(nel2p))*cosAlpha(nel2p) + &
                          (n3y - ycent(nel2p))*cosBeta(nel2p)  + &
                          (n1z - zcent(nel2p))*cosGamma(nel2p)

           w1NormDis(n) = (n1x - xcent(nel2p))*cosAlpha(nel2p) + &
                          (n1y - ycent(nel2p))*cosBeta(nel2p)  + &
                          (n2z - zcent(nel2p))*cosGamma(nel2p)          

           w2NormDis(n) = (n1x - xcent(nel2p))*cosAlpha(nel2p) + &
                          (n1y - ycent(nel2p))*cosBeta(nel2p)  + &
                          (n3z - zcent(nel2p))*cosGamma(nel2p)
 10      CONTINUE
        !$acc end parallel
	
         DEALLOCATE (minElemcell)
         print*, 'computeNormDistance done'
        
     END SUBROUTINE computeNormDistance


