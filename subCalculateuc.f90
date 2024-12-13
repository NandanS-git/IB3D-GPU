!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
       SUBROUTINE calculateUc
       USE global
         
       IMPLICIT NONE
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER            :: i, j, k, n, uccount1
       REAL (KIND = 8)    :: uc1, uc0, ucout         
      	
       uc1 = 0._rk
       uccount1 = 0
       i = nx-4

       !$acc parallel loop gang vector collapse (2) reduction(+: uc1, uccount1)  &
       !$acc present (cell, ut) &
       !$acc firstprivate (i)
       DO 30 k = 2, nz+1 
       DO 30 j = 2, ny+1 
       !n = i-1  + nx*(j-2)  + nx*ny*(k-2)
       IF(cell(i,j,k) .eq. 0) THEN
       uc1 = uc1 + ut(i,j,k) 
       uccount1 = uccount1 + 1
       END IF
 30    CONTINUE 
       !$acc end parallel 
  
       uc1 = uc1/uccount1
      	uc0 = 0._rk
      	ucout = dmax1(uc1,uc0)       
 	uc = ucout  
 	  
      	OPEN(98,FILE='uc.dat',ACCESS='Append',STATUS='unknown')
      	WRITE(98,*) ita, totime, uc
	CLOSE(98)       	
!************************ uc calculation end ***************************

!*********************************************************************** 
      END SUBROUTINE calculateUc
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
       SUBROUTINE calculateUcMO
       USE global
         
       IMPLICIT NONE
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER            :: i, j, k, n, uccount1
       REAL (KIND = 8)    :: uc1, uc0, ucout         
      
       !outlet1
       uc1 = 0._rk
       uccount1 = 0     
       i = nx-4

       !$acc parallel loop gang vector collapse (2) reduction(+: uc1, uccount1)  &
       !$acc present (cell, xp, yp, zp, ut) &
       !$acc firstprivate (i, y11, y12, z11, z12)
       DO 30 k = 2, nz+1 
       DO 30 j = 2, ny+1 
       !n = i-1  + nx*(j-2)  + nx*ny*(k-2)     
       IF(yp(j).gt.y11.and.yp(j).lt.y12)THEN
       IF(zp(k).gt.z11.and.zp(k).lt.z12)THEN
       IF(cell(i,j,k) .eq. 0) THEN
       uc1 = uc1 + ut(i,j,k) 
       uccount1 = uccount1 + 1
       END IF
       END IF
       END IF
 30    CONTINUE   
       !$acc end parallel

       uc1 = uc1/uccount1
      	uc0 = 0._rk
      	ucout = dmax1(uc1,uc0)       
     	uc11 = ucout   	  
    
       !outlet2
       uc1 = 0._rk
       uccount1 = 0      
       i = nx-4

       !$acc parallel loop gang vector collapse (2) reduction(+: uc1, uccount1)  &
       !$acc present (cell, xp, yp, zp, ut) &
       !$acc firstprivate (i, y21, y22, z21, z22)
       DO 40 k = 2, nz+1 
       DO 40 j = 2, ny+1 
       !n = i-1  + nx*(j-2)  + nx*ny*(k-2)      
       IF(yp(j).gt.y21.and.yp(j).lt.y22)THEN
       IF(zp(k).gt.z21.and.zp(k).lt.z22)THEN
       IF(cell(i,j,k) .eq. 0) THEN
       uc1 = uc1 + ut(i,j,k) 
       uccount1 = uccount1 + 1
       END IF
       END IF
       END IF
 40    CONTINUE
       !$acc end parallel
   
       uc1 = uc1/uccount1
      	uc0 = 0._rk
      	ucout = dmax1(uc1,uc0)       
     	uc22 = ucout   	  
    
        
       OPEN(100,FILE='uc.dat',ACCESS='Append',STATUS='unknown')
       WRITE(100,*) ita, totime, uc11,uc22
       CLOSE(100)       	
!************************ uc calculation end ***************************

!*********************************************************************** 
      END SUBROUTINE calculateUcMO
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss









