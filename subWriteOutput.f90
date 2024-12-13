!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE writeOutput
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8) 
       CHARACTER*70  filename1
       CHARACTER*70  filename2
       INTEGER  :: k, i, j   
       REAL (KIND = 8) :: u1, v1, w1      
         IF(mod(ita,5000).EQ.0)THEN        
            WRITE(filename1,1) ita
 1          FORMAT('out/fielddata.',i9.9,".dat")
            OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables= "x","y","z","u","v","w","p","totime","cellid"'	    	     
            WRITE(786,*) 'zone, ', 'i = ', nx,' j = ', ny, ' k = ', nz   
            !k = nz/2
            DO 30 k = 2, nz+1
            DO 30 j = 2, ny+1
            DO 30 i = 2, nx+1
               u1 = 0.5*(u(i,j,k)+u(i-1,j,k))
               v1 = 0.5*(v(i,j,k)+v(i,j-1,k))
               w1 = 0.5*(w(i,j,k)+w(i,j,k-1)) 
               WRITE(786,*) xp(i), yp(j), zp(k), u1, v1, w1, p(i,j,k), totime, cell(i,j,k)      
 30         CONTINUE 
            CLOSE(786)      
         ENDIF        
      END SUBROUTINE writeOutput
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE writeOutput1
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       CHARACTER*70  filename1
       CHARACTER*70  filename2
       INTEGER  :: k, i, j
       REAL (KIND = 8) :: u1, v1, w1
         IF(ita.eq.62500)THEN
            WRITE(filename1,1) ita
 1          FORMAT('out/fielddata.',i9.9,".dat")
            OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables="x","y","z","u","v","w","p","totime","cellid"'                     
            WRITE(786,*) 'zone, ', 'i = ', nx,' j = ', ny, ' k = ', nz   
            !k = nz/2
            DO 30 k = 2, nz+1
            DO 30 j = 2, ny+1
            DO 30 i = 2, nx+1
               u1 = 0.5*(u(i,j,k)+u(i-1,j,k))
               v1 = 0.5*(v(i,j,k)+v(i,j-1,k))
               w1 = 0.5*(w(i,j,k)+w(i,j,k-1)) 
               WRITE(786,*) xp(i), yp(j), zp(k), u1, v1, w1, p(i,j,k), totime, cell(i,j,k)      
 30         CONTINUE 
            CLOSE(786)      
         ENDIF        
      END SUBROUTINE writeOutput1
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

!***********************************************************************

!***********************************************************************
      SUBROUTINE writeResult     
        USE global
        implicit none 
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k    
        IF(mod(ita,1000).EQ.0)THEN   
        OPEN (1,FILE='result',FORM='formatted')     
 	 DO 30 k = 1, nz+2
        DO 30 j = 1, ny+2
        DO 30 i = 1, nx+2       
          WRITE(1,*) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), totime, ita, ita1
 30     CONTINUE 
 	 CLOSE(1) 
 	 END IF	
      END SUBROUTINE writeResult
!***********************************************************************

!***********************************************************************
     SUBROUTINE writeTagging
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8) 
       CHARACTER*70  filename1
       CHARACTER*70  filename2
       INTEGER  :: k, i, j                        
            WRITE(filename1,1) ita
 1          FORMAT('tag/tagdata.',i2.2,".dat")
            OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables= "x","y","z","cellid"'	    	     
            WRITE(786,*) 'zone, ', 'i = ', nx,' j = ', ny, ' k = ', nz     
            !k = nz/2
            DO 30 k = 2, nz+1
            DO 30 j = 2, ny+1
            DO 30 i = 2, nx+1
               WRITE(786,*) xp(i), yp(j), zp(k), cell(i,j,k)      
 30         CONTINUE 
            CLOSE(786)             
      END SUBROUTINE writeTagging
!***********************************************************************

!***********************************************************************
      SUBROUTINE computeSumData
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8) 
       INTEGER         :: i, j, k
	     DO 30 k = 1, nz+2
            DO 30 j = 1, ny+2
            DO 30 i = 1, nx+2
               u_sum(i,j,k) = u_sum(i,j,k) + u(i,j,k)
               v_sum(i,j,k) = v_sum(i,j,k) + v(i,j,k)
               w_sum(i,j,k) = w_sum(i,j,k) + w(i,j,k)			   
               p_sum(i,j,k) = p_sum(i,j,k) + p(i,j,k)			                                                            
 30         CONTINUE                 
      END SUBROUTINE computeSumData  
!***********************************************************************

!***********************************************************************
      SUBROUTINE writeComputeSumData    
        USE global
        implicit none 
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k    
        IF(mod(ita1,1000).EQ.0)THEN   
        OPEN (1,FILE='sumdata',FORM='formatted')     
 	 DO 30 k = 1, nz+2
        DO 30 j = 1, ny+2
        DO 30 i = 1, nx+2       
          WRITE(1,*) u_sum(i,j,k), v_sum(i,j,k), w_sum(i,j,k), p_sum(i,j,k), ita, ita1
 30     CONTINUE 
 	 CLOSE(1) 
 	 END IF	
      END SUBROUTINE  writeComputeSumData
!***********************************************************************

!***********************************************************************
      SUBROUTINE readComputeSumData
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  i, j, k
        OPEN (1,FILE='sumdata',FORM='formatted')
        DO 30 k = 1, nz+2
        DO 30 j = 1, ny+2
        DO 30 i = 1, nx+2
          READ(1,*) u_sum(i,j,k), v_sum(i,j,k), w_sum(i,j,k), p_sum(i,j,k), ita, ita1
 30     CONTINUE
        CLOSE(1)
      END SUBROUTINE  readComputeSumData
!***********************************************************************

!***********************************************************************
      SUBROUTINE readstressData
        USE global
        implicit none
        INTEGER, PARAMETER :: rk = selected_real_kind(8)
        INTEGER::  inode
        OPEN (1,FILE='stressdata',FORM='formatted')
        DO inode = 1, ibNodes
          READ(1,*) SUMWSS(inode), SIGNWSS(inode), ita1, ita
        END DO
        CLOSE(1)
      END SUBROUTINE readstressData
!***********************************************************************

!***********************************************************************      
      SUBROUTINE computeAvgData
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8) 
       INTEGER  :: n, k, i, j
	     DO 30 k = 1, nz+2
            DO 30 j = 1, ny+2
            DO 30 i = 1, nx+2              	   
               u_avg(i,j,k) = u_sum(i,j,k)/ita1  
               v_avg(i,j,k) = v_sum(i,j,k)/ita1    
               w_avg(i,j,k) = w_sum(i,j,k)/ita1    			   
               p_avg(i,j,k) = p_sum(i,j,k)/ita1                                                           
 30         CONTINUE                  
      END SUBROUTINE computeAvgData
!***********************************************************************

!***********************************************************************
      SUBROUTINE writeAvgoutput
       USE global
       INTEGER, PARAMETER :: rk = selected_real_kind(8) 
       CHARACTER*70  filename1
       CHARACTER*70  filename2
       INTEGER  :: k, i, j   
       REAL (KIND = 8) :: u_avg1, v_avg1, w_avg1     
            WRITE(filename1,1) ita1
 1          FORMAT('avgout/fielddata.',i9.9,".dat")
            OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables= "x","y","z","u_avg","v_avg","w_avg","p_avg","totime","cellid"'	    	     
            WRITE(786,*) 'zone, ', 'i = ', nx,' j = ', ny, ' k = ', nz
            !k = nz/2
            DO 30 k = 2, nz+1
            DO 30 j = 2, ny+1
            DO 30 i = 2, nx+1
               u_avg1 = 0.5*(u_avg(i,j,k)+u_avg(i-1,j,k))
               v_avg1 = 0.5*(v_avg(i,j,k)+v_avg(i,j-1,k))
               w_avg1 = 0.5*(w_avg(i,j,k)+w_avg(i,j,k-1))
               WRITE(786,*) xp(i), yp(j), zp(k), u_avg1, v_avg1, w_avg1, p_avg(i,j,k), totime, cell(i,j,k)      
 30         CONTINUE 
            CLOSE(786)      
      END SUBROUTINE writeAvgoutput
!***********************************************************************



