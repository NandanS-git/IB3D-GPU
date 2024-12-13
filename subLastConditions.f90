
      SUBROUTINE lastConditions     
       USE global
       implicit none 
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER::  i, j, k
        WRITE(*,*) 'Enter lastcondtitions'
        u_sum = 0._rk
        v_sum = 0._rk
        w_sum = 0._rk
        p_sum = 0._rk
        u_avg = 0._rk
        v_avg = 0._rk
        w_avg = 0._rk
        p_avg = 0._rk
        resi_u = 0._rk
        resi_v = 0._rk
        resi_w = 0._rk
        SUMWSS = 0._rk            !SUMWSS global real array(nsurf)           
        SIGNWSS = 0._rk           !SIGNWSS global real array(nsurf)                   
        ita = 0
        ita1 = 0
       
       OPEN (1, FILE='result', FORM='formatted')      
 	DO 30 k = 1, nz+2
       DO 30 j = 1, ny+2
       DO 30 i = 1, nx+2       
         READ(1,*) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), totime, ita, ita1
 30    CONTINUE 
 	CLOSE(1)
 
 	DO 40 k = 1, nz+2
       DO 40 j = 1, ny+2
       DO 40 i = 1, nx+2       
          ut(i,j,k) = u(i,j,k)
          vt(i,j,k) = v(i,j,k)
          wt(i,j,k) = w(i,j,k) 
 40    CONTINUE  

      !$acc update device(u, v, w, ut, vt, wt, p, resi_u, resi_v, resi_w)

      END SUBROUTINE lastConditions
      
      
      
      
