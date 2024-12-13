
      SUBROUTINE initialConditions     
       USE global
       implicit none 
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER::  i, j, k, n
        WRITE(*,*) 'Enter initialcondtitions'
       !cell variables 
        u = 0._rk
        ut = 0._rk
        v = 0._rk
        vt = 0._rk
        w = 0._rk
        wt = 0._rk
        p = 0._rk 
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
        totime = 0.
        DO n = 1, fluidCellCount
           i = fluidIndexPtr(n, 1)
           j = fluidIndexPtr(n, 2) 
	    k = fluidIndexPtr(n, 3) 
           u(i,j,k)  = 1. !396.33054782262406 !116.236233
           v(i,j,k)  = 0._rk
	    w(i,j,k)  = 0._rk
           ut(i,j,k) = 0._rk
           vt(i,j,k) = 0._rk  
	    wt(i,j,k) = 0._rk  
        END DO

       !$acc update device(u, v, w, ut, vt, wt, p, resi_u, resi_v, resi_w)

        print*, 'initial'
         
      END SUBROUTINE initialConditions
      
      
      
      
