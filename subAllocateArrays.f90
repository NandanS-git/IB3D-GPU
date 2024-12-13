      SUBROUTINE allocateArrays
        USE global
        IMPLICIT NONE
        ALLOCATE ( u(nx+2,ny+2,nz+2), ut(nx+2,ny+2,nz+2), v(nx+2,ny+2,nz+2), vt(nx+2,ny+2,nz+2),  &
        w(nx+2,ny+2,nz+2), wt(nx+2,ny+2,nz+2), p(nx+2,ny+2,nz+2),   &
	 u_sum(nx+2,ny+2,nz+2), v_sum(nx+2,ny+2,nz+2),  w_sum(nx+2,ny+2,nz+2), p_sum(nx+2,ny+2,nz+2), &
        u_avg(nx+2,ny+2,nz+2), v_avg(nx+2,ny+2,nz+2),  w_avg(nx+2,ny+2,nz+2), p_avg(nx+2,ny+2, nz+2), &
	 resi_u(nx+2,ny+2,nz+2), resi_v(nx+2,ny+2,nz+2), resi_w(nx+2,ny+2,nz+2) )
                 
        ALLOCATE ( cell(nx+2,ny+2,nz+2) )
        ALLOCATE ( Acx(nx,3), Acy(ny,3), Acz(nz,3) )
        ALLOCATE ( pc(nx+2,ny+2, nz+2), pco(nx+2,ny+2, nz+2) ) 
        
        ALLOCATE ( SUMWSS(ibNodes), SIGNWSS(ibNodes) )
      END SUBROUTINE allocateArrays
