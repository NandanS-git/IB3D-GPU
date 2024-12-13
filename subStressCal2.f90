!csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
       SUBROUTINE stressCal2
       USE global
       IMPLICIT NONE
       INTEGER, PARAMETER :: rk = selected_real_kind(8) 

       INTEGER:: i, j, k, ielem, inode, i_x1, i_y1, i_z1
       
       INTEGER:: i_cell, j_cell, k_cell

       REAL(KIND = 8):: diagdis, normdis, aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z

       REAL(KIND = 8):: xsurf, ysurf, zsurf, pos1_x, pos1_y, pos1_z,              & 
                        psurf, p_pos1, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2,  &
                        p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,     &
                        p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e, &
                        usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,        &
                        vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,        &
                        wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,        &
                        dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e,   & 
                        dvdz_e, dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, &
                        u_x1_z2, u_x2_z2, u_z1_x1, u_z2_x1, u_z1_x2, u_z2_x2,     &
                        v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1,     &
                        v_z1_x2, v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2,     &
                        w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, dudn_s, dvdn_s, dwdn_s

       !REAL(KIND = 8):: alen, area  
  
       !REAL(KIND = 8), DIMENSION(ibNodes):: Anode, Afnode, stressNode, TAWSS, OSI       
       
       !REAL(KIND = 8), DIMENSION(ibElems):: stressElem     
       
       !CHARACTER(LEN=70):: filename1
       
       REAL(KIND = 8):: alen, area, area_xz, area_yz, area_xy
     
       REAL(KIND = 8):: shear_x_force, shear_y_force, shear_z_force,  &
                        f_surf, f_surf_x, f_surf_y, f_surf_z

       REAL(KIND = 8):: pressureDrag, viscousDrag, viscousLift,                        &
                        pressureLift, viscousDragcoefficient, pressureDragcoefficient, &
                        viscousLiftcoefficient, PressureLiftcoefficient, area_Sx,      &
                        area_Sy, surf_area
         
       !stressElem = 0._rk
       
       viscousDrag = 0.
       pressureDrag = 0.
       viscousLift = 0.
       PressureLift = 0.
       surf_area = 0.
	area_Sx = 0.
       area_Sy = 0
       
       !$acc parallel loop gang vector reduction(+: pressureDrag, viscousDrag, viscousLift, pressureLift, area_Sx, area_Sy, surf_area)  &
       !$acc private (i_x1, i_y1, i_z1, i_cell, j_cell, k_cell, diagdis, normdis,                    &
       !$acc          aval, bval, cval, stx1, sty1, stz1, del_X, del_Y, del_Z,                       &
       !$acc          xsurf, ysurf, zsurf, pos1_x, pos1_y, pos1_z,                                   & 
       !$acc          psurf, p_pos1, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2,                       &
       !$acc          p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,                          &
       !$acc          p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e,                      &
       !$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,                             &
       !$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,                             &
       !$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,                             &
       !$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e,                        & 
       !$acc          dvdz_e, dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1,                      &
       !$acc          u_x1_z2, u_x2_z2, u_z1_x1, u_z2_x1, u_z1_x2, u_z2_x2,                          &
       !$acc          v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1,                          &
       !$acc          v_z1_x2, v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2,                          &
       !$acc          w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, dudn_s, dvdn_s, dwdn_s,                    &
       !$acc          alen, area, area_xz, area_yz, area_xy, shear_x_force, shear_y_force,           &
       !$acc          shear_z_force, f_surf, f_surf_x, f_surf_y, f_surf_z, viscousDragcoefficient,   &
       !$acc          pressureDragcoefficient,viscousLiftcoefficient, PressureLiftcoefficient)       &
       !$acc present (xcent, ycent, zcent, x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp, &
       !$acc          ut, u, vt, v, wt, w, p, deltax, deltay, deltaz, cosAlpha, cosBeta, cosGamma,     &
       !$acc          alpha3, beta3, gamma3, xnode1, ynode1, znode1, Elemcell, ucell, vcell, wcell, pcell)    &
       !$acc firstprivate (nx, ny, nz, deltat)
       DO ielem = 1, ibElems
       !IF((xcent(ielem).ge.1.0).and.(xcent(ielem).le.131.0)) THEN
!***********************interpolation points****************************
       i_cell = Elemcell(ielem,1) 
       j_cell = Elemcell(ielem,2) 
       k_cell = Elemcell(ielem,3)
       
       xsurf = xcent(ielem)
       ysurf = ycent(ielem)
       zsurf = zcent(ielem)

	del_X = x1(i_cell+1)-x1(i_cell)
       del_Y = y1(j_cell+1)-y1(j_cell)
       del_Z = z1(k_cell+1)-z1(k_cell)

       diagdis = dsqrt(del_X**2 + del_Y**2 + del_Z**2)
       normdis = diagdis
	pos1_x = xsurf + normdis*cosAlpha(ielem)
	pos1_y = ysurf + normdis*cosBeta(ielem)
	pos1_z = zsurf + normdis*cosGamma(ielem)

!**************************velocity at the surface**********************	
	usurf = 0.
       vsurf = 0.
       wsurf = 0.
       dpdn = 0.
!*******************velocity interpolation at point 2******************

!******************u velocity interpolation at point 2******************    
        i_x1 = ucell(ielem,1)  
        i_y1 = ucell(ielem,2)  
        i_z1 = ucell(ielem,3)  
       
       !interpolation along x @ z1 plane 
         u_x1_z1 = ut(i_x1-1, i_y1, i_z1)   + (ut(i_x1, i_y1, i_z1)   - ut(i_x1-1, i_y1, i_z1))  *(pos1_x - xu(i_x1))/(xu(i_x1+1) - xu(i_x1))
         u_x2_z1 = ut(i_x1-1, i_y1+1, i_z1) + (ut(i_x1, i_y1+1, i_z1) - ut(i_x1-1, i_y1+1, i_z1))*(pos1_x - xu(i_x1))/(xu(i_x1+1) - xu(i_x1))
         
         !interpolation along x @ z2 plane 
         u_x1_z2 = ut(i_x1-1, i_y1, i_z1+1)   + (ut(i_x1, i_y1, i_z1+1)   - ut(i_x1-1, i_y1, i_z1+1))  *(pos1_x - xu(i_x1))/(xu(i_x1+1) - xu(i_x1))
         u_x2_z2 = ut(i_x1-1, i_y1+1, i_z1+1) + (ut(i_x1, i_y1+1, i_z1+1) - ut(i_x1-1, i_y1+1, i_z1+1))*(pos1_x - xu(i_x1))/(xu(i_x1+1) - xu(i_x1))
         
         !This will be used for dudz calculation point 1
         u_z1 = u_x1_z1 + (u_x2_z1 - u_x1_z1)*(pos1_y - yu(i_y1))/(yu(i_y1+1)-yu(i_y1))
         u_z2 = u_x1_z2 + (u_x2_z2 - u_x1_z2)*(pos1_y - yu(i_y1))/(yu(i_y1+1)-yu(i_y1))
         
         !This will be used for dudy calculation at point 1
         u_y1 = u_x1_z1 + (u_x1_z2 - u_x1_z1)*(pos1_z - zu(i_z1))/(zu(i_z1+1)-zu(i_z1))
         u_y2 = u_x2_z1 + (u_x2_z2 - u_x2_z1)*(pos1_z - zu(i_z1))/(zu(i_z1+1)-zu(i_z1))
         
         !interpolation along z @ x1 plane
         u_z1_x1 = ut(i_x1-1, i_y1, i_z1)   + (ut(i_x1-1, i_y1, i_z1+1) - ut(i_x1-1, i_y1, i_z1))  *(pos1_z - zu(i_z1))/(zu(i_z1+1) - zu(i_z1))
         u_z2_x1 = ut(i_x1-1, i_y1+1, i_z1)   + (ut(i_x1-1, i_y1+1, i_z1+1) - ut(i_x1-1, i_y1+1, i_z1))  *(pos1_z - zu(i_z1))/(zu(i_z1+1) - zu(i_z1))
         
         !interpolation along z @ x2 plane
         u_z1_x2 = ut(i_x1, i_y1, i_z1)   + (ut(i_x1, i_y1, i_z1+1) - ut(i_x1, i_y1, i_z1))  *(pos1_z - zu(i_z1))/(zu(i_z1+1) - zu(i_z1))
         u_z2_x2 = ut(i_x1, i_y1+1, i_z1)   + (ut(i_x1, i_y1+1, i_z1+1) - ut(i_x1, i_y1+1, i_z1))  *(pos1_z - zu(i_z1))/(zu(i_z1+1) - zu(i_z1))
         
         !This will be used for dudx calculation at point 1
         u_x1 = u_z1_x1 + (u_z2_x1 - u_z1_x1)*(pos1_y - yu(i_y1))/(yu(i_y1+1)-yu(i_y1))
         u_x2 = u_z1_x2 + (u_z2_x2 - u_z1_x2)*(pos1_y - yu(i_y1))/(yu(i_y1+1)-yu(i_y1))
         
         u_pos1 = u_x1 + (u_x2 - u_x1)*(pos1_x - xu(i_x1))/(xu(i_x1+1)-xu(i_x1)) 
         
         h2 = dabs(xu(i_x1+1) - pos1_x)
         h1 = dabs(xu(i_x1)   - pos1_x)
         dudx_e = (h1**2*u_x2 - h2**2*u_x1 + (h2**2- h1**2)*u_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(yu(i_y1+1) - pos1_y)
         h1 = dabs(yu(i_y1)   - pos1_y)
         dudy_e = (h1**2*u_y2 - h2**2*u_y1 + (h2**2- h1**2)*u_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(zu(i_z1+1) - pos1_z)
         h1 = dabs(zu(i_z1)   - pos1_z)
         dudz_e = (h1**2*u_z2 - h2**2*u_z1 + (h2**2- h1**2)*u_pos1)/(h1*h2*(h1+h2)+1e-16)

         dudn_e = dudx_e*cosAlpha(ielem) + dudy_e*cosBeta(ielem) + dudz_e*cosGamma(ielem)         

         dudn_s = (2./normdis)*(u_pos1 - usurf) - dudn_e
         
!******************v velocity interpolation in point 2******************	 
        i_x1 = vcell(ielem,1)  
        i_y1 = vcell(ielem,2)  
        i_z1 = vcell(ielem,3) 
       
       !interpolation along x @ z1 plane 
         v_x1_z1 = vt(i_x1, i_y1-1, i_z1)   + (vt(i_x1+1, i_y1-1, i_z1) - vt(i_x1, i_y1-1, i_z1))*(pos1_x - xv(i_x1))/(xv(i_x1+1) - xv(i_x1))
         v_x2_z1 = vt(i_x1, i_y1, i_z1)     + (vt(i_x1+1, i_y1, i_z1)   - vt(i_x1, i_y1, i_z1))*(pos1_x - xv(i_x1))/(xv(i_x1+1) - xv(i_x1))
         
         !interpolation along x @ z2 plane 
         v_x1_z2 = vt(i_x1, i_y1-1, i_z1+1)   + (vt(i_x1+1, i_y1-1, i_z1+1) - vt(i_x1, i_y1-1, i_z1+1))*(pos1_x - xv(i_x1))/(xv(i_x1+1) - xv(i_x1))
         v_x2_z2 = vt(i_x1, i_y1, i_z1+1)     + (vt(i_x1+1, i_y1, i_z1+1)   - vt(i_x1, i_y1, i_z1+1))*(pos1_x - xv(i_x1))/(xv(i_x1+1) - xv(i_x1))
         
         !This will be used for dvdz calculation point 1
         v_z1 = v_x1_z1 + (v_x2_z1 - v_x1_z1)*(pos1_y - yv(i_y1))/(yv(i_y1+1)-yv(i_y1))
         v_z2 = v_x1_z2 + (v_x2_z2 - v_x1_z2)*(pos1_y - yv(i_y1))/(yv(i_y1+1)-yv(i_y1))
         
         !This will be used for dvdy calculation at point 1
         v_y1 = v_x1_z1 + (v_x1_z2 - v_x1_z1)*(pos1_z - zv(i_z1))/(zv(i_z1+1)-zv(i_z1))
         v_y2 = v_x2_z1 + (v_x2_z2 - v_x2_z1)*(pos1_z - zv(i_z1))/(zv(i_z1+1)-zv(i_z1))
         
         !interpolation along z @ x1 plane
         v_z1_x1 = vt(i_x1, i_y1-1, i_z1)   + (vt(i_x1, i_y1-1, i_z1+1) - vt(i_x1, i_y1-1, i_z1))  *(pos1_z - zv(i_z1))/(zv(i_z1+1) - zv(i_z1))
         v_z2_x1 = vt(i_x1, i_y1, i_z1)   + (vt(i_x1, i_y1, i_z1+1) - vt(i_x1, i_y1, i_z1))  *(pos1_z - zv(i_z1))/(zv(i_z1+1) - zv(i_z1))
         
         !interpolation along z @ x2 plane
         v_z1_x2 = vt(i_x1+1, i_y1-1, i_z1)   + (vt(i_x1+1, i_y1-1, i_z1+1) - vt(i_x1+1, i_y1-1, i_z1))  *(pos1_z - zv(i_z1))/(zv(i_z1+1) - zv(i_z1))
         v_z2_x2 = vt(i_x1+1, i_y1, i_z1)   + (vt(i_x1+1, i_y1, i_z1+1) - vt(i_x1+1, i_y1, i_z1))  *(pos1_z - zv(i_z1))/(zv(i_z1+1) - zv(i_z1))
         
         !This will be used for dvdx calculation at point 1
         v_x1 = v_z1_x1 + (v_z2_x1 - v_z1_x1)*(pos1_y - yv(i_y1))/(yv(i_y1+1)-yv(i_y1))
         v_x2 = v_z1_x2 + (v_z2_x2 - v_z1_x2)*(pos1_y - yv(i_y1))/(yv(i_y1+1)-yv(i_y1))
         
         v_pos1 = v_x1 + (v_x2 - v_x1)*(pos1_x - xv(i_x1))/(xv(i_x1+1)-xv(i_x1))
         
         h2 = dabs(xv(i_x1+1) - pos1_x)
         h1 = dabs(xv(i_x1)   - pos1_x)
         dvdx_e = (h1**2*v_x2 - h2**2*v_x1 + (h2**2- h1**2)*v_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(yv(i_y1+1) - pos1_y)
         h1 = dabs(yv(i_y1)   - pos1_y)
         dvdy_e = (h1**2*v_y2 - h2**2*v_y1 + (h2**2- h1**2)*v_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(zv(i_z1+1) - pos1_z)
         h1 = dabs(zv(i_z1)   - pos1_z)
         dvdz_e = (h1**2*v_z2 - h2**2*v_z1 + (h2**2- h1**2)*v_pos1)/(h1*h2*(h1+h2)+1e-16)
         

         dvdn_e = dvdx_e*cosAlpha(ielem) + dvdy_e*cosBeta(ielem) +  dvdz_e*cosGamma(ielem)

	  dvdn_s = (2./normdis)*(v_pos1 - vsurf) - dvdn_e
	  
!******************w velocity interpolation in point 2******************
        i_x1 = wcell(ielem,1)  
        i_y1 = wcell(ielem,2)  
        i_z1 = wcell(ielem,3) 
       
       !interpolation along x @ z1 plane 
         w_x1_z1 = wt(i_x1, i_y1, i_z1-1)   + (wt(i_x1+1, i_y1, i_z1-1) - wt(i_x1, i_y1, i_z1-1))*(pos1_x - xw(i_x1))/(xw(i_x1+1) - xw(i_x1))
         w_x2_z1 = wt(i_x1, i_y1+1, i_z1-1)   + (wt(i_x1+1, i_y1+1, i_z1-1) - wt(i_x1, i_y1+1, i_z1-1))*(pos1_x - xw(i_x1))/(xw(i_x1+1) - xw(i_x1))
         !interpolation along x @ z2 plane 
         w_x1_z2 = wt(i_x1, i_y1, i_z1)   + (wt(i_x1+1, i_y1, i_z1) - wt(i_x1, i_y1, i_z1))*(pos1_x - xw(i_x1))/(xw(i_x1+1) - xw(i_x1))
         w_x2_z2 = wt(i_x1, i_y1+1, i_z1)   + (wt(i_x1+1, i_y1+1, i_z1) - wt(i_x1, i_y1+1, i_z1))*(pos1_x - xw(i_x1))/(xw(i_x1+1) - xw(i_x1))
         
         !This will be used for dwdz calculation point 1
         w_z1 = w_x1_z1 + (w_x2_z1 - w_x1_z1)*(pos1_y - yw(i_y1))/(yw(i_y1+1)-yw(i_y1))
         w_z2 = w_x1_z2 + (w_x2_z2 - w_x1_z2)*(pos1_y - yw(i_y1))/(yw(i_y1+1)-yw(i_y1))
         
         !This will be used for dwdy calculation at point 1
         w_y1 = w_x1_z1 + (w_x1_z2 - w_x1_z1)*(pos1_z - zw(i_z1))/(zw(i_z1+1)-zw(i_z1))
         w_y2 = w_x2_z1 + (w_x2_z2 - w_x2_z1)*(pos1_z - zw(i_z1))/(zw(i_z1+1)-zw(i_z1))
         
         !interpolation along z @ x1 plane
         w_z1_x1 = wt(i_x1, i_y1, i_z1-1)   + (wt(i_x1, i_y1, i_z1) - wt(i_x1, i_y1, i_z1-1))  *(pos1_z - zw(i_z1))/(zw(i_z1+1) - zw(i_z1))
         w_z2_x1 = wt(i_x1, i_y1+1, i_z1-1)   + (wt(i_x1, i_y1+1, i_z1) - wt(i_x1, i_y1+1, i_z1-1))  *(pos1_z - zw(i_z1))/(zw(i_z1+1) - zw(i_z1))
         
         !interpolation along z @ x2 plane
         w_z1_x2 = wt(i_x1+1, i_y1, i_z1-1)   + (wt(i_x1+1, i_y1, i_z1) - wt(i_x1+1, i_y1, i_z1-1))  *(pos1_z - zw(i_z1))/(zw(i_z1+1) - zw(i_z1))
         w_z2_x2 = wt(i_x1+1, i_y1+1, i_z1-1)   + (wt(i_x1+1, i_y1+1, i_z1) - wt(i_x1+1, i_y1+1, i_z1-1))  *(pos1_z - zw(i_z1))/(zw(i_z1+1) - zw(i_z1))
         
         !This will be used for dwdx calculation at point 1
         w_x1 = w_z1_x1 + (w_z2_x1 - w_z1_x1)*(pos1_y - yw(i_y1))/(yw(i_y1+1)-yw(i_y1))
         w_x2 = w_z1_x2 + (w_z2_x2 - w_z1_x2)*(pos1_y - yw(i_y1))/(yw(i_y1+1)-yw(i_y1))
         
         w_pos1 = w_x1 + (w_x2 - w_x1)*(pos1_x - xw(i_x1))/(xw(i_x1+1)-xw(i_x1))
         
         h2 = dabs(xw(i_x1+1) - pos1_x)
         h1 = dabs(xw(i_x1)   - pos1_x)
         dwdx_e = (h1**2*w_x2 - h2**2*w_x1 + (h2**2- h1**2)*w_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(yw(i_y1+1) - pos1_y)
         h1 = dabs(yw(i_y1)   - pos1_y)
         dwdy_e = (h1**2*w_y2 - h2**2*w_y1 + (h2**2- h1**2)*w_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(zw(i_z1+1) - pos1_z)
         h1 = dabs(zw(i_z1)   - pos1_z)
         dwdz_e = (h1**2*w_z2 - h2**2*w_z1 + (h2**2- h1**2)*w_pos1)/(h1*h2*(h1+h2)+1e-16)
         

         dwdn_e = dwdx_e*cosAlpha(ielem) + dwdy_e*cosBeta(ielem) +  dwdz_e*cosGamma(ielem)
       
         dwdn_s = (2./normdis)*(w_pos1 - wsurf) - dwdn_e
         
!***********************calculate area of the elements******************    
       alen = sqrt(alpha3(ielem)**2.+beta3(ielem)**2.+ gamma3(ielem)**2.) 
       area = alen/2.
       area_yz = 0.5*abs(alpha3(ielem))
       area_xz = 0.5*abs(beta3(ielem))
       area_xy = 0.5*abs(gamma3(ielem))
       !areaElem(ielem) = area

!*********non-dimensional viscous stress & force calculation************
       stx1 = dudn_s - (dudn_s*cosAlpha(ielem) + dvdn_s*cosBeta(ielem) + dwdn_s*cosGamma(ielem))*cosAlpha(ielem)

       sty1 = dvdn_s - (dudn_s*cosAlpha(ielem) + dvdn_s*cosBeta(ielem) + dwdn_s*cosGamma(ielem))*cosBeta(ielem)

       stz1 = dwdn_s - (dudn_s*cosAlpha(ielem) + dvdn_s*cosBeta(ielem) + dwdn_s*cosGamma(ielem))*cosGamma(ielem)

       shear_x_force = (1.0/re)*stx1*area

       shear_y_force = (1.0/re)*sty1*area
       
       shear_z_force = (1.0/re)*stz1*area
       
       !stressElem(ielem) = 0.0035*stx1
       
!************************presssure interpolation************************

!*******************pressure interpolation at point 2******************* 
       i_x1 = pcell(ielem,1)
       i_y1 = pcell(ielem,2)
       i_z1 = pcell(ielem,3)
	
        !interpolation along x  @ z1 plane
         p_x1_z1 = p(i_x1, i_y1, i_z1)   + (p(i_x1+1, i_y1, i_z1)   - p(i_x1, i_y1, i_z1))  *(pos1_x - xp(i_x1))/(xp(i_x1+1) - xp(i_x1))
         p_x2_z1 = p(i_x1, i_y1+1, i_z1) + (p(i_x1+1, i_y1+1, i_z1) - p(i_x1, i_y1+1, i_z1))*(pos1_x - xp(i_x1))/(xp(i_x1+1) - xp(i_x1))
               
         !interpolation along x  @ z2 plane
         p_x1_z2 = p(i_x1, i_y1, i_z1+1)   + (p(i_x1+1, i_y1, i_z1+1)   - p(i_x1, i_y1, i_z1+1))  *(pos1_x - xp(i_x1))/(xp(i_x1+1) - xp(i_x1))
         p_x2_z2 = p(i_x1, i_y1+1, i_z1+1) + (p(i_x1+1, i_y1+1, i_z1+1) - p(i_x1, i_y1+1, i_z1+1))*(pos1_x - xp(i_x1))/(xp(i_x1+1) - xp(i_x1))
                 
         !This will be used for dpdz calculation point 1
         p_z1 = p_x1_z1 + (p_x2_z1 - p_x1_z1)*(pos1_y - yp(i_y1))/(yp(i_y1+1)-yp(i_y1))
         p_z2 = p_x1_z2 + (p_x2_z2 - p_x1_z2)*(pos1_y - yp(i_y1))/(yp(i_y1+1)-yp(i_y1))
         
         !This will be used for dpdy calculation at point 1
         p_y1 = p_x1_z1 + (p_x1_z2 - p_x1_z1)*(pos1_z - zp(i_z1))/(zp(i_z1+1)-zp(i_z1))
         p_y2 = p_x2_z1 + (p_x2_z2 - p_x2_z1)*(pos1_z - zp(i_z1))/(zp(i_z1+1)-zp(i_z1))
                  
         !interpolation along z @ x1 plane
         p_z1_x1 = p(i_x1, i_y1, i_z1)   + (p(i_x1, i_y1, i_z1+1) - p(i_x1, i_y1, i_z1))  *(pos1_z - zp(i_z1))/(zp(i_z1+1) - zp(i_z1))
         p_z2_x1 = p(i_x1, i_y1+1, i_z1)   + (p(i_x1, i_y1+1, i_z1+1) - p(i_x1, i_y1+1, i_z1))  *(pos1_z - zp(i_z1))/(zp(i_z1+1) - zp(i_z1))
         
         !interpolation along z @ x2 plane
         p_z1_x2 = p(i_x1+1, i_y1, i_z1)   + (p(i_x1+1, i_y1, i_z1+1) - p(i_x1+1, i_y1, i_z1))  *(pos1_z - zp(i_z1))/(zp(i_z1+1) - zp(i_z1))
         p_z2_x2 = p(i_x1+1, i_y1+1, i_z1)   + (p(i_x1+1, i_y1+1, i_z1+1) - p(i_x1+1, i_y1+1, i_z1))  *(pos1_z - zp(i_z1))/(zp(i_z1+1) - zp(i_z1))
         
         !This will be used for dpdx calculation at point 1
         p_x1 = p_z1_x1 + (p_z2_x1 - p_z1_x1)*(pos1_y - yp(i_y1))/(yp(i_y1+1)-yp(i_y1))
         p_x2 = p_z1_x2 + (p_z2_x2 - p_z1_x2)*(pos1_y - yp(i_y1))/(yp(i_y1+1)-yp(i_y1))
         
         p_pos1 = p_x1 + (p_x2 - p_x1)*(pos1_x - xp(i_x1))/(xp(i_x1+1)-xp(i_x1))
         
         h2 = dabs(xp(i_x1+1) - pos1_x)
         h1 = dabs(xp(i_x1)   - pos1_x)
         dpdx_e = (h1**2*p_x2 - h2**2*p_x1 + (h2**2- h1**2)*p_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(yp(i_y1+1) - pos1_y)
         h1 = dabs(yp(i_y1)   - pos1_y)
         dpdy_e = (h1**2*p_y2 - h2**2*p_y1 + (h2**2- h1**2)*p_pos1)/(h1*h2*(h1+h2)+1e-16)
         
         h2 = dabs(zp(i_z1+1) - pos1_z)
         h1 = dabs(zp(i_z1)   - pos1_z)
         dpdz_e = (h1**2*p_z2 - h2**2*p_z1 + (h2**2- h1**2)*p_pos1)/(h1*h2*(h1+h2)+1e-16)

         dpdn_e = dpdx_e*cosAlpha(ielem) + dpdy_e*cosBeta(ielem) + dpdz_e*cosGamma(ielem)        
 
         bval = dpdn  !dpdn=-dudt
         aval = (dpdn_e - dpdn)/(2*diagdis)
         cval = p_pos1 - (dpdn_e + dpdn)*diagdis*.5

         psurf = cval

         f_surf = cval*area 

         f_surf_x = -f_surf*cosAlpha(ielem)

         f_surf_y = -f_surf*cosBeta(ielem)

	  f_surf_z = -f_surf*cosGamma(ielem)
       
!***********************drag calculation********************************
	  viscousDrag = viscousDrag + shear_x_force
         pressureDrag = pressureDrag + f_surf_x
         viscousLift = viscousLift + shear_y_force
         PressureLift = PressureLift + f_surf_y
         surf_area = surf_area + area
	  area_Sx = area_Sx + area_xz
         area_Sy = area_Sy + area_yz
!***********************************************************************          
	!END IF
	END DO
	!$acc end parallel
!***********************************************************************
      
       area_Sx = 0.5*area_Sx
       area_Sy = 0.5*area_Sy
!*********************drag coefficient calculation**********************
	viscousDragcoefficient = 2.*(viscousDrag/area_Sx)
	
       pressureDragcoefficient = 2.*(pressureDrag/area_Sy)	
       
       viscousLiftcoefficient = 2.*(viscousLift/area_Sy) 
       
       pressureLiftcoefficient = 2.*(pressureLift/area_Sx)
       
!*********************drag file writing*********************************
       OPEN(899,file='dragcoff.dat',Access='Append',status='unknown')
       WRITE(899,*) viscousDragcoefficient, pressureDragcoefficient, totime
       CLOSE(899)
       OPEN(999,file='liftcoff.dat',Access='Append',status='unknown')
       WRITE(999,*) viscousLiftcoefficient,  pressureLiftcoefficient, totime
       CLOSE(999)
      
	END SUBROUTINE stressCal2








