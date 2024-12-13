
SUBROUTINE pressureForcing1
      USE global
      IMPLICIT NONE
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (KIND = 8) :: diagCell, n1, pos1_x, pos1_y, pos1_z, pt1, &
                         aval, bval, cval, p_pos1, sur2nodeDis, dpdn, p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, &
                         p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1, p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e
                         
      dpdn = 0._rk
!$acc parallel loop gang vector                                                                                          &
!$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1, aval, bval, cval, p_pos1, sur2nodeDis, dpdn,                   &
!$acc           p_x1, p_x2, p_y1, p_y2, p_z1, p_z2, p_x1_z1, p_x2_z1, p_x1_z2, p_x2_z2, p_z1_x1, p_z2_x1,                &
!$acc           p_z1_x2, p_z2_x2, h1, h2, dpdn_e, dpdx_e, dpdy_e, dpdz_e, k, j, i, il, jl, kl, i_x1, i_y1, i_z1)         &
!$acc present (interceptedIndexPtr, pNormDis, deltax, deltay, deltaz, xp, yp, zp, cosAlpha, cosBeta, cosGamma, nelp, p)  &
!$acc firstprivate (nx, ny,nz) 
      DO n = 1, ibCellCount
         
         !dpdn = (-(v_curr-v_prev)/deltat)*cosBeta(nelp(n))
         
         i = interceptedIndexPtr(n, 1) 
         j = interceptedIndexPtr(n, 2) 
         k = interceptedIndexPtr(n, 3) 
         
         sur2nodeDis = pNormDis(n)

         pt1 = 1.5_rk*dsqrt(deltax(i)**2 + deltay(j)**2 + deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5

         !coordinates of three points from interceptd cell pressure node
         pos1_x = xp(i) + pt1*cosAlpha(nelp(n))	  
         pos1_y = yp(j) + pt1*cosBeta(nelp(n))
         pos1_z = zp(k) + pt1*cosGamma(nelp(n))
         
         !$acc loop seq         
         DO il = 1, nx+1         
            if(pos1_x.ge.xp(il).and.pos1_x.lt.xp(il+1)) i_x1 = il	   
         END DO       
         !$acc loop seq    
         DO jl = 1, ny+1       
            if(pos1_y.ge.yp(jl).and.pos1_y.lt.yp(jl+1)) i_y1 = jl	     
         END DO
         !$acc loop seq
         DO kl = 1, nz+1         
            if(pos1_z.ge.zp(kl).and.pos1_z.lt.zp(kl+1)) i_z1 = kl
         END DO
 
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

         dpdn_e = dpdx_e*cosAlpha(nelp(n)) + dpdy_e*cosBeta(nelp(n)) + dpdz_e*cosGamma(nelp(n))        

         n1 = pt1 + sur2nodeDis
       
         bval = dpdn
         aval = (dpdn_e - dpdn)/(2*n1)
         cvaL = p_pos1 - (dpdn_e + dpdn)*n1*.5

         p(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval 
      ENDDO
!$acc end parallel

      print*, "Leaving PressureForcing"
END SUBROUTINE pressureForcing1
!***********************************************************************

!***********************************************************************
SUBROUTINE velocityForcing1
      USE global
      IMPLICIT NONE
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: n, k, j, i, il, jl, kl, i_x1, i_y1, i_z1
      REAL (KIND = 8) :: diagCell, n1, pos1_x, pos1_y, pos1_z, pt1,  &
                         aval, bval, cval, sur2nodeDis, h1, h2, &
                         usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2, &
                         vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2, &
                         wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2, &
                         dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e, & 
                         dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1, &
                         u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
                         v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2

!$acc parallel loop gang vector         &                                                                    
!$acc private (diagCell, n1, pos1_x, pos1_y, pos1_z, pt1,           &
!$acc          aval, bval, cval, sur2nodeDis, h1, h2,               &
!$acc          usurf, u_pos1, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2,   &
!$acc          vsurf, v_pos1, v_x1, v_x2, v_y1, v_y2, v_z1, v_z2,   &
!$acc          wsurf, w_pos1, w_x1, w_x2, w_y1, w_y2, w_z1, w_z2,   &
!$acc          dudn_e, dudx_e, dudy_e, dudz_e, dvdn_e, dvdx_e, dvdy_e, dvdz_e,  & 
!$acc          dwdn_e, dwdx_e, dwdy_e, dwdz_e, u_x1_z1, u_x2_z1, u_x1_z2, u_x2_z2, u_z1_x1,  &
!$acc          u_z2_x1, u_z1_x2, u_z2_x2, v_x1_z1, v_x2_z1, v_x1_z2, v_x2_z2, v_z1_x1, v_z2_x1, v_z1_x2, &
!$acc          v_z2_x2, w_x1_z1, w_x2_z1, w_x1_z2, w_x2_z2, w_z1_x1, w_z2_x1, w_z1_x2, w_z2_x2, k, j, i, il, jl, kl, i_x1, i_y1, i_z1) &                                 
!$acc present (interceptedIndexPtr, deltax, deltay, deltaz, cosAlpha, cosBeta, cosGamma,    &
!$acc           ut, u2NormDis, u1NormDis, nelu2, nelu1, xu, yu, zu,      &
!$acc           vt, v2NormDis, v1NormDis, nelv2, nelv1, xv, yv, zv,      &
!$acc           wt, w2NormDis, w1NormDis, nelw2, nelw1, xw, yw, zw)      &
!$acc firstprivate (nx, ny,nz)    
      DO n = 1, ibCellCount
                  
         i = interceptedIndexPtr(n, 1) 
         j = interceptedIndexPtr(n, 2) 
         k = interceptedIndexPtr(n, 3) 

!***********************U(i,j,k)****************************************
         !usurf = u_curr
         usurf = 0._rk
         sur2nodeDis = u2NormDis(n)

         pt1 = 1.5_rk*dsqrt(deltax(i)**2 + deltay(j)**2 + deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5
         
         !coordinates of three points from interceptd cell pressure node
         pos1_x = xu(i+1) + pt1*cosAlpha(nelu2(n))
         pos1_y = yu(j)   + pt1*cosBeta(nelu2(n))
         pos1_z = zu(k)   + pt1*cosGamma(nelu2(n))
        
         !$acc loop seq
         DO il = 1, nx+2         
            if(pos1_x.ge.xu(il).and.pos1_x.lt.xu(il+1)) i_x1 = il
         END DO        
         !$acc loop seq
         DO jl = 1, ny+1         
            if(pos1_y.ge.yu(jl).and.pos1_y.lt.yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, nz+1        
            if(pos1_z.ge.zu(kl).and.pos1_z.lt.zu(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         IF(i_x1.EQ.nx+2) i_x1 = nx+1
         !IF(i_y1.EQ.ny+2) i_y1 = ny+1
         !IF(i_z1.EQ.nz+2) i_z1 = nz+1

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

         dudn_e = dudx_e*cosAlpha(nelu2(n)) + dudy_e*cosBeta(nelu2(n)) + dudz_e*cosGamma(nelu2(n))         

         n1 = pt1 + sur2nodeDis
         
         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         ut(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval  
!***********************************************************************  
               
!******************************U(i-1,j,k)*******************************       
         !usurf = u_curr
         usurf = 0._rk
         sur2nodeDis = u1NormDis(n)

         pt1 = 1.5_rk*dsqrt(deltax(i)**2 + deltay(j)**2 + deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5
         
         !coordinates of three points from interceptd cell pressure node
         pos1_x = xu(i) + pt1*cosAlpha(nelu1(n))
         pos1_y = yu(j) + pt1*cosBeta(nelu1(n))
         pos1_z = zu(k) + pt1*cosGamma(nelu1(n))
         
         !$acc loop seq
         DO il = 1, nx+1        
            if(pos1_x.ge.xu(il).and.pos1_x.lt.xu(il+1)) i_x1 = il
         END DO        
         !$acc loop seq
         DO jl = 1, ny+1         
            if(pos1_y.ge.yu(jl).and.pos1_y.lt.yu(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, nz+1         
            if(pos1_z.ge.zu(kl).and.pos1_z.lt.zu(kl+1)) i_z1 = kl
         END DO

         IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.nx+2) i_x1 = nx+1
         !IF(i_y1.EQ.ny+2) i_y1 = ny+1
         !IF(i_z1.EQ.nz+2) i_z1 = nz+1
 
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

         dudn_e = dudx_e*cosAlpha(nelu1(n)) + dudy_e*cosBeta(nelu1(n)) + dudz_e*cosGamma(nelu1(n))            

         n1 = pt1 + sur2nodeDis
         
         cval = usurf
         bval = 2./n1*(u_pos1 - usurf) - dudn_e
         avaL = dudn_e/n1 - (u_pos1 - usurf)/n1**2
         ut(i-1,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval                          
!***********************************************************************
         
!**************************V(i,j,k)*************************************         
         !vsurf = v_curr
         vsurf = 0._rk
         sur2nodeDis = v2NormDis(n)
          
         pt1 = 1.5_rk*dsqrt(deltax(i)**2 + deltay(j)**2 + deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5
	   
         !coordinates of three points from interceptd cell pressure node
         pos1_x = xv(i) + pt1*cosAlpha(nelv2(n))
         pos1_y = yv(j+1) + pt1*cosBeta(nelv2(n))
         pos1_z = zv(k) + pt1*cosGamma(nelv2(n))
         
         !$acc loop seq
         DO il = 1, nx+1         
            if(pos1_x.ge.xv(il).and.pos1_x.lt.xv(il+1)) i_x1 = il
         END DO       
         !$acc loop seq
         DO jl = 1, ny+2         
            if(pos1_y.ge.yv(jl).and.pos1_y.lt.yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, nz+1         
            if(pos1_z.ge.zv(kl).and.pos1_z.lt.zv(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.nx+2) i_x1 = nx+1
         IF(i_y1.EQ.ny+2) i_y1 = ny+1
         !IF(i_z1.EQ.nz+2) i_z1 = nz+1

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
         

         dvdn_e = dvdx_e*cosAlpha(nelv2(n)) + dvdy_e*cosBeta(nelv2(n)) +  dvdz_e*cosGamma(nelv2(n))
         
         n1 = pt1 + sur2nodeDis
         
         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         vt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval 
!*********************************************************************** 
       
!**************************V(i,j-1,k)*************************************         
         !vsurf = v_curr
         vsurf = 0._rk
         sur2nodeDis = v1NormDis(n)
          
         pt1 = 1.5_rk*dsqrt(deltax(i)**2 + deltay(j)**2 + deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5
	   
         !coordinates of three points from interceptd cell pressure node
         pos1_x = xv(i) + pt1*cosAlpha(nelv1(n))
         pos1_y = yv(j) + pt1*cosBeta(nelv1(n))
         pos1_z = zv(k) + pt1*cosGamma(nelv1(n))
         
         !$acc loop seq
         DO il = 1, nx+1         
            if(pos1_x.ge.xv(il).and.pos1_x.lt.xv(il+1)) i_x1 = il
         END DO       
         !$acc loop seq
         DO jl = 1, ny+1         
            if(pos1_y.ge.yv(jl).and.pos1_y.lt.yv(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, nz+1        
            if(pos1_z.ge.zv(kl).and.pos1_z.lt.zv(kl+1)) i_z1 = kl
         END DO
         
         !IF(i_x1.EQ.1) i_x1 = 2
         IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.nx+2) i_x1 = nx+1
         !IF(i_y1.EQ.ny+2) i_y1 = ny+1
         !IF(i_z1.EQ.nz+2) i_z1 = nz+1

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
         

         dvdn_e = dvdx_e*cosAlpha(nelv1(n)) + dvdy_e*cosBeta(nelv1(n)) +  dvdz_e*cosGamma(nelv1(n))
         
         n1 = pt1 + sur2nodeDis
         
         cval = vsurf
         bval = 2./n1*(v_pos1 - vsurf) - dvdn_e
         avaL = dvdn_e/n1 - (v_pos1 - vsurf)/n1**2
         vt(i,j-1,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval 
!*********************************************************************** 
              
!**************************W(i,j,k)*************************************         
         !wsurf = w_curr
         wsurf = 0._rk
         sur2nodeDis = w2NormDis(n)
          
         pt1 = 1.5_rk*dsqrt(deltax(i)**2 + deltay(j)**2 + deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5
	   
         !coordinates of three points from interceptd cell pressure node
         pos1_x = xw(i) + pt1*cosAlpha(nelw2(n))
         pos1_y = yw(j) + pt1*cosBeta(nelw2(n))
         pos1_z = zw(k+1) + pt1*cosGamma(nelw2(n))
         
         !$acc loop seq
         DO il = 1, nx+1         
            if(pos1_x.ge.xw(il).and.pos1_x.lt.xw(il+1)) i_x1 = il
         END DO       
         !$acc loop seq
         DO jl = 1, ny+1         
            if(pos1_y.ge.yw(jl).and.pos1_y.lt.yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, nz+2         
            if(pos1_z.ge.zw(kl).and.pos1_z.lt.zw(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         !IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.nx+2) i_x1 = nx+1
         !IF(i_y1.EQ.ny+2) i_y1 = ny+1
         IF(i_z1.EQ.nz+2) i_z1 = nz+1
 
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
         

         dwdn_e = dwdx_e*cosAlpha(nelw2(n)) + dwdy_e*cosBeta(nelw2(n)) +  dwdz_e*cosGamma(nelw2(n))
         
         n1 = pt1 + sur2nodeDis
         
         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         wt(i,j,k) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval 
!*********************************************************************** 
         
!**************************W(i,j,k-1)*************************************         
         !wsurf = w_curr
         wsurf = 0._rk
         sur2nodeDis = w1NormDis(n)
          
         pt1 = 1.5_rk*dsqrt(deltax(i)**2 + deltay(j)**2 + deltaz(k)**2) + (dabs(sur2nodeDis)-sur2nodeDis)*0.5
	   
         !coordinates of three points from interceptd cell pressure node
         pos1_x = xw(i) + pt1*cosAlpha(nelw1(n))
         pos1_y = yw(j) + pt1*cosBeta(nelw1(n))
         pos1_z = zw(k) + pt1*cosGamma(nelw1(n))
         
         !$acc loop seq
         DO il = 1, nx+1         
            if(pos1_x.ge.xw(il).and.pos1_x.lt.xw(il+1)) i_x1 = il
         END DO       
         !$acc loop seq
         DO jl = 1, ny+1         
            if(pos1_y.ge.yw(jl).and.pos1_y.lt.yw(jl+1)) i_y1 = jl
         END DO
         !$acc loop seq
         DO kl = 1, nz+1         
            if(pos1_z.ge.zw(kl).and.pos1_z.lt.zw(kl+1)) i_z1 = kl
         END DO

         !IF(i_x1.EQ.1) i_x1 = 2
         !IF(i_y1.EQ.1) i_y1 = 2
         IF(i_z1.EQ.1) i_z1 = 2
         !IF(i_x1.EQ.nx+2) i_x1 = nx+1
         !IF(i_y1.EQ.ny+2) i_y1 = ny+1
         !IF(i_z1.EQ.nz+2) i_z1 = nz+1

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
         

         dwdn_e = dwdx_e*cosAlpha(nelw1(n)) + dwdy_e*cosBeta(nelw1(n)) +  dwdz_e*cosGamma(nelw1(n))
         
         n1 = pt1 + sur2nodeDis
         
         cval = wsurf
         bval = 2./n1*(w_pos1 - wsurf) - dwdn_e
         avaL = dwdn_e/n1 - (w_pos1 - wsurf)/n1**2
         wt(i,j,k-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval  
!***********************************************************************  
      ENDDO
!$acc end parallel

      print*, "Leaving VelocityForcing"
END SUBROUTINE velocityForcing1
