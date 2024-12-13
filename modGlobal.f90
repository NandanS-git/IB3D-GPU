       MODULE global
       IMPLICIT NONE
       CHARACTER (LEN = 128) :: line
       INTEGER               :: istart, id1, ibCellCount, fluidCellCount, solidCellCount, redCellCount, blackCellCount
       INTEGER ( KIND = 8)   :: nx, ny, nz, itamax, pcItaMax,      &
                                ita, nIterPcor, nc, ita1, ital,    &
                                sumIterPc,itaSola, totIterPc, inor, &
				    k_startSearch, k_endSearch, &
				    j_startSearch, j_endSearch, &
				    i_startSearch, i_endSearch
       REAL (KIND =8)        :: lx, ly,  dt_order,   & 
                                deltx2, delty2,  deltz2, &
                                a0, freq, &
                                u0, v0, w0, p0,  &  
                                eps1, epsDiv, re, rev, divmax,   &
                                fx, fy,  alpha, pct, pct1, &
                                deltat, totime, dfinish, dstart, solverTime, pi , msTime, mindx, al, uc
       
       !!!variables for Orlanski multiple outlet
       REAL (KIND = 8)       :: y11, y12, z11, z12, y21, y22, z21, z22, uc11, uc22

                         
       REAL (KIND =8)        :: u_init, u_final, v_init, v_final, w_init, w_final, xmove, ymove, &
                                zmove, Total_Force_X, Total_Force_Y, Total_Force_Z                               
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:,:,:)        :: cell         
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:,:,:)    :: minElemcell
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:,:)      :: Elemcell, ucell, vcell, wcell, pcell
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)           :: SUMWSS, SIGNWSS, areaElem                                              
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)       :: deltax, deltay, deltaz, x1, y1, z1,             &
	                                                     xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp, &
                                                            ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
                                                            ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
                                                            ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
                                                            ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
							           ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
                                                            ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww, &
                                                            ca1_uv, ca2_uv, ca3_uv, ca4_uv, ca5_uv, ca6_uv, & 
                                                            ck1_uv, ck2_uv, ck3_uv, ck4_uv, ck5_uv, ck6_uv, &
                                                            ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, & 
                                                            ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
                                                            ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
                                                            ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
							           ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
                                                            ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
							           ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
                                                            ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
                                                            ca1_wv, ca2_wv, ca3_wv, ca4_wv, ca5_wv, ca6_wv, &
                                                            ck1_wv, ck2_wv, ck3_wv, ck4_wv, ck5_wv, ck6_wv
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:, :)    :: Ac, Acx, Acy, Acz
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:, :, :) :: b, u, ut, u_sum, u_avg,   &
                                                            v, vt, v_sum, v_avg,      &
                                                            w, wt, w_sum, w_avg,      &
                                                            p, p_sum, p_avg, pc, pco, &
                                                            resi_u, resi_v, resi_w
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:, :) :: interceptedIndexPtr, fluidIndexPtr, solidIndexPtr,     &
                                                            redCellIndexPtr, blackCellIndexPtr, nodeId, fluidInterceptedIndexPtr
       ! ibm variables
       INTEGER  :: bcType
       REAL (KIND = 8)    :: xShift, yShift, zshift, aoa, piv_pt, var_surf, thetaDot, thetaDDot, piv_x, piv_y, piv_z, aoa1
       INTEGER (KIND=8)   :: surGeoPoints, ibElems, ibNodes 
       INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:)   :: ibElP1, ibElP2, ibElP3, nelp, nelu1, nelu2, nelv1, nelv2, nelw1, nelw2
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: xnode, ynode,  znode, xnode1, ynode1, znode1, xcent, ycent, zcent, &
                                                      cosAlpha, cosBeta, cosGamma, alpha3, beta3, gamma3, bcSurf, &
                                                      pNormDis, u1NormDis, u2NormDis , v1NormDis, v2NormDis, w1NormDis, w2NormDis


       !!$acc commands..
       !$acc declare create (deltax, deltay, deltaz)
       !$acc declare create (x1, y1, z1, xu, yu, zu, xv, yv, zv, xw, yw, zw, xp, yp, zp)
       !$acc declare create (b, Acx, Acy, Acz)
       !$acc declare create (pc, pco)
       !$acc declare create (xnode, ynode, znode, ibElP1, ibElP2, ibElP3)
       !$acc declare create (xnode1, ynode1, znode1)
       !$acc declare create (xcent, ycent, zcent)
       !$acc declare create (cosAlpha, cosBeta, cosGamma)
       !$acc declare create (alpha3, beta3, gamma3)
       !$acc declare create (cell)
       !$acc declare create (minElemcell)
       !$acc declare create (fluidIndexPtr, redCellIndexPtr, blackCellIndexPtr)
       !$acc declare create (interceptedIndexPtr, solidIndexPtr)
       !$acc declare create (pNormDis, u1NormDis, u2NormDis, v1NormDis, &
       !$acc                           v2NormDis, w1NormDis, w2NormDis)
       !$acc declare create (nelp, nelu1, nelu2, nelv1, nelv2, nelw1, nelw2)
       !$acc declare create (u, ut, v, vt, w, wt, p)
       !$acc declare create (resi_u, resi_v, resi_w)
       !$acc declare create (                                &
       !$acc ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
       !$acc ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
       !$acc ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
       !$acc ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
       !$acc ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
       !$acc ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww, &
       !$acc ca1_uv, ca2_uv, ca3_uv, ca4_uv, ca5_uv, ca6_uv, &
       !$acc ck1_uv, ck2_uv, ck3_uv, ck4_uv, ck5_uv, ck6_uv, &
       !$acc ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
       !$acc ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
       !$acc ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
       !$acc ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
       !$acc ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
       !$acc ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
       !$acc ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
       !$acc ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
       !$acc ca1_wv, ca2_wv, ca3_wv, ca4_wv, ca5_wv, ca6_wv, &
       !$acc ck1_wv, ck2_wv, ck3_wv, ck4_wv, ck5_wv, ck6_wv)         
       !$acc declare create (Elemcell, ucell, vcell, wcell, pcell)      
       END MODULE global
                
                
                
                
