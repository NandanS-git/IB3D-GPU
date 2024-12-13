!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	SUBROUTINE non_uni_coeff
       USE global

	!write(*,*)'entered non_uni_coeff'
	IMPLICIT NONE
	INTEGER           :: i, j, k
	REAL (KIND = 8)   :: f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8,                             &
                         s51, s52, s53, s54, s55, s56, s57, s61, s62, s63, s64, s65, s66, s67,  &
                         s71, s72, s73, s74, s75, s76, s81, s82, s83, s84, s85, s86, s91, s92,  &
                         s93, s94, s95, s96, s97, s101, s102, s103, s104, s105, s106, s107,     &
                         ak_1, ak_2, ak_3, ak_4, ak_5, ak_6, ak_7, theta_1, theta_2, theta_3,   &
                         tmp_dx1, tmp_dx2, tmp_dx3, tmp_dx4, tmp_dy1, tmp_dy2, tmp_dy3,         &
                         tmp_dy4, tmp_dz1, tmp_dz2, tmp_dz3, tmp_dz4
                         
       ALLOCATE(ca1_uu(nx+2), ca2_uu(nx+2), ca3_uu(nx+2), ca4_uu(nx+2), ca5_uu(nx+2), ca6_uu(nx+2), &
	ck1_uu(nx+2), ck2_uu(nx+2), ck3_uu(nx+2), ck4_uu(nx+2), ck5_uu(nx+2), ck6_uu(nx+2), &
	ca1_vv(ny+2), ca2_vv(ny+2), ca3_vv(ny+2), ca4_vv(ny+2), ca5_vv(ny+2), ca6_vv(ny+2), &
	ck1_vv(ny+2), ck2_vv(ny+2), ck3_vv(ny+2), ck4_vv(ny+2), ck5_vv(ny+2), ck6_vv(ny+2), &
	ca1_ww(nz+2), ca2_ww(nz+2), ca3_ww(nz+2), ca4_ww(nz+2), ca5_ww(nz+2), ca6_ww(nz+2), &
	ck1_ww(nz+2), ck2_ww(nz+2), ck3_ww(nz+2), ck4_ww(nz+2), ck5_ww(nz+2), ck6_ww(nz+2), &
	ca1_uv(nx+2), ca2_uv(nx+2), ca3_uv(nx+2), ca4_uv(nx+2), ca5_uv(nx+2), ca6_uv(nx+2), &
	ck1_uv(nx+2), ck2_uv(nx+2), ck3_uv(nx+2), ck4_uv(nx+2), ck5_uv(nx+2), ck6_uv(nx+2), &
	ca1_uw(nx+2), ca2_uw(nx+2), ca3_uw(nx+2), ca4_uw(nx+2), ca5_uw(nx+2), ca6_uw(nx+2), &
	ck1_uw(nx+2), ck2_uw(nx+2), ck3_uw(nx+2), ck4_uw(nx+2), ck5_uw(nx+2), ck6_uw(nx+2), &
	ca1_vu(ny+2), ca2_vu(ny+2), ca3_vu(ny+2), ca4_vu(ny+2), ca5_vu(ny+2), ca6_vu(ny+2), &
	ck1_vu(ny+2), ck2_vu(ny+2), ck3_vu(ny+2), ck4_vu(ny+2), ck5_vu(ny+2), ck6_vu(ny+2), &
	ca1_vw(ny+2), ca2_vw(ny+2), ca3_vw(ny+2), ca4_vw(ny+2), ca5_vw(ny+2), ca6_vw(ny+2), &
	ck1_vw(ny+2), ck2_vw(ny+2), ck3_vw(ny+2), ck4_vw(ny+2), ck5_vw(ny+2), ck6_vw(ny+2), &
	ca1_wu(nz+2), ca2_wu(nz+2), ca3_wu(nz+2), ca4_wu(nz+2), ca5_wu(nz+2), ca6_wu(nz+2), &
	ck1_wu(nz+2), ck2_wu(nz+2), ck3_wu(nz+2), ck4_wu(nz+2), ck5_wu(nz+2), ck6_wu(nz+2), &
	ca1_wv(nz+2), ca2_wv(nz+2), ca3_wv(nz+2), ca4_wv(nz+2), ca5_wv(nz+2), ca6_wv(nz+2), &
	ck1_wv(nz+2), ck2_wv(nz+2), ck3_wv(nz+2), ck4_wv(nz+2), ck5_wv(nz+2), ck6_wv(nz+2))
	
	do i=2,nx
	theta_1=deltax(i)/deltax(i-1)
	theta_2=deltax(i+1)/deltax(i)
	theta_3=deltax(i+2)/deltax(i+1)

	f_1=theta_3
	f_2=2.0*theta_3+theta_3**2.0
	f_3=3.0*theta_3+3.0*theta_3**2.0+theta_3**3.0
	f_4=4.0*theta_3+6.0*theta_3**2.0+4.0*theta_3**3.0+theta_3**4.0
	f_5=1.0/(theta_1*theta_2)
	f_6=(2.0*theta_1+1.0)/((theta_1**2.0)*(theta_2**2.0))
	f_7=(3.0*theta_1**2.0+3.0*theta_1+1.0)/((theta_1**3.0)*(theta_2**3.0))
       f_8=(4.0*theta_1**3.0+6.0*theta_1**2.0+4.0*theta_1+1.0)/((theta_1**4.0)*(theta_2**4.0))

	s51=-1.0/(theta_2**4.0)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0)-f_2/(theta_2**4.0))
	s57=-(f_3/(theta_2**4.0)+f_4/(theta_2**3.0))

	s61=-1.0
	s62=f_8
	s63=f_8
	s64=-1.0
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0
	s72=(1.0+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0/(theta_2**4.0)
	s82=(f_8+(1.0/theta_2**4.0))
	s83=f_8
	s84=(f_5/(theta_2**4.0)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0)-f_6/(theta_2**4.0))
	s86=(f_7/(theta_2**4.0)-f_8/(theta_2**3.0))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	ca1_uu(i)=(s97*s101-s107*s91)
	ca2_uu(i)=(s97*s102+s107*s93)
	ca3_uu(i)=(s95*s107-s105*s97)
	ca4_uu(i)=(s94*s107+s104*s97)
	ca5_uu(i)=(s97*s103-s107*s92)
	ca6_uu(i)=(s97*s106-s107*s96)

	ak_1=(1.0+2.0*theta_1)*(theta_3+theta_3**2.0)*theta_2
	ak_2=(1.0+theta_1)*(2.0*theta_3+theta_3**2.0)
	ak_3=(1.0+theta_1)
	ak_4=(1.0+theta_1)*((1.0+theta_3)**2.0)
	ak_5=((1.0+theta_1)**2.0)*(theta_3+theta_3**2.0)*theta_2
	ak_6=(theta_1**2.0)*theta_2*(theta_3+theta_3**2.0)
	ak_7=(1.0+theta_1)*(theta_3+theta_3**2.0)*theta_2

	ck1_uu(i)=ak_3
	ck2_uu(i)=-ak_4
	ck3_uu(i)=(ak_1+ak_2)
	ck4_uu(i)=-ak_5
	ck5_uu(i)=ak_6
	ck6_uu(i)=ak_7
	enddo

	do j=2,ny
	theta_1=deltay(j)/deltay(j-1)
	theta_2=deltay(j+1)/deltay(j)
	theta_3=deltay(j+2)/deltay(j+1)

	f_1=theta_3
	f_2=2.0*theta_3+theta_3**2.0
	f_3=3.0*theta_3+3.0*theta_3**2.0+theta_3**3.0
	f_4=4.0*theta_3+6.0*theta_3**2.0+4.0*theta_3**3.0+theta_3**4.0
	f_5=1.0/(theta_1*theta_2)
	f_6=(2.0*theta_1+1.0)/((theta_1**2.0)*(theta_2**2.0))
	f_7=(3.0*theta_1**2.0+3.0*theta_1+1.0)/((theta_1**3.0)*(theta_2**3.0))
       f_8=(4.0*theta_1**3.0+6.0*theta_1**2.0+4.0*theta_1+1.0)/((theta_1**4.0)*(theta_2**4.0))

	s51=-1.0/(theta_2**4.0)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0)-f_2/(theta_2**4.0))
	s57=-(f_3/(theta_2**4.0)+f_4/(theta_2**3.0))

	s61=-1.0
	s62=f_8
	s63=f_8
	s64=-1.0
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0
	s72=(1.0+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0/(theta_2**4.0)
	s82=(f_8+(1.0/theta_2**4.0))
	s83=f_8
	s84=(f_5/(theta_2**4.0)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0)-f_6/(theta_2**4.0))
	s86=(f_7/(theta_2**4.0)-f_8/(theta_2**3.0))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	ca1_vv(j)=(s97*s101-s107*s91)
	ca2_vv(j)=(s97*s102+s107*s93)
	ca3_vv(j)=(s95*s107-s105*s97)
	ca4_vv(j)=(s94*s107+s104*s97)
	ca5_vv(j)=(s97*s103-s107*s92)
	ca6_vv(j)=(s97*s106-s107*s96)

	ak_1=(1.0+2.0*theta_1)*(theta_3+theta_3**2.0)*theta_2
	ak_2=(1.0+theta_1)*(2.0*theta_3+theta_3**2.0)
	ak_3=(1.0+theta_1)
	ak_4=(1.0+theta_1)*((1.0+theta_3)**2.0)
	ak_5=((1.0+theta_1)**2.0)*(theta_3+theta_3**2.0)*theta_2
	ak_6=(theta_1**2.0)*theta_2*(theta_3+theta_3**2.0)
	ak_7=(1.0+theta_1)*(theta_3+theta_3**2.0)*theta_2

	ck1_vv(j)=ak_3
	ck2_vv(j)=-ak_4
	ck3_vv(j)=(ak_1+ak_2)
	ck4_vv(j)=-ak_5
	ck5_vv(j)=ak_6
	ck6_vv(j)=ak_7
	enddo

	do k=2,nz
	theta_1=deltaz(k)/deltaz(k-1)
	theta_2=deltaz(k+1)/deltaz(k)
	theta_3=deltaz(k+2)/deltaz(k+1)

	f_1=theta_3
	f_2=2.0*theta_3+theta_3**2.0
	f_3=3.0*theta_3+3.0*theta_3**2.0+theta_3**3.0
	f_4=4.0*theta_3+6.0*theta_3**2.0+4.0*theta_3**3.0+theta_3**4.0
	f_5=1.0/(theta_1*theta_2)
	f_6=(2.0*theta_1+1.0)/((theta_1**2.0)*(theta_2**2.0))
	f_7=(3.0*theta_1**2.0+3.0*theta_1+1.0)/((theta_1**3.0)*(theta_2**3.0))
       f_8=(4.0*theta_1**3.0+6.0*theta_1**2.0+4.0*theta_1+1.0)/((theta_1**4.0)*(theta_2**4.0))

	s51=-1.0/(theta_2**4.0)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0)-f_2/(theta_2**4.0))
	s57=-(f_3/(theta_2**4.0)+f_4/(theta_2**3.0))

	s61=-1.0
	s62=f_8
	s63=f_8
	s64=-1.0
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0
	s72=(1.0+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0/(theta_2**4.0)
	s82=(f_8+(1.0/theta_2**4.0))
	s83=f_8
	s84=(f_5/(theta_2**4.0)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0)-f_6/(theta_2**4.0))
	s86=(f_7/(theta_2**4.0)-f_8/(theta_2**3.0))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	ca1_ww(k)=(s97*s101-s107*s91)
	ca2_ww(k)=(s97*s102+s107*s93)
	ca3_ww(k)=(s95*s107-s105*s97)
	ca4_ww(k)=(s94*s107+s104*s97)
	ca5_ww(k)=(s97*s103-s107*s92)
	ca6_ww(k)=(s97*s106-s107*s96)

	ak_1=(1.0+2.0*theta_1)*(theta_3+theta_3**2.0)*theta_2
	ak_2=(1.0+theta_1)*(2.0*theta_3+theta_3**2.0)
	ak_3=(1.0+theta_1)
	ak_4=(1.0+theta_1)*((1.0+theta_3)**2.0)
	ak_5=((1.0+theta_1)**2.0)*(theta_3+theta_3**2.0)*theta_2
	ak_6=(theta_1**2.0)*theta_2*(theta_3+theta_3**2.0)
	ak_7=(1.0+theta_1)*(theta_3+theta_3**2.0)*theta_2

	ck1_ww(k)=ak_3
	ck2_ww(k)=-ak_4
	ck3_ww(k)=(ak_1+ak_2)
	ck4_ww(k)=-ak_5
	ck5_ww(k)=ak_6
	ck6_ww(k)=ak_7
	enddo

	do i=2,nx
	if(i.eq.2)then
	tmp_dx1=deltax(i-1)
	else
	tmp_dx1=0.5*(deltax(i-1)+deltax(i-2))
	endif
	tmp_dx2=0.5*(deltax(i)+deltax(i-1))
	tmp_dx3=0.5*(deltax(i)+deltax(i+1))
	tmp_dx4=0.5*(deltax(i+1)+deltax(i+2))

	theta_1=tmp_dx2/tmp_dx1
	theta_2=tmp_dx3/tmp_dx2
	theta_3=tmp_dx4/tmp_dx3

	f_1=theta_3
	f_2=2.0*theta_3+theta_3**2.0
	f_3=3.0*theta_3+3.0*theta_3**2.0+theta_3**3.0
	f_4=4.0*theta_3+6.0*theta_3**2.0+4.0*theta_3**3.0+theta_3**4.0
	f_5=1.0/(theta_1*theta_2)
	f_6=(2.0*theta_1+1.0)/((theta_1**2.0)*(theta_2**2.0))
	f_7=(3.0*theta_1**2.0+3.0*theta_1+1.0)/((theta_1**3.0)*(theta_2**3.0))
       f_8=(4.0*theta_1**3.0+6.0*theta_1**2.0+4.0*theta_1+1.0)/((theta_1**4.0)*(theta_2**4.0))

	s51=-1.0/(theta_2**4.0)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0)-f_2/(theta_2**4.0))
	s57=-(f_3/(theta_2**4.0)+f_4/(theta_2**3.0))

	s61=-1.0
	s62=f_8
	s63=f_8
	s64=-1.0
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0
	s72=(1.0+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0/(theta_2**4.0)
	s82=(f_8+(1.0/theta_2**4.0))
	s83=f_8
	s84=(f_5/(theta_2**4.0)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0)-f_6/(theta_2**4.0))
	s86=(f_7/(theta_2**4.0)-f_8/(theta_2**3.0))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	ca1_uv(i)=(s97*s101-s107*s91)
	ca1_uw(i)=ca1_uv(i)
	ca2_uv(i)=(s97*s102+s107*s93)
	ca2_uw(i)=ca2_uv(i)
	ca3_uv(i)=(s95*s107-s105*s97)
	ca3_uw(i)=ca3_uv(i)
	ca4_uv(i)=(s94*s107+s104*s97)
	ca4_uw(i)=ca4_uv(i)
	ca5_uv(i)=(s97*s103-s107*s92)
	ca5_uw(i)=ca5_uv(i)
	ca6_uv(i)=(s97*s106-s107*s96)
	ca6_uw(i)=ca6_uv(i)

	ak_1=(1.0+2.0*theta_1)*(theta_3+theta_3**2.0)*theta_2
	ak_2=(1.0+theta_1)*(2.0*theta_3+theta_3**2.0)
	ak_3=(1.0+theta_1)
	ak_4=(1.0+theta_1)*((1.0+theta_3)**2.0)
	ak_5=((1.0+theta_1)**2.0)*(theta_3+theta_3**2.0)*theta_2
	ak_6=(theta_1**2.0)*theta_2*(theta_3+theta_3**2.0)
	ak_7=(1.0+theta_1)*(theta_3+theta_3**2.0)*theta_2

	ck1_uv(i)=ak_3
	ck1_uw(i)=ck1_uv(i)
	ck2_uv(i)=-ak_4
	ck2_uw(i)=ck2_uv(i)
	ck3_uv(i)=(ak_1+ak_2)
	ck3_uw(i)=ck3_uv(i)
	ck4_uv(i)=-ak_5
	ck4_uw(i)=ck4_uv(i)
	ck5_uv(i)=ak_6
	ck5_uw(i)=ck5_uv(i)
	ck6_uv(i)=ak_7
	ck6_uw(i)=ck6_uv(i)
	enddo

	do j=2,ny
	if(j.eq.2)then
	tmp_dy1=deltay(j-1)
	else
	tmp_dy1=0.5*(deltay(j-1)+deltay(j-2))
	endif
	tmp_dy2=0.5*(deltay(j)+deltay(j-1))
	tmp_dy3=0.5*(deltay(j)+deltay(j+1))
	tmp_dy4=0.5*(deltay(j+1)+deltay(j+2))

	theta_1=tmp_dy2/tmp_dy1
	theta_2=tmp_dy3/tmp_dy2
	theta_3=tmp_dy4/tmp_dy3

	f_1=theta_3
	f_2=2.0*theta_3+theta_3**2.0
	f_3=3.0*theta_3+3.0*theta_3**2.0+theta_3**3.0
	f_4=4.0*theta_3+6.0*theta_3**2.0+4.0*theta_3**3.0+theta_3**4.0
	f_5=1.0/(theta_1*theta_2)
	f_6=(2.0*theta_1+1.0)/((theta_1**2.0)*(theta_2**2.0))
	f_7=(3.0*theta_1**2.0+3.0*theta_1+1.0)/((theta_1**3.0)*(theta_2**3.0))
       f_8=(4.0*theta_1**3.0+6.0*theta_1**2.0+4.0*theta_1+1.0)/((theta_1**4.0)*(theta_2**4.0))

	s51=-1.0/(theta_2**4.0)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0)-f_2/(theta_2**4.0))
	s57=-(f_3/(theta_2**4.0)+f_4/(theta_2**3.0))

	s61=-1.0
	s62=f_8
	s63=f_8
	s64=-1.0
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0
	s72=(1.0+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0/(theta_2**4.0)
	s82=(f_8+(1.0/theta_2**4.0))
	s83=f_8
	s84=(f_5/(theta_2**4.0)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0)-f_6/(theta_2**4.0))
	s86=(f_7/(theta_2**4.0)-f_8/(theta_2**3.0))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	ca1_vu(j)=(s97*s101-s107*s91)
	ca1_vw(j)=ca1_vu(j)
	ca2_vu(j)=(s97*s102+s107*s93)
	ca2_vw(j)=ca2_vu(j)
	ca3_vu(j)=(s95*s107-s105*s97)
	ca3_vw(j)=ca3_vu(j)
	ca4_vu(j)=(s94*s107+s104*s97)
	ca4_vw(j)=ca4_vu(j)
	ca5_vu(j)=(s97*s103-s107*s92)
	ca5_vw(j)=ca5_vu(j)
	ca6_vu(j)=(s97*s106-s107*s96)
	ca6_vw(j)=ca6_vu(j)

	ak_1=(1.0+2.0*theta_1)*(theta_3+theta_3**2.0)*theta_2
	ak_2=(1.0+theta_1)*(2.0*theta_3+theta_3**2.0)
	ak_3=(1.0+theta_1)
	ak_4=(1.0+theta_1)*((1.0+theta_3)**2.0)
	ak_5=((1.0+theta_1)**2.0)*(theta_3+theta_3**2.0)*theta_2
	ak_6=(theta_1**2.0)*theta_2*(theta_3+theta_3**2.0)
	ak_7=(1.0+theta_1)*(theta_3+theta_3**2.0)*theta_2

	ck1_vu(j)=ak_3
	ck1_vw(j)=ck1_vu(j)
	ck2_vu(j)=-ak_4
	ck2_vw(j)=ck2_vu(j)
	ck3_vu(j)=(ak_1+ak_2)
	ck3_vw(j)=ck3_vu(j)
	ck4_vu(j)=-ak_5
	ck4_vw(j)=ck4_vu(j)
	ck5_vu(j)=ak_6
	ck5_vw(j)=ck5_vu(j)
	ck6_vu(j)=ak_7
	ck6_vw(j)=ck6_vu(j)
	enddo

	do k=2,nz
	if(k.eq.2)then
	tmp_dz1=deltaz(k-1)
	else
	tmp_dz1=0.5*(deltaz(k-1)+deltaz(k-2))
	endif
	tmp_dz2=0.5*(deltaz(k)+deltaz(k-1))
	tmp_dz3=0.5*(deltaz(k)+deltaz(k+1))
	tmp_dz4=0.5*(deltaz(k+1)+deltaz(k+2))

	theta_1=tmp_dz2/tmp_dz1
	theta_2=tmp_dz3/tmp_dz2
	theta_3=tmp_dz4/tmp_dz3

	f_1=theta_3
	f_2=2.0*theta_3+theta_3**2.0
	f_3=3.0*theta_3+3.0*theta_3**2.0+theta_3**3.0
	f_4=4.0*theta_3+6.0*theta_3**2.0+4.0*theta_3**3.0+theta_3**4.0
	f_5=1.0/(theta_1*theta_2)
	f_6=(2.0*theta_1+1.0)/((theta_1**2.0)*(theta_2**2.0))
	f_7=(3.0*theta_1**2.0+3.0*theta_1+1.0)/((theta_1**3.0)*(theta_2**3.0)) 	
	f_8=(4.0*theta_1**3.0+6.0*theta_1**2.0+4.0*theta_1+1.0)/((theta_1**4.0)*(theta_2**4.0))

	s51=-1.0/(theta_2**4.0)
	s52=f_4
	s53=f_4
	s54=s51
	s55=-(f_1/(theta_2**4.0)+f_4/theta_2)
	s56=(f_4/(theta_2**2.0)-f_2/(theta_2**4.0))
	s57=-(f_3/(theta_2**4.0)+f_4/(theta_2**3.0))

	s61=-1.0
	s62=f_8
	s63=f_8
	s64=-1.0
	s65=(f_5+f_8)
	s66=(f_8-f_6)
	s67=(f_7+f_8)

	s71=-1.0
	s72=(1.0+f_4)
	s73=f_4
	s74=(f_4-f_1)
	s75=(f_4-f_2)
	s76=(f_4-f_3)

	s81=-1.0/(theta_2**4.0)
	s82=(f_8+(1.0/theta_2**4.0))
	s83=f_8
	s84=(f_5/(theta_2**4.0)-f_8/theta_2)
	s85=(f_8/(theta_2**2.0)-f_6/(theta_2**4.0))
	s86=(f_7/(theta_2**4.0)-f_8/(theta_2**3.0))

	s91=-s66*s51
	s92=s56*s61
	s93=-(s54*s66+s56*s62)
	s94=(s56*s64+s66*s52)
	s95=(s56*s63-s66*s53)
	s96=(s56*s65-s66*s55)
	s97=(s56*s67-s66*s57)

	s101=-s85*s71
	s102=-s85*s72
	s103=s75*s81
	s104=s75*s82
	s105=(s75*s83-s85*s73)
	s106=(s75*s84-s85*s74)
	s107=(s75*s86-s85*s76)

	ca1_wu(k)=(s97*s101-s107*s91)
	ca1_wv(k)=ca1_wu(k)
	ca2_wu(k)=(s97*s102+s107*s93)
	ca2_wv(k)=ca2_wu(k)
	ca3_wu(k)=(s95*s107-s105*s97)
	ca3_wv(k)=ca3_wu(k)
	ca4_wu(k)=(s94*s107+s104*s97)
	ca4_wv(k)=ca4_wu(k)
	ca5_wu(k)=(s97*s103-s107*s92)
	ca5_wv(k)=ca5_wu(k)
	ca6_wu(k)=(s97*s106-s107*s96)
	ca6_wv(k)=ca6_wu(k)

	ak_1=(1.0+2.0*theta_1)*(theta_3+theta_3**2.0)*theta_2
	ak_2=(1.0+theta_1)*(2.0*theta_3+theta_3**2.0)
	ak_3=(1.0+theta_1)
	ak_4=(1.0+theta_1)*((1.0+theta_3)**2.0)
	ak_5=((1.0+theta_1)**2.0)*(theta_3+theta_3**2.0)*theta_2
	ak_6=(theta_1**2.0)*theta_2*(theta_3+theta_3**2.0)
	ak_7=(1.0+theta_1)*(theta_3+theta_3**2.0)*theta_2

	ck1_wu(k)=ak_3
	ck1_wv(k)=ck1_wu(k)
	ck2_wu(k)=-ak_4
	ck2_wv(k)=ck2_wu(k)
	ck3_wu(k)=(ak_1+ak_2)
	ck3_wv(k)=ck3_wu(k)
	ck4_wu(k)=-ak_5
	ck4_wv(k)=ck4_wu(k)
	ck5_wu(k)=ak_6
	ck5_wv(k)=ck5_wu(k)
	ck6_wu(k)=ak_7
	ck6_wv(k)=ck6_wu(k)
	enddo

       !$acc update device (                                 &
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
	
	write(*,*)'leaving non_uni_coeff'

	END SUBROUTINE non_uni_coeff

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE nsMomentum
!***********************************************************************
!     navier-stokes equations for constant properties                                              
!***********************************************************************
         USE global
!***********************************************************************
      !write(*,*)'has entered nseqcp'   
	  IMPLICIT NONE
	  INTEGER          :: i, j, k, n 
	  REAL (KIND = 8)  :: dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,            &
	                      u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
			        w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
			        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
			        dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          & 
			        r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
			        dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
			        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,            &
			        dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
			        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
			        duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw		  
			  
	  al = 1.

!$acc parallel loop gang vector  &
!$acc private (i, j, k, dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,   &
!$acc	        u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
!$acc		 w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
!$acc	        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
!$acc		 dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          & 
!$acc		 r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
!$acc	        dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
!$acc	        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,            &
!$acc	        dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
!$acc	        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
!$acc		 duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw) &
!$acc present (fluidIndexPtr, deltax, deltay, deltaz, u, v, w,             &
!$acc          ut, vt, wt, p, cell, x1, y1, z1, xp, yp, zp,                &
!$acc          xu, yu, zu, xv, yv, zv, xw, yw, zw, resi_u, resi_v, resi_w, &
!$acc          ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
!$acc          ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
!$acc          ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
!$acc          ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
!$acc          ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
!$acc          ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
!$acc          ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
!$acc          ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
!$acc          ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
!$acc          ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
!$acc          ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
!$acc          ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
!$acc          ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
!$acc          ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
!$acc          ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
!$acc          ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
!$acc          ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
!$acc          ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww) &
!$acc firstprivate(nx, ny, nz, rev, deltat, al)  
	DO n = 1, fluidCellCount
	
       i = fluidIndexPtr(n, 1)
       j = fluidIndexPtr(n, 2)
       k = fluidIndexPtr(n, 3) 
 
	dxr=deltax(i+1)
	dx=deltax(i)
	dxl=deltax(i-1)
	dye=deltay(j+1)
	dy=deltay(j)
	dyw=deltay(j-1)
	dzt=deltaz(k+1)
	dz=deltaz(k)
	dzb=deltaz(k-1)
!cccccccccccccccccccccccccc  diff-u     ccccccccccccccccccccccccccccccccc
!     duu / dx
       u1a = u(i-1,j,k) + u(i,j,k)
       u22 = u(i-1,j,k) - u(i,j,k)
       u3  = u(i,j,k)   + u(i+1,j,k)
       u4  = u(i,j,k)   - u(i+1,j,k)

!     duv  / dy
       u5 = u(i,j-1,k)  + u(i,j,k)
       u6 = u(i,j-1,k)  - u(i,j,k)
       u7 = u(i,j,k)    + u(i,j+1,k)
       u8 = u(i,j,k)    - u(i,j+1,k)

!     dwu  /  dz
       u9 = u(i,j,k-1)  + u(i,j,k)
       u10= u(i,j,k-1)  - u(i,j,k)
       u11= u(i,j,k)    + u(i,j,k+1)
       u12= u(i,j,k)    - u(i,j,k+1)

!     duv/dx
       u13 = u(i-1,j,k) + u(i-1,j+1,k)
       u14 = u7

!     duw / dx
       u15 = u(i-1,j,k) + u(i-1,j,k+1)
       u16 = u11
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccc    diff -v    ccccccccccccccccccccccccccccccc
!     dvu / dx
       v1a = v(i,j-1,k)  + v(i+1,j-1,k)     
       v22 = v(i,j,k)    + v(i+1,j,k) 

!     duv / dx
       v3 = v(i-1,j,k)   + v(i,j,k)
       v4 = v(i-1,j,k)   - v(i,j,k)
       v5 = v(i,j,k)     + v(i+1,j,k)
       v6 = v(i,j,k)     - v(i+1,j,k)

!     dvv / dy
       v7 = v(i,j-1,k)   + v(i,j,k)
       v8 = v(i,j-1,k)   - v(i,j,k)
       v9 = v(i,j,k)     + v(i,j+1,k)
       v10= v(i,j,k)     - v(i,j+1,k)  

!     dwu / dz
       v11 = v(i,j,k-1)  + v(i,j,k)
       v12 = v(i,j,k-1)  - v(i,j,k) 
       v13 = v(i,j,k)    + v(i,j,k+1)
       v14 = v(i,j,k)    - v(i,j,k+1)    

!     dvw / dy
       v15 = v(i,j-1,k)  + v(i,j-1,k+1)
       v16 = v13
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccc  diff - w  cccccccccccccccccccccccccccccccc
!     dwu / dz
       w1a = w(i,j,k-1)  + w(i+1,j,k-1)
       w22 = w(i,j,k)    + w(i+1,j,k)

!     dwv / dz     
       w3 = w(i,j,k-1)   + w(i,j+1,k-1)
       w4 = w(i,j,k)     + w(i,j+1,k)

!     duw / dx   
       w5 = w(i-1,j,k)   + w(i,j,k)
       w6 = w(i-1,j,k)   - w(i,j,k)
       w7 = w(i,j,k)     + w(i+1,j,k)
       w8 = w(i,j,k)     - w(i+1,j,k)

!     dvw / dy
       w9 = w(i,j-1,k)   + w(i,j,k)
       w10 = w(i,j-1,k)  - w(i,j,k)
       w11 = w(i,j,k)    + w(i,j+1,k)
       w12 = w(i,j,k)    - w(i,j+1,k)

!     dww / dz 
       w13 = w(i,j,k-1)  + w(i,j,k)
       w14 = w(i,j,k-1)  - w(i,j,k)
       w15 = w(i,j,k)    + w(i,j,k+1)
       w16 = w(i,j,k)    - w(i,j,k+1)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       dpdx = (p(i,j,k) - p(i+1,j,k))/(0.5*(dxr+dx))
       dpdy = (p(i,j,k) - p(i,j+1,k))/(0.5*(dye+dy))
       dpdz = (p(i,j,k) - p(i,j,k+1))/(0.5*(dzt+dz))
       
    	dx2xr = dx + dxr
       dx2xl = dx + dxl
	dy2ye = dy + dye
	dy2yw = dy + dyw
	dz2zt = dz + dzt
	dz2zb = dz + dzb
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c*********************** U - Momentum **********************************
	vu_e=v(i+1,j,k)+(dxr/(dx2xr))*(v(i,j,k)-v(i+1,j,k))
	vu_w=v(i+1,j-1,k)+(dxr/(dx2xr))*(v(i,j-1,k)-v(i+1,j-1,k))
	v_in_um=0.5*(vu_e+vu_w)
		 
	wu_n=w(i+1,j,k)+(dxr/(dx2xr))*(w(i,j,k)-w(i+1,j,k))
	wu_s=w(i+1,j,k-1)+(dxr/(dx2xr))*(w(i,j,k-1)-w(i+1,j,k-1))
	w_in_um=0.5*(wu_n+wu_s)
	
!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
       if(i.ne.2.and.i.lt.nx.and.j.ne.2.and.j.lt.ny+1.and.    &
       k.ne.2.and.k.lt.nz+1.and.cell(i+1,j,k).ne.2.and.       &
       cell(i-1,j,k).ne.2.and.cell(i,j+1,k).ne.2.and.         &
       cell(i,j-1,k).ne.2.and.cell(i,j,k+1).ne.2.and.         &
       cell(i,j,k-1).ne.2) then
 
	ddy=0.5*(deltay(j)+deltay(j-1))
       ddye=0.5*(deltay(j)+deltay(j+1))
       ddz=0.5*(deltaz(k)+deltaz(k-1))
       ddzr=0.5*(deltaz(k)+deltaz(k+1))

       duutdx=u(i,j,k)*(ca1_uu(i)*u(i+2,j,k)+ca2_uu(i)*u(i+1,j,k) &
       +ca3_uu(i)*u(i,j,k)+ca4_uu(i)*u(i-1,j,k)+ca5_uu(i)* &
       u(i-2,j,k))/(ca6_uu(i)*deltax(i+1))+dabs(u(i,j,k))* &
       (ck1_uu(i)*u(i+2,j,k)+ck2_uu(i)*u(i+1,j,k) &
       +ck3_uu(i)*u(i,j,k)+ck4_uu(i)*u(i-1,j,k)+ck5_uu(i)* &
       u(i-2,j,k))/(2.0*ck6_uu(i)*deltax(i))

       dvutdy=v_in_um*(ca1_vu(j)*u(i,j+2,k)+ca2_vu(j)*u(i,j+1,k) &
       +ca3_vu(j)*u(i,j,k)+ca4_vu(j)*u(i,j-1,k)+ca5_vu(j)* &
       u(i,j-2,k))/(ca6_vu(j)*ddye)+dabs(v_in_um)* &
       (ck1_vu(j)*u(i,j+2,k)+ck2_vu(j)*u(i,j+1,k) &
       +ck3_vu(j)*u(i,j,k)+ck4_vu(j)*u(i,j-1,k)+ck5_vu(j)* &
       u(i,j-2,k))/(2.0*ck6_vu(j)*ddy) 

       dwutdz=w_in_um*(ca1_wu(k)*u(i,j,k+2)+ca2_wu(k)*u(i,j,k+1) &
       +ca3_wu(k)*u(i,j,k)+ca4_wu(k)*u(i,j,k-1)+ca5_wu(k)* &
       u(i,j,k-2))/(ca6_wu(k)*ddzr)+dabs(w_in_um)* &
       (ck1_wu(k)*u(i,j,k+2)+ck2_wu(k)*u(i,j,k+1) &
       +ck3_wu(k)*u(i,j,k)+ck4_wu(k)*u(i,j,k-1)+ck5_wu(k)* &
       u(i,j,k-2))/(2.0*ck6_wu(k)*ddz)
     
     	duudx=duutdx
       dvudy=dvutdy
       dwudz=dwutdz
       
	else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
	r1x=deltax(i)/deltax(i+1)
       r1yn=0.5*(deltay(j)+deltay(j-1))
       r1yd=0.5*(deltay(j)+deltay(j+1))
       r1y=r1yn/r1yd
	r1zn=0.5*(deltaz(k)+deltaz(k-1))
	r1zd=0.5*(deltaz(k)+deltaz(k+1))
	r1z=r1zn/r1zd

	duuwdx=(u(i,j,k)/deltax(i+1))*((-1.0/(r1x*(r1x+1.0)))* &
     	u(i-1,j,k)-(1.0-1.0/r1x)*u(i,j,k)+   &
     	(r1x/(r1x+1.0))*u(i+1,j,k))-al*dabs(u(i,j,k))* &
     	(1.0/(2.0*deltax(i+1)))*(2.0/(r1x*(r1x+1.0))* &
     	u(i-1,j,k)-2.0/r1x*u(i,j,k)+2.0/(r1x+1.0)* &
       u(i+1,j,k))

       dvuwdy=(v_in_um/r1yd)*((-1.0/(r1y*(r1y+1.0)))* &
       u(i,j-1,k)-(1.0-1.0/r1y)*u(i,j,k)+ &
       (r1y/(r1y+1.0))*u(i,j+1,k))-al*dabs(v_in_um)* &
       (1.0/(2.0*r1yd))*(2.0/(r1y*(r1y+1.0))* &
       u(i,j-1,k)-2.0/r1y*u(i,j,k)+2.0/(r1y+1.0)* &
       u(i,j+1,k))

	dwuwdz=(w_in_um/r1zd)*((-1.0/(r1z*(r1z+1.0)))* &
     	u(i,j,k-1)-(1.0-1.0/r1z)*u(i,j,k)+ &
    	(r1z/(r1z+1.0))*u(i,j,k+1))-al*dabs(w_in_um)* &
     	(1.0/(2.0*r1zd))*(2.0/(r1z*(r1z+1.0))* &
     	u(i,j,k-1)-2.0/r1z*u(i,j,k)+2.0/(r1z+1.0)* &
      	u(i,j,k+1))
 
       duudx=duuwdx
       dvudy=dvuwdy
       dwudz=dwuwdz
        
       endif
   
!ccccccccccccccccccccccccc  grad of u part cccccccccccccccccccccccccccccc
	!d2udx2=(2.0/dx2xr)*((-u4/dxr)+(u22/dx))
       !d2udy2=(2.0/dy)*((-u8/dy2ye)+(u6/dy2yw))
	!d2udz2=(2.0/dz)*((-u12/dz2zt)+(u10/dz2zb))
	
	d2udx2=(2.0/dx2xr)*((-u4/dxr)+(u22/dx))
	
       d2udy2=(2.0/(0.5*(dy2ye+dy2yw)))*((-u8/(0.5*dy2ye))+ &
       (u6/(0.5*dy2yw)))
     
	d2udz2=(2.0/(0.5*(dz2zt+dz2zb)))*((-u12/(0.5*dz2zt))+ &
       (u10/(0.5*dz2zb)))

       xtt2=rev*(d2udx2+d2udy2+d2udz2)
       
       residu=(-duudx-dvudy-dwudz+xtt2)

       !if(ita.eq.1)then
       ut(i,j,k)=u(i,j,k)+deltat*(residu+dpdx)
       !else
       !ut(i,j,k)=u(i,j,k)+deltat*(0.5*(3.0*residu-resi_u(i,j,k))+dpdx)
       !endif

       resi_u(i,j,k)=residu
!***********************************************************************

!c*********************** V - Momentum *********************************
	uv_e=u(i,j+1,k)+(dye/(dy2ye))*(u(i,j,k)-u(i,j+1,k))
	uv_w=u(i-1,j+1,k)+(dye/(dy2ye))*(u(i-1,j,k)-u(i-1,j+1,k))
	u_in_vm=0.5*(uv_e+uv_w)
	
	wv_n=w(i,j+1,k)+(dye/(dy2ye))*(w(i,j,k)-w(i,j+1,k))
	wv_s=w(i,j+1,k-1)+(dye/(dy2ye))*(w(i,j,k-1)-w(i,j+1,k-1))
	w_in_vm=0.5*(wv_n+wv_s)

!cccccccccccccc---Third Order Upwinding ----cccccccccccccccccccccccccccc
      if(i.ne.2.and.i.lt.nx+1.and.j.ne.2.and.j.lt.ny.and. &
      k.ne.2.and.k.lt.nz+1.and.cell(i+1,j,k).ne.2.and.    &
      cell(i-1,j,k).ne.2.and.cell(i,j+1,k).ne.2.and.      &
      cell(i,j-1,k).ne.2.and.cell(i,j,k+1).ne.2.and.      &
      cell(i,j,k-1).ne.2) then
      
	ddx=0.5*(deltax(i)+deltax(i-1))
       ddxr=0.5*(deltax(i)+deltax(i+1))
       ddz=0.5*(deltaz(k)+deltaz(k-1))
       ddzr=0.5*(deltaz(k)+deltaz(k+1))

       duvtdx=u_in_vm*(ca1_uw(i)*v(i+2,j,k)+ca2_uw(i)*v(i+1,j,k) &
       +ca3_uw(i)*v(i,j,k)+ca4_uw(i)*v(i-1,j,k)+ca5_uw(i)* &
       v(i-2,j,k))/(ca6_uw(i)*ddxr)+dabs(u_in_vm)* &
       (ck1_uw(i)*v(i+2,j,k)+ck2_uw(i)*v(i+1,j,k) &
       +ck3_uw(i)*v(i,j,k)+ck4_uw(i)*v(i-1,j,k)+ck5_uw(i)* &
       v(i-2,j,k))/(2.0*ck6_uw(i)*ddx)

       dvvtdy=v(i,j,k)*(ca1_vv(j)*v(i,j+2,k)+ca2_vv(j)*v(i,j+1,k) &
       +ca3_vv(j)*v(i,j,k)+ca4_vv(j)*v(i,j-1,k)+ca5_vv(j)* &
       v(i,j-2,k))/(ca6_vv(j)*deltay(j+1))+dabs(v(i,j,k))* &
       (ck1_vv(j)*v(i,j+2,k)+ck2_vv(j)*v(i,j+1,k) &
       +ck3_vv(j)*v(i,j,k)+ck4_vv(j)*v(i,j-1,k)+ck5_vv(j)* &
       v(i,j-2,k))/(2.0*ck6_vv(j)*deltay(j))

       dwvtdz=w_in_vm*(ca1_wu(k)*v(i,j,k+2)+ca2_wu(k)*v(i,j,k+1) &
       +ca3_wu(k)*v(i,j,k)+ca4_wu(k)*v(i,j,k-1)+ca5_wu(k)* &
       v(i,j,k-2))/(ca6_wu(k)*ddzr)+dabs(w_in_vm)* &
       (ck1_wu(k)*v(i,j,k+2)+ck2_wu(k)*v(i,j,k+1) &
       +ck3_wu(k)*v(i,j,k)+ck4_wu(k)*v(i,j,k-1)+ck5_wu(k)* &
       v(i,j,k-2))/(2.0*ck6_wu(k)*ddz)
     
       duvdx=duvtdx
       dvvdy=dvvtdy
       dwvdz=dwvtdz
       
       else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
	r1xn=0.5*(deltax(i)+deltax(i-1))
	r1xd=0.5*(deltax(i)+deltax(i+1))
	r1x=r1xn/r1xd
       r1y=deltay(j)/deltay(j+1)
	r1zn=0.5*(deltaz(k)+deltaz(k-1))
	r1zd=0.5*(deltaz(k)+deltaz(k+1))
	r1z=r1zn/r1zd

	duvwdx=(u_in_vm/r1xd)*((-1.0/(r1x*(r1x+1.0)))* &
     	v(i-1,j,k)-(1.0-1.0/r1x)*v(i,j,k)+ &
     	(r1x/(r1x+1.0))*v(i+1,j,k))-al*dabs(u_in_vm)* &
     	(1.0/(2.0*r1xd))*(2.0/(r1x*(r1x+1.0))* &
     	v(i-1,j,k)-2.0/r1x*v(i,j,k)+2.0/(r1x+1.0)* &
     	v(i+1,j,k))

       dvvwdy=(v(i,j,k)/deltay(j+1))*((-1.0/(r1y*(r1y+1.0)))* &
       v(i,j-1,k)-(1.0-1.0/r1y)*v(i,j,k)+ &
       (r1y/(r1y+1.0))*v(i,j+1,k))-al*dabs(v(i,j,k))* &
       (1.0/(2.0*deltay(j+1)))*(2.0/(r1y*(r1y+1.0))* &
       v(i,j-1,k)-2.0/r1y*v(i,j,k)+2.0/(r1y+1.0)* &
       v(i,j+1,k))

	dwvwdz=(w_in_vm/r1zd)*((-1.0/(r1z*(r1z+1.0)))* &
     	v(i,j,k-1)-(1.0-1.0/r1z)*v(i,j,k)+ &
     	(r1z/(r1z+1.0))*v(i,j,k+1))-al*dabs(w_in_vm)* &
     	(1.0/(2.0*r1zd))*(2.0/(r1z*(r1z+1.0))* &
     	v(i,j,k-1)-2.0/r1z*v(i,j,k)+2.0/(r1z+1.0)* &
     	v(i,j,k+1))	 

       duvdx=duvwdx
       dvvdy=dvvwdy
       dwvdz=dwvwdz
       
       endif

!ccccccccccccccccccccccccc  grad of v part cccccccccccccccccccccccccccccc
	!d2vdx2=(2.0/dx)*((-v6/dx2xr)+(v4/dx2xl))
       !d2vdy2=(2.0/dy2ye)*((-v10/dye)+(v8/dy))
	!d2vdz2=(2.0/dz)*((-v14/dz2zt)+(v12/dz2zb))
	
	d2vdx2=(2.0/(0.5*(dx2xr+dx2xl)))*((-v6/(0.5*dx2xr))+ &
	(v4/(0.5*dx2xl)))  
     
       d2vdy2=(2.0/dy2ye)*((-v10/dye)+(v8/dy))
       
	d2vdz2=(2.0/(0.5*(dz2zt+dz2zb)))*((-v14/(0.5*dz2zt))+ &
	(v12/(0.5*dz2zb)))
		
       ytt2=rev*(d2vdx2+d2vdy2+d2vdz2)

       residv=(-duvdx-dvvdy-dwvdz+ytt2)

       !if(ita.eq.1 .or. irest.eq.1)then
       vt(i,j,k)=v(i,j,k)+deltat*(residv+dpdy)
       !else
       !vt(i,j,k)=v(i,j,k)+deltat*(0.5*(3.0*residv-resi_v(i,j,k))+dpdy)
       !endif

       resi_v(i,j,k)=residv
!c***********************************************************************

!c*********************** W - Momentum **********************************
	uw_e=u(i,j,k+1)+(dzt/(dz2zt))*(u(i,j,k)-u(i,j,k+1))
	uw_w=u(i-1,j,k+1)+(dzt/(dz2zt))*(u(i-1,j,k)-u(i-1,j,k+1))
	u_in_wm=0.5*(uw_e+uw_w)
	
	vw_n=v(i,j,k+1)+(dzt/(dz2zt))*(v(i,j,k)-v(i,j,k+1))
	vw_s=v(i,j-1,k+1)+(dzt/(dz2zt))*(v(i,j-1,k)-v(i,j-1,k+1))
	v_in_wm=0.5*(vw_n+vw_s)	

!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
      if(i.ne.2.and.i.lt.nx+1.and.j.ne.2.and.j.lt.ny+1.and. &
      k.ne.2.and.k.lt.nz.and.cell(i+1,j,k).ne.2.and.        &
      cell(i-1,j,k).ne.2.and.cell(i,j+1,k).ne.2.and.        &
      cell(i,j-1,k).ne.2.and.cell(i,j,k+1).ne.2.and.        &
      cell(i,j,k-1).ne.2) then

	ddx=0.5*(deltax(i)+deltax(i-1))
       ddxr=0.5*(deltax(i)+deltax(i+1))
       ddy=0.5*(deltay(j)+deltay(j-1))
       ddye=0.5*(deltay(j)+deltay(j+1))

       duwtdx=u_in_wm*(ca1_uw(i)*w(i+2,j,k)+ca2_uw(i)*w(i+1,j,k) &
       +ca3_uw(i)*w(i,j,k)+ca4_uw(i)*w(i-1,j,k)+ca5_uw(i)* &
       w(i-2,j,k))/(ca6_uw(i)*ddxr)+dabs(u_in_wm)* &
       (ck1_uw(i)*w(i+2,j,k)+ck2_uw(i)*w(i+1,j,k) &
       +ck3_uw(i)*w(i,j,k)+ck4_uw(i)*w(i-1,j,k)+ck5_uw(i)* &
       w(i-2,j,k))/(2.0*ck6_uw(i)*ddx)

       dvwtdy=v_in_wm*(ca1_vw(j)*w(i,j+2,k)+ca2_vw(j)*w(i,j+1,k) &
       +ca3_vw(j)*w(i,j,k)+ca4_vw(j)*w(i,j-1,k)+ca5_vw(j)* &
       w(i,j-2,k))/(ca6_vw(j)*ddye)+dabs(v_in_wm)* &
       (ck1_vw(j)*w(i,j+2,k)+ck2_vw(j)*w(i,j+1,k) &
       +ck3_vw(j)*w(i,j,k)+ck4_vw(j)*w(i,j-1,k)+ck5_vw(j)* &
       w(i,j-2,k))/(2.0*ck6_vw(j)*ddy)

       dwwtdz=w(i,j,k)*(ca1_ww(k)*w(i,j,k+2)+ca2_ww(k)*w(i,j,k+1) &
       +ca3_ww(k)*w(i,j,k)+ca4_ww(k)*w(i,j,k-1)+ca5_ww(k)* &
       w(i,j,k-2))/(ca6_ww(k)*deltaz(k+1))+dabs(w(i,j,k))* &
       (ck1_ww(k)*w(i,j,k+2)+ck2_ww(k)*w(i,j,k+1) &
       +ck3_ww(k)*w(i,j,k)+ck4_ww(k)*w(i,j,k-1)+ck5_ww(k)* &
       w(i,j,k-2))/(2.0*ck6_ww(k)*deltaz(k))
     
     	duwdx=duwtdx
       dvwdy=dvwtdy
       dwwdz=dwwtdz

	else
!cccccccccccccc---First Order Upwinding ----cccccccccccccccccccccccccccc
	r1xn=0.5*(deltax(i)+deltax(i-1))
	r1xd=0.5*(deltax(i)+deltax(i+1))
	r1x=r1xn/r1xd
       r1yn=0.5*(deltay(j)+deltay(j-1))
       r1yd=0.5*(deltay(j)+deltay(j+1))
       r1y=r1yn/r1yd
	r1z=deltaz(k)/deltaz(k+1)

	duwwdx=(u_in_wm/r1xd)*((-1.0/(r1x*(r1x+1.0)))* &
     	w(i-1,j,k)-(1.0-1.0/r1x)*w(i,j,k)+ &
     	(r1x/(r1x+1.0))*w(i+1,j,k))-al*dabs(u_in_wm)* &
     	(1.0/(2.0*r1xd))*(2.0/(r1x*(r1x+1.0))* &
     	w(i-1,j,k)-2.0/r1x*w(i,j,k)+2.0/(r1x+1.0)* &
     	w(i+1,j,k))

       dvwwdy=(v_in_wm/r1yd)*((-1.0/(r1y*(r1y+1.0)))* &
       w(i,j-1,k)-(1.0-1.0/r1y)*w(i,j,k)+ &
       (r1y/(r1y+1.0))*w(i,j+1,k))-al*dabs(v_in_wm)* &
       (1.0/(2.0*r1yd))*(2.0/(r1y*(r1y+1.0))* &
       w(i,j-1,k)-2.0/r1y*w(i,j,k)+2.0/(r1y+1.0)* &
       w(i,j+1,k))

	dwwwdz=(w(i,j,k)/deltaz(k+1))*((-1.0/(r1z*(r1z+1.0)))* &
     	w(i,j,k-1)-(1.0-1.0/r1z)*w(i,j,k)+ &
     	(r1z/(r1z+1.0))*w(i,j,k+1))-al*dabs(w(i,j,k))* &
     	(1.0/(2.0*deltaz(k+1)))*(2.0/(r1z*(r1z+1.0))* &
     	w(i,j,k-1)-2.0/r1z*w(i,j,k)+2.0/(r1z+1.0)* &
     	w(i,j,k+1))
	
      
       duwdx=duwwdx
       dvwdy=dvwwdy
       dwwdz=dwwwdz
       
       endif

!ccccccccccccccccccccccccc  grad of w part cccccccccccccccccccccccccccccc
	!d2wdx2=(2.0/dx)*((-w8/dx2xr)+(w6/dx2xl))
	!d2wdy2=(2.0/dy)*((-w12/dy2ye)+(w10/dy2yw)) 
	!d2wdz2=(2.0/dz2zt)*((-w16/dzt)+(w14/dz))
	
	d2wdx2=(2.0/(0.5*(dx2xr+dx2xl)))*((-W8/(0.5*dx2xr))+ &
       (W6/(0.5*dx2xl))) 
	
	d2wdy2=(2.0/(0.5*(dy2ye+dy2yw)))*((-W12/(0.5*dy2ye))+ &
       (W10/(0.5*dy2yw)))
	
	d2wdz2=(2.0/dz2zt)*((-w16/dzt)+(w14/dz))

       ztt2=rev*(d2wdx2+d2wdy2+d2wdz2)

       residw =(-duwdx-dvwdy-dwwdz+ztt2)

       !if(ita.eq.1 .or. irest.eq.1)then
       wt(i,j,k)=w(i,j,k)+deltat*(residw+dpdz)
       !else
       !wt(i,j,k)=w(i,j,k)+deltat*(0.5*(3.0*residw-resi_w(i,j,k))+dpdz)
       !endif

       resi_w(i,j,k)=residw
!***********************************************************************
	ENDDO
!$acc end parallel
       !write(6,*) 'leaving nseqcp '
!***********************************************************************
      END SUBROUTINE nsMomentum
!csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss



!csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE nsMomentumAB
!***********************************************************************
!     navier-stokes equations for constant properties                                              
!***********************************************************************
         USE global
!***********************************************************************
      !write(*,*)'has entered nseqcp'   
	  IMPLICIT NONE
	  INTEGER          :: i, j, k, n 
	  REAL (KIND = 8)  :: dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,            &
	                      u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
			        w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
			        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
			        dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          & 
			        r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
			        dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
			        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,            &
			        dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
			        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
			        duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw		  
			  
	al = 1.

!$acc parallel loop gang vector  &
!$acc private (i, j, k, dpdx,dpdy,dpdz,u1a,u22,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,   &
!$acc	        u15,u16,v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,           &
!$acc		 w1a,w22,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,dx2xr,             &
!$acc	        dx2xl,dy2ye,dy2yw,dxr,dx,dxl,dye,dy,dyw,dzt,dz,dzb,                         &
!$acc		 dz2zt,dz2zb,ddx,ddxr,ddy,ddye,ddz,ddzr,r1x,r1y,r1z,r1xn,r1xd,r1yn,          & 
!$acc		 r1yd,r1zn,r1zd,wu_n,wu_s,w_in_um,vu_e,vu_w,v_in_um,duutdx,dvutdy,           &
!$acc	        dwutdz,duuwdx,dvuwdy,dwuwdz,duudx,dvudy,dwudz,d2udx2,d2udy2,d2udz2,         &
!$acc	        uv_e,uv_w,u_in_vm,wv_n,wv_s,w_in_vm,duvtdx,dvvtdy,dwvtdz,duvwdx,            &
!$acc	        dvvwdy,dwvwdz,duvdx,dvvdy,dwvdz,d2vdx2,d2vdy2,d2vdz2,uw_e,uw_w,             &
!$acc	        u_in_wm,vw_n,vw_s,v_in_wm,duwtdx,dvwtdy,dwwtdz,duwwdx,dvwwdy,dwwwdz,        &
!$acc		 duwdx,dvwdy,dwwdz,d2wdx2,d2wdy2,d2wdz2,xtt2,residu,ytt2,residv,ztt2,residw) &
!$acc present (fluidIndexPtr, deltax, deltay, deltaz, u, v, w,             &
!$acc          ut, vt, wt, p, cell, x1, y1, z1, xp, yp, zp,                &
!$acc          xu, yu, zu, xv, yv, zv, xw, yw, zw, resi_u, resi_v, resi_w, &
!$acc          ca1_uu, ca2_uu, ca3_uu, ca4_uu, ca5_uu, ca6_uu, &
!$acc          ck1_uu, ck2_uu, ck3_uu, ck4_uu, ck5_uu, ck6_uu, &
!$acc          ca1_vu, ca2_vu, ca3_vu, ca4_vu, ca5_vu, ca6_vu, &
!$acc          ck1_vu, ck2_vu, ck3_vu, ck4_vu, ck5_vu, ck6_vu, &
!$acc          ca1_wu, ca2_wu, ca3_wu, ca4_wu, ca5_wu, ca6_wu, &
!$acc          ck1_wu, ck2_wu, ck3_wu, ck4_wu, ck5_wu, ck6_wu, &
!$acc          ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
!$acc          ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
!$acc          ca1_vv, ca2_vv, ca3_vv, ca4_vv, ca5_vv, ca6_vv, &
!$acc          ck1_vv, ck2_vv, ck3_vv, ck4_vv, ck5_vv, ck6_vv, &
!$acc          ca1_uw, ca2_uw, ca3_uw, ca4_uw, ca5_uw, ca6_uw, &
!$acc          ck1_uw, ck2_uw, ck3_uw, ck4_uw, ck5_uw, ck6_uw, &
!$acc          ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
!$acc          ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
!$acc          ca1_vw, ca2_vw, ca3_vw, ca4_vw, ca5_vw, ca6_vw, &
!$acc          ck1_vw, ck2_vw, ck3_vw, ck4_vw, ck5_vw, ck6_vw, &
!$acc          ca1_ww, ca2_ww, ca3_ww, ca4_ww, ca5_ww, ca6_ww, &
!$acc          ck1_ww, ck2_ww, ck3_ww, ck4_ww, ck5_ww, ck6_ww) &
!$acc firstprivate(nx, ny, nz, rev, deltat, al)         
	DO n = 1, fluidCellCount
	
       i = fluidIndexPtr(n, 1)
       j = fluidIndexPtr(n, 2)
       k = fluidIndexPtr(n, 3) 
    
	dxr=deltax(i+1)
	dx=deltax(i)
	dxl=deltax(i-1)
	dye=deltay(j+1)
	dy=deltay(j)
	dyw=deltay(j-1)
	dzt=deltaz(k+1)
	dz=deltaz(k)
	dzb=deltaz(k-1)
!cccccccccccccccccccccccccc  diff-u     ccccccccccccccccccccccccccccccccc
!     duu / dx
       u1a = u(i-1,j,k) + u(i,j,k)
       u22 = u(i-1,j,k) - u(i,j,k)
       u3  = u(i,j,k)   + u(i+1,j,k)
       u4  = u(i,j,k)   - u(i+1,j,k)

!     duv  / dy
       u5 = u(i,j-1,k)  + u(i,j,k)
       u6 = u(i,j-1,k)  - u(i,j,k)
       u7 = u(i,j,k)    + u(i,j+1,k)
       u8 = u(i,j,k)    - u(i,j+1,k)

!     dwu  /  dz
       u9 = u(i,j,k-1)  + u(i,j,k)
       u10= u(i,j,k-1)  - u(i,j,k)
       u11= u(i,j,k)    + u(i,j,k+1)
       u12= u(i,j,k)    - u(i,j,k+1)

!     duv/dx
       u13 = u(i-1,j,k) + u(i-1,j+1,k)
       u14 = u7

!     duw / dx
       u15 = u(i-1,j,k) + u(i-1,j,k+1)
       u16 = u11
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccc    diff -v    ccccccccccccccccccccccccccccccc
!     dvu / dx
       v1a = v(i,j-1,k)  + v(i+1,j-1,k)     
       v22 = v(i,j,k)    + v(i+1,j,k) 

!     duv / dx
       v3 = v(i-1,j,k)   + v(i,j,k)
       v4 = v(i-1,j,k)   - v(i,j,k)
       v5 = v(i,j,k)     + v(i+1,j,k)
       v6 = v(i,j,k)     - v(i+1,j,k)

!     dvv / dy
       v7 = v(i,j-1,k)   + v(i,j,k)
       v8 = v(i,j-1,k)   - v(i,j,k)
       v9 = v(i,j,k)     + v(i,j+1,k)
       v10= v(i,j,k)     - v(i,j+1,k)  

!     dwu / dz
       v11 = v(i,j,k-1)  + v(i,j,k)
       v12 = v(i,j,k-1)  - v(i,j,k) 
       v13 = v(i,j,k)    + v(i,j,k+1)
       v14 = v(i,j,k)    - v(i,j,k+1)    

!     dvw / dy
       v15 = v(i,j-1,k)  + v(i,j-1,k+1)
       v16 = v13
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccc  diff - w  cccccccccccccccccccccccccccccccc
!     dwu / dz
       w1a = w(i,j,k-1)  + w(i+1,j,k-1)
       w22 = w(i,j,k)    + w(i+1,j,k)

!     dwv / dz     
       w3 = w(i,j,k-1)   + w(i,j+1,k-1)
       w4 = w(i,j,k)     + w(i,j+1,k)

!     duw / dx   
       w5 = w(i-1,j,k)   + w(i,j,k)
       w6 = w(i-1,j,k)   - w(i,j,k)
       w7 = w(i,j,k)     + w(i+1,j,k)
       w8 = w(i,j,k)     - w(i+1,j,k)

!     dvw / dy
       w9 = w(i,j-1,k)   + w(i,j,k)
       w10 = w(i,j-1,k)  - w(i,j,k)
       w11 = w(i,j,k)    + w(i,j+1,k)
       w12 = w(i,j,k)    - w(i,j+1,k)

!     dww / dz 
       w13 = w(i,j,k-1)  + w(i,j,k)
       w14 = w(i,j,k-1)  - w(i,j,k)
       w15 = w(i,j,k)    + w(i,j,k+1)
       w16 = w(i,j,k)    - w(i,j,k+1)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       dpdx = (p(i,j,k) - p(i+1,j,k))/(0.5*(dxr+dx))
       dpdy = (p(i,j,k) - p(i,j+1,k))/(0.5*(dye+dy))
       dpdz = (p(i,j,k) - p(i,j,k+1))/(0.5*(dzt+dz))
       
    	dx2xr = dx + dxr
       dx2xl = dx + dxl
	dy2ye = dy + dye
	dy2yw = dy + dyw
	dz2zt = dz + dzt
	dz2zb = dz + dzb
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*********************** U - Momentum **********************************
	vu_e=v(i+1,j,k)+(dxr/(dx2xr))*(v(i,j,k)-v(i+1,j,k))
	vu_w=v(i+1,j-1,k)+(dxr/(dx2xr))*(v(i,j-1,k)-v(i+1,j-1,k))
	v_in_um=0.5*(vu_e+vu_w)
		 
	wu_n=w(i+1,j,k)+(dxr/(dx2xr))*(w(i,j,k)-w(i+1,j,k))
	wu_s=w(i+1,j,k-1)+(dxr/(dx2xr))*(w(i,j,k-1)-w(i+1,j,k-1))
	w_in_um=0.5*(wu_n+wu_s)
	
!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
       if(i.ne.2.and.i.lt.nx.and.j.ne.2.and.j.lt.ny+1.and.    &
       k.ne.2.and.k.lt.nz+1.and.cell(i+1,j,k).ne.2.and.       &
       cell(i-1,j,k).ne.2.and.cell(i,j+1,k).ne.2.and.         &
       cell(i,j-1,k).ne.2.and.cell(i,j,k+1).ne.2.and.         &
       cell(i,j,k-1).ne.2) then
 
	ddy=0.5*(deltay(j)+deltay(j-1))
       ddye=0.5*(deltay(j)+deltay(j+1))
       ddz=0.5*(deltaz(k)+deltaz(k-1))
       ddzr=0.5*(deltaz(k)+deltaz(k+1))

       duutdx=u(i,j,k)*(ca1_uu(i)*u(i+2,j,k)+ca2_uu(i)*u(i+1,j,k) &
       +ca3_uu(i)*u(i,j,k)+ca4_uu(i)*u(i-1,j,k)+ca5_uu(i)* &
       u(i-2,j,k))/(ca6_uu(i)*deltax(i+1))+dabs(u(i,j,k))* &
       (ck1_uu(i)*u(i+2,j,k)+ck2_uu(i)*u(i+1,j,k) &
       +ck3_uu(i)*u(i,j,k)+ck4_uu(i)*u(i-1,j,k)+ck5_uu(i)* &
       u(i-2,j,k))/(2.0*ck6_uu(i)*deltax(i))

       dvutdy=v_in_um*(ca1_vu(j)*u(i,j+2,k)+ca2_vu(j)*u(i,j+1,k) &
       +ca3_vu(j)*u(i,j,k)+ca4_vu(j)*u(i,j-1,k)+ca5_vu(j)* &
       u(i,j-2,k))/(ca6_vu(j)*ddye)+dabs(v_in_um)* &
       (ck1_vu(j)*u(i,j+2,k)+ck2_vu(j)*u(i,j+1,k) &
       +ck3_vu(j)*u(i,j,k)+ck4_vu(j)*u(i,j-1,k)+ck5_vu(j)* &
       u(i,j-2,k))/(2.0*ck6_vu(j)*ddy) 

       dwutdz=w_in_um*(ca1_wu(k)*u(i,j,k+2)+ca2_wu(k)*u(i,j,k+1) &
       +ca3_wu(k)*u(i,j,k)+ca4_wu(k)*u(i,j,k-1)+ca5_wu(k)* &
       u(i,j,k-2))/(ca6_wu(k)*ddzr)+dabs(w_in_um)* &
       (ck1_wu(k)*u(i,j,k+2)+ck2_wu(k)*u(i,j,k+1) &
       +ck3_wu(k)*u(i,j,k)+ck4_wu(k)*u(i,j,k-1)+ck5_wu(k)* &
       u(i,j,k-2))/(2.0*ck6_wu(k)*ddz)
     
     	duudx=duutdx
       dvudy=dvutdy
       dwudz=dwutdz
       
	else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
	r1x=deltax(i)/deltax(i+1)
       r1yn=0.5*(deltay(j)+deltay(j-1))
       r1yd=0.5*(deltay(j)+deltay(j+1))
       r1y=r1yn/r1yd
	r1zn=0.5*(deltaz(k)+deltaz(k-1))
	r1zd=0.5*(deltaz(k)+deltaz(k+1))
	r1z=r1zn/r1zd

	duuwdx=(u(i,j,k)/deltax(i+1))*((-1.0/(r1x*(r1x+1.0)))* &
     	u(i-1,j,k)-(1.0-1.0/r1x)*u(i,j,k)+   &
     	(r1x/(r1x+1.0))*u(i+1,j,k))-al*dabs(u(i,j,k))* &
     	(1.0/(2.0*deltax(i+1)))*(2.0/(r1x*(r1x+1.0))* &
     	u(i-1,j,k)-2.0/r1x*u(i,j,k)+2.0/(r1x+1.0)* &
       u(i+1,j,k))

       dvuwdy=(v_in_um/r1yd)*((-1.0/(r1y*(r1y+1.0)))* &
       u(i,j-1,k)-(1.0-1.0/r1y)*u(i,j,k)+ &
       (r1y/(r1y+1.0))*u(i,j+1,k))-al*dabs(v_in_um)* &
       (1.0/(2.0*r1yd))*(2.0/(r1y*(r1y+1.0))* &
       u(i,j-1,k)-2.0/r1y*u(i,j,k)+2.0/(r1y+1.0)* &
       u(i,j+1,k))

	dwuwdz=(w_in_um/r1zd)*((-1.0/(r1z*(r1z+1.0)))* &
     	u(i,j,k-1)-(1.0-1.0/r1z)*u(i,j,k)+ &
    	(r1z/(r1z+1.0))*u(i,j,k+1))-al*dabs(w_in_um)* &
     	(1.0/(2.0*r1zd))*(2.0/(r1z*(r1z+1.0))* &
     	u(i,j,k-1)-2.0/r1z*u(i,j,k)+2.0/(r1z+1.0)* &
      	u(i,j,k+1))
 
       duudx=duuwdx
       dvudy=dvuwdy
       dwudz=dwuwdz
        
       endif
   
!ccccccccccccccccccccccccc  grad of u part cccccccccccccccccccccccccccccc
	!d2udx2=(2.0/dx2xr)*((-u4/dxr)+(u22/dx))
       !d2udy2=(2.0/dy)*((-u8/dy2ye)+(u6/dy2yw))
	!d2udz2=(2.0/dz)*((-u12/dz2zt)+(u10/dz2zb))
	
	d2udx2=(2.0/dx2xr)*((-u4/dxr)+(u22/dx))
	
       d2udy2=(2.0/(0.5*(dy2ye+dy2yw)))*((-u8/(0.5*dy2ye))+ &
       (u6/(0.5*dy2yw)))
     
	d2udz2=(2.0/(0.5*(dz2zt+dz2zb)))*((-u12/(0.5*dz2zt))+ &
       (u10/(0.5*dz2zb)))

       xtt2=rev*(d2udx2+d2udy2+d2udz2)
       
       residu=(-duudx-dvudy-dwudz+xtt2)

       !if(ita.eq.1)then
       !ut(i,j,k)=u(i,j,k)+deltat*(residu+dpdx)
       !else
       ut(i,j,k)=u(i,j,k)+deltat*(0.5*(3.0*residu-resi_u(i,j,k))+dpdx)
       !endif

       resi_u(i,j,k)=residu
!c***********************************************************************

!c*********************** V - Momentum *********************************
	uv_e=u(i,j+1,k)+(dye/(dy2ye))*(u(i,j,k)-u(i,j+1,k))
	uv_w=u(i-1,j+1,k)+(dye/(dy2ye))*(u(i-1,j,k)-u(i-1,j+1,k))
	u_in_vm=0.5*(uv_e+uv_w)
	
	wv_n=w(i,j+1,k)+(dye/(dy2ye))*(w(i,j,k)-w(i,j+1,k))
	wv_s=w(i,j+1,k-1)+(dye/(dy2ye))*(w(i,j,k-1)-w(i,j+1,k-1))
	w_in_vm=0.5*(wv_n+wv_s)

!cccccccccccccc---Third Order Upwinding ----cccccccccccccccccccccccccccc
      if(i.ne.2.and.i.lt.nx+1.and.j.ne.2.and.j.lt.ny.and. &
      k.ne.2.and.k.lt.nz+1.and.cell(i+1,j,k).ne.2.and.    &
      cell(i-1,j,k).ne.2.and.cell(i,j+1,k).ne.2.and.      &
      cell(i,j-1,k).ne.2.and.cell(i,j,k+1).ne.2.and.      &
      cell(i,j,k-1).ne.2) then
      
	ddx=0.5*(deltax(i)+deltax(i-1))
       ddxr=0.5*(deltax(i)+deltax(i+1))
       ddz=0.5*(deltaz(k)+deltaz(k-1))
       ddzr=0.5*(deltaz(k)+deltaz(k+1))

       duvtdx=u_in_vm*(ca1_uw(i)*v(i+2,j,k)+ca2_uw(i)*v(i+1,j,k) &
       +ca3_uw(i)*v(i,j,k)+ca4_uw(i)*v(i-1,j,k)+ca5_uw(i)* &
       v(i-2,j,k))/(ca6_uw(i)*ddxr)+dabs(u_in_vm)* &
       (ck1_uw(i)*v(i+2,j,k)+ck2_uw(i)*v(i+1,j,k) &
       +ck3_uw(i)*v(i,j,k)+ck4_uw(i)*v(i-1,j,k)+ck5_uw(i)* &
       v(i-2,j,k))/(2.0*ck6_uw(i)*ddx)

       dvvtdy=v(i,j,k)*(ca1_vv(j)*v(i,j+2,k)+ca2_vv(j)*v(i,j+1,k) &
       +ca3_vv(j)*v(i,j,k)+ca4_vv(j)*v(i,j-1,k)+ca5_vv(j)* &
       v(i,j-2,k))/(ca6_vv(j)*deltay(j+1))+dabs(v(i,j,k))* &
       (ck1_vv(j)*v(i,j+2,k)+ck2_vv(j)*v(i,j+1,k) &
       +ck3_vv(j)*v(i,j,k)+ck4_vv(j)*v(i,j-1,k)+ck5_vv(j)* &
       v(i,j-2,k))/(2.0*ck6_vv(j)*deltay(j))

       dwvtdz=w_in_vm*(ca1_wu(k)*v(i,j,k+2)+ca2_wu(k)*v(i,j,k+1) &
       +ca3_wu(k)*v(i,j,k)+ca4_wu(k)*v(i,j,k-1)+ca5_wu(k)* &
       v(i,j,k-2))/(ca6_wu(k)*ddzr)+dabs(w_in_vm)* &
       (ck1_wu(k)*v(i,j,k+2)+ck2_wu(k)*v(i,j,k+1) &
       +ck3_wu(k)*v(i,j,k)+ck4_wu(k)*v(i,j,k-1)+ck5_wu(k)* &
       v(i,j,k-2))/(2.0*ck6_wu(k)*ddz)
     
       duvdx=duvtdx
       dvvdy=dvvtdy
       dwvdz=dwvtdz
       
       else
!cccccccccccccc---First Order Upwinding ----ccccccccccccccccccccccccccccc
	r1xn=0.5*(deltax(i)+deltax(i-1))
	r1xd=0.5*(deltax(i)+deltax(i+1))
	r1x=r1xn/r1xd
       r1y=deltay(j)/deltay(j+1)
	r1zn=0.5*(deltaz(k)+deltaz(k-1))
	r1zd=0.5*(deltaz(k)+deltaz(k+1))
	r1z=r1zn/r1zd

	duvwdx=(u_in_vm/r1xd)*((-1.0/(r1x*(r1x+1.0)))* &
     	v(i-1,j,k)-(1.0-1.0/r1x)*v(i,j,k)+ &
     	(r1x/(r1x+1.0))*v(i+1,j,k))-al*dabs(u_in_vm)* &
     	(1.0/(2.0*r1xd))*(2.0/(r1x*(r1x+1.0))* &
     	v(i-1,j,k)-2.0/r1x*v(i,j,k)+2.0/(r1x+1.0)* &
     	v(i+1,j,k))

       dvvwdy=(v(i,j,k)/deltay(j+1))*((-1.0/(r1y*(r1y+1.0)))* &
       v(i,j-1,k)-(1.0-1.0/r1y)*v(i,j,k)+ &
       (r1y/(r1y+1.0))*v(i,j+1,k))-al*dabs(v(i,j,k))* &
       (1.0/(2.0*deltay(j+1)))*(2.0/(r1y*(r1y+1.0))* &
       v(i,j-1,k)-2.0/r1y*v(i,j,k)+2.0/(r1y+1.0)* &
       v(i,j+1,k))

	dwvwdz=(w_in_vm/r1zd)*((-1.0/(r1z*(r1z+1.0)))* &
     	v(i,j,k-1)-(1.0-1.0/r1z)*v(i,j,k)+ &
     	(r1z/(r1z+1.0))*v(i,j,k+1))-al*dabs(w_in_vm)* &
     	(1.0/(2.0*r1zd))*(2.0/(r1z*(r1z+1.0))* &
     	v(i,j,k-1)-2.0/r1z*v(i,j,k)+2.0/(r1z+1.0)* &
     	v(i,j,k+1))	 

       duvdx=duvwdx
       dvvdy=dvvwdy
       dwvdz=dwvwdz
       
       endif

!ccccccccccccccccccccccccc  grad of v part cccccccccccccccccccccccccccccc
	!d2vdx2=(2.0/dx)*((-v6/dx2xr)+(v4/dx2xl))
       !d2vdy2=(2.0/dy2ye)*((-v10/dye)+(v8/dy))
	!d2vdz2=(2.0/dz)*((-v14/dz2zt)+(v12/dz2zb))
	
	d2vdx2=(2.0/(0.5*(dx2xr+dx2xl)))*((-v6/(0.5*dx2xr))+ &
	(v4/(0.5*dx2xl)))  
     
       d2vdy2=(2.0/dy2ye)*((-v10/dye)+(v8/dy))
       
	d2vdz2=(2.0/(0.5*(dz2zt+dz2zb)))*((-v14/(0.5*dz2zt))+ &
	(v12/(0.5*dz2zb)))
		
       ytt2=rev*(d2vdx2+d2vdy2+d2vdz2)

       residv=(-duvdx-dvvdy-dwvdz+ytt2)

       !if(ita.eq.1 .or. irest.eq.1)then
       !vt(i,j,k)=v(i,j,k)+deltat*(residv+dpdy)
       !else
       vt(i,j,k)=v(i,j,k)+deltat*(0.5*(3.0*residv-resi_v(i,j,k))+dpdy)
       !endif

       resi_v(i,j,k)=residv
!c***********************************************************************

!c*********************** W - Momentum **********************************
	uw_e=u(i,j,k+1)+(dzt/(dz2zt))*(u(i,j,k)-u(i,j,k+1))
	uw_w=u(i-1,j,k+1)+(dzt/(dz2zt))*(u(i-1,j,k)-u(i-1,j,k+1))
	u_in_wm=0.5*(uw_e+uw_w)
	
	vw_n=v(i,j,k+1)+(dzt/(dz2zt))*(v(i,j,k)-v(i,j,k+1))
	vw_s=v(i,j-1,k+1)+(dzt/(dz2zt))*(v(i,j-1,k)-v(i,j-1,k+1))
	v_in_wm=0.5*(vw_n+vw_s)	

!cccccccccccccc---Third Order Upwinding ----ccccccccccccccccccccccccccccc
      if(i.ne.2.and.i.lt.nx+1.and.j.ne.2.and.j.lt.ny+1.and. &
      k.ne.2.and.k.lt.nz.and.cell(i+1,j,k).ne.2.and.        &
      cell(i-1,j,k).ne.2.and.cell(i,j+1,k).ne.2.and.        &
      cell(i,j-1,k).ne.2.and.cell(i,j,k+1).ne.2.and.        &
      cell(i,j,k-1).ne.2) then

	ddx=0.5*(deltax(i)+deltax(i-1))
       ddxr=0.5*(deltax(i)+deltax(i+1))
       ddy=0.5*(deltay(j)+deltay(j-1))
       ddye=0.5*(deltay(j)+deltay(j+1))

       duwtdx=u_in_wm*(ca1_uw(i)*w(i+2,j,k)+ca2_uw(i)*w(i+1,j,k) &
       +ca3_uw(i)*w(i,j,k)+ca4_uw(i)*w(i-1,j,k)+ca5_uw(i)* &
       w(i-2,j,k))/(ca6_uw(i)*ddxr)+dabs(u_in_wm)* &
       (ck1_uw(i)*w(i+2,j,k)+ck2_uw(i)*w(i+1,j,k) &
       +ck3_uw(i)*w(i,j,k)+ck4_uw(i)*w(i-1,j,k)+ck5_uw(i)* &
       w(i-2,j,k))/(2.0*ck6_uw(i)*ddx)

       dvwtdy=v_in_wm*(ca1_vw(j)*w(i,j+2,k)+ca2_vw(j)*w(i,j+1,k) &
       +ca3_vw(j)*w(i,j,k)+ca4_vw(j)*w(i,j-1,k)+ca5_vw(j)* &
       w(i,j-2,k))/(ca6_vw(j)*ddye)+dabs(v_in_wm)* &
       (ck1_vw(j)*w(i,j+2,k)+ck2_vw(j)*w(i,j+1,k) &
       +ck3_vw(j)*w(i,j,k)+ck4_vw(j)*w(i,j-1,k)+ck5_vw(j)* &
       w(i,j-2,k))/(2.0*ck6_vw(j)*ddy)

       dwwtdz=w(i,j,k)*(ca1_ww(k)*w(i,j,k+2)+ca2_ww(k)*w(i,j,k+1) &
       +ca3_ww(k)*w(i,j,k)+ca4_ww(k)*w(i,j,k-1)+ca5_ww(k)* &
       w(i,j,k-2))/(ca6_ww(k)*deltaz(k+1))+dabs(w(i,j,k))* &
       (ck1_ww(k)*w(i,j,k+2)+ck2_ww(k)*w(i,j,k+1) &
       +ck3_ww(k)*w(i,j,k)+ck4_ww(k)*w(i,j,k-1)+ck5_ww(k)* &
       w(i,j,k-2))/(2.0*ck6_ww(k)*deltaz(k))
     
     	duwdx=duwtdx
       dvwdy=dvwtdy
       dwwdz=dwwtdz

	else
!cccccccccccccc---First Order Upwinding ----cccccccccccccccccccccccccccc
	r1xn=0.5*(deltax(i)+deltax(i-1))
	r1xd=0.5*(deltax(i)+deltax(i+1))
	r1x=r1xn/r1xd
       r1yn=0.5*(deltay(j)+deltay(j-1))
       r1yd=0.5*(deltay(j)+deltay(j+1))
       r1y=r1yn/r1yd
	r1z=deltaz(k)/deltaz(k+1)

	duwwdx=(u_in_wm/r1xd)*((-1.0/(r1x*(r1x+1.0)))* &
     	w(i-1,j,k)-(1.0-1.0/r1x)*w(i,j,k)+ &
     	(r1x/(r1x+1.0))*w(i+1,j,k))-al*dabs(u_in_wm)* &
     	(1.0/(2.0*r1xd))*(2.0/(r1x*(r1x+1.0))* &
     	w(i-1,j,k)-2.0/r1x*w(i,j,k)+2.0/(r1x+1.0)* &
     	w(i+1,j,k))

       dvwwdy=(v_in_wm/r1yd)*((-1.0/(r1y*(r1y+1.0)))* &
       w(i,j-1,k)-(1.0-1.0/r1y)*w(i,j,k)+ &
       (r1y/(r1y+1.0))*w(i,j+1,k))-al*dabs(v_in_wm)* &
       (1.0/(2.0*r1yd))*(2.0/(r1y*(r1y+1.0))* &
       w(i,j-1,k)-2.0/r1y*w(i,j,k)+2.0/(r1y+1.0)* &
       w(i,j+1,k))

	dwwwdz=(w(i,j,k)/deltaz(k+1))*((-1.0/(r1z*(r1z+1.0)))* &
     	w(i,j,k-1)-(1.0-1.0/r1z)*w(i,j,k)+ &
     	(r1z/(r1z+1.0))*w(i,j,k+1))-al*dabs(w(i,j,k))* &
     	(1.0/(2.0*deltaz(k+1)))*(2.0/(r1z*(r1z+1.0))* &
     	w(i,j,k-1)-2.0/r1z*w(i,j,k)+2.0/(r1z+1.0)* &
     	w(i,j,k+1))
	
      
       duwdx=duwwdx
       dvwdy=dvwwdy
       dwwdz=dwwwdz
       
       endif

!ccccccccccccccccccccccccc  grad of w part cccccccccccccccccccccccccccccc
	!d2wdx2=(2.0/dx)*((-w8/dx2xr)+(w6/dx2xl))
	!d2wdy2=(2.0/dy)*((-w12/dy2ye)+(w10/dy2yw)) 
	!d2wdz2=(2.0/dz2zt)*((-w16/dzt)+(w14/dz))
	
	d2wdx2=(2.0/(0.5*(dx2xr+dx2xl)))*((-W8/(0.5*dx2xr))+ &
       (W6/(0.5*dx2xl))) 
	
	d2wdy2=(2.0/(0.5*(dy2ye+dy2yw)))*((-W12/(0.5*dy2ye))+ &
       (W10/(0.5*dy2yw)))
	
	d2wdz2=(2.0/dz2zt)*((-w16/dzt)+(w14/dz))

       ztt2=rev*(d2wdx2+d2wdy2+d2wdz2)

       residw =(-duwdx-dvwdy-dwwdz+ztt2)

       !if(ita.eq.1 .or. irest.eq.1)then
       !wt(i,j,k)=w(i,j,k)+deltat*(residw+dpdz)
       !else
       wt(i,j,k)=w(i,j,k)+deltat*(0.5*(3.0*residw-resi_w(i,j,k))+dpdz)
       !endif

       resi_w(i,j,k)=residw
!***********************************************************************
	ENDDO
!$acc end parallel
       !write(6,*) 'leaving nseqcp '
!***********************************************************************
      END SUBROUTINE nsMomentumAB









