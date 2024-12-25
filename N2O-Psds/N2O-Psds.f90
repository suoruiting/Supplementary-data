Program main    
  implicit none
  double precision                   :: T                  
  double precision                   :: Ps, denLs, dengs    
  double precision                   :: dl, dg              
  integer                            :: i
  
  write(*,*)"**********************************************************************"
  write(*,*) "T-P range : 182.33 - 309.52 K "
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)
   
100  write(*,*) 
  write(*,*) "Please input T(K):"
  read(*,*) T
  
 if( T < 182.33 .or. T >309.52 ) then
	write(*,*) "T is wrong, try again!"
    goto 100 
 end if
 
  write(*,*)
  write(*,*)
 
  call N2OEOS_VLE( T, Ps, denLs, dengs )
  
  dl = denLs/44.0128
  dg = dengs/44.0128
  
  write(*,*)'the results are as follows:'
  write(*,*)"**********************************************************************"
  write(*,'(a8,a12,2a15)')'  T(K)','   Ps(bar)','  dl(mol/dm3)','  dg(mol/dm3)'
  write(*,*)"______________________________________________________________________"
  write(*,'(f9.2,f9.3,2f14.3)') T, Ps, dl, dg
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)
  

  open(40,file='Psds.out.txt')
  write(40,'(a25,a25,2a25)')'  T(K)','   Ps(bar)','  dl(mol/dm3)','  dg(mol/dm3)'
  write(40,'(f25.2,f22.3,f25.3,f30.3)') T, Ps, dl, dg
  
  goto 100    
  
  
  
  stop
end program main


subroutine N2OEOS_VLE( T, Ps_cal, dl_cal, dg_cal)
implicit none
  double precision  :: T, Ps_cal, dl_cal, dg_cal
  double precision  :: dN2O1, dN2O2                
  double precision  :: dc, Tc, R
  double precision  :: tau, delta1, phird1, J1, K1, Jd1, Kd1
  double precision  :: delta2, phird2, J2, K2, Jd2, Kd2
  double precision  :: angle
  integer           :: iter


  dc = 452.01               
  Tc = 309.52                
  R = 8.314472d0 / 44.0128   

  
	call N2Osaturation( T, Ps_cal, dl_cal,dg_cal  ) 
	dN2O1 = dl_cal 
	dN2O2 = dg_cal
  
  
  delta1 = dN2O1 / dc
  delta2 = dN2O2 / dc
  tau = Tc / T 

  call N2O( tau, delta1, phird1, J1, K1, Jd1, Kd1 )
  call N2O( tau, delta2, phird2, J2, K2, Jd2, Kd2 )

  iter = 0
  do while( ( dabs(K2-K1)+dabs(J2-J1) ) >= 1.0d-8 )
    angle = Jd2*Kd1 - Jd1*Kd2
    delta1 = delta1 + 1.0/angle * ( (K2-K1)*Jd2 - (J2-J1)*Kd2 )
	delta2 = delta2 + 1.0/angle * ( (K2-K1)*Jd1 - (J2-J1)*Kd1 )
    call N2O( tau, delta1, phird1, J1, K1, Jd1, Kd1 )
    call N2O( tau, delta2, phird2, J2, K2, Jd2, Kd2 )
    iter = iter + 1
  end do

  Ps_cal = delta2*dc * R * T * ( 1.0 + delta2 * phird2 ) / 100 
  dl_cal = delta1 * dc 
  dg_cal = delta2 * dc 

  return
end subroutine N2OEOS_VLE



subroutine N2Osaturation( Ts, Ps_cal, dl_cal,dg_cal )  
  implicit none
  double precision :: Ts,t, Ps_cal, dl_cal, dg_cal    
  double precision :: a(6), b(6), c(5) 
  data a / -6.86756968137918,  1.99959566385607, -4.32957440667956, 9.14498700879084,  -17.0562465175674, 9.27917419595028 /
  data b /  1.76129126313810 , -0.901947823843190, 0.948599809065473,  -1.39709870576305, 5.32895289468383,  -133.113915278205 /
  data c / -1.85504542355836 , -3.95309455117505,  1.37240598427993, -1.40918801740346, -3.38977734035671 /
	 
  double precision, parameter :: Tc = 309.52         
  double precision, parameter :: Pc = 72.450         
  double precision, parameter :: denc = 452.01       

      t = 1.0 - Ts/Tc
	  
	  Ps_cal = Pc * dexp( Tc/Ts * ( a(1)*t + a(2)*t**1.5 + a(3)*t**2.5 + a(4)*t**3.5 + a(5)*t**4.5 + a(6)*t**5.5) )

	  dl_cal= denc * dexp( b(1)*t**(1.0/3.0) + b(2)*t**(2.0/3.0) + b(3)*t**(4.0/3.0)              &    
	                 + b(4)*t**(8.0/3.0) + b(5)*t**(16.0/3.0) + b(6)*t**(32.0/3.0) )
     
	  dg_cal = denc * dexp ( Tc/Ts * ( c(1)*t**(1.0/3.0) + c(2)*t**(5.0/6.0) + c(3)*t**(7.0/6.0) + &
                     c(4)*t**(13.0/6.0) + c(5)*t**(14.0/3.0) ) )


  return
end subroutine N2Osaturation


subroutine N2O( tau, delta, phird, J, K, Jd, Kd )
  implicit none
  double precision  :: tau, delta, phird, J, K, Jd, Kd
  double precision  :: phir, phirdd
  integer           :: i
  double precision  :: ii(14), jj(14), para(14), kk(14)
	
  data ii / 1.0, 1.0, 1.0, 3.0, 7.0, 2.0, 1.0, 1.0, 2.0, 5.0, 1.0, 1.0, 4.0, 2.0 /

  data jj / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /

  data kk / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0 /
  

 data para / 0.40209218D+00, 0.87844892D+00, -0.24241584D+01, 0.69011106D-01, 0.20412510D-03, -0.52403092D-02, -0.43148215D-02, &
             0.11046048D+00, 0.44400973D+00, -0.35551109D-02, -0.22637699D+00, -0.43690514D-03, -0.41132000D-01, -0.22998680D-01 /

!----------------------------------to calculate phir-----------------------------------------
  phir = 0.0
  
  do i =1, 6
    phir = phir + para(i) * delta**ii(i) * tau**jj(i)
  end do

  do i = 7, 14
    phir = phir + para(i)  * delta**ii(i) * tau**jj(i) * dexp( -delta**kk(i) ) 
  end do

!-----------------------------------------------------------------------------------------------

!----------------------------------to calculate phird-----------------------------------------
  phird = 0.0
  
  do i =1, 6
    phird = phird + para(i) * ii(i) * delta**(ii(i)-1) * tau**jj(i)
  end do

  do i = 7, 14
    phird = phird + para(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-1)*tau**jj(i)* ( ii(i)-kk(i)*delta**kk(i) )
  end do

!-----------------------------------------------------------------------------------------------

!------------------------------------to calculate phirdd----------------------------------------
  phirdd = 0.0
  do i = 1, 6
    phirdd = phirdd + para(i) * ii(i) * ( ii(i)-1 ) * delta**(ii(i)-2) * tau**jj(i)
  end do

  do i = 7, 14
    phirdd = phirdd + para(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-2) * tau**jj(i) *  &
	         ( ( ii(i)-kk(i)*delta**kk(i) )*( ii(i)-1-kk(i)*delta**kk(i) ) - kk(i)**2 * delta**kk(i) )
  end do
!------------------------------------------------------------------------------------------------------

  J = delta * ( 1.0 + delta*phird )
  K = delta*phird + phir + dlog(delta)
  Jd = 1.0 + 2.0*delta*phird + delta**2*phirdd
  Kd = 2.0*phird + delta*phirdd + 1.0/delta

  return
end subroutine N2O
