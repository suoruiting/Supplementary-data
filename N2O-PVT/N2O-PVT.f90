Program main
implicit none

  double precision  :: dN2O,d_N2O     
  double precision  :: T, P   
  
  
  write(*,*)"**********************************************************************"
  write(*,*) "T-P range : 183.33-525 K and pressure to 500 bar"
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)
   100 Write(*,*) "Please input T(K)and P(bar): "
  read(*,*) T, P

   if( T < 182.33 .or. T >525 .or. P < 0 .or. P > 500 ) then
	write(*,*) "Input error, please try again!"
    goto 100 
 end if
   write(*,*)
   write(*,*)
  
  call N2OEOS( T, P, dN2O )
  
  d_N2O = dN2O/44.0128
  

  write(*,*)'the results are as follows:'
  write(*,*)"**********************************************************************"
  write(*,'(a9,a12,a20)') "T(K)", "P(bar)", "d_N2O(mol/dm3)"
  write(*,*)"______________________________________________________________________"
  write(*,'(f10.2,f10.2,f16.4)') T, P, d_N2O
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)
  
  open(50,file='PVT.out.txt')
  write(50,'(a20,a23,a28)') "T(K)", "P(bar)", "dN2O(mol/dm3)"
  write(50,'(f20.2,f20.2,f28.4)') T, P, d_N2O
  
  goto 100
   stop
  end program
  


                       
subroutine N2Osaturation( Ts, Ps, denls, dengs )  
implicit none
  double precision :: Ts, Ps, denls, dengs     
  double precision :: a(6)
  double precision :: b(6)
  double precision :: c(5)
  double precision :: t
	 
  double precision, parameter :: Tc = 309.52         
  double precision, parameter :: Pc = 72.450         
  double precision, parameter :: denc = 452.01       
  
  data a / -6.86756968137918,  1.99959566385607, -4.32957440667956, 9.14498700879084,  -17.0562465175674, 9.27917419595028 /
  data b /  1.76129126313810, -0.901947823843190, 0.948599809065473, -1.39709870576305, 5.32895289468383, -133.113915278205 /
  data c / -1.85504542355836, -3.95309455117505, 1.37240598427993, -1.40918801740346, -3.38977734035671 /

	  t = 1 - Ts/Tc

	  ps = Pc * dexp( Tc/Ts * ( a(1)*t + a(2)*t**1.5 + a(3)*t**2.5 + a(4)*t**3.5 + a(5)*t**4.5 + a(6)*t**5.5) )

	  denls = denc * dexp( b(1)*t**(1.0/3.0) + b(2)*t**(2.0/3.0) + b(3)*t**(4.0/3.0)              &    
	                 + b(4)*t**(8.0/3.0) + b(5)*t**(16.0/3.0) + b(6)*t**(32.0/3.0) )
     
	  dengs = denc * dexp ( Tc/Ts * ( c(1)*t**(1.0/3.0) + c(2)*t**(5.0/6.0) + c(3)*t**(7.0/6.0) + &
	                       c(4)*t**(13.0/6.0) + c(5)*t**(14.0/3.0) ) )
  return
end subroutine N2Osaturation



Subroutine DNewt( T, P, X, X_OUT )
    Double precision, parameter  :: EPS=1.0d-8
	double precision             :: T, P
	Double precision             :: X, X1, F, DF,X_OUT

    L = 60
	Call F_N2O( T, P, X, F, DF )

10  IF(dABS(DF)+1.0.EQ.1.0)THEN
      L = 0
	  WRITE(*,20)
	  RETURN
	ENDIF

20  FORMAT(1X,'ERR')
    X1 = X-F/DF
	CALL F_N2O( T, P, X1, F, DF )
	
	IF((dABS(X1-X).GE.EPS).OR.(dABS(F).GE.EPS))THEN
	  L = L - 1
	  X = X1
	  IF(L.EQ.0)RETURN
	  GOTO 10
	ENDIF
	n_time=60-l
	X_Out = X1
	RETURN

End Subroutine DNewt


Subroutine F_N2O( T, P, delta, F, dF )
  implicit none
  double precision  :: T, P, delta, F, DF, phird, phirdd
  double precision  :: R, dc, tau
  
  R = 8.314472/44.0128       
  dc = 452.01    
  tau = 309.52/T

   call phir_N2O( tau, delta, phird, phirdd )


    F = P - ( dc * R * T * delta + dc * R * T * delta**2 * phird ) / 100

    dF = ( -dc * R * T - 2.0 * dc * R * T * delta * phird - dc * R * T * delta**2 * phirdd ) / 100

    Return

End Subroutine F_N2O


Subroutine N2OEOS( T, P, dN2O ) 
implicit none
  double precision  :: T, P, R, dc, Tc, Ps, denls, dengs, delta, deltacal
  double precision  :: dN2O, phic
  double precision  :: phir, phird, phirdd, tau
  integer           :: i
  double precision  :: im(14),jm(14),am(14),km(14)

 data am / 0.40209218D+00, 0.87844892D+00, -0.24241584D+01, 0.69011106D-01, 0.20412510D-03, -0.52403092D-02, -0.43148215D-02, &
           0.11046048D+00, 0.44400973D+00, -0.35551109D-02, -0.22637699D+00, -0.43690514D-03, -0.41132000D-01, -0.22998680D-01 /

   data im / 1.0, 1.0, 1.0, 3.0, 7.0, 2.0, 1.0, 1.0, 2.0, 5.0, 1.0, 1.0, 4.0, 2.0 /

  data jm / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /

  data km / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0 /
  
  dc = 452.01       
  Tc = 309.52       
  R = 8.314472      
  
  if( T<Tc ) then
    call N2Osaturation( T, Ps, denls, dengs ) 
  else 
    goto 200
  end if

200 if( T < Tc .and. P > Ps ) then
      dN2O = denls                            
    else
	  dN2O = P / R / T * 100/44.0128            
    end if
  
  delta = dN2O / dc                           
  tau = Tc / T
 
  Call DNewt( T, P, delta, deltacal ) 

  delta = deltacal
  dN2O = dc * delta
 
  
  return
end subroutine N2OEOS

subroutine phir_N2O( tau, delta, phird, phirdd ) 
implicit none
  double precision  :: tau, delta, phird, phirdd
  integer           :: i
  double precision  :: im(14), jm(14), am(14), km(14)

   data am / 0.40209218D+00, 0.87844892D+00, -0.24241584D+01, 0.69011106D-01, 0.20412510D-03, -0.52403092D-02, -0.43148215D-02, &
             0.11046048D+00, 0.44400973D+00, -0.35551109D-02, -0.22637699D+00, -0.43690514D-03, -0.41132000D-01, -0.22998680D-01 /

   data im / 1.0, 1.0, 1.0, 3.0, 7.0, 2.0, 1.0, 1.0, 2.0, 5.0, 1.0, 1.0, 4.0, 2.0 /

   data jm / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /

   data km / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0 /
  

!----------------------------------to calculate phird-----------------------------------------
  phird = 0.0
  
  do i =1, 6
    phird = phird + am(i) * im(i) * delta**(im(i)-1) * tau**jm(i)   
  end do

  do i = 7, 14
    phird = phird + am(i) * dexp( -delta**km(i) ) * delta**(im(i)-1)*tau**jm(i)* ( im(i)-km(i)*delta**km(i) )
  end do

!-----------------------------------------------------------------------------------------------

!------------------------------------to calculate phirdd----------------------------------------
  phirdd = 0.0
  do i = 1, 6
    phirdd = phirdd + am(i) * im(i) * ( im(i)-1 ) * delta**(im(i)-2) * tau**jm(i)
  end do

  do i = 7, 14
    phirdd = phirdd + am(i) * dexp( -delta**km(i) ) * delta**(im(i)-2) * tau**jm(i) *  &
	         ( ( im(i)-km(i)*delta**km(i) )*( im(i)-1-km(i)*delta**km(i) ) - km(i)**2 * delta**km(i) )
  end do

!------------------------------------------------------------------------------------------------------

  return
end subroutine phir_N2O

