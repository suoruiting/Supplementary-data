Module expdata
  double precision:: T, P, d_exp
  double precision:: x1, x2 
End Module expdata
    
Module Parameters         
  double precision  :: tau, dc1, dc2, tc1, tc2 ,R
End Module Parameters

 
Program main

  use expdata
  
  implicit none
  double precision :: dm
  
  write(*,*)"**********************************************************************"
  write(*,*) "T-d range : 298.15-473.15 K and pressure to 100 bar"
  write(*,*) "x1(CO2) range: 0-1; x2(N2O)= 1-x1(CO2) "
  write(*,*)"**********************************************************************"
  write(*,*) 
  write(*,*) 

  100 write(*,*)  
  write(*,*)  'Plese input T(K), P(bar), x1, x2'
  read(*,*)  T, P, x1, x2     
       

  call GET_DE(T,P,D_EXP)
  dm=D_EXP
      
  call PVT_CO2_N2O( T, P, x1, x2, dm )
  
  write(*,*) 
  write(*,*)     
  write(*,*)'the results are as follows:'
  write(*,*)"**********************************************************************"
  write(*,'(4a10,a28)') "T(K)", "P(bar)","x1", "x2","d_CO2-N2O(mol/dm3)"
  write(*,*)"______________________________________________________________________"
  write(*,'(2f10.3,2f11.2,f20.4)') T,P,x1,x2,dm
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)
  
  open(40,file='PVTx.out.txt')
  write(40,'(a10,a13,a13,a13,a28)') "T(K)", "P(bar)","x1", "x2","dCO2-N2O(mol/dm3)"
  write(40,'(2f10.3,2f12.2,f20.4)') T,P,x1,x2,dm
  

  go to 100
      
  stop
  
end Program main
      
      
Subroutine GET_DE(Ts,Ps,D_EXP)

  implicit none
  double precision :: Ts, Ps, D_EXP
  double precision :: Denls, Dengs, PPS, T
  double precision :: a(4), b(4), c(5) 
  double precision, parameter :: Tc = 304.1282      
  double precision, parameter :: Pc = 73.773         
  double precision, parameter :: denc = 467.6  
  double precision, parameter :: R= 8.314472d0 / 44.0098
  
  data a / -7.0602087, 1.9391218, -1.6463597, -3.2995634 /
  data b / 1.9245108, -0.62385555, -0.32731127, 0.39245142 /
  data c / -1.7074879, -0.82274670, -4.6008549, -10.111178, -29.742252 / 
  
      t = 1.0 - Ts/Tc
	  
	  Pps = Pc * dexp( Tc/Ts * ( a(1)*t + a(2)*t**1.5 + a(3)*t**2.0 +  a(4)*t**4.0) )

	  denls = denc * dexp( b(1)*t**(1.0/3.0) + b(2)*t**(1.0/2.0) + b(3)*t**(10.0/6.0)              &    
	                 + b(4)*t**(11.0/6.0) )
     
	  dengs = denc * dexp ( Tc/Ts * ( c(1)*t**(1.0/3.0) + c(2)*t**(1.0/2.0) + c(3)*t+ &
	                       c(4)*t**(7.0/3.0) + c(5)*t**(14.0/3.0) ) )
      if(Ts<=Tc.AND.Ps>=PPS)then
          D_EXP=Denls
      else
          D_EXP=Ps/(R*Ts)
      end if 
      return
      
end subroutine GET_DE
    
    
Subroutine PVT_CO2_N2O( T, P, x1, x2, dm )
  use parameters                
  implicit none
  double precision  :: T, P            
  double precision  :: dm                 
  double precision  :: dc                  
  double precision  :: Tc                  
  double precision  :: x1, x2
  double precision  :: delta, deltacal
  double precision  :: Para(4)
  data Para / -0.31722984D+00, 0.19725457D-01, -0.75656513D+02, 0.13796067D+01 /

  R = 8.314472  
  
  Tc1 = 304.1282  
  dc1 = 0.4676   
  
  Tc2 = 309.52   
  dc2 = 0.45201   

  Tc = x1*Tc1 + x2*Tc2  + x1**para(4)*x2*para(3)  
  
  dc = 1.0 / ( x1/(dc1*1000./44.0098) + x2/(dc2*1000./44.0128) + x1*x2*para(2)  ) 

  tau = Tc / T
  delta = dm / dc

  Call DNewt( T, P, x1, x2,  delta, deltacal  )

  delta = deltacal
  dm = dc * delta  
  
  return
end subroutine PVT_CO2_N2O
    
    
      
Subroutine DNewt( T, P, xCO2, xN2O,  x, x_OUT )
    Double precision, parameter   :: EPS=1.0d-8
	double precision              :: T, P, xCO2, xN2O
	Double precision              :: x, x1, F, DF,x_OUT
 
    L= 60
	F= 0.0
	DF=0.0
	
	Call F_CO2_N2O( T, P, xCO2, xN2O,  x, F, DF )

10  IF(dABS(DF)+1.0.EQ.1.0)THEN
      L = 0
	  WRITE(*,20)
	  RETURN
	ENDIF

20  FORMAT(1X,'ERR')
    x1 = x-F/DF
	CALL F_CO2_N2O( T, P, xCO2, xN2O,  x1, F, DF  )
	
	IF((dABS(x1-x).GE.EPS).OR.(dABS(F).GE.EPS))THEN
	  L = L - 1
	  x = x1
	  IF(L.EQ.0)RETURN
	  GOTO 10
	ENDIF
	n_time=60-l
	x_Out = x1
	RETURN

 End Subroutine DNewt
      
Subroutine F_CO2_N2O( T, P, x1, x2,  delta, F, dF  )
use parameters
  implicit none
  double precision  :: T, P, x1, x2, delta, dc
  double precision  :: phir, phir1, phir2, phir0
  double precision  :: phird, phird1, phird2, phird0
  double precision  :: phirdd, phirdd1, phirdd2, phirdd0
  double precision  ::  F, dF
  double precision  :: Para(4)
  data Para / -0.31722984D+00, 0.19725457D-01, -0.75656513D+02, 0.13796067D+01 /

   
	call phirdd_CO2( tau, delta, phir1, phird1, phirdd1 )

    call phirdd_N2O( tau, delta, phir2, phird2, phirdd2 )

	call phirdd_binary( tau, delta, phir0, phird0, phirdd0 )
	
    phir =  x1*phir1  + x2*phir2  + x1*x2*phir0*para(1)
    
	phird = x1*phird1 + x2*phird2 + x1*x2*phird0*para(1) 
	
	phirdd = x1*phirdd1 + x2*phirdd2 + x1*x2*phirdd0*para(1) 
	
	dc = 1.0 / ( x1/(dc1*1000./44.0098) + x2/(dc2*1000./44.0128)   + x1*x2*para(2)  ) 
  
    F = P - ( dc * R * T * delta + dc * R * T *delta**2 * phird ) / 100 

    dF = ( -dc * R * T - 2.0 * dc * R * T * delta * phird - dc * R * T * delta**2 * phirdd ) / 100
    
    Return

 End Subroutine F_CO2_N2O
    
    

 Subroutine phirdd_CO2( tau, delta, phir, phird, phirdd )
  implicit none
  double precision  :: tau, delta, phir, phird, phirdd
  integer           :: i
  double precision  :: ii(14),jj(14),aa(14),kk(14)

  data aa / 0.58478135D+00, 0.95777220D+00, -0.26415042D+01, 0.74935594D-01, 0.21519215D-03, -0.26666316D-01, 0.13294773D-01, &
            0.12546908D+00, 0.42015287D+00, -0.12243601D-01, -0.25984348D+00, -0.17963087D-01, -0.61607536D-01, -0.19553446D-01 /   
  
  data ii / 1, 1, 1, 3, 7, 2, 1, 1, 2, 5, 1, 1, 4, 2 /
	
  data jj / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /
	
  data kk / 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3 /
  
!----------------------------------to calculate phir-----------------------------------------
 phir = 0.0
  
  do i = 1, 6
    phir = phir + aa(i) * delta**ii(i) * tau**jj(i)
  end do
  
  do i = 7, 14
    phir = phir + aa(i) * delta**ii(i) * tau**jj(i) * dexp( -delta**kk(i) )  
  end do

!-----------------------------------------------------------------------------------------------
!----------------------------------to calculate phird-----------------------------------------
  phird = 0.0
  
  do i =1, 6
    phird = phird + aa(i) * ii(i) * delta**(ii(i)-1) * tau**jj(i)   
  end do

  do i = 7, 14
    phird = phird + aa(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-1)*tau**jj(i)* ( ii(i)-kk(i)*delta**kk(i) )
  end do

!-----------------------------------------------------------------------------------------------

!------------------------------------to calculate phirdd----------------------------------------
  phirdd = 0.0
  do i = 1, 6
    phirdd = phirdd + aa(i) * ii(i) * ( ii(i)-1 ) * delta**(ii(i)-2) * tau**jj(i)
  end do

  do i = 7, 14
    phirdd = phirdd + aa(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-2) * tau**jj(i) *  &
	         ( ( ii(i)-kk(i)*delta**kk(i) )*( ii(i)-1-kk(i)*delta**kk(i) ) - kk(i)**2 * delta**kk(i) )
  end do

!------------------------------------------------------------------------------------------------------

  return
 end subroutine phirdd_CO2

    
    
    
Subroutine phirdd_N2O( tau, delta, phir, phird, phirdd )
  implicit none
  double precision  :: tau, delta, phir, phird, phirdd
  integer           :: i
  double precision  :: ii(14),jj(14),aa(14),kk(14)
	
  data ii / 1.0, 1.0, 1.0, 3.0, 7.0, 2.0, 1.0, 1.0, 2.0, 5.0, 1.0, 1.0, 4.0, 2.0 /

  data jj / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /

  data kk / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0 /


 data aa /  0.40209218D+00, 0.87844892D+00, -0.24241584D+01, 0.69011106D-01, 0.20412510D-03,&
	        -0.52403092D-02, -0.43148215D-02, 0.11046048D+00, 0.44400973D+00, -0.35551109D-02,&
			-0.22637699D+00, -0.43690514D-03, -0.41132000D-01, -0.22998680D-01/
  

!----------------------------------to calculate phir-----------------------------------------
  phir = 0.0
  
  do i = 1, 6
    phir = phir + aa(i) * delta**ii(i) * tau**jj(i)
  end do
  
  do i = 7, 14
    phir = phir + aa(i) * delta**ii(i) * tau**jj(i) * dexp( -delta**kk(i) )  
  end do

!-----------------------------------------------------------------------------------------------
!----------------------------------to calculate phird-----------------------------------------
  phird = 0.0
  
  do i =1, 6
    phird = phird + aa(i) * ii(i) * delta**(ii(i)-1) * tau**jj(i)
  end do

  do i = 7, 14
    phird = phird + aa(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-1)*tau**jj(i)* ( ii(i)-kk(i)*delta**kk(i) )
  end do

!-----------------------------------------------------------------------------------------------

!------------------------------------to calculate phirdd----------------------------------------
  phirdd = 0.0
  do i = 1, 6
    phirdd = phirdd + aa(i) * ii(i) * ( ii(i)-1 ) * delta**(ii(i)-2) * tau**jj(i)
  end do

  do i = 7, 14
    phirdd = phirdd + aa(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-2) * tau**jj(i) *  &
	         ( ( ii(i)-kk(i)*delta**kk(i) )*( ii(i)-1-kk(i)*delta**kk(i) ) - kk(i)**2 * delta**kk(i) )
  end do
 !----------------------------------------------------------------------------------------------- 
  

  return
end subroutine phirdd_N2O

!---------------------------------------------------------------------------------------------------------

Subroutine phirdd_binary( tau, delta, phir, phird, phirdd )
  implicit none
  double precision  :: tau, delta, phir, phird, phirdd
  double precision  :: Nk(10), Ik(10), Jk(10)
  integer            :: i
  
  data Nk / -2.45476271425d-2, -2.41206117483d-1, -5.13801950309d-3, -2.39824834123d-2, 2.59772344008d-1, &
            -1.72014123104d-1,  4.29490028551d-2, -2.02108593862d-4, -3.82984234857d-3, 2.629923313540d-6 /
  data Ik / 1, 1, 1, 2, 3, 4, 5, 6, 6, 8 /
  data Jk / 2, 4, -2, 1, 4, 4, 4, 0, 4, -2 /
  
  phir = 0.0
  do i = 1, 10
    phir = phir + Nk(i) * delta**Ik(i) * tau**Jk(i)
  end do
    
  phird = 0.0
  do i = 1, 10
    phird = phird + Nk(i) * Ik(i) * delta**(Ik(i)-1) * tau**Jk(i)
  end do

  
  phirdd = 0.0
  do i = 1, 10
    phirdd = phirdd + Nk(i) * Ik(i) * (Ik(i)-1) * delta**(Ik(i)-2) * tau**Jk(i)
  end do
  

  return
end subroutine phirdd_binary