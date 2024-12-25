
Program main

  Implicit none
  
  double precision   :: P 
  double precision   :: T, x1, x2, dm
  
  write(*,*)"**********************************************************************"
  write(*,*) "T-d_CO2-N2O range : 298.15-523.15 K and density to 7.6022 mol/dm3"
  write(*,*) "x1(CO2) range: 0-1; x2(N2O)= 1-x1(CO2) "
  write(*,*)"**********************************************************************"
  write(*,*) 
  write(*,*) 
  
  
100 write(*,*)    
  write(*,*) 'Plese input T(K), x1, x2, d_CO2-N2O(mol/dm3)'
  read(*,*) T, x1, x2, dm
 
  call isochore_CO2_N2O( T, x1, x2, dm, P )
      
  write(*,*) 
  write(*,*)     
  write(*,*)'the results are as follows:'
  write(*,*)"**********************************************************************"
  write(*,'(3a10,a28,a10)') "T(K)","x1", "x2","d_CO2-N2O(mol/dm3)", "P(bar)"
  write(*,*)"______________________________________________________________________"
  write(*,'(f10.3,2f11.2,f20.4,f15.3)') T, x1, x2, dm, P
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)      
  
  
  open(40,file='isochore.out.txt')    
  write(40,'(3a15,a28,a15)') "T(K)","x1", "x2","dCO2-N2O(mol/dm3)", "P(bar)" 
  write(40,'(f15.3,2f13.2,f22.4,f32.3)') T, x1, x2, dm, P
  
  
   go to 100
    
  stop
End Program main
      

    
Subroutine isochore_CO2_N2O( T, x1, x2, dm, P)             
  implicit none
  double precision  :: tau, dc1, dc2, tc1, tc2 ,R
  double precision  :: T, P                       
  double precision  :: dm                         
  double precision  :: dc                         
  double precision  :: Tc                         
  double precision  :: x1, x2
  double precision  :: delta
  double precision  :: phird, phird1, phird2, phird0
  double precision  :: para(4)
  DATA para /-0.31722984D+00, 0.19725457D-01, -0.75656513D+02, 0.13796067D+01 /
 
  R = 8.314472    
  
  Tc1 = 304.1282   
  dc1 = 0.4676     
  
  Tc2 = 309.52     
  dc2 = 0.45201    

  Tc = x1*Tc1 + x2*Tc2  + x1**para(4)*x2*para(3)                                 
  dc = 1.0 / ( x1/(dc1*1000./44.0098) + x2/(dc2*1000./44.0128) + x1*x2*para(2)  ) 
  
  tau = Tc / T
  delta = dm / dc
  	
  
  call phird_CO2( tau, delta, phird1)

  call phird_N2O( tau, delta, phird2)

  call phird_binary( tau, delta, phird0)

    
  phird = x1*phird1 + x2*phird2 + x1*x2*phird0*para(1) 
  

  P = ( dc * R * T * delta + dc * R * T *delta**2 * phird ) / 100

  return
end subroutine isochore_CO2_N2O
    


 Subroutine phird_CO2( tau, delta, phird)
  implicit none
  double precision  :: tau, delta,  phird
  integer           :: i
  double precision  :: ii(14), jj(14), aa(14), kk(14)


  data aa / 0.58478135D+00, 0.95777220D+00, -0.26415042D+01, 0.74935594D-01, 0.21519215D-03, -0.26666316D-01, 0.13294773D-01, &
            0.12546908D+00, 0.42015287D+00, -0.12243601D-01, -0.25984348D+00, -0.17963087D-01, -0.61607536D-01, -0.19553446D-01/
  
  data ii / 1, 1, 1, 3, 7, 2, 1, 1, 2, 5, 1, 1, 4, 2 /
	
  data jj / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /
	
  data kk / 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3 /
  

!----------------------------------to calculate phird-----------------------------------------
  phird = 0.0
  
  do i =1, 6
    phird = phird + aa(i) * ii(i) * delta**(ii(i)-1) * tau**jj(i)   
  end do

  do i = 7, 14
    phird = phird + aa(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-1)*tau**jj(i)* ( ii(i)-kk(i)*delta**kk(i) )
  end do

!-----------------------------------------------------------------------------------------------

  return
 end subroutine phird_CO2

    
    
    
Subroutine phird_N2O( tau, delta,  phird)
  implicit none
  double precision  :: tau, delta,  phird
  integer           :: i
  double precision  :: ii(14),jj(14),aa(14),kk(14)
	
  data ii / 1.0, 1.0, 1.0, 3.0, 7.0, 2.0, 1.0, 1.0, 2.0, 5.0, 1.0, 1.0, 4.0, 2.0 /

  data jj / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /

  data kk / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0 /


 data aa /  0.40209218D+00, 0.87844892D+00, -0.24241584D+01, 0.69011106D-01, 0.20412510D-03,&
	        -0.52403092D-02, -0.43148215D-02, 0.11046048D+00, 0.44400973D+00, -0.35551109D-02,&
			-0.22637699D+00, -0.43690514D-03, -0.41132000D-01, -0.22998680D-01/
  

!----------------------------------to calculate phird-----------------------------------------
  phird = 0.0
  
  do i =1, 6
    phird = phird + aa(i) * ii(i) * delta**(ii(i)-1) * tau**jj(i)
  end do

  do i = 7, 14
    phird = phird + aa(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-1)*tau**jj(i)* ( ii(i)-kk(i)*delta**kk(i) )
  end do

!-----------------------------------------------------------------------------------------------


  return
end subroutine phird_N2O

!---------------------------------------------------------------------------------------------------------

Subroutine phird_binary( tau, delta,  phird )
  implicit none
  double precision  :: tau, delta,  phird
  double precision  :: Nk(10), Ik(10), Jk(10)
  integer            :: i
  
  data Nk / -2.45476271425d-2, -2.41206117483d-1, -5.13801950309d-3, -2.39824834123d-2, 2.59772344008d-1, &
            -1.72014123104d-1,  4.29490028551d-2, -2.02108593862d-4, -3.82984234857d-3, 2.629923313540d-6 /
  data Ik / 1, 1, 1, 2, 3, 4, 5, 6, 6, 8 /
  data Jk / 2, 4, -2, 1, 4, 4, 4, 0, 4, -2 /  
    
  phird = 0.0
  do i = 1, 10
    phird = phird + Nk(i) * Ik(i) * delta**(Ik(i)-1) * tau**Jk(i)
  end do


  return
end subroutine phird_binary