Module Parameters        
  Double precision  :: dc                 
  double precision  :: R
  double precision  :: Tc                     
  double precision  :: tau                   
End Module Parameters

Program main
  Implicit none
  double precision  :: T,d_N2O,P      
  
  write(*,*)"**********************************************************************"
  write(*,*) "T-d_N2O range : 183.33-525 K and density to 28.11 mol/dm3"
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)
  
  
  100  write(*,*) 
  write(*,*) "Please input T(K) and d_N2O(mol/dm3)  :"
  read(*,*) T,d_N2O

  if( T < 182.33 .or. T >525 .or. d_N2O < 0 .or. d_N2O > 28.11 ) then
  write(*,*) "Input error, please try again!"
    goto 100 
  end if
  write(*,*)
  write(*,*)
  
  call isochore_N2O( T, d_N2O, P )
  
  write(*,*)'the results are as follows:'
  write(*,*)"**********************************************************************"
  write(*,'(a9,a20,a15)') "T(K)","d_N2O(mol/dm3)","P(bar)"
  write(*,*)"______________________________________________________________________"
  write(*,'(f10.2,f13.4,f20.2)') T,d_N2O, P
  write(*,*)"**********************************************************************"
  write(*,*)
  write(*,*)
  
  open(40,file='isochore.out.txt')
  write(40,'(a20,a30,a15)') "T(K)","dN2O(mol/dm3)","P(bar)"
  write(40,'(f20.2,f25.4,f25.2)') T,d_N2O, P
  
  go to 100
  stop
end program main


Subroutine isochore_N2O( T, d, P ) 
   use parameters
   implicit none
   integer           :: i  
   double precision  :: T, d, P        
   double precision  :: phird        
   double precision  :: delta 
   double precision  :: ii(14), jj(14), kk(14)
   double precision  :: para(14)  
   
   data ii / 1.0, 1.0, 1.0, 3.0, 7.0, 2.0, 1.0, 1.0, 2.0, 5.0, 1.0, 1.0, 4.0, 2.0 /

   data jj / 1.5, 0.25, 1.25, 0.25, 0.875, 1.375, 0.0, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5 /

   data kk / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0 /
 
   data para / 0.40209218D+00, 0.87844892D+00, -0.24241584D+01, 0.69011106D-01, 0.20412510D-03, -0.52403092D-02, -0.43148215D-02, &
               0.11046048D+00, 0.44400973D+00, -0.35551109D-02, -0.22637699D+00, -0.43690514D-03, -0.41132000D-01, -0.22998680D-01 /
 
   dc = 452.01 / 44.0128           
   Tc = 309.52                     
   R = 8.314472                    
      
   delta = d / dc              
   tau = Tc / T

!----------------------------------to calculate phird-----------------------------------------
  phird = 0.0
  
  do i =1, 6
    phird = phird + para(i) * ii(i) * delta**(ii(i)-1) * tau**jj(i)
  end do

  do i = 7, 14
    phird = phird +para(i) * dexp( -delta**kk(i) ) * delta**(ii(i)-1)*tau**jj(i)* ( ii(i)-kk(i)*delta**kk(i) )
  end do

   P = ( dc * R * T * delta + dc * R * T * delta**2 * phird ) /100  

End Subroutine isochore_N2O