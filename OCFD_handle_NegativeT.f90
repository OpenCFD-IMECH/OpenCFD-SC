  subroutine handle_NegativeT(i,j,k,Num_NegT)
  use flow_data 
  implicit none 
  integer:: i,j,k,i1,j1,k1,Num_NegT
  integer,parameter:: Max_NegT=10
  ! input you code here !!!
  i1=i_offset(npx)+i-1
  j1=j_offset(npy)+j-1
  k1=k_offset(npz)+k-1 
  print*, "Negative T found !",  i1,j1,k1 
  T(i,j,k)=(T(i+1,j,k)+T(i-1,j,k)+T(i,j+1,k)+T(i,j-1,k)+ T(i,j,k+1)+T(i,j,k-1))/6.d0
  Num_NegT=Num_NegT+1
  if(Num_NegT > 10) then 
   print*, "Number of Negative Temperature Points >  Limit, STOP !"  
   stop  
  endif 
  end 
  
