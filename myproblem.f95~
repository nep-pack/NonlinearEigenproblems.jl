subroutine compute_Mder_legacy(s,n,der,I,J,F)
  real*8, intent(in) :: s
  integer*8, intent(in) :: n
  integer*8, intent(in) :: der
  integer*8, intent(out), dimension(3*n):: I
  integer*8, intent(out), dimension(3*n):: J
  real*8, intent(out), dimension(3*n):: F
  integer*8 :: p
  real*8 :: factor
  if (der==0) then
     factor=1;
  else
     factor=0;
  end if
  do p = 1, n
     I(p) = p
     J(p) = p
     F(p) = 2.0*factor;
  end do
  do p = 1, n-1
     I(n+p) = p
     J(n+p) = p+1
     F(n+p) = -1.0*factor;
     I(2*n+p) = p+1
     J(2*n+p) = p
     F(2*n+p) = -1.0*factor;
  end do
  I(2*n)=n;
  J(2*n)=1;
  if (der == 0) then
     F(2*n)=s*s*s;
  else if (der == 1) then
     F(2*n)=3*s*s;
  else if (der == 2) then
     F(2*n)=3*2*s;
  else if (der == 3) then
     F(2*n)=3*2;
  else
     F(2*n)=0;
  end if
  I(3*n)=1;
  J(3*n)=n;
  F(3*n)=-exp(s);
end subroutine compute_Mder_legacy
