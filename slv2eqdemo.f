	program main
	implicit none
	
	integer, parameter :: n=2
	double precision a(n,n), c(n,n), r(n)
	double precision y(n)
	integer i,j,m
	character ans*1
	character ch(n,n)
	data ch/'a','b','c','d'/
	
	! matrix A
10	print *, "Please enter system of equations"

	do i=1,n
		do j=1,n
			print *, "Enter the value of", ch(j,i) 
			read (*,*)a(i,j)
		end do
	end do
	
	do i=1,n
		print *, "Enter the equality values:"
		read (*,*)r(i)
	end do

	!print confirmation for system of equations
	write (*,203)
	do i=1,n
		write (*,*) (a(i,j),j=1,n),"=", (r(i))
	end do
	!300	format (2x,a,i1,a1f10.5)

	print *, "Is the system of equations setup correctly? (Y/N)"
	read *, ans
	print *,
	if (ans .eq. 'N' .or. ans .eq. 'n') go to 10

	! print a header and the original matrix
20  	write (*,200)
  	
	do i=1,n
     		write (*,201) (a(i,j),j=1,n)
  	end do

  	call inverse(a,c,n)

	! print the inverse matrix C = A^{-1} 
	write (*,202)
  	do i = 1,n
     		write (*,201)  (c(i,j),j=1,n)
  	end do

	!print the multiplication
	write (*,204)
	
	call multi(c,r,y)
	do i=1,n
		print 100, "x", i, "=", y(i)
	end do
100	format (2x,a,i1,a,f10.5)
	
200 	format (/, ' Computing Inverse matrix ',/,/, ' Matrix A')
201 	format (3f12.6)
202 	format (/, ' Inverse matrix A^{-1}')
203	format (/, 'The equations are:')
204	format (/, ' The value of x and y are: ')
	end

	

	subroutine inverse(a,c,n)
    	implicit none 
	integer n
	double precision a(n,n), c(n,n)
	double precision L(n,n), U(n,n), b(n), d(n), x(n)
	double precision coeff
	integer i, j, k

	! step 0: initialization for matrices L and U and b
	L=0.0
	U=0.0
	b=0.0

	! step 1: forward elimination
	do k=1, n-1
	   do i=k+1,n
	      coeff=a(i,k)/a(k,k)
	      L(i,k) = coeff
	      do j=k+1,n
	         a(i,j) = a(i,j)-coeff*a(k,j)
	      end do
	   end do
	end do
	
	! Step 2: prepare L and U matrices 
	! L matrix is a matrix of the elimination coefficient
	! + the diagonal elements are 1.0
	do i=1,n
	  L(i,i) = 1.0
	end do
	! U matrix is the upper triangular part of A
	do j=1,n
	  do i=1,j
	    U(i,j) = a(i,j)
	  end do
	end do
	
	! Step 3: compute columns of the inverse matrix C
	do k=1,n
	  b(k)=1.0
	  d(1) = b(1)
	! Step 3a: Solve Ld=b using the forward substitution
	  do i=2,n
	    d(i)=b(i)
	    do j=1,i-1
	      d(i) = d(i) - L(i,j)*d(j)
	    end do
	  end do
	! Step 3b: Solve Ux=d using the back substitution
	  x(n)=d(n)/U(n,n)
	  do i = n-1,1,-1
	    x(i) = d(i)
	    do j=n,i+1,-1
	      x(i)=x(i)-U(i,j)*x(j)
	    end do
	    x(i) = x(i)/u(i,i)
	  end do
	! Step 3c: fill the solutions x(n) into column k of C
	  do i=1,n
	    c(i,k) = x(i)
	  end do
	  b(k)=0.0
	end do
	end subroutine inverse

	!Subroutine for multiply inverse and results
	subroutine multi(c,r,y)
	implicit none
	integer, parameter :: n=2
	double precision r(n), y(n), c(n,n)
	integer i,j,k
	
	do i=1,n
		do j=1,n
			y(i)=0
			do k=1,n
				y(i)=y(i)+c(i,k)*r(k)
			end do
		end do
	end do
	end subroutine multi
