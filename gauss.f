	program gaussian
	implicit none

	!To solve the Gauss Law in integral form for a spherical shape
	
	real eps, r, q, pi
	real E, flux, flux2, theta, area, check, n
	parameter (pi=3.14159)
	parameter (eps = 8.85419E-12 ) !(C**2)/(N*m**2)
	
	!Section for selecting the kind of shape
	!print *, "Enter the number for the shape"
	!print *, "1) sphere, 2) circle, 3)cylindrical, 4)"
	!read *, n
	!go to (1,2,3)n
	
	print *, "What is the charge?"
	read *, q
	print *, "Whar is the radius?"
	read *, r
	
	!Calculatte the elctric field of a sphere
	E = (q/(4*pi*eps*r**2))
	

	print *, "The electric field is: ", E

	!Calculation the electric flux
	flux = E*(r**2)*(sin(pi/2)+sin(pi/2))*(2*pi-0)
	
	!Confirm Gauss' Law	
	flux2= q/eps

	!Printing results
	check =abs( ((flux - flux2)/((flux+flux2)/2))*100 )
	if (abs(check) .le. 1.0D-5) then
		write (*,200)
		write (*,201) flux, "(N*m^2)/C"		
	else	
		!print *, /, "The electric flux is:", flux, "(N*m^2)/C"
		print *, "Error in calculation of Gauss' Law..."
		!print *, "flux is: ", flux
		!print *, "equality is:", flux2
		!print *, "what was the difference:", check, "%"

	end if
	
200	format (/, 'The electric flux is:')
201	format (d20.6,a)

	end program gaussian

	
