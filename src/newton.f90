subroutine newton(f, fp, x0, x, iters, debug)

    ! Estimate the zero of f(x) using Newton's method.
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!)
    !   the number of iterations iters

    implicit none
    interface
      function f(r)
        implicit none
        double precision, intent(in) :: r
        double precision f
      end function f
      function fp(r)
        implicit none
        double precision, intent(in) :: r
        double precision fp
      end function fp
    end interface
    integer, parameter :: maxiter = 40
    double precision, parameter :: tol = 1.d-6
    double precision, intent(in) :: x0
    ! double precision:: f, fp
    logical, intent(in) :: debug
    double precision, intent(out) :: x
    integer, intent(out) :: iters

    ! Declare any local variables:
    double precision :: deltax, fx, fxprime
    integer :: k


    ! initial guess
    x = x0

    if (debug) then
        print 11, x
 11     format('Initial guess: x = ', e22.15)
        endif

    ! Newton iteration to find a zero of f(x)

    do k=1,maxiter

        ! evaluate function and its derivative:
        fx = f(x)
        fxprime = fp(x)

        if (abs(fx) < tol) then
            exit  ! jump out of do loop
            endif

        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        x = x - deltax

        if (debug) then
            print 12, k,x
 12         format('After', i3, ' iterations, x = ', e22.15)
            endif

        enddo


    if (k > maxiter) then
        ! might not have converged

        fx = f(x)
        if (abs(fx) > tol) then
            print *, '*** Warning: has not yet converged'
            endif
        endif

    ! number of iterations taken:
    iters = k-1


end subroutine newton
