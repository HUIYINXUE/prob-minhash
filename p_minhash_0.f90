subroutine p_minhash_0(W, n, m, upper, z)
    implicit none
    integer, intent(in) :: n, m
    integer(4), intent(in) :: upper
    integer(4), intent(in), dimension(1:2, 1:n) :: W
    integer(4), intent(inout), dimension(1:m) :: z
    real(8) :: inv, h
    real(8), dimension(1:m) :: q
    integer :: i, k
    real(8) :: mu
    real(8) :: ZBQLEXP

    ! f2py intent(in) :: n, m, upper, W
    ! f2py intent(in,out) :: z

    mu = 1
    ! Initialize q
    do i = 1,m
        q(m) = real(upper, 8)
    end do

    
    do i = 1,n
        inv = 1 / real(W(2, i), 8)
        call ZBQLINI(W(1, i))
        do k = 1,m
            h = inv * ZBQLEXP(mu)
            if (h .lt. q(k)) then
                q(k) = h
                z(k) = W(1, i)
            end if
        end do
    end do

    return

end subroutine p_minhash_0
