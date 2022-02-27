subroutine p_minhash_1(W, n, m, upper, z)
    implicit none
    integer, intent(in) :: n, m
    integer(4), intent(in) :: upper
    integer(4), intent(in), dimension(1:2, 1:n) :: W
    integer(4), intent(inout), dimension(1:m) :: z
    real(8) :: inv, h
    real(8), dimension(1:m) :: q
    integer :: i, j, k, d
    real(8) :: mu
    real(8) :: ZBQLEXP, ZBQLUAB

    ! f2py intent(in) :: n, m, upper, W
    ! f2py intent(in,out) :: z

    mu = 1
    ! Initialize q
    do i = 1,m
        q(m) = real(upper, 8)
    end do

    
    do d = 1,n
        inv = 1 / real(W(2, i), 8)
        call ZBQLINI(W(1, d))
        h = inv * ZBQLEXP(mu)
        do while (h .lt. MAXVAL(q))
            k = floor(ZBQLUAB(real(1, 8), real(m+1, 8)))
            if (k .eq. (m+1)) then
                k = m
            end if
            if (h .lt. q(k)) then
                q(k) = h
                z(k) = W(1, d)
                ! update qmax using Algorithm 4
                do while (h .lt. q(k))
                    q(k) = h
                    i = m + ceiling(real(k, 8)/2)
                    if (i .ge. (2*m)) then
                        exit
                    end if
                    j = IEOR((k-1), 1) + 1
                    if (q(j) .ge. q(i)) then
                        exit
                    end if
                    if (h .lt. q(j)) then
                        h = q(j)
                    end if
                    k = i
                end do 

                if (h .ge. MAXVAL(q)) then
                    exit
                end if                   
            h = h + inv * ZBQLEXP(mu)
            end if
        end do
    end do

    return

end subroutine p_minhash_1
