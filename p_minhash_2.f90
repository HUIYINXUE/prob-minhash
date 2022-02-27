subroutine p_minhash_2(W, n, m, upper, z)
    implicit none
    integer, intent(in) :: n, m
    integer(4), intent(in) :: upper
    integer(4), intent(in), dimension(1:2, 1:n) :: W
    integer(4), intent(inout), dimension(1:m) :: z
    real(8) :: inv, h
    real(8), dimension(1:m) :: q
    integer :: i, j, k, d, c
    real(8) :: mu
    real(8) :: ZBQLEXP, ZBQLUAB
    integer(4), dimension(1:m) :: v, g


    ! f2py intent(in) :: n, m, upper, W
    ! f2py intent(in,out) :: z

    mu = 1
    ! Initialize q
    do i = 1,m
        q(m) = real(upper, 8)
        v(m) = 0
        g(m) = i
    end do

    ! InitPermutationGenerator
    c = 0

    do d = 1,n
        inv = 1 / real(W(2, d), 8)
        call ZBQLINI(W(1, d))
        ! ResetPermutationGenerator
        i = 1
        c = c + 1
        h = inv * ZBQLEXP(mu)
        do while(h .lt. MAXVAL(q))
            ! GenerateNextPermutationElement(R)
            i = i + 1
            k = floor(ZBQLUAB(real(0, 8), real(m-i+1, 8)))
            if (k .eq. m-i+1) then
                k = m-i
            end if
            j = i + k
            if (v(j) .eq. c) then
                k = g(j)
            else
                k = j
            end if
            if (v(i) .eq. c) then
                g(j) = g(i)
            else
                g(j) = i
            end if
            v(j) = c
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
            end if
            i = i + 1
            h = h + inv * (real(m, 8)/(m-i+1)) * ZBQLEXP(mu)
        end do

    end do
    return

end subroutine p_minhash_2