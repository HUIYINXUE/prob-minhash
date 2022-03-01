subroutine p_minhash_4(W, n, m, upper, z)
    implicit none
    integer, intent(in) :: n, m
    integer(4), intent(in) :: upper
    integer(4), intent(in), dimension(1:2, 1:n) :: W
    integer(4), intent(inout), dimension(1:m) :: z
    real(8) :: inv, h, tmp0
    real(8), dimension(1:m) :: q
    real(8), dimension(1:m-1) :: lambda, gamma
    integer :: i, j, k, d, c0
    real(8) :: ZBQLUAB, ZEXPT01, ZBQLEXP
    integer(4), dimension(1:m) :: v, g
    integer :: size
    integer(4),allocatable :: seed(:)


    ! f2py intent(in) :: n, m, upper, W
    ! f2py intent(in,out) :: z

    if (m .lt. 2) then
        print *,'m must be great than 1 !'
        return
    end if

    
    ! Initialize q
    do i = 1,m-1
        q(i) = real(upper, 8)
        v(i) = 0
        g(i) = i
        lambda(i) = LOG(1.0D0 + 1.0D0 / real(m-i, 8))
        gamma(i) = LOG(1.0D0 + real(i, 8) / real(m-i, 8)) / lambda(1)
    end do
    q(m) = real(upper, 8)
    v(m) = 0
    g(m) = m
    ! InitPermutationGenerator
    c0 = 0
    size = 16

    do d=1,n
        inv = 1 / real(W(2, d), 8)
        ! set seed
        if(allocated(seed)) then
            deallocate(seed)
        endif
        call random_seed(size=size)
        allocate(seed(size))
        seed = W(1, d)
        call random_seed(put=seed)
        ! Reset PermutationGenerator
        i = 1
        c0 = c0 + 1
        h = inv * ZEXPT01(lambda(1))
        do while (h .lt. MAXVAL(q))
            ! GenerateNextPermutationElement(R)
            i = i + 1
            k = floor(ZBQLUAB(real(0, 8), real(m-i+1, 8)))
            j = i + k
            if (v(j) .eq. c0) then
                k = g(j)
            else
                k = j
            end if
            if (v(i) .eq. c0) then
                g(j) = g(i)
            else
                g(j) = i
            end if
            v(j) = c0


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
            end if

            i = i + 1
            if ((inv * gamma(i-1)) .ge. MAXVAL(q)) then
                exit
            end if
            if (i .lt. m) then
                tmp0 = (gamma(i) - gamma(i-1)) * ZEXPT01(lambda(i))
                h = inv * (gamma(i-1) + tmp0)
            else
                h = inv * (ZBQLEXP(1.0D0) / lambda(1) + gamma(m-1))
                if (h .lt. MAXVAL(q)) then
                    ! GenerateNextPermutationElement(R)
                    i = i + 1
                    k = floor(ZBQLUAB(real(0, 8), real(m-i+1, 8)))
                    j = i + k
                    if (v(j) .eq. c0) then
                        k = g(j)
                    else
                        k = j
                    end if
                    if (v(i) .eq. c0) then
                        g(j) = g(i)
                    else
                        g(j) = i
                    end if
                    v(j) = c0

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
                end if
                exit
            end if
                

        end do


    end do


return
end subroutine p_minhash_4