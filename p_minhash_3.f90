subroutine p_minhash_3(W, n, m, upper, z)
    implicit none
    integer, intent(in) :: n, m
    integer(4), intent(in) :: upper
    integer(4), intent(in), dimension(1:2, 1:n) :: W
    integer(4), intent(inout), dimension(1:m) :: z
    real(8) :: inv, h
    real(8), dimension(1:m) :: q
    integer :: i, j, k, d
    real(8) :: lambda
    real(8) :: ZBQLUAB, ZEXPT01
    integer :: size
    integer(4),allocatable :: seed(:)


    ! f2py intent(in) :: n, m, upper, W
    ! f2py intent(in,out) :: z

    if (m .lt. 2) then
        print *,'m must be great than 1 !'
        return
    end if

    lambda = LOG(1.0D0 + 1.0D0/(m-1.0D0))
    ! Initialize q
    do i = 1,m
        q(i) = real(upper, 8)
    end do

    size = 16
    do d = 1,n
        inv = 1 / real(W(2, d), 8)
        ! set seed
        if(allocated(seed)) then
            deallocate(seed)
        endif
        call random_seed(size=size)
        allocate(seed(size))
        seed = W(1, d)
        call random_seed(put=seed)
        h = inv * ZEXPT01(lambda)
        i =1
        do while (h .lt. MAXVAL(q))
            k = floor(ZBQLUAB(real(1, 8), real(m+1, 8)))
            if (k .eq. (m+1)) then
                k = m
            end if
            if (h .lt. q(k)) then
                q(k) = h
                z(k) = W(1, d)
                ! update qmax using Algorithm4
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
            h = inv * (i - 1.0D0)
            if (h .ge. MAXVAL(q)) then
                exit
            end if
            h = h + inv * ZEXPT01(lambda)

        end do


    end do


    return

end subroutine p_minhash_3