MODULE SORT
CONTAINS


! Source:   https://github.com/certik/fortran-utils
!
! Copyright (c) 2012 Ondřej Čertík

! Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!The above copyright notice and this permission notice shall be included in
!all copies or substantial portions of the Software.
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
!THE SOFTWARE.


function rargsort(a) result(b)
! Returns the indices that would sort an array.
!
! Arguments
! ---------
!
real(kind = 8), intent(in):: a(:)   ! array of numbers
integer :: b(size(a))         ! indices into the array 'a' that sort it
!
! Example
! -------
!
! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

integer :: N                           ! number of numbers/vectors
integer :: i,imin                      ! indices: i, i of smallest
integer :: temp1                       ! temporary
real(kind = 8) :: temp2
real(kind = 8) :: a2(size(a))
a2 = a
N=size(a)
do i = 1, N
b(i) = i
end do
do i = 1, N-1
! find ith smallest in 'a'
imin = minloc(a2(i:),1) + i - 1
! swap to position i in 'a' and 'b', if not already there
if (imin /= i) then
temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
temp1 = b(i); b(i) = b(imin); b(imin) = temp1
end if
end do
end function


END MODULE SORT
