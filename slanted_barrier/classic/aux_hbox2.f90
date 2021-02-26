module aux_module_hbox

  use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
  use, intrinsic :: iso_fortran_env, only: real32
  implicit none
contains

 !
        function dist(A,B) result(distance)
          implicit none
          real(kind=8) :: A(2),B(2)
          real(kind=8) :: distance
          distance = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
        end function dist

      SUBROUTINE KB07AD(COUNT,N,INDEX)
!  TO BE SORTED.

!     .. Scalar Arguments ..
      INTEGER N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION COUNT(*)
      INTEGER INDEX(*)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION AV,X
      INTEGER I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M
      INTEGER MLOOP
!     ..
!     .. Local Arrays ..
      INTEGER MARK(100)

!     ..
!     .. Executable Statements ..
!  SET INDEX ARRAY TO ORIGINAL ORDER .
      DO 10 I = 1,N
        INDEX(I) = I
   10 CONTINUE
!  CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED .
      IF (N.EQ.1) GO TO 200
      IF (N.GE.1) GO TO 30
      WRITE (6,FMT=20)

   20 FORMAT (/,/,/,20X,' ***KB07AD***NO NUMBERS TO BE SORTED ** ',&
      'RETURN TO CALLING PROGRAM')

      GO TO 200
!  'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
!  THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
   30 M = 12
!  SET UP INITIAL VALUES.
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
!  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE .
      IFKA = IF - IS
       IF ((IFKA+1).GT.M) GO TO 70
!********* FINAL SORTING ***
!  ( A SIMPLE BUBBLE SORT )
       IS1 = IS + 1
       DO 60 J = IS1,IF
         I = J
   40     IF (COUNT(I-1).LT.COUNT(I)) GO TO 60
         IF (COUNT(I-1).GT.COUNT(I)) GO TO 50
         IF (INDEX(I-1).LT.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          IF (I.GT.IS) GO TO 40
   60   CONTINUE
        LA = LA - 2
        GO TO 170
!             *******  QUICKSORT  ********
!  SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
!  THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S
!  HIGHEST ADDRESS.
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
!  THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
!  OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
!  OF X .
        K = 1
        IFK = IF
!  WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
!  INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS
!  NECESSARY, UNTIL THEY MEET .
        DO 110 I = IS,IF
          IF (X.GT.COUNT(I)) GO TO 110
          IF (X.LT.COUNT(I)) GO TO 80
          IF (INTEST.GT.INDEX(I)) GO TO 110
   80     IF (I.GE.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          DO 100 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).GT.X) GO TO 100
            IF (COUNT(IFK).LT.X) GO TO 90
            IF (INTEST.LE.INDEX(IFK)) GO TO 100
   90       IF (I.GE.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110

  100     CONTINUE
          GO TO 120

  110   CONTINUE
!  RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
!  WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
!  2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL
!  TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED
!  INDEPENDENTLY .
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140

  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
!  STORE THE LONGER SUBDIVISION IN WORKSPACE.
  140   IF ((IP-IS).GT. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160

  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
!  FIND THE LENGTH OF THE SHORTER SUBDIVISION.
  160   LNGTH = IF - IS
        IF (LNGTH.LE.0) GO TO 180
!  IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE .
        LA = LA + 2
        GO TO 190

  170   IF (LA.LE.0) GO TO 200
!  OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 CONTINUE
  200 RETURN

      END SUBROUTINE KB07AD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

FUNCTION area_polygon(x, y) RESULT(fn_val)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-04  Time: 12:24:06

IMPLICIT NONE

REAL(8), INTENT(IN)     :: x(:)
REAL(8), INTENT(IN)     :: y(:)
INTEGER  :: nb
REAL(8)                 :: fn_val, v1(2),v2(2)

!*****************************************************************

!   GIVEN A SEQUENCE OF NB POINTS (X(I),Y(I)),  polyarea COMPUTES THE AREA
! BOUNDED BY THE CLOSED POLYGONAL CURVE WHICH PASSES THROUGH THE POINTS IN
! THE ORDER THAT THEY ARE INDEXED.  THE FINAL POINT OF THE CURVE IS ASSUMED
! TO BE THE FIRST POINT GIVEN.  THEREFORE, IT NEED NOT BE LISTED AT THE END
! OF X AND Y.  THE CURVE IS NOT REQUIRED TO BE SIMPLE.  e.g. It may cross over
! itself.

!*****************************************************************

INTEGER  :: i, n, nm1
REAL     :: a

nb = size(x)
n = nb
a = 0.d0
do i=1,nb-1
  v1 = (/x(i),y(i)/)
  v2 = (/x(i+1),y(i+1)/)
  a = a+v1(1)*v2(2) - v2(1)*v1(2)
end do
fn_val = abs(a/2.d0)
end function

!
   function find_centroid(cycl) result(cent)
     implicit none
     real(8):: cycl(:,:)
     real(8) :: c_x,c_y,cent(2), det, tempdet
     integer :: n,j,i
     j = 1
     n = size(cycl,2)
     c_x = 0.d0
     c_y = 0.d0
     det = 0.d0
     do i=1,n
       if (i .eq. n) then
         j = 1
       else
         j = i + 1
       end if
       tempdet = cycl(1,i) * cycl(2,j) - cycl(1,j)*cycl(2,i)
       det = det + tempdet
       c_x = c_x + (cycl(1,i)+cycl(1,j))*tempdet
       c_y = c_y + (cycl(2,i)+cycl(2,j))*tempdet
     end do
     c_x = c_x/(3*det)
     c_y = c_y/(3*det)
     cent = (/c_x,c_y/)
   end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   subroutine find_ind(cycl,rough_guess,xe,ye,mx,my,mbc,right_ind,centroid)
    ! find the index of the cell that cycle-polygon covers
    !
    ! INPUT :: cycl (list), list of N vertices (list) of the (h-box fragment) polygon, first element same with last
    !          rough_guess (list), list of integers representing index close to actual index
    !          xe, ye (array), physical grid nodes
    ! OUTPUT :: right_ind (list), list of integers representing actual index of the cycle-polygon
    implicit none
    integer :: N,rough_guess(2),mbc,mx,my,right_ind(2),rot,i,right_ind_x,right_ind_y
    real(kind=8) :: cycl(:,:),xe(-mbc:mx+mbc),ye(-mbc:my+mbc),centroid(2)
    real(kind=8) :: A, c_x,c_y, Ver3(2),Ver2(2),Ver1(2),vec2(2),vec1(2)
    real(kind=8) :: V1(2), V2(2),x_guess,y_guess

    N = size(cycl,2)
    ! # centroid coordinate calculation:
    A = area_polygon(cycl(1,:),cycl(2,:))  ! switch the dmensions atrribute of lists... for all affected parts
    c_x = 0.d0
    c_y = 0.d0
    ! determine counter or clock:
    centroid = find_centroid(cycl)
    c_x = centroid(1)
    c_y = centroid(2)
    ! rough guess and adjusting
    x_guess = xe(rough_guess(1))
    y_guess = ye(rough_guess(2))

    ! search:
    if (c_x .le. x_guess) then
        do while (c_x .le. x_guess)
            rough_guess(1) = rough_guess(1)-1
            ! print*, rough_guess(1)
            x_guess = xe(max(rough_guess(1),-mbc))
        end do
        right_ind_x = rough_guess(1)
    else if (c_x .ge. x_guess) then
        do while (c_x .ge. x_guess)
            rough_guess(1) = rough_guess(1) + 1
            x_guess = xe(min(rough_guess(1),mx+mbc))
        end do
        right_ind_x = rough_guess(1)-1
    end if
    if (c_y .le. y_guess) then
        do while (c_y .le. y_guess)
            rough_guess(2) = rough_guess(2)-1
            y_guess = ye(max(rough_guess(2),-mbc))
        end do
        right_ind_y = rough_guess(2)
    else if (c_y .ge. y_guess) then
        do while (c_y .ge. y_guess)
            rough_guess(2) = rough_guess(2) + 1
            y_guess = ye(min(rough_guess(2),my+mbc))
        end do
        right_ind_y = rough_guess(2)-1
    end if
    right_ind(1) = right_ind_x+1
    right_ind(2) = right_ind_y+1
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ray_intersects_seg(P,A1,B1,intersects,on_line)
    ! checks if right going ray starting from P intersects line segment A1-B1,
    ! if point P is on the line, then returns true for "on_line" bool (needed to show such point will not be
    ! counted as inside polygon)
    !
    ! INPUT :: P (list), coordinate of starting point of ray
    !         A1, B1 (list,list), coordinates of line segment's endpoints
    ! OUTPUT :: intersects (bool), true if intersects or false if not
    !           on_line (bool), true if point P is on A1-B1 or false if not
    implicit none
    integer, allocatable :: index_ord(:)
    real(kind=8) :: P(2),A1(2),B1(2), points(2,2),A(2),B(2),tol
    real(kind=8) :: m_red,m_blue
    logical :: on_line, intersects
    tol = 1d-8
    on_line = .false.
    intersects = .false.
    points(1:2,1) = A1
    points(1:2,2) = B1
    allocate(index_ord(2))
    call KB07AD(points(2,1:2),2,index_ord)
    points(1,1:2) = points(1,index_ord)
    A = points(1:2,1)
    B = points(1:2,2)
    if (abs(P(2)-A(2)).le.1d-8 .or. abs(P(2)-B(2)).le.1d-8) then
        P(2) = P(2) + 1d-8
    end if
    if (P(2) .le. A(2) .or. P(2) .gt. B(2)) then
      GOTO 99
    else if (P(1) .ge. max(A(1),B(1))) then
      GOTO 99
    else
        if (P(1) .lt. min(A(1),B(1))) then
            intersects = .true.
            GOTO 99
        else
            if (abs(A(1)-B(1)) .gt. tol) then
                m_red = (B(2)-A(2))/(B(1)-A(1))
            else
                m_red = huge(0)
            end if
            if (abs(A(1)-P(1)) .gt. tol) then
                m_blue =(P(2)-A(2))/(P(1)-A(1))
            else
                m_blue = huge(0)
            end if
            if (abs(m_blue-m_red) .lt. tol) then
                on_line = .true.
               GOTO 99
            end if
            if (m_blue .gt. m_red) then
                intersects = .true.
                GOTO 99
            else
                intersects = .false.
                GOTO 99
            end if
          end if
      end if
  99 continue
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function in_or_out(P,sides,N) result(code)
    ! finds whether point P is inside or outside of polygon made by sides
    !
    ! INPUT :: P (list), coordinate of point to check whether in polygon
    !         sides (list), collection of *N* pairs of coordinates that form a polygon
    ! OUTPUT :: code (integer), 1 if inside or -1 if outside
    implicit none
    integer :: count,code,i,N
    real(kind=8) :: P(2),A1(2),B1(2)
    real(kind=8) :: sides(N,2,2)
    logical :: intersects, on_line
    count = 0
    code = 1
    do i=1,N
      call ray_intersects_seg(P,sides(i,1:2,1),sides(i,1:2,2),intersects,on_line)
        if (on_line) then
          return
        end if
        if (intersects) then
            count = count + 1
        end if
    end do
    if (modulo(count, 2) .ne.  0) then
        GOTO 99
    else
        code = code * -1
        GOTO 99
    end if
  99 continue
    end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine line_crossing(m1,b1,m2,b2,x_star,y_star,code)
      implicit none
      real(8) :: m1,m2,b1,b2,x_star,y_star
      integer :: code

      ! input: slopes m1, m2, and y-ints b1, b2
      ! output: intersection x_star,y_star, and integer code: 1, found intersection, 2, parallel

      if (abs(m1-m2).le.1d-8) then
        code = 2
        x_star =0.d0
        y_star =0.d0
        return
      else
        x_star = (b2-b1)/(m1-m2)
        y_star = m1*x_star + b1
        code =1
      end if
    end subroutine

     subroutine edge_in_region(edge,region,N,true)
       implicit none
       integer :: N,i
       real(8) :: edge(2,2), region(N,2,2)
       logical :: true

       do i = 1,N
         if (dist(edge(1:2,1),region(i,1:2,1)).le.1d-8 .and. dist(edge(1:2,2),region(i,1:2,2)).le.1d-8) then
           true = .true.
           exit
         else
           true = .false.

         end if
       end do
     end subroutine

     function vert_in_edge(vert,edge) result(true)
       implicit none
       real(8) :: vert(2), edge(2,2)
       logical :: true

       if (dist(vert,edge(1:2,1)).le.1d-8 .or. dist(vert,edge(1:2,2)).le.1d-8) then
         true = .true.
       else
         true = .false.
       end if
     end function

     function ind_in_inds(ind,inds) result(truea)
       implicit none
       integer :: ind(2),inds(:,:), i,N
       logical :: truea
       N = size(inds,2)
       truea=.false.
       do i=1,N
         if (ind(1)==inds(1,i) .and. ind(2)==inds(2,i)) then
           truea = .true.
         end if
       end do
     end function


    function orientation(p,q,r) result(o)
      ! input : verts p,q,r
      ! output : integer o (0 if colinear, 1 if CW, 2 if CCW)
      implicit none
      real(8) :: p(2),q(2),r(2)
      integer :: o
      ! local ::
      real(8) :: val

      val = (q(2)-p(2))*(r(1)-q(1)) - (q(1)-p(1))*(r(2)-q(2))
      if (val .gt. 0.d0) then
        o = 1
      else if (val.lt.0.d0) then
        o = 2
      else
        o = 0
      end if
    end function

    function onSegment(p,q,r) result(yes)
      ! input : verts p,q,r
      ! output:  logical yes (true if q lies on pr, false if not)
      implicit none
      real(8) :: p(2),q(2),r(2)
      logical :: yes
      if (q(1).le.max(p(1),r(1)) .and. q(1).ge.min(p(1),r(1)) .and. q(2).le.max(p(2),r(2)) .and. q(2).ge.min(p(2),r(2))) then
        yes = .true.
      else
        yes = .false.
      end if
    end function

    function doIntersect(p1,q1,p2,q2) result(yes)
      ! input : verts p1,q1,p2,q2, the same numbered ones form a segment
      ! output : logical yes, true if two line segs intersect and false if not
      implicit none
      real(8) :: p1(2),q1(2),p2(2),q2(2)
      logical :: yes
      ! local ::
      integer :: o1,o2,o3,o4

      o1 = orientation(p1,q1,p2)
      o2 = orientation(p1,q1,q2)
      o3 = orientation(p2,q2,p1)
      o4 = orientation(p2,q2,q1)

      yes = .false. ! default
      if ((o1 /= o2) .and. (o3 /= o4)) then
        yes = .true.
      end if


      if ((o1 == 0) .and. onSegment(p1,p2,q1)) then
        yes = .true.
      else if ((o2 == 0) .and. onSegment(p1, q2, q1)) then
       yes = .true.
     else if ((o3 == 0) .and. onSegment(p2, p1, q2)) then
       yes = .true.
     else if ((o4 == 0) .and. onSegment(p2, q1, q2)) then
       yes = .true.
     end if

    end function

    subroutine lineseg_intersections(p1,q1,p2,q2,yes,p_star)
      ! input: verts p1,q1,p2,q2
      ! output: logical yes if intersect , and vert p_star, the intersection
      implicit none
      real(8) :: p1(2),q1(2),p2(2),q2(2),p_star(2)
      logical :: yes
      ! local
      real(8) :: m1,b1,m2,b2,x_star,y_star
      integer :: code

      yes = doIntersect(p1,q1,p2,q2)
      if (abs(min(p1(1),q1(1))-min(p2(1),q2(1))).le.1d-8 .or. abs(max(p1(1),q1(1)) - max(p2(1),q2(1))).le.1d-8) then
        yes = .false.
      else if (abs(min(p1(2),q1(2))-min(p2(2),q2(2))).le.1d-8 .or. abs(max(p1(2),q1(2))-max(p2(2),q2(2))).le.1d-8) then
        yes = .false.
      end if
      if (yes) then
        if (abs(q1(1)-p1(1)).le.1d-8) then
          x_star = p1(1)
          y_star = (q2(2)-p2(2))/(q2(1)-p2(1)) * x_star + (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
        else if (abs(q2(1)-p2(1)).le.1d-8) then
          x_star = p2(1)
          y_star = (q1(2)-p1(2))/(q1(1)-p1(1)) * x_star + (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
        else
          m1 = (q1(2)-p1(2))/(q1(1)-p1(1))
          m2 = (q2(2)-p2(2))/(q2(1)-p2(1))
          b1 = (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
          b2 = (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
          call line_crossing(m1,b1,m2,b2,x_star,y_star,code)
        end if
      else if (.not. yes) then
        return
      end if
    end subroutine

    subroutine lineseg_intersections_endpoint(p1,q1,p2,q2,yes,p_star)
      ! input: verts p1,q1,p2,q2
      ! output: logical yes if intersect , and vert p_star, the intersection
      implicit none
      real(8) :: p1(2),q1(2),p2(2),q2(2),p_star(2)
      logical :: yes
      ! local
      real(8) :: m1,b1,m2,b2,x_star,y_star
      integer :: code

      yes = doIntersect(p1,q1,p2,q2)
      if (yes) then
        if (abs(q1(1)-p1(1)).le.1d-8) then
          x_star = p1(1)
          y_star = (q2(2)-p2(2))/(q2(1)-p2(1)) * x_star + (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
        else if (abs(q2(1)-p2(1)).le.1d-8) then
          x_star = p2(1)
          y_star = (q1(2)-p1(2))/(q1(1)-p1(1)) * x_star + (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
        else
          m1 = (q1(2)-p1(2))/(q1(1)-p1(1))
          m2 = (q2(2)-p2(2))/(q2(1)-p2(1))
          b1 = (q1(2)-(q1(2)-p1(2))*q1(1)/(q1(1)-p1(1)))
          b2 = (q2(2)-(q2(2)-p2(2))*q2(1)/(q2(1)-p2(1)))
          call line_crossing(m1,b1,m2,b2,x_star,y_star,code)
        end if
      else if (.not. yes) then
        return
      end if
    end subroutine


      function same_side_as_small_cell(up,centroid,m,b) result(code)
        implicit none
        logical :: up
        real(8) :: centroid(2),m,b
        integer :: code

        if (up) then
          if ((centroid(2) .gt. m*centroid(1) + b)) then
            code = 1
          else
            code = -1
          end if
        else
          if ((centroid(2) .lt. m*centroid(1) + b)) then
            code = 1
          else
            code = -1
          end if
        end if
      end function

      function reflect(P,m,b) result(P_prime)
        ! point P to be reflected; mx+b line of equation to be reflected across; P_prime the reflection
        implicit none
        real(8) :: P(2),m,b,P_prime(2)
        real(8) :: m_s,t,L(2),delx,dely
        ! check if point P is actually on the line, then return just P back
        if (abs(P(2)-m*P(1)-b).le.1d-7) then
          P_prime = P
          return
        end if
        ! all other cases:
        m_s = -1/m
        t= P(2) - m_s*P(1)
        L(1) = (b-t)/(m_s-m)
        L(2) = L(1)*m_s + t
        delx = (L(1)-P(1))
        dely = (L(2)-P(2))
        P_prime(1) = P(1) + 2*delx
        P_prime(2) = P(2) + 2*dely
      end function

      function colinear_check(P1,P2,P3) result(yes)
        ! check if points P1,P2,P3 are colinear, returns yes=.true. if so, and if not, yes=.false.
        implicit none
        real(8) :: P1(2),P2(2),P3(2)
        logical :: yes
        integer :: i
        real(8) :: A, x1,x2,x3,y1,y2,y3

        yes = .false.
        x1 = P1(1)
        y1 = P1(2)
        x2 = P2(1)
        y2 = P2(2)
        x3 = P3(1)
        y3 = P3(2)

        A = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
        if (abs(A).le.1d-7) then
          yes = .true.
        end if
      end function

      subroutine int_in_list(int,list,N,yes)
        implicit none
        integer :: int, N, list(N), i
        logical :: yes
        yes = .false.
        do i = 1,N
          if (int == list(i)) then
            yes = .true.
          end if
        end do
      end subroutine

      function vert_in_list(vert, list) result(yes_or_no)
        implicit none
        real(8) :: vert(2)
        real(8) :: list(:,:)
        integer :: N ,i
        logical :: yes_or_no
        N = size(list,2)
        yes_or_no = .false.
        do i = 1,N
          if (dist(vert,list(:,i)).lt.1d-7) then
            yes_or_no = .true.
            return
          end if
        end do
        return
      end function

end module aux_module_hbox
