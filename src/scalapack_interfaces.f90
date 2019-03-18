Module scalapack_interfaces
 
  Use numbers_module, Only : wp

  Interface 
 
     Pure Function numroc( n, nb, iproc, isrcproc, nprocs  ) Result( loc_n )
       Implicit None
       Integer               :: loc_n
       Integer, Intent( In ) :: n
       Integer, Intent( In ) :: nb
       Integer, Intent( In ) :: iproc
       Integer, Intent( In ) :: isrcproc
       Integer, Intent( In ) :: nprocs
     End Function numroc

     Subroutine pdgemm( ta, tb, m, n, k, alpha, A, ia, ja, descA, B, ib, jb, descB, beta, C, ic, jc, descC )
       Import :: wp
       Implicit None
       Character                   , Intent( In    ) :: ta
       Character                   , Intent( In    ) :: tb
       Integer                     , Intent( In    ) :: m
       Integer                     , Intent( In    ) :: n
       Integer                     , Intent( In    ) :: k
       Real( wp )                  , Intent( In    ) :: alpha
       Real( wp ), Dimension(  *  ), Intent( In    ) :: A
       Integer                     , Intent( In    ) :: ia
       Integer                     , Intent( In    ) :: ja
       Integer,    Dimension( 1:9 ), Intent( In    ) :: descA
       Real( wp ), Dimension(  *  ), Intent( In    ) :: B
       Integer                     , Intent( In    ) :: ib
       Integer                     , Intent( In    ) :: jb
       Integer,    Dimension( 1:9 ), Intent( In    ) :: descB
       Real( wp )                  , Intent( In    ) :: beta
       Real( wp ), Dimension(  *  ), Intent( InOut ) :: C
       Integer                     , Intent( In    ) :: ic
       Integer                     , Intent( In    ) :: jc
       Integer,    Dimension( 1:9 ), Intent( In    ) :: descC
     End Subroutine pdgemm

     Subroutine pzgemm( ta, tb, m, n, k, alpha, A, ia, ja, descA, B, ib, jb, descB, beta, C, ic, jc, descC )
       Import :: wp
       Implicit None
       Character                      , Intent( In    ) :: ta
       Character                      , Intent( In    ) :: tb
       Integer                        , Intent( In    ) :: m
       Integer                        , Intent( In    ) :: n
       Integer                        , Intent( In    ) :: k
       Complex( wp )                  , Intent( In    ) :: alpha
       Complex( wp ), Dimension(  *  ), Intent( In    ) :: A
       Integer                        , Intent( In    ) :: ia
       Integer                        , Intent( In    ) :: ja
       Integer,       Dimension( 1:9 ), Intent( In    ) :: descA
       Complex( wp ), Dimension(  *  ), Intent( In    ) :: B
       Integer                        , Intent( In    ) :: ib
       Integer                        , Intent( In    ) :: jb
       Integer,       Dimension( 1:9 ), Intent( In    ) :: descB
       Complex( wp )                  , Intent( In    ) :: beta
       Complex( wp ), Dimension(  *  ), Intent( InOut ) :: C
       Integer                        , Intent( In    ) :: ic
       Integer                        , Intent( In    ) :: jc
       Integer,       Dimension( 1:9 ), Intent( In    ) :: descC
     End Subroutine pzgemm

!!$     Subroutine pdsyevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, work, lwork, iwork, liwork, info )
!!$       Import :: wp
!!$       Implicit None
!!$       Character                    , Intent( In    ) :: jobz
!!$       Character                    , Intent( In    ) :: uplo
!!$       Integer                      , Intent( In    ) :: n
!!$       Real( wp ), Dimension( :, : ), Intent( InOut ) :: a
!!$       Integer                      , Intent( In    ) :: ia
!!$       Integer                      , Intent( In    ) :: ja
!!$       Integer   , Dimension( :    ), Intent( In    ) :: desca
!!$       Real( wp ), Dimension( :    ), Intent(   Out ) :: w
!!$       Real( wp ), Dimension( :, : ), Intent(   Out ) :: z
!!$       Integer                      , Intent( In    ) :: iz
!!$       Integer                      , Intent( In    ) :: jz
!!$       Integer   , Dimension( :    ), Intent( In    ) :: descz
!!$       Real( wp ), Dimension( :    ), Intent(   Out ) :: work
!!$       Integer                      , Intent( In    ) :: lwork
!!$       Integer   , Dimension( :    ), Intent(   Out ) :: iwork
!!$       Integer   ,                    Intent( In    ) :: liwork
!!$       Integer   ,                    Intent( In    ) :: info
!!$     End Subroutine pdsyevd

!!$     Subroutine pzheevd( jobz, uplo, n, a, ia, ja, desca, w, z, iz, jz, descz, cwork, lcwork, rwork, lrwork, iwork, liwork, info )
!!$       Import :: wp
!!$       Implicit None
!!$       Character                       , Intent( In    ) :: jobz
!!$       Character                       , Intent( In    ) :: uplo
!!$       Integer                         , Intent( In    ) :: n
!!$       Complex( wp ), Dimension( :, : ), Intent( InOut ) :: a
!!$       Integer                         , Intent( In    ) :: ia
!!$       Integer                         , Intent( In    ) :: ja
!!$       Integer      , Dimension( :    ), Intent( In    ) :: desca
!!$       Real   ( wp ), Dimension( :    ), Intent(   Out ) :: w
!!$       Complex( wp ), Dimension( :, : ), Intent(   Out ) :: z
!!$       Integer                         , Intent( In    ) :: iz
!!$       Integer                         , Intent( In    ) :: jz
!!$       Integer      , Dimension( :    ), Intent( In    ) :: descz
!!$       Complex( wp ), Dimension( :    ), Intent(   Out ) :: cwork
!!$       Integer                         , Intent( In    ) :: lcwork
!!$       Real( wp )   , Dimension( :    ), Intent(   Out ) :: rwork
!!$       Integer                         , Intent( In    ) :: lrwork
!!$       Integer      , Dimension( :    ), Intent(   Out ) :: iwork
!!$       Integer      ,                    Intent( In    ) :: liwork
!!$       Integer      ,                    Intent( In    ) :: info
!!$     End Subroutine pzheevd
 
  End Interface

End Module scalapack_interfaces
