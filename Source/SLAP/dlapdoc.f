*DECK DLPDOC
      SUBROUTINE DLPDOC
C***BEGIN PROLOGUE  DLPDOC
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2B, Z
C***KEYWORDS  LIBRARY=SLATEC(DLAP),
C             TYPE=DOUBLE PRECISION(DLPDOC-A),
C             Sparse Iterative Methods, Iterative Improvement,
C             Preconditioned Conjugate Gradient, Normal Equations,
C             BiConjugate Gradient, BiConjugate Gradient Squared, 
C             Orthomin, Generalize Minimum Residual, Documentation
C***AUTHOR  Seager, Mark. K., (LLNL)
C             User Systems Division
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550
C             (FTS) 543-3141, (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Sparse Linear Algebra Package Version 2.0
C            Routines to  solve large sparse symmetric and nonsymmetric 
C            positive definite linear systems, Ax = b, using precondit-
C            ioned iterative methods.
C
C            This package was originally derived from a set of iterative
C            routines written by Anne Greenbaum, "Routines for Solving
C            Large Sparse Linear Systems", Lawrence Livermore Nat.
C            Laboratory, Livermore Computing Center, January 1986 
C            Tentacle, pp 15-21.
C***DESCRIPTION
C                                 The
C                    Sparse Linear Algebra Package
C                      Double Precision Routines
C
C                @@@@@@@  @            @@@    @@@@@@@@
C               @       @ @           @   @   @       @
C               @         @          @     @  @       @
C                @@@@@@@  @         @       @ @@@@@@@@
C                       @ @         @@@@@@@@@ @
C               @       @ @         @       @ @        
C                @@@@@@@  @@@@@@@@@ @       @ @        
C                                        
C      @       @                            @@@@@@@        @@@@@
C      @       @                           @       @      @    @@ 
C      @       @  @@@@@@@  @ @@                    @     @    @  @
C      @       @ @       @ @@  @             @@@@@@      @   @   @
C       @     @  @@@@@@@@@ @                @            @  @    @
C        @   @   @         @               @         @@@  @@    @ 
C         @@@     @@@@@@@  @               @@@@@@@@@ @@@   @@@@@  
C
C     
C    =================================================================
C    ========================== Introduction =========================
C    =================================================================
C    This document  contains the specifications for  the  SLAP Version
C    2.0 package, a Fortran 77  package  for  the  solution  of  large
C    sparse   linear systems, Ax  =  b,  via  preconditioned iterative
C    methods.   Included in  this  package are "core"  routines  to do
C    Iterative   Refinement  (Jacobi's  method),  Conjugate  Gradient,
C    Conjugate Gradient on the normal equations, AA'y = b,  (where x =
C    A'y and  A' denotes the  transpose of   A), BiConjugate Gradient,
C    BiConjugate  Gradient  Squared, Orthomin and  Generalized Minimum
C    Residual Iteration.    These "core" routines   do  not  require a
C    "fixed"   data  structure   for storing  the   matrix  A  and the
C    preconditioning   matrix  M.   The  user  is free  to  choose any
C    structure that facilitates  efficient solution  of the problem at
C    hand.  The drawback  to this approach  is that the user must also
C    supply at least two routines  (MATVEC and MSOLVE,  say).   MATVEC
C    must calculate, y = Ax, given x and the user's data structure for
C    A.  MSOLVE must solve,  r = Mz, for z (*NOT*  r) given r  and the
C    user's data  structure for  M (or its  inverse).  The user should
C    choose   M  so that MA   is  approximately the   identity and the
C    solution step r = Mz is "easy" to  solve.  For some of the "core"
C    routines (Orthomin,  BiConjugate Gradient and  Conjugate Gradient
C    on the  normal equations)   the user must  also  supply  a matrix
C    transpose times   vector  routine  (MTTVEC,  say)  and (possibly,
C    depending    on the "core"  method)   a  routine  that solves the
C    transpose  of   the   preconditioning    step     (MTSOLV,  say).
C    Specifically, MTTVEC is a routine which calculates y = A'x, given
C    x and the user's data structure for A (A' is the transpose of A).
C    MTSOLV is a routine which solves the system r = M'z for z given r
C    and the user's data structure for M.
C    
C    This process of writing the matrix vector operations  can be time
C    consuming and error  prone.  To alleviate  these problems we have
C    written drivers   for  the  "core" methods  that  assume the user
C    supplies one of two specific data structures (SLAP Triad and SLAP
C    Column format), see  below.  Utilizing these  data structures  we
C    have augmented   each  "core" method  with   two preconditioners:
C    Diagonal  Scaling and Incomplete Factorization.  Diagonal scaling
C    is easy to implement, vectorizes very  well and for problems that
C    are  not too  ill-conditioned  reduces the  number  of iterations
C    enough   to warrant its use.  On   the other  hand, an Incomplete
C    factorization  (Incomplete  Cholesky for  symmetric systems   and
C    Incomplete LU for nonsymmetric  systems) may  take much longer to
C    calculate, but it reduces the iteration count (for most problems)
C    significantly.  Our implementations  of IC and ILU  vectorize for
C    machines with hardware gather scatter, but the vector lengths can
C    be quite short if  the  number  of nonzeros   in a column is  not
C    large.
C    
C    =================================================================
C    ==================== Supplied Data Structures ===================
C    =================================================================
C    The following describes the data   structures supplied  with  the
C    package: SLAP Triad and Column formats.
C    
C    ====================== S L A P Triad format =====================
C    In the SLAP Triad format only the non-zeros are stored.  They may
C    appear in *ANY* order.  The user supplies three  arrays of length
C    NELT, where NELT  is the   number of  non-zeros  in the   matrix:
C    (IA(NELT),  JA(NELT), A(NELT)).  If  the matrix is symmetric then
C    one need only store the lower triangle (including  the  diagonal)
C    and NELT would be the corresponding  number  of non-zeros stored.
C    For each non-zero the user puts the row and column  index of that
C    matrix  element   in the  IA  and JA  arrays.  The  value  of the
C    non-zero matrix element is placed  in  the corresponding location
C    of  the A array.   This  is an extremely  easy  data structure to
C    generate.  On the other hand, it is not very  efficient on vector
C    computers for the iterative  solution of  linear systems.  Hence,
C    SLAP changes this input data structure to  the SLAP Column format
C    for the iteration (but does not change it back).
C    
C    Here  is an example   of  the  SLAP  Triad storage  format  for a
C    nonsymmetric 5x5 Matrix.  NELT=11.   Recall that the  entries may
C    appear in any order.
C     
C     5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                           1  2  3  4  5  6  7  8  9 10 11
C    |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C    |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C    | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C    | 0  0  0 44  0|
C    |51  0 53  0 55|
C    
C       =================== S L A P Column format ==================
C       This routine requires  that the  matrix  A be  stored in the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C       
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C       1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C    =================================================================
C    ======================= Naming Conventions ======================
C    =================================================================
C    SLAP  iterative  methods,    matrix vector    and  preconditioner
C    calculation  routines   follow a naming   convention  which, when
C    understood, allows one to determine the iterative method and data
C    structure(s) used.  The  subroutine  naming convention  takes the
C    following form:
C                          P[S][M]DESC
C    where P stands  for the precision (or  data type)  of the routine
C    and is required in   all names, S denotes   whether  or  not  the
C    routine requires the SLAP Triad or Column format (it  does if the
C    second letter  of the name   is S and   does not otherwise),  the
C    optional  M   stands for the   type of preconditioner used  (only
C    appears in drivers for  "core" routines) and  DESC is some number
C    of letters describing the method  or purpose of the  routine.  In
C    this  incarnation of SLAP  only single  precision  data types are
C    supported (no double precision or complex data type routines have
C    been written).   Hence,  all routines start with an   S,  boring.
C    The  brackets around  S  and M  designate  that these  fields are
C    optional.
C    
C    The  following is a  list  of  the  "DESC"  fields for  iterative
C    methods  and   their meaning: BCG:   BiConjugate   Gradient;  CG:
C    Conjugate  Gradient; CGN:   Conjugate  Gradient   on   the Normal
C    equations; CGS, CS: biConjugate Gradient Squared; GMRES, GMR, GM:
C    Generalized Minimum  RESidual; IR, R:  Iterative Refinement; JAC:
C    JACobi's method; GS: Gauss-Seidel; OMN, OM: Orthomin.
C    
C    Here are some examples of the routines:
C    1) DBCG: Double  precision BiConjugate Gradient  "core"  routine.
C    On can deduce that this is a "core" routine,  because the S and M
C    fields are   missing  and BiConjugate  Gradient  is an  iterative
C    method. 2) DSDBCG: Dingle precision, SLAP data structure BCG with
C    Diagonal scaling.  3)   DSLUBC: Double    precision,  BCG    with
C    incomplete  LU factorization  as   the preconditioning.  4)  DCG:
C    Double precision Conjugate  Gradient  "core" routine.  5)  DSDCG:
C    Double precision,  SLAP data  structure  Conjugate  Gradient with
C    Diagonal scaling.  6)  DSICCG:  Double   precision,    SLAP  data
C    structure   Conjugate    Gradient  with    Incomplete    Cholesky
C    factorization preconditioning.
C    
C    
C    =================================================================
C    ===================== USER CALLABLE ROUTINES ====================
C    =================================================================
C    The following is a list of  the "user callable" SLAP routines and
C    their one  line description.  The headers denote  the  file names
C    where the routines can be found, as distributed for Unix systems.
C    
C    ============================= DBCG.F ============================
C    DBCG: Preconditioned BiConjugate Gradient Sparse Ax=b solver.
C    DSDBCG: Diagonally Scaled BiConjugate Gradient Sparse Ax=b solver.
C    DSLUBC: Incomplete LU BiConjugate Gradient Sparse Ax=b solver.
C    
C    ============================= DCG.F =============================
C    DCG: Preconditioned Conjugate Gradient iterative Ax=b solver.
C    DSDCG: Diagonally Scaled Conjugate Gradient Sparse Ax=b Solver.
C    DSICCG: Incomplete Cholesky Conjugate Gradient Sparse Ax=b Solver.
C    
C    ============================= DCGN.F ============================
C    DCGN: Preconditioned CG Sparse Ax=b Solver for Normal Equations.
C    DSDCGN: Diagonally Scaled CG Sparse Ax=b Solver for Normal Eqn's.
C    DSLUCN: Incomplete LU CG Sparse Ax=b Solver for Normal Equations.
C    
C    ============================= DCGS.F ============================
C    DCGS: Preconditioned BiConjugate Gradient Sparse Ax=b solver.
C    DSDCGS: Diagonally Scaled CGS Sparse Ax=b Solver.
C    DSLUCS: Incomplete LU BiConjugate Gradient Sparse Ax=b solver.
C    
C    ============================ DGMRES.F ===========================
C    DGMRES: Preconditioned GMRES iterative sparse Ax=b solver.
C    DSDGMR: Diagonally scaled GMRES iterative sparse Ax=b solver.
C    DSLUGM: Incomplete LU GMRES iterative sparse Ax=b solver.
C    
C    ============================= DIR.F =============================
C    DIR: Preconditioned Iterative Refinement sparse Ax = b solver.
C    DSJAC: Jacobi's method iterative sparse Ax = b solver.
C    DSGS: Gauss-Seidel method iterative sparse Ax = b solver.
C    DSILUR: Incomplete LU Iterative Refinement sparse Ax = b solver.
C    
C    ============================ DOMN.F =============================
C    DOMN: Preconditioned Orthomin Sparse Iterative Ax=b Solver.
C    DSDOMN: Diagonally Scaled Orthomin Sparse Iterative Ax=b Solver.
C    DSLUOM: Incomplete LU Orthomin Sparse Iterative Ax=b Solver.
C    
C    ============================ DMSET.F ============================
C    DSDS: Diagonal Scaling Preconditioner SLAP Set Up.
C    DSDSCL: Diagonally Scales/Unscales a SLAP Column Matrix.
C    DSD2S: Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up.
C    DS2LT: Lower Triangle Preconditioner SLAP Set Up.
C    DSICS: Incomplete Cholesky Decomp. Preconditioner SLAP Set Up.
C    DSILUS: Incomplete LU Decomposition Preconditioner SLAP Set Up.
C    
C    =========================== DMVOPS.F ============================
C       Most of the incomplete  factorization  (LL' and LDU) solvers
C       in this  file require an  intermediate routine  to translate
C       from the SLAP MSOLVE(N, R, Z, NELT, IA,  JA, A, ISYM, RWORK,
C       IWORK) calling  convention to the calling  sequence required
C       by  the solve routine.   This generally  is  accomplished by
C       fishing out pointers to the preconditioner (stored in RWORK)
C       from the  IWORK  array and then making a call to the routine
C       that actually does the backsolve.
C
C    DSMV: SLAP Column Format Sparse Matrix Vector Product.
C    DSMTV: SLAP Column Format Sparse Matrix (transpose) Vector Prod.
C    DSDI: Diagonal Matrix Vector Multiply.
C    DSLI: SLAP MSOLVE for Lower Triangle Matrix (set up for SSLI2).
C    DSLI2: Lower Triangle Matrix Backsolve.
C    DSLLTI: SLAP MSOLVE for LDL' (IC) Fact. (set up for SLLTI2).
C    DLLTI2: Back solve routine for LDL' Factorization.
C    DSLUI: SLAP MSOLVE for LDU Factorization (set up for SSLUI2).
C    DSLUI2: SLAP Back solve for LDU Factorization.
C    DSLUTI: SLAP MTSOLV for LDU Factorization (set up for SSLUI4).
C    DSLUI4: SLAP back solve for LDU Factorization.
C    DSMMTI: SLAP MSOLVE for LDU Fact of Normal Eq (set up for SSMMI2).
C    DSMMI2: SLAP Back solve for LDU Factorization of Normal Equations.
C    
C    ========================== SLAPUTIL.F ==========================
C    DCHKW: SLAP WORK/IWORK Array Bounds Checker.
C    QS2I1D: Quick Sort Integer array coupled with int and double 
C            precision arrays.
C    DS2Y: SLAP Triad to SLAP Column Format Converter.
C    DCPPLT: Printer Plot of SLAP Column Format Matrix.
C    DTOUT: Write out SLAP Triad Format Linear System.
C    DTIN: Read in SLAP Triad Format Linear System.
C    
C    =================================================================
C    ====================== Which Method To Use ======================
C    =================================================================
C                          BACKGROUND  
C    In solving a large sparse linear system Ax = b using an iterative
C    method, it   is  not necessary to actually   store  the matrix A.
C    Rather, what is needed is a procedure  for multiplying the matrix
C    A times a given vector y to obtain the matrix-vector product, Ay.
C    SLAP has been written to take advantage of this fact.  The higher
C    level routines in the package require storage only of the nonzero
C    elements of   A (and  their  positions), and  even this   can  be
C    avoided, if the  user  writes his own subroutine for  multiplying
C    the matrix times a vector  and   calls the lower-level  iterative
C    routines in the package.
C    
C    If  the matrix A is ill-conditioned,  then most iterative methods
C    will be slow to converge (if they converge  at all!).  To improve
C    the  convergence  rate,  one  may use  a "matrix  splitting," or,
C    "preconditioning matrix," say, M.  It is then necessary to solve,
C    at each iteration, a linear system  with coefficient matrix M.  A
C    good preconditioner  M should have  two  properties: (1) M should
C    "approximate" A, in the sense that the  matrix inv(M)*A  (or some
C    variant  thereof) is better conditioned  than the original matrix
C    A; and  (2) linear  systems with coefficient  matrix M should  be
C    much easier  to solve  than  the original system with coefficient
C    matrix   A.   Preconditioning routines  in the   SLAP package are
C    separate from the  iterative   routines,  so   that any of    the
C    preconditioners provided in the package,   or one that the   user
C    codes himself, can be used with any of the iterative routines.
C    
C                        CHOICE OF PRECONDITIONER  
C    If you  willing   to live with   either the SLAP Triad or  Column
C    matrix data structure  you  can then  choose one  of two types of
C    preconditioners   to   use:   diagonal  scaling    or  incomplete
C    factorization.  To  choose   between these two   methods requires
C    knowing  something  about the computer you're going  to run these
C    codes on  and how well incomplete factorization  approximates the
C    inverse of your matrix.
C    
C    Let's   suppose you have   a scalar  machine.   Then,  unless the
C    incomplete factorization is very,  very poor this  is *GENERALLY*
C    the method to choose.  It  will reduce the  number of  iterations
C    significantly and is not all  that expensive  to compute.  So  if
C    you have just one  linear system to solve  and  "just want to get
C    the job  done" then try  incomplete factorization first.   If you
C    are thinking of integrating some SLAP  iterative method into your
C    favorite   "production  code" then  try incomplete  factorization
C    first,  but  also check  to see that  diagonal  scaling is indeed
C    slower for a large sample of test problems.
C    
C    Let's now  suppose  you have  a  vector  computer  with  hardware
C    gather/scatter support (Cray X-MP, Y-MP, SCS-40 or Cyber 205, ETA
C    10, ETA Piper  or Convex C-1,  etc.).   Then  it's much harder to
C    choose  between the  two  methods.   The  versions  of incomplete
C    factorization in SLAP do in fact vectorize, but have short vector
C    lengths and the factorization step is relatively  more expensive.
C    Hence,  for  most problems (i.e.,  unless  your  problem  is  ill
C    conditioned,  sic!)  diagonal  scaling is  faster,  with its very
C    fast    set up  time    and  vectorized  (with   long    vectors)
C    preconditioning step (even though  it  may take more iterations).
C    If you have several systems (or  right hand sides) to  solve that
C    can  utilize  the  same  preconditioner  then the   cost   of the
C    incomplete factorization can   be  amortized over these  several
C    solutions.  This situation gives more advantage to the incomplete
C    factorization methods.  If  you have  a  vector  machine  without
C    hardware  gather/scatter (Cray  1,  Cray  2  &  Cray 3) then  the
C    advantages for incomplete factorization are even less.
C    
C    If you're trying to shoehorn SLAP into your  favorite "production
C    code" and can not easily generate either the SLAP Triad or Column
C    format  then  you are  left  to   your  own  devices in terms  of
C    preconditioning.  Also,  you may  find that the   preconditioners
C    supplied with SLAP are not sufficient  for your problem.  In this
C    situation we would  recommend  that you   talk  with a  numerical
C    analyst  versed in   iterative   methods   about   writing  other
C    preconditioning  subroutines (e.g.,  polynomial  preconditioning,
C    shifted incomplete factorization,  SOR  or SSOR  iteration).  You
C    can always "roll your own"  by using the "core" iterative methods
C    and supplying your own MSOLVE and MATVEC (and possibly MTSOLV and
C    MTTVEC) routines.  If you do develop a new preconditioner for the
C    SLAP data structure send the code to us  (if you can do that with
C    no strings attached!, i.e.  copyright restrictions) and we'll add
C    it to the package!
C    
C                          SYMMETRIC SYSTEMS
C    If your matrix is symmetric then you would want to use one of the
C    symmetric system  solvers.    If  your  system  is  also positive
C    definite,   (Ax,x) (Ax dot  product  with x) is  positive for all
C    non-zero  vectors x,  then use   Conjugate Gradient (SCG,  SSDCG,
C    SSICSG).  If you're  not sure it's SPD   (symmetric and  Positive
C    Definite)  then try SCG anyway and  if it works, fine.  If you're
C    sure your matrix is not  positive definite  then you  may want to
C    try the iterative refinement   methods  (SIR)  or the  GMRES code
C    (SGMRES) if SIR converges too slowly.
C    
C                         NONSYMMETRIC SYSTEMS
C    This   is currently  an  area  of  active research  in  numerical
C    analysis  and   there   are   new  strategies  being   developed.
C    Consequently take the following advice with a grain of salt.   If
C    you matrix is positive definite, (Ax,x)  (Ax  dot product  with x
C    is positive for all non-zero  vectors x), then you can use any of
C    the    methods   for   nonsymmetric   systems (Orthomin,   GMRES,
C    BiConjugate Gradient, BiConjugate Gradient  Squared and Conjugate
C    Gradient applied to the normal equations).  If your system is not
C    too ill conditioned then try  BiConjugate Gradient Squared (BCGS)
C    or GMRES (SGMRES).  Both  of  these methods converge very quickly
C    and do  not require A'  or M' ('  denotes transpose) information.
C    SGMRES  does require  some  additional storage,  though.  If  the
C    system is very  ill conditioned  or   nearly positive  indefinite
C    ((Ax,x) is positive,  but may be  very small),  then GMRES should
C    be the first choice,  but try the  other  methods  if you have to
C    fine tune  the solution process for a  "production code".  If you
C    have a great preconditioner for the normal  equations (i.e., M is
C    an approximation to the inverse of AA' rather than  just  A) then
C    this is not a bad route to travel.  Old wisdom would say that the
C    normal equations are a disaster  (since it squares the  condition
C    number of the system and SCG convergence is linked to this number
C    of    infamy), but   some     preconditioners    (like incomplete
C    factorization) can reduce the condition number back below that of
C    the original system.
C***REFERENCES  1. M. Seager, ``A SLAP for the Masses,'' Lawrence 
C                 Livermore Nat. Laboratory Technical Report, 
C                 UCRL-100267, December 1988.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DLPDOC
C
C     This is a *DUMMY* subroutine and should be called only for
C     the futility of it all.  Aren't *SLATEC* conventions FUN!!!
C
C***FIRST EXECUTABLE STATEMENT  DLPDOC
      RETURN
C------------- LAST LINE OF DLPDOC FOLLOWS -----------------------------
      END
