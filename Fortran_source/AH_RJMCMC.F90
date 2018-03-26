PROGRAM Main_Age_Uncertain_RJMCMC
!
!  This program uses a RJMCMC algorithm to sample the model: comprising unknown ages and the data fit.
! It then marginalises to find the posterior of all variables of interest.

! Phil Livermore / Alex Fournier, Apr 2016
! 
! Jan 2017: adapted to conform to same inputfile as Python code
!

USE AGE_HYPERPARAMETER_RJMCMC


IMPLICIT NONE
INTEGER, PARAMETER :: MAX_DATA  = 1000

CHARACTER(500) :: ARG, Data_file_name, Header, WRITE_MODEL_FILE_NAME
CHARACTER(500) :: inputline, junk, outputs_directory
INTEGER ::  NARG
CHARACTER(1) :: AGE_DISTRIBUTION

INTEGER :: I, K, BURN_IN, NSAMPLE, K_INIT, K_MAX, K_MIN, K_MAX_ARRAYBOUND, discretise_size, show, thin, num,j, &
           NUM_DATA, IOS, K_MAX_ARRAY_BOUND, s, ind, NBINS, I_MODEL, FREQ_WRITE_MODELS, FREQ_WRITE_JOINT_DISTRIB, SAMPLE_INDEX_JOINT_DISTRIBUTION
REAL( KIND = 8) :: I_MAX, I_MIN, sigma_move, sigma_change_value, sigma_birth, int_j,&
                   X_MIN, X_MAX, age_frac

CHARACTER(100), Allocatable :: LINE_READ(:)
INTEGER, ALLOCATABLE :: SEED(:)

INTEGER :: stratification(1:MAX_DATA)
INTEGER ::  age_col, d_age_col, F_col, dF_col, Strat_col

REAL( KIND = 8) :: age(1:MAX_DATA), delta_age(1:MAX_DATA),  intensity(1:MAX_DATA), delta_intensity(1:MAX_DATA)
REAL( KIND = 8), ALLOCATABLE :: X(:)

LOGICAL :: CALC_CREDIBLE
Real ( KIND = 8) :: credible

TYPE (RETURN_INFO_STRUCTURE) :: RETURN_INFO

CALC_CREDIBLE = .TRUE.

stratification(:) = 0 !default to no stratification.


! Load parameter file
NARG = iargc()
      IF(NARG .gt. 0) THEN
      CALL getarg(1,ARG)
      ENDIF


      IF(ARG(1:1) .eq. '-' .OR. ARG(1:4) .eq. 'help' .OR. &
      ARG(1:4) .eq. 'HELP' .OR.  ARG(1:1) .eq. '?' .OR. &
      NARG .eq. 0 ) THEN
         PRINT*, 'ERROR: FIRST ARGUMENT MUST BE INPUT FILE'
         STOP
      ENDIF


      OPEN(30, FILE = TRIM(ARG), STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT = IOS)

       IF( IOS  .NE. 0) THEN
       PRINT*, 'CANNOT OPEN INPUT FILE'
        STOP
        ENDIF

DO
READ(30,'(A)',END=101) inputline
IF ( to_upper(inputline(1:len('Data_file'))) == to_upper('Data_file') ) read(inputline(len('Data_file')+2:),'(A)') Data_file_name

IF ( to_upper(inputline(1:len('File_format'))) == to_upper('File_format') ) read(inputline(len('File_format')+2:),*) age_col, d_age_col, F_col, dF_col, Strat_col
IF ( to_upper(inputline(1:len('Intensity_prior'))) == to_upper('Intensity_prior') ) read(inputline(len('Intensity_prior')+2:),*) I_min, I_max
IF ( to_upper(inputline(1:len('Burn_in'))) == to_upper('Burn_in') ) read(inputline(len('Burn_in')+2:),*) burn_in
IF ( to_upper(inputline(1:len('Nsamples'))) == to_upper('Nsamples') ) read(inputline(len('Burn_in')+2:),*) NSAMPLE
IF ( to_upper(inputline(1:len('model_discretisation'))) == to_upper('model_discretisation') ) read(inputline(len('model_discretisation')+2:),*) discretise_size
IF ( to_upper(inputline(1:len('Chain_parameters'))) == to_upper('Chain_parameters') ) read(inputline(len('Chain_parameters')+2:),*) show, thin
IF ( to_upper(inputline(1:len('Age_distribution'))) == to_upper('Age_distribution') ) read(inputline(len('Age_distribution')+2:),*) Age_distribution
IF ( to_upper(inputline(1:len('running_mode'))) == to_upper('running_mode') ) read(inputline(len('running_mode')+2:),*) RUNNING_MODE
IF ( to_upper(inputline(1:len('Age_bounds'))) == to_upper('Age_bounds') ) read(inputline(len('Age_bounds')+2:),*) X_min, X_max
IF ( to_upper(inputline(1:len('Sigmas'))) == to_upper('Sigmas') ) read(inputline(len('Sigmas')+2:),*) sigma_move, sigma_change_value, sigma_birth
IF ( to_upper(inputline(1:len('Age_frac'))) == to_upper('Age_frac') ) read(inputline(len('Age_frac')+2:),*) age_frac
IF ( to_upper(inputline(1:len('Num_change_points'))) == to_upper('Num_change_points') ) read(inputline(len('Num_change_points')+2:),*) K_MIN, K_MAX
IF ( to_upper(inputline(1:len('Nbins'))) == to_upper('Nbins') ) read(inputline(len('Nbins')+2:),*) NBINS
IF ( to_upper(inputline(1:len('output_model'))) == to_upper('output_model') ) read(inputline(len('output_model')+2:),*) WRITE_MODEL_FILE_NAME, FREQ_WRITE_MODELS
IF ( to_upper(inputline(1:len('output_joint_distribution_freq'))) == to_upper('output_joint_distribution_freq') ) read(inputline(len('output_joint_distribution_freq')+2:),*) FREQ_WRITE_JOINT_DISTRIB
IF ( to_upper(inputline(1:len('Credible'))) == to_upper('Credible') ) read(inputline(len('Credible')+2:),*) Credible
IF ( to_upper(inputline(1:len('Outputs_directory'))) == to_upper('Outputs_directory') ) read(inputline(len('Outputs_directory')+2:),'(A)') outputs_directory

ENDDO
101     CONTINUE

CLOSE( 30 )


IF( to_upper(AGE_DISTRIBUTION(1:1)) .NE. 'U' .AND. to_upper(AGE_DISTRIBUTION(1:1)) .NE. 'N') THEN
PRINT*, 'INVALID AGE DISTRIBUTION: N or U'
PRINT*, 'CODE READ ', to_upper(AGE_DISTRIBUTION(1:1)), ' FROM INPUT FILE'
STOP
ENDIF

PRINT*, 'DATA file name: ', TRIM(Data_file_name)
OPEN(30, FILE = Data_file_name, STATUS = 'OLD', FORM = 'FORMATTED', &
                             IOSTAT = IOS, ACTION = 'READ')
    IF( IOS .NE. 0) THEN
    PRINT*, 'ERROR IN OPENING FILE ', TRIM(Data_file_name)
    STOP
    ENDIF

    ALLOCATE( LINE_READ(1: max(age_col, d_age_col, F_col, dF_col, Strat_col)+1) )

    I=1
    DO
    READ(30,'(A)', END = 1001) inputline
    IF( inputline(1:1) == '#') CYCLE

    READ(inputline,*) LINE_READ(:)
    READ(LINE_READ(age_col+1),*) age(i)
    READ(LINE_READ(d_age_col+1),*) delta_age(i)
    READ(LINE_READ(F_col+1),*) intensity(i)
    READ(LINE_READ(dF_col+1),*) delta_intensity(i)


    IF(Strat_col > 0) THEN
    READ(LINE_READ(Strat_col+1),*) stratification(i)

    ENDIF

    i = i + 1
    IF( i > MAX_DATA ) THEN
     PRINT*, 'NUMBER OF DATA EXCEEDS HARDWIRED MAX', MAX_DATA
     STOP
    ENDIF

    ENDDO
1001    CONTINUE
    CLOSE(30)
    NUM_DATA = i-1



k_max_array_bound = k_max + 1
NUM=ceiling((nsample-burn_in)*(100-credible)/200.0_8/thin) ! number of collected samples for credible intervals

ALLOCATE( X(1: discretise_size) )
ALLOCATE( RETURN_INFO%AV(1:discretise_size), RETURN_INFO%BEST(1:discretise_size),  &
          RETURN_INFO%INF(1:discretise_size), RETURN_INFO%sup(1:discretise_size), &
          RETURN_INFO%change_points(nsample * k_max_array_bound), RETURN_INFO%n_changepoint_hist(k_max_array_bound), &
          RETURN_INFO%convergence(nsample * k_max_array_bound), &
          RETURN_INFO%MARGINAL_AGES(1:NBINS,1:NUM_DATA)    )

ALLOCATE( RETURN_INFO%MEDIAN(1:discretise_size), RETURN_INFO%MODE(1:discretise_size), &
          RETURN_INFO%MARGINAL_DENSITY_INTENSITY(1:discretise_size,1:NBINS) )


      PRINT*, 'READ IN ', NUM_DATA, 'DATA ITEMS'
      PRINT*, 'DATA AGE RANGE IS ', MINVAL(AGE(1:NUM_DATA)), ' TO ', MAXVAL( AGE(1:NUM_DATA))


CALL RANDOM_SEED(SIZE = i)
ALLOCATE( SEED(1:i) )
SEED(:) = 1
CALL RANDOM_SEED(PUT = SEED)


IF( X_MIN .eq. X_MAX) THEN
X_MIN = MINVAL( AGE(1:NUM_DATA) ) - X_MIN
X_MAX = MAXVAL( AGE(1:NUM_DATA) ) + X_MAX
ENDIF

! Check to see that X_MIN and X_MAX enclose the data set:
IF( to_upper(AGE_DISTRIBUTION(1:1)) == 'U') THEN
IF( X_MIN .GT. MINVAL(AGE(1:NUM_DATA)-delta_age(1:NUM_DATA)) .OR. X_MAX .LT. MAXVAL(AGE(1:NUM_DATA)+delta_age(1:NUM_DATA)) ) THEN
PRINT*, 'INCREASE RANGE OF END POINT AGES AS RANGE DOES NOT '
PRINT*, 'SPAN ALL DATA POINTS (WITH UNIFORM ERRORS)'
PRINT*, 'RANGE NEEDS TO INCLUDE ', MINVAL(AGE(1:NUM_DATA)-delta_age(1:NUM_DATA)),' : ',MAXVAL(AGE(1:NUM_DATA)+delta_age(1:NUM_DATA))
STOP
ENDIF
ELSE !assumed normal distribution
IF( X_MIN .GT. MINVAL(AGE(1:NUM_DATA)-2.0_8 * delta_age(1:NUM_DATA)) .OR. X_MAX .LT. MAXVAL(AGE(1:NUM_DATA)+2.0_8 * delta_age(1:NUM_DATA)) ) THEN
PRINT*, 'INCREASE RANGE OF END POINT AGES AS RANGE DOES NOT '
PRINT*, 'SPAN ALL DATA POINTS (ASSUMING A RANGE OF 2 X STANDARD DEVIATION)'
PRINT*, 'RANGE NEEDS TO INCLUDE ', MINVAL(AGE(1:NUM_DATA)-2*delta_age(1:NUM_DATA)),' : ',MAXVAL(AGE(1:NUM_DATA)+2*delta_age(1:NUM_DATA))
STOP
ENDIF

ENDIF

PRINT*, 'MODEL AGE RANGE IS ', X_MIN, ' : ', X_MAX


    DO I=1, discretise_size
    X(I) = X_MIN + REAL(I-1, KIND = 8)/REAL(discretise_size-1, KIND = 8) * (X_MAX - X_MIN)
    ENDDO

IF( MAXVAL( stratification(1:NUM_DATA)) .eq. 0) THEN
PRINT*, 'NO STRATIFICATION CONSTRAINTS'
ELSE
PRINT*, 'STRATIFICATION CONSTRAINTS BEING USED'
ENDIF

CALC_CREDIBLE = .TRUE.


call system('mkdir -p '//TRIM(Outputs_directory))

! copy input file
call system('cp '//TRIM(ARG)//' '//TRIM(Outputs_directory)//'/input_file')
CALL RJMCMC(burn_in, NUM_DATA, age(1:NUM_DATA), delta_age(1:NUM_DATA), intensity(1:NUM_DATA), delta_intensity(1:NUM_DATA), stratification(1: NUM_DATA), NSAMPLE, I_MIN, I_MAX, X_MIN, X_MAX, K_MIN, K_MAX, SIGMA_MOVE, sigma_change_value, sigma_birth, age_frac, discretise_size, SHOW, THIN, NBINS, RETURN_INFO, CALC_CREDIBLE, FREQ_WRITE_MODELS, WRITE_MODEL_FILE_NAME, FREQ_WRITE_JOINT_DISTRIB,  AGE_DISTRIBUTION, credible, Outputs_directory)



OPEN(24, FILE = TRIM(Outputs_directory)//'/data.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, NUM_DATA
WRITE(24,'(4(F10.3,X),X,I)') AGE(i), delta_age(i), intensity(I), delta_intensity(I), stratification(i)
enddo
CLOSE(24)

OPEN(23, FILE = TRIM(Outputs_directory)//'/k_histogram.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, k_max
WRITE(23,*) i,RETURN_INFO%n_changepoint_hist(i)
enddo
CLOSE(23)

OPEN(24, FILE = TRIM(Outputs_directory)//'/average.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i),RETURN_INFO%av(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/best_fit.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i),RETURN_INFO%best(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/mode.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i),RETURN_INFO%mode(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/median.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i),RETURN_INFO%median(i)
enddo
CLOSE(24)



OPEN(24, FILE = TRIM(Outputs_directory)//'/misfit.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, nsample, 100
WRITE(24,'(i,X,F10.3)') i, RETURN_INFO%convergence(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/credible_upper.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i), RETURN_INFO%sup(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/credible_lower.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i), RETURN_INFO%inf(i)
enddo
CLOSE(24)

OPEN(23, FILE = TRIM(Outputs_directory)//'/changepoints.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, RETURN_INFO%MAX_NUMBER_CHANGE_POINTS_HISTORY
WRITE(23,'(F10.3)') RETURN_INFO%change_points(i)
enddo
CLOSE(23)

OPEN(23, FILE = TRIM(Outputs_directory)//'/age_marginals.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
WRITE(23,*) NUM_DATA, NBINS
DO i=1, NUM_DATA
do j=1, NBINS
IF( AGE_DISTRIBUTION(1:1) == 'U' .OR. AGE_DISTRIBUTION(1:1) == 'u') THEN
WRITE(23,*) age(i) - delta_age(i) + REAL(j-1, KIND = 8)/nbins * 2.0_8 * delta_age(i),  RETURN_INFO%MARGINAL_AGES(j,i)
ELSE
WRITE(23,*) age(i) - 2*delta_age(i) + REAL(j-1, KIND = 8)/nbins * 4.0_8 * delta_age(i),  RETURN_INFO%MARGINAL_AGES(j,i)
ENDIF
enddo
ENDDO
CLOSE(23)


OPEN(23, FILE = TRIM(Outputs_directory)//'/intensity_density.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
WRITE(23,*) discretise_size, NBINS
DO i=1, discretise_size
do j=1, NBINS
WRITE(23,*) REAL(X(i), KIND = 4), REAL((REAL(j-1, KIND = 8)+0.5_8)/NBINS * (I_MAX-I_MIN) + I_MIN,KIND =4),&
            REAL(RETURN_INFO%MARGINAL_DENSITY_INTENSITY(i,j), KIND = 4)
enddo
ENDDO
CLOSE(23)


STOP



END PROGRAM

