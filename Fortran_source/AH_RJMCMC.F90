PROGRAM Main_Age_Uncertain_RJMCMC
!
! This program uses the AH-RJMCMC algorithm to sample the posterior distribution of intensity and data ages.

! Authors: Phil Livermore, Alex Fournier, Thomas Bodin, March 27 2018
!
!
use iso_fortran_env, only : error_unit, output_unit
USE AGE_HYPERPARAMETER_RJMCMC


IMPLICIT NONE
INTEGER, PARAMETER :: MAX_DATA  = 1000

CHARACTER(len=500) :: ARG, Data_file_name, Header, WRITE_MODEL_FILE_NAME
CHARACTER(len=500) :: inputline, junk, outputs_directory
INTEGER ::  NARG
CHARACTER(len=1) :: AGE_distribution(1:MAX_DATA), AGE_DISTRIBUTION_TYPE, STRATIFICATION_INDEX(1:MAX_DATA)

INTEGER :: I, K, BURN_IN, NSAMPLE, K_INIT, K_MAX, K_MIN, K_MAX_ARRAYBOUND, discretise_size, show, thin, num,j, &
NUM_DATA, IOS, K_MAX_ARRAY_BOUND, s, ind, NBINS, I_MODEL, FREQ_WRITE_MODELS, FREQ_WRITE_JOINT_DISTRIB, SAMPLE_INDEX_JOINT_DISTRIBUTION
INTEGER :: input_random_seed

REAL( KIND = 8) :: I_MAX, I_MIN, sigma_move, sigma_change_value, sigma_age, sigma_birth, int_j,&
X_MIN, X_MAX, age_frac

INTEGER :: AGE_INDEX(1 : MAX_DATA), NUM_AGE_PARAMETERS
CHARACTER(len=100), Allocatable :: LINE_READ(:)
INTEGER, ALLOCATABLE :: SEED(:)

INTEGER :: stratification(1:MAX_DATA)
INTEGER ::  age_col, d_age_col, F_col, dF_col, distribution_col, id_col, type_col, strat_col

REAL( KIND = 8) :: age(1:MAX_DATA), delta_age(1:MAX_DATA),  intensity(1:MAX_DATA), delta_intensity(1:MAX_DATA)
CHARACTER(len=20) :: ID(1:MAX_DATA)
CHARACTER(1) :: Data_type(1: MAX_DATA), Data_type_specified
CHARACTER(len=10) :: stratification_read_line(1: MAX_DATA)

REAL( KIND = 8), ALLOCATABLE :: X(:)

LOGICAL :: CALC_CREDIBLE, MULTIPLE_STRATIFIED_DATA
Real ( KIND = 8) :: credible

TYPE (RETURN_INFO_STRUCTURE) :: RETURN_INFO

CALC_CREDIBLE = .TRUE.

stratification(:) = 0 !default to no stratification.


! Load parameter file
! af we should use the get_command_argument (2003 standard)
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

!default
input_random_seed = 1

DO
READ(30,'(A)',END=101) inputline
IF ( to_upper(inputline(1:len('Data_file'))) == to_upper('Data_file') ) read(inputline(len('Data_file')+2:),'(A)') Data_file_name
IF ( to_upper(inputline(1:len('Age_distribution'))) == to_upper('Age_distribution') ) read(inputline(len('Age_distribution')+2:),*) AGE_DISTRIBUTION_TYPE
IF ( to_upper(inputline(1:len('Data_type'))) == to_upper('Data_type') ) read(inputline(len('Data_type')+2:),*) Data_type_specified
IF ( to_upper(inputline(1:len('File_format'))) == to_upper('File_format') ) read(inputline(len('File_format')+2:),*) id_col, age_col, d_age_col, F_col, dF_col,  type_col, distribution_col, strat_col
IF ( to_upper(inputline(1:len('Intensity_prior'))) == to_upper('Intensity_prior') ) read(inputline(len('Intensity_prior')+2:),*) I_min, I_max
IF ( to_upper(inputline(1:len('Burn_in'))) == to_upper('Burn_in') ) read(inputline(len('Burn_in')+2:),*) burn_in
IF ( to_upper(inputline(1:len('Nsamples'))) == to_upper('Nsamples') ) read(inputline(len('Burn_in')+2:),*) NSAMPLE
IF ( to_upper(inputline(1:len('model_discretisation'))) == to_upper('model_discretisation') ) read(inputline(len('model_discretisation')+2:),*) discretise_size
IF ( to_upper(inputline(1:len('Chain_parameters'))) == to_upper('Chain_parameters') ) read(inputline(len('Chain_parameters')+2:),*) show, thin
!IF ( to_upper(inputline(1:len('Age_distribution'))) == to_upper('Age_distribution') ) read(inputline(len('Age_distribution')+2:),*) Age_distribution
IF ( to_upper(inputline(1:len('running_mode'))) == to_upper('running_mode') ) read(inputline(len('running_mode')+2:),*) RUNNING_MODE
IF ( to_upper(inputline(1:len('Age_bounds'))) == to_upper('Age_bounds') ) read(inputline(len('Age_bounds')+2:),*) X_min, X_max
IF ( to_upper(inputline(1:len('Sigmas'))) == to_upper('Sigmas') ) read(inputline(len('Sigmas')+2:),*) sigma_move, sigma_change_value, sigma_birth, sigma_age
IF ( to_upper(inputline(1:len('Age_frac'))) == to_upper('Age_frac') ) read(inputline(len('Age_frac')+2:),*) age_frac
IF ( to_upper(inputline(1:len('Num_change_points'))) == to_upper('Num_change_points') ) read(inputline(len('Num_change_points')+2:),*) K_MIN, K_MAX
IF ( to_upper(inputline(1:len('Nbins'))) == to_upper('Nbins') ) read(inputline(len('Nbins')+2:),*) NBINS
IF ( to_upper(inputline(1:len('output_model'))) == to_upper('output_model') ) read(inputline(len('output_model')+2:),*) WRITE_MODEL_FILE_NAME, FREQ_WRITE_MODELS
IF ( to_upper(inputline(1:len('output_joint_distribution_freq'))) == to_upper('output_joint_distribution_freq') ) read(inputline(len('output_joint_distribution_freq')+2:),*) FREQ_WRITE_JOINT_DISTRIB
IF ( to_upper(inputline(1:len('Credible'))) == to_upper('Credible') ) read(inputline(len('Credible')+2:),*) Credible
IF ( to_upper(inputline(1:len('Outputs_directory'))) == to_upper('Outputs_directory') ) read(inputline(len('Outputs_directory')+2:),'(A)') outputs_directory
IF ( to_upper(inputline(1:len('SEED'))) == to_upper('SEED') ) read(inputline(len('SEED')+2:),*) input_random_seed
ENDDO
101     CONTINUE

CLOSE( 30 )


PRINT*, 'DATA file name: ', TRIM(Data_file_name)
OPEN(30, FILE = Data_file_name, STATUS = 'OLD', FORM = 'FORMATTED', &
IOSTAT = IOS, ACTION = 'READ')
IF( IOS .NE. 0) THEN
PRINT*, 'ERROR IN OPENING FILE ', TRIM(Data_file_name)
STOP
ENDIF

ALLOCATE( LINE_READ(1: max(id_col, age_col, d_age_col, F_col, dF_col, distribution_col, type_col, strat_col)+1) )

write(OUTPUT_unit, fmt="(a,2x,i4)") "strat_col = ", strat_col

I=1
DO
READ(30,'(A)', END = 1001) inputline

! remove leading white spaces - needed to sense a '#' character as the first non-space character
inputline = ADJUSTL(TRIM(inputline))

IF( inputline(1:1) == '#') CYCLE

READ(inputline,*) LINE_READ(:)

READ(LINE_READ(id_col+1),*) id(i)
READ(LINE_READ(age_col+1),*) age(i)
READ(LINE_READ(d_age_col+1),*) delta_age(i)
READ(LINE_READ(F_col+1),*) intensity(i)
READ(LINE_READ(dF_col+1),*) delta_intensity(i)

IF( type_col .GE. 0) READ(LINE_READ(type_col+1),*) data_type(i)
IF( distribution_col .GE. 0) READ(LINE_READ(distribution_col+1),*) AGE_distribution(i)
IF( strat_col .NE. -1) READ(LINE_READ(strat_col+1),*) stratification_read_line(i)


i = i + 1
IF( i > MAX_DATA ) THEN
PRINT*, 'NUMBER OF DATA EXCEEDS HARDWIRED MAX', MAX_DATA
STOP
ENDIF

ENDDO
1001    CONTINUE
CLOSE(30)
NUM_DATA = i-1

! # Data type:  brick (B), baked clay (C), slag (S), pottery (P), other (O).
! If brick then data are grouped together; otherwise the type is ignored.



k_max_array_bound = k_max + 1
NUM=ceiling((nsample-burn_in)*(100-credible)/200.0_8/thin) ! number of collected samples for credible intervals

ALLOCATE( X(1: discretise_size) )
ALLOCATE( RETURN_INFO%AV(1:discretise_size), RETURN_INFO%BEST(1:discretise_size),  &
RETURN_INFO%INF(1:discretise_size), RETURN_INFO%sup(1:discretise_size), &
!         RETURN_INFO%change_points(nsample * k_max_array_bound), & af to save memory
RETURN_INFO%n_changepoint_hist(k_max_array_bound), &
RETURN_INFO%convergence(nsample),  & ! af
RETURN_INFO%MARGINAL_AGES(1:NBINS,1:NUM_DATA)    )


call system('mkdir -p '//TRIM(Outputs_directory))
open(newunit = RETURN_INFO%output_changepoints_unit, file=TRIM(Outputs_directory)//'/changepoints.dat', &
status="replace", form="formatted")

ALLOCATE( RETURN_INFO%MEDIAN(1:discretise_size), RETURN_INFO%MODE(1:discretise_size), &
RETURN_INFO%MARGINAL_DENSITY_INTENSITY(1:discretise_size,1:NBINS) )

ALLOCATE( RETURN_INFO%INF_DFDT(1:discretise_size), RETURN_INFO%sup_dFDT(1:discretise_size) )
ALLOCATE( RETURN_INFO%AV_DFDT(1:discretise_size), RETURN_INFO%MEDIAN_DFDT(1:discretise_size) )
ALLOCATE( RETURN_INFO%MODE_DFDT(1:discretise_size) )




CALL RANDOM_SEED(SIZE = i)
ALLOCATE( SEED(1:i) )
SEED(:) = input_random_seed
CALL RANDOM_SEED(PUT = SEED)
PRINT*, 'INITIALISING RANDON NUMBERS USING SEED ', input_random_seed

IF( X_MIN .eq. X_MAX) THEN
X_MIN = MINVAL( AGE(1:NUM_DATA) ) - X_MIN
X_MAX = MAXVAL( AGE(1:NUM_DATA) ) + X_MAX
ENDIF

! Check to see that X_MIN and X_MAX enclose the data set:

IF( to_upper(AGE_distribution(i)) == 'U') THEN
IF( X_MIN .GT. MINVAL(AGE(1:NUM_DATA)-delta_age(1:NUM_DATA)) .OR. X_MAX .LT. MAXVAL(AGE(1:NUM_DATA)+delta_age(1:NUM_DATA) ) ) THEN
PRINT*, 'INCREASE RANGE OF END POINT AGES'
WRITE(6,'(A,F8.1,A,F8.1)') 'RANGE NEEDS TO SPAN AT LEAST ', MINVAL(AGE(1:NUM_DATA)-delta_age(1:NUM_DATA)),' : ', MAXVAL(AGE(1:NUM_DATA)+delta_age(1:NUM_DATA))
STOP
ENDIF
!
ELSE !assumed normal distribution
IF( X_MIN .GT. MINVAL(AGE(1:NUM_DATA)-2.0_8 * delta_age(1:NUM_DATA)) .OR. X_MAX .LT. MAXVAL(AGE(1:NUM_DATA)+2.0_8 * delta_age(1:NUM_DATA) ) ) THEN
PRINT*, 'INCREASE RANGE OF END POINT AGES'
WRITE(6,*) 'ALL DATA DOES NOT LIE WITHIN 2 S.D. OF THE END POINTS'
WRITE(6,'(A,F8.1,A,F8.1)') 'RANGE NEEDS TO SPAN AT LEAST ', MINVAL(AGE(1:NUM_DATA)-2.0_8 * delta_age(1:NUM_DATA)),' : ', MAXVAL(AGE(1:NUM_DATA)+2.0_8 * delta_age(1:NUM_DATA))
STOP
ENDIF
ENDIF


PRINT*, 'READ IN ', NUM_DATA, 'DATA ITEMS'
PRINT*, 'DATA AGE RANGE IS ', MINVAL(AGE(1:NUM_DATA)), ' TO ', MAXVAL( AGE(1:NUM_DATA))
PRINT*, 'MODEL AGE RANGE IS ', X_MIN, ' : ', X_MAX

IF( type_col < 0 ) THEN
data_type(:) = 'O'   ! Set to O(ther): never referenced.
PRINT*, 'NO SPECIFICATION OF DATA TYPE: SETTING GLOBAL TYPE TO (O)THER'
ENDIF

IF( distribution_col .EQ. -1) THEN
AGE_DISTRIBUTION(:) = 'U'
PRINT*, 'SETTING GLOBAL AGE DISTRIBUTION TO UNIFORM'
ELSEIF( distribution_col .EQ. -2) THEN
AGE_DISTRIBUTION(:) = 'N'
PRINT*, 'SETTING GLOBAL AGE DISTRIBUTION TO NORMAL'
ELSEIF( distribution_col < 0) THEN
PRINT*, 'ERROR: UNKNOWN VALUE FOR PARAMETER DESCRIBING AGE DISTRIBUTION ', distribution_col
STOP
ENDIF




DO I=1, discretise_size
X(I) = X_MIN + REAL(I-1, KIND = 8)/REAL(discretise_size-1, KIND = 8) * (X_MAX - X_MIN)
ENDDO

!****************
! STRATIFICATION
!****************
! Interpret stratification information
! First, see if the user has specified different datasets using 1a, 2a, 3a; 1b, 2b; etc notation, or simply 1,2,3,4 etc.
! We can tell these apart by looking for the first non-zero instance of stratification_read_line(:) and seeing whether it ends with an 'a'.
MULTIPLE_STRATIFIED_DATA = .FALSE.
DO i = 1, NUM_DATA
IF(index(stratification_read_line(i), 'a') .NE. 0) MULTIPLE_STRATIFIED_DATA = .TRUE.
ENDDO
!PRINT*, MULTIPLE_STRATIFIED_DATA

IF( MULTIPLE_STRATIFIED_DATA ) THEN ! read in
! separate into integer and character if non-zero
STRATIFICATION(:) = 0
STRATIFICATION_INDEX(:) = ' '
DO i = 1, NUM_DATA
IF( stratification_read_line(i)(1:1) .NE. '0') THEN
J = LEN(TRIM(stratification_read_line(i)))
READ( stratification_read_line(i)(1:j-1), fmt="(i1)") STRATIFICATION(i)
READ( stratification_read_line(i)(j:j),'(A)') STRATIFICATION_INDEX(i)
ENDIF
ENDDO
ELSE  !either no stratification, or only a single dataset is specified without the 'a' notation. In either case, set all elements of stratification_index to 'a'
if ( strat_col .NE. -1) then
write(ERROR_UNIT,fmt="(a)") 'here'
write(ERROR_UNIT,*) stratification_read_line(1:num_data)
DO I = 1, NUM_DATA
!af
write(output_unit,*) stratification_read_line(i)
READ( stratification_read_line(i), *) STRATIFICATION(I)
ENDDO
stratification_index(1:NUM_DATA) = 'a'
end if
ENDIF


IF( MAXVAL( stratification(1:NUM_DATA)) .eq. 0) THEN
PRINT*, 'NO STRATIFICATION CONSTRAINTS'
ELSE
PRINT*, 'STRATIFICATION CONSTRAINTS BEING USED'
ENDIF

!DO i=1,NUM_DATA
!PRINT*, STRATIFICATION(i), STRATIFICATION_INDEX(i)
!ENDDO

! If data_type is used, then gather together data of type 'B' with common ID. Assumes grouped data is listed sequentially.
! Grouped samples are moved collectively by the MCMC algorithm when changing a sample date.
! Otherwise treat as independent samples

IF( type_col < 0) THEN  ! no data type set - no grouping needed.
DO i=1, NUM_DATA
AGE_INDEX(i) = i
ENDDO
NUM_AGE_PARAMETERS = NUM_DATA
ELSE
J = 1
AGE_INDEX(J) = 1; J = J + 1
DO i=2, NUM_DATA

IF( ID(I) == ID(I-1) .AND. DATA_TYPE(I) == 'B' .AND. DATA_TYPE(I-1) == 'B') THEN
CYCLE
ELSE
! NEW AGE PARAMETER
AGE_INDEX(J) = I; J = J + 1
ENDIF
ENDDO
NUM_AGE_PARAMETERS = J-1
OPEN(13, FILE = 'AGE_PARAMETERS', FORM = 'FORMATTED', STATUS = 'REPLACE')
DO i = 1,NUM_AGE_PARAMETERS
WRITE(13,*) I, AGE_INDEX(I), AGE_INDEX(I+1)
ENDDO
CLOSE(13)

!PRINT*, AGE_INDEX(1:NUM_AGE_PARAMETERS )

! Check to make sure the data associated with single age parameter has been listed sequentially, otherwise it will be split into multiple age indices.
DO I=1,NUM_AGE_PARAMETERS
IF  (DATA_TYPE(AGE_INDEX(I)) == 'B') THEN
DO J = I+1,NUM_AGE_PARAMETERS
IF( DATA_TYPE(AGE_INDEX(J)) == 'B'  .AND. ID(AGE_INDEX(I)) == ID(AGE_INDEX(J)) ) THEN
write(unit=error_unit, fmt=*) 'FOUND DISJOINT OCCURENCES OF SITE ID ', ID(I), ' WITH TYPE B', AGE_INDEX(I), AGE_INDEX(J)
STOP
ENDIF
ENDDO
ENDIF
ENDDO

ENDIF
PRINT*, 'NUMBER OF AGE HYPERPARAMETERS : ', NUM_AGE_PARAMETERS

CALC_CREDIBLE = .TRUE.


!call system('mkdir -p '//TRIM(Outputs_directory))

! copy input file
call system('cp '//TRIM(ARG)//' '//TRIM(Outputs_directory)//'/input_file')
CALL RJMCMC(burn_in, NUM_DATA, age(1:NUM_DATA), delta_age(1:NUM_DATA), intensity(1:NUM_DATA), delta_intensity(1:NUM_DATA), stratification(1: NUM_DATA), STRATIFICATION_INDEX(1:NUM_DATA), AGE_DISTRIBUTION(1:NUM_DATA), AGE_INDEX(1:NUM_AGE_PARAMETERS), NSAMPLE, I_MIN, I_MAX, X_MIN, X_MAX, K_MIN, K_MAX, SIGMA_MOVE, sigma_change_value, sigma_birth, sigma_age, age_frac, discretise_size, SHOW, THIN, NBINS, RETURN_INFO, CALC_CREDIBLE, FREQ_WRITE_MODELS, WRITE_MODEL_FILE_NAME, FREQ_WRITE_JOINT_DISTRIB,   credible, Outputs_directory)



OPEN(24, FILE = TRIM(Outputs_directory)//'/data.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, NUM_DATA
WRITE(24, fmt="(4(F10.3,X),X,I3)") AGE(i), delta_age(i), intensity(I), delta_intensity(I), stratification(i)
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
WRITE(24,fmt="(i9,X,ES12.3)") i, RETURN_INFO%convergence(i)
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

!OPEN(23, FILE = TRIM(Outputs_directory)//'/changepoints.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
!DO i=1, RETURN_INFO%MAX_NUMBER_CHANGE_POINTS_HISTORY
!WRITE(23,'(F10.3)') RETURN_INFO%change_points(i)
!enddo
CLOSE(unit=RETURN_info%output_changepoints_unit)

OPEN(23, FILE = TRIM(Outputs_directory)//'/age_marginals.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
WRITE(23,*) NUM_DATA, NBINS
DO i=1, NUM_DATA
do j=1, NBINS
IF( AGE_DISTRIBUTION(i) == 'U' .OR. AGE_DISTRIBUTION(i) == 'u') THEN
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

OPEN(24, FILE = TRIM(Outputs_directory)//'/credible_dFdt_upper.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i), RETURN_INFO%sup_dFdt(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/credible_dFdt_lower.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i), RETURN_INFO%inf_dFdt(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/average_dFdt.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i),RETURN_INFO%av_DFDT(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/mode_dFdt.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i),RETURN_INFO%mode_DFDT(i)
enddo
CLOSE(24)

OPEN(24, FILE = TRIM(Outputs_directory)//'/median_dFdt.dat', STATUS = 'REPLACE', FORM = 'FORMATTED')
DO i=1, discretise_size
WRITE(24,'(F10.3,X,F10.3)') x(i),RETURN_INFO%median_DFDT(i)
enddo
CLOSE(24)



STOP



END PROGRAM
