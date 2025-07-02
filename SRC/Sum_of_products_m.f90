!==================================================================================================
!==================================================================================================
! This file is part of MolecCav.
!
!==================================================================================================
! MIT License
!
! Copyright (c) 2025 Mathéo Segaud
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!==================================================================================================
! README :
! The only module related to general HO that the others modules will need to call in a "USE".  
! Initialize_operator_ND : reads the namelist and initialize the type, then constructs the operat-
! or using parameters of the HO1D_para object from the so called derived type.
! Append_operator_ND     : add an HO operator to a already initialized object of Operator_ND_t type.
! Write_operator_ND      : display values of the type in the output
! Deallocate_operator_ND : deallocate all tables of the type
! The module to initialize the HO by reading its parameters from the namelist.  
! Read_HO1D_parameters  : reads the namelist and initialize the type.
! Write_HO1D_parameters : displays values of the type in the output.
!==================================================================================================
!==================================================================================================
MODULE Sum_of_products_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT, real64
  USE QDUtil_m                                                                 ! gives Rkind=real64; out_unit=OUTPUT_UNIT; INPUT_UNIT=in_unit; EYE=i and other numbers; TO_LOWERCASE; TO_UPPERCASE;... We thereby use ZERO instead of 0.0_real64
  USE Operator_ND_m
  IMPLICIT NONE


  TYPE                               :: Sum_of_products_t                      ! N_mat, N_cav canNOT be part of the derived type because they are not specific of one operator_ND, but parameters of the whole system/calculation. All OpND will have the same. (therefore only one namelist per mode is needed, and not one per mode and oOND)
    integer                          :: N_products = 0                         ! the number of terms in the sum of products operator
    TYPE(Operator_ND_t), allocatable :: tab_opnd(:)                            ! the list of the products terms in the sum, which one being an OpND
    real(kind=Rkind),    allocatable :: tab_coeffs(:)
  END TYPE


  PRIVATE

  PUBLIC Sum_of_products_t, Initialize_totH, Initialize_dipmomt, Initialize_sop, Action, Get, Write, Dealloc

    INTERFACE Initialize_totH
      MODULE PROCEDURE MolecCav_Initialize_total_hamiltonian
     END INTERFACE
    INTERFACE Initialize_dipmomt
      MODULE PROCEDURE MolecCav_Initialize_dipole_moment
     END INTERFACE
    INTERFACE Initialize_sop !/!\/!\/!\ NEVER HAVE BEEN TESTED /!\/!\/!\
      MODULE PROCEDURE MolecCav_Initialize_sum_of_products
    END INTERFACE
    INTERFACE Read_pdt !/!\/!\/!\ NEVER HAVE BEEN TESTED /!\/!\/!\
      MODULE PROCEDURE MolecCav_Read_products
    END INTERFACE
  INTERFACE Action
    MODULE PROCEDURE MolecCav_Action_SOP_R1_real, MolecCav_Action_SOP_R1_complex
  END INTERFACE
   INTERFACE Get
     MODULE PROCEDURE MolecCav_Get_SOP_parameter_integer
   END INTERFACE
   INTERFACE Write
     MODULE PROCEDURE MolecCav_Write_sum_of_products
   END INTERFACE
   INTERFACE Dealloc
     MODULE PROCEDURE MolecCav_Deallocate_sum_of_products
   END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_total_hamiltonian(TotH, nio, Dense, Verbose, Debug) ! no need for N_mat and N_cav explicitly : they are SIZE(Mat_op and Cav_op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE
  
    TYPE(Sum_of_products_t), intent(inout) :: TotH
    integer,                 intent(in)    :: nio
    logical, optional,       intent(in)    :: Dense                                                                        ! cf. comments in HO1D_parameters_m
    integer, optional,       intent(in)    :: Verbose                                                                      ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                        ! cf. comments in HO1D_parameters_m

    integer                                :: N_mat, N_cav, i_product, i_mat, i_cav
    real(kind=Rkind), allocatable          :: tab_coeffs(:)
    character(len=:), allocatable          :: Mat_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_mat>', ex : 'hamiltonian, Identity'. Not case sensitive, ' ' <=> \otimes
    character(len=:), allocatable          :: Cav_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_cav>', ex : 'hamiltonian'.   Not case sensitive, ' ' <=> \otimes. This exemple means OpND = H_mat1\otimesI_mat2\otimesH_cav
    real(kind=Rkind)                       :: Cavw, Cavlambda, Matlambda
    integer                                :: err_io
    logical                                :: Dense_local                                                                  ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                                :: Verbose_local                                                                ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                                :: Debug_local
    
    NAMELIST /TOTAL_HAMILTONIAN/ N_mat, N_cav

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE TOTAL HAMILTONIA&
                                              &N OBJECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_total_hamiltonian :"
      WRITE(out_unit,*) "The <<TotH>> argument :"
      CALL Write(TotH)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_total_hamiltonian"
      FLUSH(out_unit)
    END IF
    
    !----------------------Initialization to default values--------------------
    N_mat = 0
    N_cav = 0

    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** READING THE NUMBER OF MODES *************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = TOTAL_HAMILTONIAN, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = TOTAL_HAMILTONIAN)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
      !------------------------------Check reading error-------------------------
    IF(err_io /= 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '##### Error in MolecCav_Initialize_total_hamiltonian (err_io/=0) #####'
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '####################### err_io = ', err_io, '######################'
      STOP '######################### Check basis data ########################'
    END IF

    IF (N_mat == 0 .AND. N_cav == 0) THEN
      WRITE(out_unit,*) "### The number of modes of the system CANNOT be 0 (what are are you going to stu&
      &dy if there is nothing to ???). Please check the data file '.nml'"
      STOP "### The number of modes of the system CANNOT be 0 (what are are you going to study if there i&
      &s nothing ???). Please check the data file '.nml'"
    END IF
    
    !---------------Construction of the table of coefficients of the sum-----------
    TotH%N_products = N_mat + N_cav + N_mat * N_cav 

    ALLOCATE(TotH%tab_opnd(TotH%N_products))
    ALLOCATE(TotH%tab_coeffs(TotH%N_products))
    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "TotH%N_products : "//TO_string(TotH%N_products)
    TotH%tab_coeffs = 1

    !---------------Construction of the table of OpND composing the sum of products operator-----------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

      !##### matter hamiltonians #####
    ALLOCATE(character(len=    11 + 10*(N_mat-1)    ) :: Mat_operators)                                                   ! 11 ("hamiltonian") + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    ALLOCATE(character(len=MAX(8  + 10*(N_cav-1), 0)) :: Cav_operators)                                                   ! 8  ("identity")    + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    Mat_operators = ""
    Cav_operators = ""
    IF (N_mat > 0) Mat_operators = "Hamiltonian"//REPEAT(", Identity", N_mat-1)
    IF (N_cav > 0) Cav_operators = "Identity"//REPEAT(", Identity", N_cav-1)
    IF (N_mat > 0) CALL Initialize(TotH%tab_opnd(1), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)

    DO i_product = 2, N_mat
      Mat_operators = "Identity"//REPEAT(", Identity", i_product-2)//", Hamiltonian"//REPEAT(", Identity", N_mat-i_product)
      CALL Initialize(TotH%tab_opnd(i_product), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)
    END DO
    DEALLOCATE(Mat_operators); DEALLOCATE(Cav_operators)

      !##### cavity hamiltonians #####
    ALLOCATE(character(len=MAX(8  + 10*(N_mat-1), 0)) :: Mat_operators)                                                   ! 8  ("identity")    + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    ALLOCATE(character(len=11 + 10*(N_cav-1))         :: Cav_operators)                                                   ! 11 ("hamiltonian") + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    Mat_operators = ""
    Cav_operators = ""
    IF (N_mat > 0) Mat_operators = "Identity"//REPEAT(", Identity", N_mat-1)
    IF (N_cav > 0) Cav_operators = "Hamiltonian"//REPEAT(", Identity", N_cav-1)
    CALL Initialize(TotH%tab_opnd(N_mat+1), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)

    DO i_product = 2, N_cav
      Cav_operators = "Identity"//REPEAT(", Identity", i_product-2)//", Hamiltonian"//REPEAT(", Identity", N_cav-i_product)
      CALL Initialize(TotH%tab_opnd(N_mat+i_product), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)
    END DO
    DEALLOCATE(Mat_operators); DEALLOCATE(Cav_operators)

      !##### coupling terms #####
    IF (N_mat > 0) ALLOCATE(character(len=7 + 10*(N_mat-1)) :: Mat_operators)                                                   ! 11 ("hamiltonian") + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    IF (N_cav > 0) ALLOCATE(character(len=8 + 10*(N_cav-1)) :: Cav_operators)                                                   ! 8  ("identity")    + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    i_product = N_mat + N_cav + 1
    DO i_cav = 1, N_cav ! if ncav or nmat == 0, the program does not even enter in the loops : make perfectly sens since there is then no couplings !
      DO i_mat = 1, N_mat
        IF (i_mat == 1) THEN
          Mat_operators = "DipMomt"//REPEAT(", Identity", N_mat-1)
        ELSE 
          Mat_operators = "Identity"//REPEAT(", Identity", i_mat-2)//", DipMomt"//REPEAT(", Identity", N_mat-i_mat)
        END IF 
        IF (i_cav == 1) THEN
          Cav_operators = "Position"//REPEAT(", Identity", N_cav-1)
        ELSE 
          Cav_operators = "Identity"//REPEAT(", Identity", i_cav-2)//", Position"//REPEAT(", Identity", N_cav-i_cav)
        END IF 
        CALL Get(Cavw, "w", "Cavity", i_cav)
        CALL Get(Cavlambda, "lambda", "Cavity", i_cav)
        CALL Get(Matlambda, "lambda", "Matter", i_mat)

        IF (Debug_local) THEN
          WRITE(out_unit,*)
          WRITE(out_unit,*) "--- System parameters for the coupling term between the "//TO_string(i_mat)//"^{th} matter mode and &
          &the "//TO_string(i_cav)//"^{th} cavity mode :"
          WRITE(out_unit,*) "Cavw          = "//TO_string(Cavw)
          WRITE(out_unit,*) "Cavlambda     = "//TO_string(Cavlambda)
          WRITE(out_unit,*) "Matlambda     = "//TO_string(Matlambda)
          WRITE(out_unit,*) "Mat_operators = "//Mat_operators
          WRITE(out_unit,*) "Cav_operators = "//Cav_operators
        END IF

        CALL Initialize(TotH%tab_opnd(i_product), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)
        TotH%tab_coeffs(i_product) = Cavlambda * Matlambda * Cavw
        i_product = i_product + 1
      END DO
    END DO
    IF (ALLOCATED(Mat_operators)) DEALLOCATE(Mat_operators)
    IF (ALLOCATED(Cav_operators)) DEALLOCATE(Cav_operators)

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Sum of products constructed by MolecCav_Initialize_total_hamiltonian--------------"
      CALL Write(TotH)
      WRITE(out_unit,*) "------------End Sum of products constructed by MolecCav_Initialize_total_hamiltonian------------"
    END IF
    
    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "--------------------------------------------------TOTAL HAMILTONIAN OPERATOR INITIALIZED-&
    &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_total_hamiltonian


  SUBROUTINE MolecCav_Initialize_dipole_moment(DipMomt, nio, Dense, Verbose, Debug) ! no need for N_mat and N_cav explicitly : they are SIZE(Mat_op and Cav_op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE
  
    TYPE(Sum_of_products_t), intent(inout) :: DipMomt
    integer,             intent(in)        :: nio
    logical, optional,   intent(in)        :: Dense                                                                        ! cf. comments in HO1D_parameters_m
    integer, optional,   intent(in)        :: Verbose                                                                      ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)        :: Debug                                                                        ! cf. comments in HO1D_parameters_m

    integer                                :: N_mat, N_cav, i_product, i_mat, i_cav
    real(kind=Rkind), allocatable          :: tab_coeffs(:)
    character(len=:), allocatable          :: Mat_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_mat>', ex : 'hamiltonian, Identity'. Not case sensitive, ' ' <=> \otimes
    character(len=:), allocatable          :: Cav_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_cav>', ex : 'hamiltonian'.   Not case sensitive, ' ' <=> \otimes. This exemple means OpND = H_mat1\otimesI_mat2\otimesH_cav
    real(kind=Rkind)                       :: Cavw, Cavlambda, Matlambda
    integer                                :: err_io
    logical                                :: Dense_local                                                                  ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                                :: Verbose_local                                                                ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                                :: Debug_local
    
    NAMELIST /DIPOLE_MOMENT/ N_mat, N_cav

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE DIPOLE MOMENT OBJECT---&
                                              &----------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_dipole_moment :"
      WRITE(out_unit,*) "The <<DipMomt>> argument :"
      CALL Write(DipMomt)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_dipole_moment"
      FLUSH(out_unit)
    END IF
    
    !----------------------Initialization to default values--------------------
    N_mat = 0
    N_cav = 0

    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** READING THE NUMBER OF MODES *************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = DIPOLE_MOMENT, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = DIPOLE_MOMENT)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
      !------------------------------Check reading error-------------------------
    IF(err_io /= 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '###### Error in MolecCav_Initialize_dipole_moment (err_io/=0) #####'
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '####################### err_io = ', err_io, '######################'
      STOP '######################### Check basis data ########################'
    END IF

    IF (N_mat == 0) THEN
      WRITE(out_unit,*) "### The number of matter modes CANNOT be 0 to construct the global matter dipole moment ! Please check t&
      &he data file '.nml'"
      STOP "### The number of matter modes CANNOT be 0 to construct the global matter dipole moment ! Please check the data file &
      &'.nml'"
    END IF
    
    !---------------Construction of the table of coefficients of the sum-----------
    DipMomt%N_products = N_mat 

    ALLOCATE(DipMomt%tab_opnd(DipMomt%N_products))
    ALLOCATE(DipMomt%tab_coeffs(DipMomt%N_products))
    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "DipMomt%N_products : "//TO_string(DipMomt%N_products)
    DipMomt%tab_coeffs = 1

    !---------------Construction of the table of OpND composing the sum of products operator-----------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

      !##### matter modes #####
    ALLOCATE(character(len=    7 + 10*(N_mat-1)    ) :: Mat_operators)                                                   ! 11 ("hamiltonian") + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    ALLOCATE(character(len=MAX(8 + 10*(N_cav-1), 0)) :: Cav_operators)                                                   ! 8  ("identity")    + (N_mat-1)*8 ("identity") + (N_mat-1)*2 (", ")
    Mat_operators = ""
    Cav_operators = ""
    Mat_operators                = "DipMomt"//REPEAT(", Identity", N_mat-1)
    IF (N_cav > 0) Cav_operators = "Identity"//REPEAT(", Identity", N_cav-1)
    CALL Initialize(DipMomt%tab_opnd(1), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)

    DO i_product = 2, N_mat
      Mat_operators = "Identity"//REPEAT(", Identity", i_product-2)//", DipMomt"//REPEAT(", Identity", N_mat-i_product)
      CALL Initialize(DipMomt%tab_opnd(i_product), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)
    END DO
    DEALLOCATE(Mat_operators); DEALLOCATE(Cav_operators)

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Sum of products constructed by MolecCav_Initialize_dipole_moment--------------"
      CALL Write(DipMomt)
      WRITE(out_unit,*) "------------End Sum of products constructed by MolecCav_Initialize_dipole_moment------------"
    END IF
    
    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "--------------------------------------------------DIPOLE MOMENT OPERATOR INITIALIZED-&
    &------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_dipole_moment


  SUBROUTINE MolecCav_Initialize_sum_of_products(SumProduct, nio, Dense, Verbose, Debug) ! no need for N_mat and N_cav explicitly : they are SIZE(Mat_op and Cav_op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE
  
    TYPE(Sum_of_products_t), intent(inout) :: SumProduct
    integer,                 intent(in)    :: nio
    logical, optional,       intent(in)    :: Dense                                                                        ! cf. comments in HO1D_parameters_m
    integer, optional,       intent(in)    :: Verbose                                                                      ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                        ! cf. comments in HO1D_parameters_m

    integer                                :: N_products, i_product
    real(kind=Rkind), allocatable          :: tab_coeffs(:)
    character(len=:), allocatable          :: Mat_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_mat>', ex : 'hamiltonian, Identity'. Not case sensitive, ' ' <=> \otimes
    character(len=:), allocatable          :: Cav_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_cav>', ex : 'hamiltonian'.   Not case sensitive, ' ' <=> \otimes. This exemple means OpND = H_mat1\otimesI_mat2\otimesH_cav
    integer                                :: err_io_1, err_io_2
    logical                                :: Dense_local                                                                  ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                                :: Verbose_local                                                                ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                                :: Debug_local
    
    NAMELIST /NUMBER_OF_PRODUCTS/ N_products
    NAMELIST /COEFFS_SUM_OF_PRODUCTS/ tab_coeffs

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE SUM OF PRODUCTS OPERATOR&
                                              & OBJECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_sum_of_products :"
      WRITE(out_unit,*) "The <<SumProduct>> argument :"
      CALL Write(SumProduct)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_sum_of_products"
      FLUSH(out_unit)
    END IF
    
    !----------------------Initialization to default values--------------------
    N_products = 0

    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** READING THE NUMBER OF OPND **************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = NUMBER_OF_PRODUCTS, iostat = err_io_1)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = NUMBER_OF_PRODUCTS)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
      !------------------------------Check reading error-------------------------
    IF(err_io_1 /= 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '##### Error in MolecCav_Initialize_sum_of_products (err_io_1/=0) #####'
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '####################### err_io_1 = ', err_io_1, '######################'
      STOP '######################### Check basis data ########################'
    END IF

    IF (N_products == 0) THEN
      WRITE(out_unit,*) "### The number of ND Operators in the sum of products operator CANNOT be 0 (what are are you going to stu&
      &dy if there is no operator ???). Please check the data file '.nml'"
      STOP "### The number of ND Operators in the sum of products operator CANNOT be 0 (what are are you going to study if there i&
      &s no operator ???). Please check the data file '.nml'"
    END IF
    
    !---------------Construction of the table of coefficients of the sum-----------
    SumProduct%N_products = N_products 
    ALLOCATE(SumProduct%tab_coeffs(SumProduct%N_products))
    ALLOCATE(tab_coeffs(SumProduct%N_products))
    tab_coeffs = 0

    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '********************** READING THE COEFFICIENTS OF THE SUM *********************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = COEFFS_SUM_OF_PRODUCTS, iostat = err_io_2)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = COEFFS_SUM_OF_PRODUCTS)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
      !------------------------------Check reading error-------------------------
    IF(err_io_2 /= 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '##### Error in MolecCav_Initialize_sum_of_products (err_io_2/=0) #####'
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '####################### err_io_2 = ', err_io_2, '######################'
      STOP '######################### Check basis data ########################'
    END IF

    SumProduct%tab_coeffs = tab_coeffs

    !---------------Construction of the table of OpND composing the sum of products operator-----------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    ALLOCATE(SumProduct%tab_opnd(SumProduct%N_products))

    DO i_product = 1, N_products
      CALL Read_pdt(Mat_operators, Cav_operators, nio, Verbose=Verbose_local, Debug=Debug_local)
      CALL Initialize(SumProduct%tab_opnd(i_product), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)
      DEALLOCATE(Mat_operators); DEALLOCATE(Cav_operators)
    END DO

    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** SUM OF PRODUCTS CONSTRUCTED *************************'
    WRITE(out_unit,*) '********************************************************************************'

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Sum of products constructed by MolecCav_Initialize_sum_of_products--------------"
      CALL Write(SumProduct)
      WRITE(out_unit,*) "------------End Sum of products constructed by MolecCav_Initialize_sum_of_products------------"
    END IF
    
    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "--------------------------------------------------SUM OF PRODUCTS OPERATOR INITIAL&
    &IZED-------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_sum_of_products


  SUBROUTINE MolecCav_Read_products(Mat_operators, Cav_operators, nio, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE
  
    character(len=:), allocatable, intent(inout) :: Mat_operators
    character(len=:), allocatable, intent(inout) :: Cav_operators
    integer,                       intent(in)    :: nio
    integer, optional,             intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical, optional,             intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    character(len=200)                           :: Mat_operators_local
    character(len=200)                           :: Cav_operators_local
    integer                                      :: err_io
    integer                                      :: Verbose_local                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                                      :: Debug_local

    NAMELIST /PRODUCTS_OF_OP1D/ Mat_operators_local, Cav_operators_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "-------------------------------------------------READING ONE PRODUCTS OF THE SUM OF&
                                              & PRODUCTS OPERATOR-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Read_products :"
      WRITE(out_unit,*) "Is the <<Mat_operators>> argument already allocated ?"//TO_string(ALLOCATED(Mat_operators))
      WRITE(out_unit,*) "Is the <<Cav_operators>> argument already allocated ?"//TO_string(ALLOCATED(Cav_operators))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Read_products"
      FLUSH(out_unit)
    END IF
    
    IF (ALLOCATED(Mat_operators) .OR. ALLOCATED(Cav_operators)) THEN
      WRITE(out_unit,*) "### Please mind that Mat_operator and Cav_operator have not to be already allocated when call the Read_p&
      &roduct procedure."
      STOP "### Mat_operator and Cav_operator already allocated at the Read_products call."
    END IF 

    !----------------------Initialization to default values--------------------
    Mat_operators_local = ""
    Cav_operators_local = ""

    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** READING THE PRODUCTS_OF_OP1D *************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = PRODUCTS_OF_OP1D, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = PRODUCTS_OF_OP1D)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
    !------------------------------Check reading error-------------------------
    IF(err_io < 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '########## Error in Read_products (err_io/=0) ##########'
      WRITE(out_unit,*) '#######################################################'
      WRITE(out_unit,*) '################# err_io = ', err_io, '################'
      STOP '################### Check basis data ##################'
    END IF
    
    !---------------Construction of the two strings-----------
    ALLOCATE(character(len=LEN_TRIM(Mat_operators_local)) :: Mat_operators)                                                   ! /!\ strings cannot be allocated the exact same way as tables ! /!\
    ALLOCATE(character(len=LEN_TRIM(Cav_operators_local)) :: Cav_operators)                                                   ! /!\ strings cannot be allocated the exact same way as tables ! /!\

    Mat_operators = TO_lowercase(TRIM(Mat_operators_local))
    Cav_operators = TO_lowercase(TRIM(Cav_operators_local))

    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '********************************** PRODUCTS READ ********************************'
    WRITE(out_unit,*) '********************************************************************************'

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Op1D products read by MolecCav_MolecCav_Read_products--------------"
      WRITE(out_unit,*) "Mat_operator : "//Mat_operators
      WRITE(out_unit,*) "Cav_operator : "//Cav_operators
      WRITE(out_unit,*) "------------End Op1D products read by MolecCav_MolecCav_Read_products------------"
    END IF

  END SUBROUTINE MolecCav_Read_products


  SUBROUTINE MolecCav_Action_SOP_R1_real(Op_psi, SumProduct, Psi, Verbose, Debug) ! Psi is ND AND R1
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE

    real(kind=Rkind),        intent(inout) :: Op_psi(:)
    TYPE(Sum_of_products_t), intent(in)    :: SumProduct
    real(kind=Rkind),        intent(in)    :: Psi(:)
    integer, optional,       intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                                :: i_product
    real(kind=Rkind), allocatable          :: Op_psi_local(:)
    integer                                :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "---------------------------------------COMPUTING ACTION OF THE SOP OVER &
                                              &THE R1 ND WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_SOP_R1_real :"
      WRITE(out_unit,*) "The <<SumProduct>> argument :"
      CALL Write(SumProduct)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_SOP_R1_real"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! THE DIMENSIONS OF EACH 1D MATMUL WILL BE TESTED IN THE ACTIONS CODED IN ELEM_OP_M !

    IF (SIZE(SumProduct%tab_opnd) /= SIZE(SumProduct%tab_coeffs)) THEN
      WRITE(out_unit,*) "The number of products in the sum should match the number of coefficients of the SumProduct object !"
      WRITE(out_unit,*) "SumProduct%tab_opnd   = "//TO_string(SIZE(SumProduct%tab_opnd))
      WRITE(out_unit,*) "SumProduct%tab_coeffs = "//TO_string(SIZE(SumProduct%tab_coeffs))
      STOP "The number of products in the sum should match the number of coefficients of the SumProduct object !"
    END IF 

    !----------------------------Computation---------------------------------- 
    ALLOCATE(Op_psi_local(SIZE(Op_psi)))
    Op_psi       = ZERO
    Op_psi_local = ZERO

    DO i_product = 1, SIZE(SumProduct%tab_opnd)
      CALL Action(Op_psi_local, SumProduct%tab_opnd(i_product), Psi, Verbose=Verbose, Debug=Debug)
      Op_psi = Op_psi + Op_psi_local * SumProduct%tab_coeffs(i_product)
    END DO

    !--------------Conclusion----------------
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the ND Operator on the Psi statevector operand, computed &
                        &by MolecCav_Action_SOP_R1_real :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by MolecCav_Action_SOP_R1_real"
    END IF
  
    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "----------------------------------------ACTION OF THE ND OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_SOP_R1_real

  
  SUBROUTINE MolecCav_Action_SOP_R1_complex(Op_psi, SumProduct, Psi, Verbose, Debug) ! Psi is ND AND R1
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64 
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE

    complex(kind=Rkind),     intent(inout) :: Op_psi(:)
    TYPE(Sum_of_products_t), intent(in)    :: SumProduct
    complex(kind=Rkind),     intent(in)    :: Psi(:)
    integer, optional,       intent(in)    :: Verbose                                                                              ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                                ! cf. comments in HO1D_parameters_m

    integer                                :: i_product
    complex(kind=Rkind), allocatable       :: Op_psi_local(:)
    integer                                :: Verbose_local                                                                   ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "---------------------------------------COMPUTING ACTION OF THE SOP OVER &
                                              &THE R1 ND WF---------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Action_SOP_R1_complex :"
      WRITE(out_unit,*) "The <<SumProduct>> argument :"
      CALL Write(SumProduct)
      WRITE(out_unit,*) "The <<Psi>> argument : "
      CALL Write_Vec(Psi, out_unit, 1, info="Psi")
      WRITE(out_unit,*) "The size of its vector : "//TO_string(Size(Psi))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Action_SOP_R1_complex"
      FLUSH(out_unit)
    END IF
    
    !-----------------------------------------------------Checking dimensions----------------------------------------------------
    ! THE DIMENSIONS OF EACH 1D MATMUL WILL BE TESTED IN THE ACTIONS CODED IN ELEM_OP_M !

    IF (SIZE(SumProduct%tab_opnd) /= SIZE(SumProduct%tab_coeffs)) THEN
      WRITE(out_unit,*) "The number of products in the sum should match the number of coefficients of the SumProduct object !"
      WRITE(out_unit,*) "SumProduct%tab_opnd   = "//TO_string(SIZE(SumProduct%tab_opnd))
      WRITE(out_unit,*) "SumProduct%tab_coeffs = "//TO_string(SIZE(SumProduct%tab_coeffs))
      STOP "The number of products in the sum should match the number of coefficients of the SumProduct object !"
    END IF 

    !----------------------------Computation---------------------------------- 
    ALLOCATE(Op_psi_local(SIZE(Op_psi)))
    Op_psi       = ZERO
    Op_psi_local = ZERO

    DO i_product = 1, SIZE(SumProduct%tab_opnd)
      CALL Action(Op_psi_local, SumProduct%tab_opnd(i_product), Psi, Verbose=Verbose, Debug=Debug)
      Op_psi = Op_psi + Op_psi_local * SumProduct%tab_coeffs(i_product)
    END DO

    !--------------Conclusion----------------
    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Resulting statevector from the action of the Sum of products operator on the Psi statevector operand&
                        &, computed by MolecCav_Action_SOP_R1_complex :"
      CALL Write_Vec(Op_psi, out_unit, 1, info="Op_Psi")
      WRITE(out_unit,*) "--- End resulting statevector computed by MolecCav_Action_SOP_R1_complex"
    END IF
  
    IF (Debug_local) WRITE(out_unit,*) 
    IF (Debug_local) WRITE(out_unit,*) "----------------------------------------ACTION OF THE ND OPERATOR OVER THE R1 WF&
                                              & COMPUTED---------------------------------------"; FLUSH(out_unit)
  
  END SUBROUTINE MolecCav_Action_SOP_R1_complex

  
  SUBROUTINE MolecCav_Get_SOP_parameter_integer(Parameter_value, SumProduct, Parameter_name)
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE 

    integer,                 intent(inout) :: Parameter_value                                                                            ! the current values of the indexes for each dimension
    TYPE(Sum_of_products_t), intent(in)    :: SumProduct
    character(len=*),        intent(in)    :: Parameter_name

    IF (TO_lowercase(TRIM(Parameter_name)) == "n_product") THEN
      Parameter_value = SumProduct%N_products
    ELSE 
      WRITE(out_unit,*) "Parameter_name not recognized at MolecCav_Get_SOP_parameter_integer."
      STOP "Parameter_name not recognized at MolecCav_Get_SOP_parameter_integer"
    END IF

  END SUBROUTINE MolecCav_Get_SOP_parameter_integer


  SUBROUTINE MolecCav_Write_sum_of_products(SumProduct)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE 
    
    TYPE(Sum_of_products_t), intent(in) :: SumProduct

    integer                             :: i_product

    WRITE(out_unit,*) "_________________________________The Sum of products Operator object__________________________________"
    WRITE(out_unit,*) "|The number of products in the sum (SumProduct%N_products)                     | "//TO_string(SumProduct&
    &%N_products)
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    FLUSH(out_unit)
    IF (ALLOCATED(SumProduct%tab_opnd)) THEN
      WRITE(out_unit,*) "|The list of products (ND_operators_t) in the sum (SumProduct%tab_opnd) :      |"
      DO i_product = 1, SIZE(SumProduct%tab_opnd)
        WRITE(out_unit,*) "|"//TO_string(i_product)//"^{th} TERM OF THE SUM :"
        CALL Write(SumProduct%tab_opnd(i_product))
      END DO
    ELSE 
      WRITE(out_unit,*) "| The sum of products is not allocated (SumProduct%tab_opnd)                   |"
    END IF
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    FLUSH(out_unit)
    IF (ALLOCATED(SumProduct%tab_coeffs)) THEN
      WRITE(out_unit,*) "|The associated list of coefficients (SumProduct%tab_coeffs) :                 |"
      CALL Write_Vec(SumProduct%tab_coeffs, out_unit, SIZE(SumProduct%tab_coeffs), info="SumProduct%tab_coeffs")
    ELSE 
      WRITE(out_unit,*) "| The associated list of coefficients is not allocated (SumProduct%tab_coeffs) |"
    END IF
    WRITE(out_unit,*) "|______________________________________________________________________________|______________________"
    FLUSH(out_unit)
    WRITE(out_unit,*) "|_____________________________________End ND Operator object___________________|"
    FLUSH(out_unit)

  END SUBROUTINE MolecCav_Write_sum_of_products


  SUBROUTINE MolecCav_Deallocate_sum_of_products(SumProduct, Dealloc_all, Verbose, Debug)
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE 

    TYPE(Sum_of_products_t), intent(inout) :: SumProduct
    logical, optional,       intent(in)    :: Dealloc_all
    integer, optional,       intent(in)    :: Verbose                                                                                 ! cf. comments in HO1D_parameters_m
    logical, optional,       intent(in)    :: Debug                                                                                   ! cf. comments in HO1D_parameters_m

    integer                                :: i_product                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Dealloc_all_local
    integer                                :: Verbose_local                                                                      ! goes from 25 (= 0 verbose) to 29 (= maximum verbose) at this layer
    logical                                :: Debug_local

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Dealloc_all)) THEN; Dealloc_all_local = Dealloc_all
    ELSE; Dealloc_all_local = .FALSE.; END IF 
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Debug_local) THEN
      WRITE(out_unit,*) "--- The SumProduct to be deallocated :"
      CALL Write(SumProduct)
      WRITE(out_unit,*) "--- End SumProduct to be deallocated"
    END IF 

    !-----------------------------Deallocating the HO1D operator object----------------------------
    IF (Debug_local) WRITE(out_unit,*)
    IF (Debug_local) WRITE(out_unit,*) "-----------------------------------------------Deallocating the SumProduct obje&
    &----------------------------------------------"
  
    SumProduct%N_products = 0
    IF (ALLOCATED(SumProduct%tab_opnd)) THEN 
      DO i_product = 1, SIZE(SumProduct%tab_opnd)
        CALL Dealloc(SumProduct%tab_opnd(i_product), Dealloc_all=Dealloc_all_local, Verbose=Verbose_local, Debug=Debug_local)
      END DO
      DEALLOCATE(SumProduct%tab_opnd)
    END IF
    IF (ALLOCATED(SumProduct%tab_coeffs)) DEALLOCATE(SumProduct%tab_coeffs)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- The SumProduct object after having been deallocated :"
      CALL Write(SumProduct)
      WRITE(out_unit,*) "--- End dellocating SumProduct"
    END IF

  END SUBROUTINE MolecCav_Deallocate_sum_of_products


END MODULE
