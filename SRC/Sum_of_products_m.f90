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
    integer                          :: N_products = 0                         ! the number of terms in the sum of product operator
    TYPE(Operator_ND_t), allocatable :: tab_opnd(:)                            ! the list of the product terms in the sum, which one being an OpND
  END TYPE


  PRIVATE

  PUBLIC Sum_of_products_t!, Initialize, Action, Get, Write, Dealloc

!   INTERFACE Initialize
!     MODULE PROCEDURE MolecCav_Initialize_operator_ND
!   END INTERFACE
!   INTERFACE Initialize_tabs_ops
!     MODULE PROCEDURE MolecCav_Initialize_tabs_operators
!   END INTERFACE
!   INTERFACE Action
!     MODULE PROCEDURE MolecCav_Action_operator_ND_R1_real, MolecCav_Action_operator_ND_R1_complex
!   END INTERFACE
!   INTERFACE Get
!     MODULE PROCEDURE MolecCav_Get_OpND_parameter_integer, MolecCav_Get_OpND_parameter_real
!   END INTERFACE
!   INTERFACE Write
!     MODULE PROCEDURE MolecCav_Write_operator_ND
!   END INTERFACE
!   INTERFACE Dealloc
!     MODULE PROCEDURE MolecCav_Deallocate_operator_ND
!   END INTERFACE
    

  CONTAINS


  SUBROUTINE MolecCav_Initialize_Sum_of_product(SumProduct, nio, Dense, Verbose, Debug) ! no need for N_mat and N_cav explicitly : they are SIZE(Mat_op and Cav_op)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE Operator_ND_m
    IMPLICIT NONE
  
    TYPE(Sum_of_products_t), intent(inout) :: SumProduct
    integer,             intent(in)        :: nio
    logical, optional,   intent(in)        :: Dense                                                                        ! cf. comments in HO1D_parameters_m
    integer, optional,   intent(in)        :: Verbose                                                                      ! cf. comments in HO1D_parameters_m
    logical, optional,   intent(in)        :: Debug                                                                        ! cf. comments in HO1D_parameters_m

    integer                                :: N_products, i_product
    character(len=:), allocatable          :: Mat_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_mat>', ex : 'hamiltonian, Identity'. Not case sensitive, ' ' <=> \otimes
    character(len=:), allocatable          :: Cav_operators ! syntax : '<op_mode_1>, <op_mode_2>, ..., <op_mode_N_cav>', ex : 'hamiltonian'.   Not case sensitive, ' ' <=> \otimes. This exemple means OpND = H_mat1\otimesI_mat2\otimesH_cav
    integer                                :: err_io
    logical                                :: Dense_local                                                                  ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                                :: Verbose_local                                                                ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                                :: Debug_local
    
    NAMELIST /Sum_of_products/ N_products

    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE SUM OF PRODUCT O&
                                              &PERATOR OBJECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_Sum_of_product :"
      WRITE(out_unit,*) "The <<SumProduct>> argument :"
      ! CALL Write(SumProduct)
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_Sum_of_product"
      FLUSH(out_unit)
    END IF
    
    !----------------------Initialization to default values--------------------
    N_products = 0

    !------------------------------Reading of the nml--------------------------
    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** READING THE NUMBER OF OPND **************************'
    WRITE(out_unit,*) '********************************************************************************'
    
    READ(nio, nml = Sum_of_products, iostat = err_io)                                     ! assign the values read in the nml to the declared list of parameters

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "-----------------------The namelist parameters are read as----------------------"
      WRITE(out_unit, nml = Sum_of_products)
      WRITE(out_unit,*) "-------------------------End of the namelist parameters-------------------------"
    END IF
    
      !------------------------------Check reading error-------------------------
    IF(err_io /= 0) THEN
      WRITE(out_unit,*) ''
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '##### Error in MolecCav_Initialize_Sum_of_product (err_io/=0) #####'
      WRITE(out_unit,*) '###################################################################'
      WRITE(out_unit,*) '####################### err_io = ', err_io, '######################'
      STOP '######################### Check basis data ########################'
    END IF

    IF (N_products == 0) THEN
      WRITE(out_unit,*) "### The number of ND Operators in the sum of product operator CANNOT be 0 (what are are you going to stu&
      &dy if there is no operator ???). Please check the data file '.nml'"
      STOP "### The number of ND Operators in the sum of product operator CANNOT be 0 (what are are you going to study if there i&
      &s no operator ???). Please check the data file '.nml'"
    END IF
    
    !---------------Construction of the table of OpND composing the sum of product operator-----------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    SumProduct%N_products = N_products 
    ALLOCATE(SumProduct%tab_opnd(SumProduct%N_products))

    DO i_product = 1, N_products
      CALL Read_product(Mat_operators, Cav_operators, nio, Verbose=Verbose_local, Debug=Debug_local)
      CALL Initialize(SumProduct%tab_opnd(i_product), Mat_operators, Cav_operators, nio, Dense_local, Verbose_local, Debug_local)
      DEALLOCATE(Mat_operators); DEALLOCATE(Cav_operators)
    END DO

    WRITE(out_unit,*) 
    WRITE(out_unit,*) '********************************************************************************'
    WRITE(out_unit,*) '************************** SUM OF PRODUCTS CONSTRUCTED *************************'
    WRITE(out_unit,*) '********************************************************************************'

    IF (Debug) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--------------Sum of products constructed by MolecCav_Initialize_Sum_of_product--------------"
      ! CALL Write(SumProduct)
      WRITE(out_unit,*) "------------End Sum of products constructed by MolecCav_Initialize_Sum_of_product------------"
    END IF
    
    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "--------------------------------------------------SUM OF PRODUCTS OPERATOR INITIAL&
    &IZED-------------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_Sum_of_product


END MODULE
