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
    integer                          :: N_products                             ! the number of terms in the sum of product operator
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

    ! integer                            :: i_mode, i_op
    logical                            :: Dense_local                                                                  ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    integer                            :: Verbose_local                                                                ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                            :: Debug_local
    
    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) THEN; Verbose_local = Verbose
    ELSE; Verbose_local = 20; END IF 
    IF (PRESENT(Debug))   THEN; Debug_local   = Debug
    ELSE; Debug_local = .FALSE.; END IF

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE OPERATOR_ND OBJE&
                                              &CT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_Sum_of_product :"
      WRITE(out_unit,*) "The <<OpND>> argument :"
      CALL Write(OpND)
      WRITE(out_unit,*) "The <<Mat_operators>>  argument :"//Mat_operators
      WRITE(out_unit,*) "The <<Cav_operators>>  argument :"//Cav_operators
      IF (PRESENT(Dense)) WRITE(out_unit,*) "The <<Dense>> argument : "//TO_string(Dense)
      WRITE(out_unit,*) "Are the module's <<tab_mat/cav_ops>> allocated ? "//TO_string(ALLOCATED(tab_mat_ops))//TO_string(ALLOCAT&
      &ED(tab_cav_ops))
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_Sum_of_product"
      FLUSH(out_unit)
    END IF
    
    !------------------------------------------Initializing the procedure------------------------------------------
    IF (PRESENT(Dense)) THEN; Dense_local = Dense
    ELSE; Dense_local = .FALSE.; END IF

    IF (LEN_TRIM(Mat_operators)==0) THEN
      N_mat = 0
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
      WRITE(out_unit,*) "                          The code is now used without any matter mode to compute : cavity alone "
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
    ELSE
      N_mat = 1
      DO i_mode = 1, LEN_TRIM(Mat_operators)
        IF (Mat_operators(i_mode:i_mode)==',') N_mat = N_mat + 1
      END DO 
    END IF 
    IF (LEN_TRIM(Cav_operators)==0) THEN
      N_cav = 0
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
      WRITE(out_unit,*) "                          The code is now used without any cavity mode to compute : matter alone "
      WRITE(out_unit,*) "########################## WARNING ########################## WARNING ########################## WARNING #&
                      &########################"
    ELSE
      N_cav = 1
      DO i_mode = 1, LEN_TRIM(Cav_operators)
        IF (Cav_operators(i_mode:i_mode)==',') N_cav = N_cav + 1
      END DO 
    END IF 
    
    IF (ALLOCATED(tab_mat_ops) .AND. ALLOCATED(tab_cav_ops)) THEN
      IF (SIZE(tab_mat_ops)/=N_mat) THEN
        WRITE(out_unit,*) "### The number of declared matter modes is not consistent with the one previously declared at last cal&
        &l of Initialize_operator_ND."
        WRITE(out_unit,*) "    Please check arguments to keep consistency in the system's size along the simulation."
        STOP "### The number of declared matter modes is not consistent with the one previously declared."
      ELSE IF (SIZE(tab_cav_ops)/=N_cav) THEN
        WRITE(out_unit,*) "### The number of declared cavity modes is not consistent with the one previously declared at last cal&
        &l of Initialize_operator_ND."
        WRITE(out_unit,*) "    Please check arguments to keep consistency in the system's size along the simulation."
        STOP "### The number of declared cavity modes is not consistent with the one previously declared."
      END IF
    END IF

    IF (.NOT. ALLOCATED(tab_mat_ops) .OR. .NOT. ALLOCATED(tab_cav_ops)) THEN
      CALL Initialize_tabs_ops(N_mat=N_mat, N_cav=N_cav, Dense=Dense_local, Verbose=Verbose_local, Debug=Debug_local)
    END IF
    
    !--------------------------------------Constructing the OpND = parsing the Mat/Cav_operators strings-------------------------------------
    ALLOCATE(Mat_operators_local(N_mat))
    ALLOCATE(Cav_operators_local(N_cav))
    ALLOCATE(OpND%tab_indexes_mat_op(N_mat))
    ALLOCATE(OpND%tab_indexes_cav_op(N_cav))

    READ(unit=Mat_operators, fmt=*) Mat_operators_local ! /!\ neither trimmed nor lowercased so far /!\
    READ(unit=Cav_operators, fmt=*) Cav_operators_local

    DO i_mode = 1, N_mat
      DO i_op = 0, tab_mat_ops(i_mode)%Nb_op-1
        IF (TO_lowercase(TRIM(Mat_operators_local(i_mode))) == tab_mat_ops(i_mode)%Tab_op(i_op)%Operator_type) THEN
          OpND%tab_indexes_mat_op(i_mode) = i_op
          EXIT 
        ELSE IF (i_op == tab_mat_ops(i_mode)%Nb_op-1) THEN
          WRITE(out_unit,*) "### No operator name recognized. Please check arguments of MolecCav_Initialize_operator_ND"
          STOP "### No operator name recognized in MolecCav_Initialize_operator_ND"
        END IF 
      END DO
    END DO 

    DO i_mode = 1, N_cav
      DO i_op = 0, tab_cav_ops(i_mode)%Nb_op-1
        IF (TO_lowercase(TRIM(Cav_operators_local(i_mode))) == tab_cav_ops(i_mode)%Tab_op(i_op)%Operator_type) THEN
          OpND%tab_indexes_cav_op(i_mode) = i_op
          EXIT 
        ELSE IF (i_op == tab_cav_ops(i_mode)%Nb_op-1) THEN
          WRITE(out_unit,*) "### No operator name recognized. Please check arguments of MolecCav_Initialize_operator_ND"
          STOP "### No operator name recognized in MolecCav_Initialize_operator_ND"
        END IF 
      END DO
    END DO

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "--------------------------------------------------OPERATOR_ND OBJECT INITIALIZED--&
    &-----------------------------------------------"; FLUSH(out_unit)

  END SUBROUTINE MolecCav_Initialize_Sum_of_product


END MODULE
