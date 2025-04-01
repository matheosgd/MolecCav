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
! Initialize_quantum_HO1D : reads the namelist and initialize the type, then constructs the operat-
! or using parameters of the HO1D_para object from the so called derived type.
! Append_quantum_HO1D     : add an HO operator to a already initialized object of Quantum_HO1D_t type.
! Write_quantum_HO1D      : display values of the type in the output
! Deallocate_quantum_HO1D : deallocate all tables of the type
!==================================================================================================
!==================================================================================================
MODULE Quantum_HO1D_m
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
  USE QDUtil_m
  USE HO1D_parameters_m
  USE HO1D_operator_m
  USE Actions_HO1D_operators_R1_m
  IMPLICIT NONE


  TYPE, EXTENDS(HO1D_parameters_t) :: Quantum_HO1D_t
    TYPE(HO1D_operator_t)          :: Hamiltonian  
    TYPE(HO1D_operator_t)          :: Position  
    TYPE(HO1D_operator_t)          :: NbQuanta  
  END TYPE


  PRIVATE

  PUBLIC Quantum_HO1D_t, Initialize_quantum_HO1D, Append_quantum_HO1D, Write_quantum_HO1D, Deallocate_quantum_HO1D

  INTERFACE Initialize_quantum_HO1D
    MODULE PROCEDURE MolecCav_Initialize_quantum_HO1D
  END INTERFACE
  INTERFACE Append_quantum_HO1D
    MODULE PROCEDURE MolecCav_Append_quantum_HO1D
  END INTERFACE
  INTERFACE Write_quantum_HO1D
    MODULE PROCEDURE MolecCav_Write_quantum_HO1D
  END INTERFACE
  INTERFACE Deallocate_quantum_HO1D
    MODULE PROCEDURE MolecCav_Deallocate_quantum_HO1D
  END INTERFACE
    

  CONTAINS


  SUBROUTINE Initialize_quantum_HO1D(HO1D, nio, Which_op, Verbose, Debug)
    !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64
    USE QDUtil_m
    USE HO1D_parameters_m
    USE HO1D_operator_m
    IMPLICIT NONE
  
    TYPE(Quantum_HO1D_t),       intent(inout) :: HO1D
    integer,                    intent(in)    :: nio                                                                             ! cf. comments in HO1D_parameters_m
    character(len=3), optional, intent(in)    :: Which_op                                                                        ! a string of three characters with the initial of the operators to construct. ex : "H.." constructs only the Hamiltonian, ".xN" constructs the position and Nb of quanta operators.
    integer, optional,          intent(in)    :: Verbose                                                                         ! cf. comments in HO1D_parameters_m
    logical, optional,          intent(in)    :: Debug                                                                           ! cf. comments in HO1D_parameters_m

    character                                 :: Operator_1, Operator_2, Operator_3
    integer                                   :: Verbose_local = 20                                                              ! goes from 20 (= 0 verbose) to 24 (= maximum verbose) at this layer
    logical                                   :: Debug_local   = .FALSE.
    
    !------------------------------------------------------Debugging options-----------------------------------------------------
    IF (PRESENT(Verbose)) Verbose_local = Verbose
    IF (PRESENT(Debug))   Debug_local   = Debug

    IF (Verbose_local > 20) WRITE(out_unit,*) 
    IF (Verbose_local > 20) WRITE(out_unit,*) "-------------------------------------------------INITIALIZING THE QUANTUM HO1D OBJ&
                                              &ECT-------------------------------------------------"; FLUSH(out_unit)

    IF (Debug_local) THEN
      WRITE(out_unit,*)
      WRITE(out_unit,*) "--- Arguments of MolecCav_Initialize_quantum_HO1D :"
      WRITE(out_unit,*) "The <<HO1D>> argument :"
      CALL Write_quantum_HO1D(HO1D)
      IF (PRESENT(Which_op)) THEN
        WRITE(out_unit,*) "The argument listing the operators to build within the quantum HO1D : "//Which_op
      ELSE 
        WRITE(out_unit,*) "No argument listing the operators to build has been provided. All will be"
      END IF
      WRITE(out_unit,*) "--- End arguments of MolecCav_Initialize_quantum_HO1D"
      FLUSH(out_unit)
    END IF
    
    !----------------------------------------------Selecting the operators to build----------------------------------------------
    IF (PRESENT(Which_op)) THEN
      Operator_1 = TO_lowercase(Which_op(1))
      Operator_2 = TO_lowercase(Which_op(2))
      Operator_3 = TO_lowercase(Which_op(3))
      IF (Debug_local) THEN
        WRITE(out_unit,*) "The argument listing the operators to build has been parsed into :" 
        WRITE(out_unit,*) "Operator_1 = "//Operator_1//"; Operator_2 = "//Operator_2//"; Operator_3 = "//Operator_3
      END IF 
    END IF 

    CALL READ_HO1D_parameters(HO1D_para=HO1D%HO1D_parameters_t, nio=nio, Verbose=Verbose_local, Debug=Debug_local)
    CALL Initialize_HO1D_operator(Operator=HO1D%Hamiltonian, Operator_type="Hamiltonian", HO1D_para=HO1D%HO1D_parameters_t, Verbo&
                                 &se=Verbose_local, Debug=Debug_local)
    CALL Initialize_HO1D_operator(Operator=HO1D%Position,    Operator_type="Position",    HO1D_para=HO1D%HO1D_parameters_t, Verbo&
                                 &se=Verbose_local, Debug=Debug_local)
    CALL Initialize_HO1D_operator(Operator=HO1D%Nb_Quanta,   Operator_type="NbQuanta",    HO1D_para=HO1D%HO1D_parameters_t, Verbo&
                                 &se=Verbose_local)

  END SUBROUTINE Initialize_quantum_HO1D


  SUBROUTINE Append_quantum_HO1D(HO1D, Operator_type)
    CALL Initialize_HO1D_operator(HO1D%<the chosen operator>, Operator_type="Hamiltonian" for instance, HO1D=HO1D%HO1D_parameters_t)  
    HO1D%<the chosen operator> = Operator
  END SUBROUTINE Append_quantum_HO1D

  SUBROUTINE Write_quantum_HO1D(HO1D)  
    display values of the type in the output  
    CALL Write_HO1D_parameters(HO1D%HO1D_parameter_t)
    CALL Write_HO1D_operator(HO1D%Hamiltonian)
    CALL Write_HO1D_operator(HO1D%Position)
    CALL Write_HO1D_operator(HO1D%NbQuanta)
  END SUBROUTINE Write_quantum_HO1D  

  SUBROUTINE Deallocate_quantum_HO1D(HO1D)  
    deallocate all tables of the type  
    CALL deallocate_HO1D_operator(HO1D%Hamiltonian)
    CALL deallocate_HO1D_operator(HO1D%Position)
    CALL deallocate_HO1D_operator(HO1D%NbQuanta)
  END SUBROUTINE Deallocate_quantum_HO1D  


END MODULE
