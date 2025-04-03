# MolecCav : Structure proposition  
## Contents  
- [Current structure](#Current-structure)  
- [Proposed structure](#Proposed-structure)  
- [What changed in the new structure ?](#What-changed-in-the-new-structure-?)  
    - [```Quantum_HO1D_m```](#Quantum_HO1D_m)  
    - [```HO1D_parameters_m```](#HO1D_parameters_m)  
    - [```HO1D_operator_m```](#HO1D_operator_m)  
    - [```Actions_HO1D_operators_R1_m```](#Actions_HO1D_operators_R1_m)  
    - [```Cavity_mode_m```](#Cavity_mode_m)  
    - [```Action_cavity_operators```](#Action_cavity_operators)  
    - [```Matter_mode_m```](#Matter_mode_m)  
    - [```Action_matter_operators```](#Action_matter_operators)  
    - [```Total_hamiltonian```](#Total_hamiltonian)  
    - [```Action_total_hamiltonian```](#Action_total_hamiltonian)  
    - [```Initialize_total_ham_matrix_1p1D_m```](#Initialize_total_ham_matrix_1p1D_m)  
    - [```Transition_intensities```](#Transition_intensities)  
- [```Layered structure```](#Layered-structure)  


# Current structure  
```.
MolecCav
├── APP
│
├── DATA
│
├── Ext_Lib
│
├── OBJ
│
├── OUT
│
├── SRC
│   ├── Algebra_m.f90
│   ├── Cavity_mode_m.f90
│   ├── Mapping_m.f90
│   ├── ND_indexes_m.f90
│   ├── Elem_op_m.f90
│   ├── Operator_2D_m.f90
│   ├── Psi_analysis_m.f90
│   └── Total_hamiltonian_m.f90
│
├── TESTS
│
├── .gitignore 
├── makefile
├── MolecCav_Manual.pdf
├── README.md
└── run.sh
```  


# Proposed structure  
Are written in lowercasethe SRC/ modules that are not affected by the change, and in UPPERCASE are written the changed modules.  
Rq on the nomenclature :  
- A "_R1" can be added at the end of the name of a subroutine or a module related to an operator. This means that the procedure is is designed to compute the action of this operator upon a TENSOR OF RANK 1, regardless of the number of dimensions this tensor can be associated with.
- On the contrary, the addition of "$\_$1p1D" or "$\_ \langle n \rangle$D" after the name of an object means that THIS object (and not the tensor on which this object takes action) is meant to describe a system of the corresponding number of dimensions, regardless of the rank of any tensor used.

```.
MolecCav
├── APP
│
├── DATA
│
├── Ext_Lib
│
├── OBJ
│
├── OUT
│
├── SRC
│   ├── algebra_m.f90
│   ├── nd_indexes_m.f90
│   ├── psi_analysis_m.f90
│   │
│   ├── QUANTUM_HO1D_m.f90
│   ├── HO1D_PARAMETERS_m.f90
│   ├── HO1D_OPERATOR_m.f90
│   ├── ACTIONS_HO1D_OPERATORS_R1_m.f90
│   │
│   ├── CAVITY_MODE_m.f90
│   ├── ACTION_CAVITY_OPERATORS_m.f90
│   │
│   ├── MATTER_MODE_m.f90
│   ├── ACTION_MATTER_OPERATORS_m.f90
│   │
│   ├── TOTAL_HAMILTONIAN_m.f90
│   ├── ACTIONS_TOTAL_HAMILTONIAN_m.f90
│   ├── INITIALIZE_TOTAL_HAMILTONIAN_1P1D_m.f90
│   └── TRANSITION_INTENSITIES_m.f90
│
├── TESTS
│
├── .gitignore 
├── makefile
├── MolecCav_Manual.pdf
├── README.md
└── run.sh
```  


# What changed in the new structure ?  

## Quantum_HO1D_m  
The only module related to general HO that the others modules will need to call in a "USE".  

```fortran
MODULE Quantum_HO1D_m  
  USE HO1D_parameters_m  
  USE HO1D_operator_m  
  USE Actions_HO1D_operators_R1_m  
  
  TYPE, EXTENDS(HO1D_parameters_t) :: Quantum_HO1D_t  
    TYPE(HO1D_operator_t)          :: Hamiltonian  
    TYPE(HO1D_operator_t)          :: Position  
    TYPE(HO1D_operator_t)          :: NbQuanta  
  END TYPE  

  CONTAINS  

  SUBROUTINE Initialize_quantum_HO1D(HO1D, nio, arguments to construct the operators as before + & argument to choose which operator you was to construct in case you don't want all of them)  
    CALL READ_HO1D_parameters(HO1D%HO1D_parameters_t, nio)  
    CALL Initialize_HO1D_operator(Operator=HO1D%Hamiltonian, Operator_type="Hamiltonian", HO1D_para=HO1D%HO1D_parameters_t)  
    CALL Initialize_HO1D_operator(Operator=HO1D%Position, Operator_type="Position", HO1D_para=HO1D%HO1D_parameters_t)  
    CALL Initialize_HO1D_operator(Operator=HO1D%Nb_Quanta, Operator_type="NbQuanta", HO1D_para=HO1D%HO1D_parameters_t)  
  END SUB  

  SUBROUTINE Append_quantum_HO1D(HO1D, Operator_type="Hamiltonian")
    add an HO operator to a already initialized object of Quantum_HO1D_t type
    CALL Initialize_HO1D_operator(HO1D%<the chosen operator>, Operator_type="Hamiltonian" for instance, HO1D=HO1D%HO1D_parameters_t)  
  END SUB

  SUBROUTINE Write_quantum_HO1D(HO1D)  
    display values of the type in the output  
    CALL Write_HO1D_parameters(HO1D%HO1D_parameter_t)
    CALL Write_HO1D_operator(HO1D%Hamiltonian)
    CALL Write_HO1D_operator(HO1D%Position)
    CALL Write_HO1D_operator(HO1D%NbQuanta)
  END SUB  

  SUBROUTINE Deallocate_quantum_HO1D(HO1D)  
    deallocate all tables of the type  
    CALL deallocate_HO1D_operator(HO1D%Hamiltonian)
    CALL deallocate_HO1D_operator(HO1D%Position)
    CALL deallocate_HO1D_operator(HO1D%NbQuanta)
  END SUB  

```  


## HO1D_parameters_m  
The module to initialize the HO by reading its parameters from the namelist.    

```fortran  
MODULE HO1D_parameters_m  
  
  TYPE      :: HO1D_parameters_t  
    integer :: Nb      ! basis set size  
    real()  :: w       ! frequency  
    real()  :: m       ! mass  
    real()  :: eq_pos  ! equilibrium position  
  END TYPE  

  CONTAINS  

  SUBROUTINE Read_HO1D_parameters(HO1D_para, nio)  
    reads the namelist and initialize the type as before  
  END SUB  

  SUBROUTINE Write_HO1D_parameters(HO1D_para)  
    displays values of the type in the output  
  END SUB  
```


## HO1D_operator_m  
The module to initialize the operators related to a HO.  

```fortran  
MODULE HO1D_operator_m  
  
  TYPE                               :: HO1D_operator_t  
    character(len=:),    allocatable :: Operator_type  
    logical                          :: Dense = .FALSE.  
    integer                          :: Upper_bandwidth = 0  
    integer                          :: Lower_bandwidth = 0  
    real(kind=Rkind),    allocatable :: Dense_val(:,:)  
    real(kind=Rkind),    allocatable :: Diag_val(:)  
    real(kind=Rkind),    allocatable :: Band_val(:,:)  
  END TYPE  

  CONTAINS  

  SUBROUTINE Initialize_HO1D_operator(Operator, Operator_type, Dense, HO1D_para, Debug)  
    CALL Initialize_H_HO1D(Hamiltonian, HO1D_para) or  
    CALL Initialize_x_HO1D(PosOp, HO1D_para) or  
    CALL Initialize_N_HO1D(NbQuanta, HO1D_para)  
    constructs the operator as before but using parameters of the HO1D_para object from the above derived type instead of the old Cavity_mode_t. Calls the next three procedures. 
  END SUB

  SUBROUTINE Initialize_H_HO1D(Hamiltonian=Operator, HO1D_para)  
    constructs the H operator as before but using parameters of the HO1D_para object from the above derived type instead of the old Cavity_mode_t.  
  END SUB

  SUBROUTINE Initialize_x_HO1D(PosOp=Operator, HO1D_para)  
    constructs the x operator as before but using parameters of the HO1D_para object from the above derived type instead of the old Cavity_mode_t.  
  END SUB

  SUBROUTINE Initialize_N_HO1D(NbQuanta=Operator, HO1D_para)  
    constructs the N operator as before but using parameters of the HO1D_para object from the above derived type instead of the old Cavity_mode_t.  
  END SUB  

  SUBROUTINE Write_HO1D_operator(Operator)  
    displays values of the type in the output
  END SUB  

  SUBROUTINE Deallocate_HO1D_operator(Operator)  
    deallocates all tables of the type  
  END SUB  
```


## Actions_HO1D_operators_R1_m  
The module that accounts for the actions of the operators related to a HO over an any statevector (described by RANK-1 tensors) of this HO.  

```fortran  
MODULE Actions_HO1D_operator_R1_m  
  
  CONTAINS  

  SUBROUTINE Action_HO1D_operator_R1(Op_psi, Operator, Psi)  
    computes as before the resulting vector Op_psi(RANK-1 tensor) from the action of the operator of the HO on the state vector Psi(RANK-1 tensor) written in the Eigenbasis of the H of the HO1D.  
    CALL Action_dense_HO1D_operator_R1(Op_psi, Operator, Psi) or  
    CALL Action_diag_HO1D_operator_R1(Op_psi, Operator, Psi)  or  
    CALL Action_band_HO1D_operator_R1(Op_psi, Operator, Psi)  
  END SUB  

  SUBROUTINE Action_dense_HO1D_operator_R1(Op_psi, Operator, Psi) ! as before  
  SUBROUTINE Action_diag_HO1D_operator_R1(Op_psi, Operator, Psi)  ! as before  
  SUBROUTINE Action_band_HO1D_operator_R1(Op_psi, Operator, Psi)  ! as before  
  
  SUBROUTINE Average_value_HO1D_operator_R1(Value, Operator, Psi) ! as before  

  /!\ Duplicate each for Real and Complex state vectors (then generic procedures) /!\  
```


## Cavity_mode_m  
This module uses the Quantum_HO1D procedures and provides an overlayer or physics over the general 1D QHO to describe one cavity mode (1D subsystem then). It has "USE Action_cav_op_m" so only this module needs to be USEd in the main program.  

```fortran  
MODULE Cavity_mode_m  
  USE Quantum_HO1D_m  
  USE Action_cavity_operators_m  

  TYPE, EXTENDS(Quantum_HO1D_t) :: Cavity_mode_t  
    real(kind=Rkind),           :: lambda         ! coupling strenght parameter  
  END TYPE  

  CONTAINS  

  SUBROUTINE Read_cavity_mode(Mode%lambda, nio)  
    reads the namelist to get the lambda parameter  
    Mode%lambda = lambda  
  END SUB  

  SUBROUTINE Initialize_cavity_mode_from_scratch(Mode, and same arguments as Initialize_Quantum_HO1D)  
    if no Quantum HO has been initialized yet you can initialize both at the same time.  
    CALL Initialize_Quantum_HO1D(Mode%Quantum_HO1D, nio, arguments to construct the operators)  
    CALL Read_cavity_mode(Mode%lambda, nio)  
  END SUB  

  SUBROUTINE Initialize_cavity_mode_from_QHO1D(Mode, HO1D)  
    if a Quantum HO has already been initialized and you don't want/can't read the nml a second time. May be useful if you want to use the same Quantum HO to build several different mode (maybe).  
    Mode%Quantum_HO1D_t = HO1D  
    CALL Read_cavity_mode(Mode%lambda, nio)  
  END SUB  

  SUBROUTINE Write_cavity_mode(Mode)  
    displays values of the type in the output  
    CALL Write_quantum_HO1D(Mode%Quantum_HO1D)  
    WRITE Mode%lambda  
  END SUB  

  SUBROUTINE Deallocate_cavity_mode(Mode)  
    deallocates all tables of the type  
    CALL deallocate_quantum_HO1D(Mode%Quantum_HO1D)  
  END SUB  
```


## Action_cavity_operators_m  
The module that accounts for the actions of the operators related to a cavity mode over an any statevector (described by RANK-1 tensors) (of this HO or of multi-dim systems).  

```fortran  
MODULE Action_cavity_operators_m  

  CONTAINS  

  SUBROUTINE Action_cavity_mode_operators_1D_R1(Op_psi, Operator, Psi)  
    CALL Action_HO1D_operator_R1(Op_psi, Operator, Psi)  
  END SUB  
  ! Possible de mettre directement dans une interface même si vient pas du même module ?  

  SUBROUTINE Action_cavity_mode_operators_ND_R1(Op_psi, Operator, Psi, arguments pour manipuler ND_indexes)  
    loops with ND_indexes to take action over the aimed dimension of the Psi  
      CALL Action_HO1D_operator_R1(Op_psi, Operator, Psi)  
    end loops  
  END SUB  
```


## Matter_mode_m  
This module uses the Quantum_HO1D procedures and provides an overlayer or physics over the general 1D QHO to describe one matter mode (1D subsystem then). It has "USE Action_mat_op_m" so only this module needs to be USEd in the main program.  

```fortran  
MODULE Matter_mode_m  
  USE Quantum_HO1D_m  
  USE Action_matter_operators_m  

  TYPE, EXTENDS(Quantum_HO1D_t) :: Matter_mode_t  
    TYPE(HO1D_operator_t)       :: MatDipMomt     ! the dipole moment 1D operator  
    real()                      :: CoeffDipMomt   ! in the framework of our linear   behaviour assumption : the variation intentity of the dipole moment with the DOF/matter mode  
    real(kind=Rkind)            :: lambda         ! coupling strenght parameter  
  END TYPE  

  CONTAINS  

  SUBROUTINE Read_cavity_mode(Mode, nio)  
    reads the namelist to get the lambda parameter  
    Mode%CoeffDipMomt = CoeffDipMomt  
    Mode%lambda       = lambda  
  END SUB  

  SUBROUTINE Initialize_matter_mode_from_scratch(Mode and same arguments as Initialize_Quantum_HO1D)  
    if no Quantum HO has been initialized yet you can initialize both at the same time  
    CALL Initialize_Quantum_HO1D(HO1D, nio, arguments to construct the operators)  
    CALL Read_cavity_mode(Mode, nio)  
    CALL Initialize_HO1D_operator(Operator=Mode%MatDipMomt, Operator_type="Position", HO1D_para=Mode%HO1D_parameters_t)  
  END SUB  

  SUBROUTINE Initialize_cavity_mode_from_QHO1D(Mode, HO1D)  
    if a Quantum HO has already been initialized and you don't want/can't read the nml a second time. May be useful if you want to use the same Quantum HO to build several different mode (maybe).  
    Mode%Quantum_HO1D_t = HO1D  
    CALL Read_cavity_mode(Mode%lambda, nio)  
  END SUB  

  SUBROUTINE Write_matter_mode(Mode)  
    displays values of the type in the output  
    CALL Write_quantum_HO1D(Mode%Quantum_HO1D)  
    WRITE Mode%lambda etc  
  END SUB  

  SUBROUTINE Deallocate_matter_mode(Mode)  
    deallocates all tables of the type  
    CALL deallocate_quantum_HO1D(Mode%Quantum_HO1D)  
    DEALLOCATE(MatDipMomt)  
  END SUB  
```


## Action_matter_operators_m  
The module that accounts for the actions of the operators related to a matter mode over an any statevector (described by RANK-1 tensors) (of this HO or of multi-dim systems).  
  
```fortran  
MODULE Action_matter_operators_m  

  CONTAINS  

  SUBROUTINE Action_matter_mode_operators_1D_R1(Op_psi, Operator, Psi)  
    CALL Action_HO1D_operator_R1(Op_psi, Operator, Psi)  
  END SUB  
  ! Possible de mettre directement dans une interface même si vient pas du même module ?  
  
  SUBROUTINE Action_cavity_mode_operators_ND_R1(Op_psi, Operator, Psi, arguments pour manipuler ND_indexes)  
    loops with ND_indexes to take action over the aimed dimension of the Psi  
      CALL Action_HO1D_operator_R1(Op_psi, Operator, Psi)  
    end loops  
  END SUB  
```


## Total_hamiltonian_m  
This module uses the previous procedures for Cavity modes and Matter modes to describe the coupled MpND (M matter dimensions plus N Cavity dimensions but for now 1p1D) system. It is designed to work over Rank-1 statevectors that can represent a multi-dimension statefunction, using ND_indexes procedures. It gather modules to compute specific parts of the calculations, using "USE $\langle These modules\langle$_m" so only this module "summary module" needs to be USEd in the main program.  
  
```fortran  
MODULE Total_hamilonian_m  
  USE Action_total_hamiltonian_m  
  USE Initialize_total_ham_matrix_1p1D_m  
  USE Transition_intensities_m  
```

## Action_total_hamiltonian_R1_m  
A part of the Total_hamiltonian_m module that accounts for the action of the total hamiltonian over RANK-1 statevectors of any number of dimension (MpND but for now 1p1D).  
  
```fortran  
MODULE Action_total_hamiltonian_R1_m  

  CONTAINS  
  
  SUBROUTINE Action_total_hamiltonian_R1(TotH_psi, CavMode, MatMode, Psi, Debug)  
    CALL Action_matter_1p1D_R1(Psi_1, Psi_2, Psi_3, MatMode, Psi, Debug_local)  
    CALL Action_cavity_1p1D_R1(TotH_psi, CavMode, Psi_1, Psi_2, Psi_3, Debug_local)  
  END SUB  

  SUBROUTINE Action_matter_1p1D_R1(Psi_1, Psi_2, Psi_3, MatMode, Psi, Debug_local)      ! as before but over Rank-1 statevectors  
  SUBROUTINE Action_cavity_1p1D_R1(TotH_psi, CavMode, Psi_1, Psi_2, Psi_3, Debug_local) ! as before but over Rank-1 statevectors  

  SUBROUTINE Average_value_total_hamiltonian(Value, TotH, Psi)  
    as before but souldn't use the matrix representation of the total Hamiltonian here.  
  END SUB  
```

## Initialize_total_ham_matrix_1p1D_m  
A part of the Total_hamiltonian_m module that accounts for the construction of the total hamiltonian matrix over its eigenbasis for a 1p1D system (MpND but for now 1p1D).  
  
```fortran  
MODULE Initialize_total_ham_matrix_1p1D_m  

  CONTAINS  
  
  SUBROUTINE Initialize_total_ham_matrix_1p1D(TotH, CavMode, MatMode, Debug)  
    as before
  END SUB  
```


## Transition_intensities_m  
A part of the Total_hamiltonian_m module that accounts for the computation of the transition intensities matrix.  
  
```fortran  
MODULE Transition_intensities_m  

  CONTAINS  
  
  SUBROUTINE Transition_intensity(Intensity, InitPsi, MatDipMomt, FinPsi, Rank_size1, Rank_size2)   
    ! /!\ FOR NOW DESIGNED FOR 1p1D
    as before
  END SUB  

  SUBROUTINE Initialize_transition_matrix(Intensities, MatDipMomt, REigvec, Energy_threshold, REigval, Nb_states, Debug)   
    ! /!\ FOR NOW DESIGNED FOR 1p1D
    as before
  END SUB

  SUBROUTINE MolecCav_Compute_transition_matrix(Intensities, MatDipMomt, REigvec, Debug)   
    ! /!\ FOR NOW DESIGNED FOR 1p1D
    as before
  END SUB
```

# Layered structure
```.
App_MolecCav (Layer 0)
├─ Total_hamiltonian_m ── Actions_total_hamiltonian_m ────────── Matter_mode_m ─────────── Quantum_HO1D_m ── HO1D_parameters_m ───────── Algebra_m
├─ Psi_analysis_m      ├─ Initialize_total_hamiltonian_1P1D_m ├─ Action_matter_operators_m                ├─ HO1D_operator_m
│                      └─ Transition_intensities_m            ├─ Cavity_mode_m                            └─ Actions_HO1D_operators_R1_m
│                                                             ├─ Action_cavity_oprators_m
│                                                             └─ ND_indexes_m
│
└─ Layer 1 ─────────────└ Layer 2 ─────────────────────────────└ Layer 3 ────────────────└ Layer 4 ────────└ Layer 5 ──────────────────└ Layer 6
```  
