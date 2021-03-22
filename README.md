# Nth Order Analytical Time Derivatives of Inverse Dynamics in Recursive and Closed Forms

This repository contains the MATLAB code for computing nth order analytical time derivatives of Inverse Dynamics in recursive and closed forms using body fixed representation of the twists. This repository accompanies the ICRA 2021 paper: Shivesh Kumar, Andreas Mueller, Nth Order Analytical Time Derivatives of Inverse Dynamics in Recursive and Closed Forms. IEEE ICRA 2021, Xi'an, China.

# Usage instructions:

# Main scripts
* Panda_NthOrder_InvDyn_BodyFixed_ClosedForm.m: nth order time derivatives in closed form which uses robot parameters from Franka Emika Panda robot
* Panda_NthOrder_InvDyn_BodyFixed_RecursiveForm.m: nth order time derivatives in recursive form which uses robot parameters from Franka Emika Panda robot

# Helper functions
* NthOrder_ClosedFormInvDyn_BodyFixed.m: Function to compute Nth order inverse dynamics in closed form
* NthOrder_RecursiveInvDyn_BodyFixed.m: Function to compute Nth order inverse dynamics in recursive form
* SE3Exp.m: Function to compute exponential mapping for SE(3)
* SO3Exp.m: Function to compute exponential mapping for SO(3)
* SE3Inv.m: Function to compute analytical inverse of exponential mapping for SE(3)
* SE3AdjMatrix.m: Function to compute (6x6) Adjoint Matrix for SE(3)
* SE3adjMatrix.m: Function to compute (6x6) adjoint Matrix for SE(3) - also known as spatial cross product in the literature
* SE3AdjInvMatrix.m: Function to compute Inverse of (6x6) Adjoint Matrix for SE(3)
* MassMatrixMixedData.m: Function to build mass-inertia matrix in SE(3) from mass, inertia and center of mass information

# 2nd order inverse dynamics:
* Panda_InvDyn_BodyFixed_ClosedForm.m: 2nd order time derivatives in closed form which uses robot parameters from Franka Emika Panda robot
* Panda_InvDyn_BodyFixed.m: 2nd order time derivatives in recursive form which uses robot parameters from Franka Emika Panda robot
* InvDyn_BodyFixed.m: Function to compute 2nd order inverse dynamics in recursive form
* ClosedFormInvDyn_BodyFixed.m: Function to compute 2nd order inverse dynamics in closed form

# Citation
Shivesh Kumar, Andreas Mueller, Nth Order Analytical Time Derivatives of Inverse Dynamics in Recursive and Closed Forms, In: IEEE International Conference on Robotics and Automation (ICRA) 2021, Xi'an, China 2021. 
```
@inproceedings{2021_Kumar_and_Mueller_NthOrderEOMDerivatives_ICRA,
author={Kumar, Shivesh
and Mueller, Andreas},
title={Nth Order Analytical Time Derivatives of Inverse Dynamics in Recursive and Closed Forms},
Booktitle = {IEEE International Conference on Robotics and Automation (ICRA) 2021},
Organization={IEEE},
}
```

