# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
physica_add_test(BatchedGEMM MatrixProduct/BatchedGEMM.cpp)
physica_add_test(GEMM MatrixProduct/GEMM.cpp)
physica_add_test(GEMV MatrixProduct/GEMV.cpp)
physica_add_test(GEVM MatrixProduct/GEVM.cpp)
physica_add_test(Kronecker MatrixProduct/Kronecker.cpp)
if(${PHYSICA_CUDA})
    physica_add_test(BatchedGEMM_cuda MatrixProduct/BatchedGEMM.cu)
    physica_add_test(GEMM_cuda MatrixProduct/GEMM.cu)
    physica_add_test(GEMV_cuda MatrixProduct/GEMV.cu)
    physica_add_test(Kronecker_cuda MatrixProduct/Kronecker.cu)
endif()

physica_add_test(Balance Balance.cpp)
physica_add_test(Det Det.cpp)
physica_add_test(Inverse Inverse.cpp)
physica_add_test(LValueMatrix LValueMatrix.cpp)
physica_add_test(MatrixExpr MatrixExpr.cpp)
physica_add_test(MatrixTrig MatrixTrig.cpp)
physica_add_test(RValueMatrix RValueMatrix.cpp)
if(${PHYSICA_CUDA})
    physica_add_test(LValueMatrix_cuda LValueMatrix.cu)
    physica_add_test(MatrixTrig_cuda MatrixTrig.cu)
    physica_add_test(RValueMatrix_cuda RValueMatrix.cu)
endif()
