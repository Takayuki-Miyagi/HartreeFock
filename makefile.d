obj/half_precision_floating_points.o : src/half_precision_floating_points.f90 
obj/ThreeBodyInteraction.o : src/ThreeBodyInteraction.F90 obj/TwoBodyOperator.o obj/ThreeBodyModelSpace.o obj/SingleParticleState.o obj/myfort.o 
obj/DefineOperators.o : src/DefineOperators.F90 obj/myfort.o 
obj/HFMBPT.o : src/HFMBPT.F90 obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/myfort.o 
obj/ThreeBodyModelSpace.o : src/ThreeBodyModelSpace.F90 obj/SingleParticleState.o obj/myfort.o 
obj/ModelSpace.o : src/ModelSpace.F90 obj/ThreeBodyModelSpace.o obj/TwoBodyModelSpace.o obj/OneBodyModelSpace.o obj/SingleParticleState.o obj/myfort.o 
obj/ThreeBodyMonInteraction.o : src/ThreeBodyMonInteraction.F90 obj/SingleParticleState.o obj/myfort.o 
obj/OneBodyOperator.o : src/OneBodyOperator.F90 obj/DefineOperators.o obj/OneBodyModelSpace.o obj/myfort.o 
obj/SingleParticleState.o : src/SingleParticleState.F90 obj/myfort.o 
obj/ThreeBodyOperator.o : src/ThreeBodyOperator.F90 obj/ThreeBodyInteraction.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ThreeBodyModelSpace.o obj/myfort.o 
obj/TwoBodyOperator.o : src/TwoBodyOperator.F90 obj/DefineOperators.o obj/OneBodyOperator.o obj/TwoBodyModelSpace.o obj/myfort.o 
obj/OneBodyModelSpace.o : src/OneBodyModelSpace.F90 obj/SingleParticleState.o 
obj/Operators.o : src/Operators.F90 obj/DefineOperators.o obj/ThreeBodyNO2BInteraction.o obj/ThreeBodyMonInteraction.o obj/ThreeBodyInteraction.o obj/ThreeBodyOperator.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ModelSpace.o obj/myfort.o 
obj/HartreeFock.o : src/HartreeFock.F90 obj/ThreeBodyNO2BInteraction.o obj/ThreeBodyMonInteraction.o obj/Operators.o obj/myfort.o 
obj/ThreeBodyNO2BInteraction.o : src/ThreeBodyNO2BInteraction.F90 obj/TwoBodyOperator.o obj/half_precision_floating_points.o obj/SingleParticleState.o obj/myfort.o 
obj/TwoBodyModelSpace.o : src/TwoBodyModelSpace.F90 obj/SingleParticleState.o obj/myfort.o 
obj/profiler.o : submodule/myfort/src/profiler.f90 obj/general.o 
obj/mat_vec_single.o : submodule/myfort/src/mat_vec_single.f90 obj/matrix_single.o obj/vector_single.o obj/general.o 
obj/wave_functions.o : submodule/myfort/src/wave_functions.f90 obj/general.o obj/physics_constants.o obj/functions_from_c.o 
obj/vector_double.o : submodule/myfort/src/vector_double.f90 obj/general.o 
obj/mat_vec_complex.o : submodule/myfort/src/mat_vec_complex.f90 obj/matrix_complex.o obj/vector_complex.o obj/general.o 
obj/sngl_dble_cmplx.o : submodule/myfort/src/sngl_dble_cmplx.f90 obj/matrix_complex.o obj/vector_complex.o obj/matrix_double.o obj/vector_double.o obj/matrix_single.o obj/vector_single.o 
obj/matrix_double.o : submodule/myfort/src/matrix_double.f90 obj/vector_double.o obj/general.o 
obj/vector_single.o : submodule/myfort/src/vector_single.f90 obj/general.o 
obj/angular_momentum_couplings.o : submodule/myfort/src/angular_momentum_couplings.f90 obj/functions_from_c.o obj/general.o 
obj/vector_complex.o : submodule/myfort/src/vector_complex.f90 obj/general.o 
obj/mat_vec_double.o : submodule/myfort/src/mat_vec_double.f90 obj/matrix_double.o obj/vector_double.o obj/general.o 
obj/matrix_single.o : submodule/myfort/src/matrix_single.f90 obj/vector_single.o obj/general.o 
obj/store_couplings.o : submodule/myfort/src/store_couplings.f90 obj/angular_momentum_couplings.o obj/functions_from_c.o obj/profiler.o obj/general.o 
obj/linear_algebra.o : submodule/myfort/src/linear_algebra.f90 obj/mat_vec_complex.o obj/mat_vec_double.o obj/mat_vec_single.o obj/matrix_complex.o obj/matrix_double.o obj/matrix_single.o obj/vector_complex.o obj/vector_double.o obj/vector_single.o obj/sngl_dble_cmplx.o obj/general.o 
obj/physics_constants.o : submodule/myfort/src/physics_constants.f90 obj/general.o 
obj/iteration_methods.o : submodule/myfort/src/iteration_methods.f90 obj/linear_algebra.o 
obj/functions_from_c.o : submodule/myfort/src/functions_from_c.f90 
obj/matrix_complex.o : submodule/myfort/src/matrix_complex.f90 obj/vector_complex.o obj/general.o 
obj/myfort.o : submodule/myfort/src/myfort.f90 obj/iteration_methods.o obj/linear_algebra.o obj/wave_functions.o obj/profiler.o obj/physics_constants.o obj/general.o obj/store_couplings.o obj/angular_momentum_couplings.o obj/functions_from_c.o 
obj/general.o : submodule/myfort/src/general.f90 
obj/WriteOperator.o : main/WriteOperator.F90 obj/Operators.o obj/HFInput.o obj/myfort.o 
obj/Atomic.o : main/Atomic.F90 obj/HFMBPT.o obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o obj/myfort.o 
obj/HFInput.o : main/HFInput.F90 obj/myfort.o 
obj/HFMain.o : main/HFMain.F90 obj/Atomic.o obj/WriteOperator.o obj/HFMBPT.o obj/HartreeFock.o obj/ThreeBodyMonInteraction.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o obj/myfort.o 
