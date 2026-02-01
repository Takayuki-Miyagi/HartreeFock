obj/half_precision_floating_points.o : src/half_precision_floating_points.f90 
obj/DefineOperators.o : src/DefineOperators.F90 obj/myfort.o 
obj/HFMBPT.o : src/HFMBPT.F90 obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/myfort.o 
obj/HartreeFock.o : src/HartreeFock.F90 obj/ThreeBodyNO2BInteraction.o obj/ThreeBodyMonInteraction.o obj/Operators.o obj/myfort.o 
obj/ModelSpace.o : src/ModelSpace.F90 obj/ThreeBodyModelSpace.o obj/TwoBodyModelSpace.o obj/OneBodyModelSpace.o obj/SingleParticleState.o obj/myfort.o 
obj/OneBodyModelSpace.o : src/OneBodyModelSpace.F90 obj/SingleParticleState.o 
obj/OneBodyOperator.o : src/OneBodyOperator.F90 obj/DefineOperators.o obj/OneBodyModelSpace.o obj/myfort.o 
obj/Operators.o : src/Operators.F90 obj/DefineOperators.o obj/ThreeBodyNO2BInteraction.o obj/ThreeBodyMonInteraction.o obj/ThreeBodyInteraction.o obj/ThreeBodyOperator.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ModelSpace.o obj/myfort.o 
obj/SingleParticleState.o : src/SingleParticleState.F90 obj/myfort.o 
obj/ThreeBodyInteraction.o : src/ThreeBodyInteraction.F90 obj/TwoBodyOperator.o obj/ThreeBodyModelSpace.o obj/SingleParticleState.o obj/myfort.o 
obj/ThreeBodyModelSpace.o : src/ThreeBodyModelSpace.F90 obj/SingleParticleState.o obj/myfort.o 
obj/ThreeBodyMonInteraction.o : src/ThreeBodyMonInteraction.F90 obj/SingleParticleState.o obj/myfort.o 
obj/ThreeBodyNO2BInteraction.o : src/ThreeBodyNO2BInteraction.F90 obj/TwoBodyOperator.o obj/half_precision_floating_points.o obj/SingleParticleState.o obj/myfort.o 
obj/ThreeBodyOperator.o : src/ThreeBodyOperator.F90 obj/ThreeBodyInteraction.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ThreeBodyModelSpace.o obj/myfort.o 
obj/TwoBodyModelSpace.o : src/TwoBodyModelSpace.F90 obj/SingleParticleState.o obj/myfort.o 
obj/TwoBodyOperator.o : src/TwoBodyOperator.F90 obj/DefineOperators.o obj/OneBodyOperator.o obj/TwoBodyModelSpace.o obj/myfort.o 
obj/angular_momentum_couplings.o : submodule/myfort/src/angular_momentum_couplings.f90 obj/functions_from_c.o 
obj/functions_from_c.o : submodule/myfort/src/functions_from_c.f90 
obj/general.o : submodule/myfort/src/general.f90 
obj/iteration_methods.o : submodule/myfort/src/iteration_methods.f90 obj/linear_algebra.o 
obj/linear_algebra.o : submodule/myfort/src/linear_algebra.f90 obj/matrix_definitions.o obj/vector_definitions.o 
obj/matrix_definitions.o : submodule/myfort/src/matrix_definitions.f90 obj/vector_definitions.o 
obj/myfort.o : submodule/myfort/src/myfort.f90 obj/renormalization.o obj/iteration_methods.o obj/linear_algebra.o obj/wave_functions.o obj/profiler.o obj/physics_constants.o obj/general.o obj/store_couplings.o obj/angular_momentum_couplings.o obj/functions_from_c.o 
obj/physics_constants.o : submodule/myfort/src/physics_constants.f90 
obj/profiler.o : submodule/myfort/src/profiler.f90 obj/general.o 
obj/renormalization.o : submodule/myfort/src/renormalization.f90 obj/profiler.o obj/linear_algebra.o 
obj/store_couplings.o : submodule/myfort/src/store_couplings.f90 obj/angular_momentum_couplings.o obj/functions_from_c.o obj/profiler.o 
obj/vector_definitions.o : submodule/myfort/src/vector_definitions.f90 obj/general.o 
obj/wave_functions.o : submodule/myfort/src/wave_functions.f90 obj/physics_constants.o obj/functions_from_c.o 
obj/Atomic.o : main/Atomic.F90 obj/HFMBPT.o obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o obj/myfort.o 
obj/BasisTransform.o : main/BasisTransform.F90 obj/WriteOperator.o obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o obj/myfort.o 
obj/HFInput.o : main/HFInput.F90 obj/myfort.o 
obj/HFMain.o : main/HFMain.F90 obj/Atomic.o obj/WriteOperator.o obj/HFMBPT.o obj/HartreeFock.o obj/ThreeBodyMonInteraction.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o obj/myfort.o 
obj/WriteOperator.o : main/WriteOperator.F90 obj/Operators.o obj/HFInput.o obj/myfort.o 
obj/dvode_f90_m.o : submodule/myfort/src/dvode/dvode_f90_m.f90
obj/renormalization.o : submodule/myfort/src/renormalization.f90 obj/dvode_f90_m.o obj/linear_algebra.o obj/profiler.o
