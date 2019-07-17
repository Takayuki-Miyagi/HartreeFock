obj/ClassSys.o : src/ClassSys.f90 
obj/MyLibrary.o : src/MyLibrary.f90 obj/ClassSys.o 
obj/DefineOperators.o : src/DefineOperators.F90 obj/MyLibrary.o 
obj/HFMBPT.o : src/HFMBPT.F90 obj/MyLibrary.o obj/Profiler.o obj/HartreeFock.o obj/StoreCouplings.o obj/Operators.o obj/ModelSpace.o 
obj/HartreeFock.o : src/HartreeFock.F90 obj/MyLibrary.o obj/Profiler.o obj/ThreeBodyNO2BInteraction.o obj/ThreeBodyMonInteraction.o obj/Operators.o obj/LinAlgLib.o 
obj/MPIFunction.o : src/MPIFunction.F90 
obj/ModelSpace.o : src/ModelSpace.F90 obj/MyLibrary.o obj/Profiler.o obj/ClassSys.o obj/ThreeBodyModelSpace.o obj/TwoBodyModelSpace.o obj/OneBodyModelSpace.o obj/SingleParticleState.o 
obj/OneBodyModelSpace.o : src/OneBodyModelSpace.F90 obj/SingleParticleState.o 
obj/OneBodyOperator.o : src/OneBodyOperator.F90 obj/ClassSys.o obj/DefineOperators.o obj/MyLibrary.o obj/OneBodyModelSpace.o obj/LinAlgLib.o 
obj/Operators.o : src/Operators.F90 obj/Profiler.o obj/DefineOperators.o obj/ThreeBodyNO2BInteraction.o obj/ThreeBodyMonInteraction.o obj/ThreeBodyInteraction.o obj/ThreeBodyOperator.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ModelSpace.o 
obj/Profiler.o : src/Profiler.F90 obj/MPIFunction.o obj/ClassSys.o 
obj/SingleParticleState.o : src/SingleParticleState.F90 obj/ClassSys.o 
obj/StoreCouplings.o : src/StoreCouplings.F90 obj/MyLibrary.o obj/Profiler.o 
obj/ThreeBodyInteraction.o : src/ThreeBodyInteraction.F90 obj/TwoBodyOperator.o obj/LinAlgLib.o obj/ClassSys.o obj/MyLibrary.o obj/ThreeBodyModelSpace.o obj/SingleParticleState.o obj/Profiler.o 
obj/ThreeBodyModelSpace.o : src/ThreeBodyModelSpace.F90 obj/MyLibrary.o obj/LinAlgLib.o obj/SingleParticleState.o 
obj/ThreeBodyMonInteraction.o : src/ThreeBodyMonInteraction.F90 obj/ClassSys.o obj/MyLibrary.o obj/SingleParticleState.o obj/Profiler.o 
obj/ThreeBodyNO2BInteraction.o : src/ThreeBodyNO2BInteraction.F90 obj/ClassSys.o obj/MyLibrary.o obj/SingleParticleState.o obj/Profiler.o 
obj/ThreeBodyOperator.o : src/ThreeBodyOperator.F90 obj/ClassSys.o obj/ThreeBodyInteraction.o obj/MyLibrary.o obj/TwoBodyOperator.o obj/OneBodyOperator.o obj/ThreeBodyModelSpace.o obj/LinAlgLib.o 
obj/TwoBodyModelSpace.o : src/TwoBodyModelSpace.F90 obj/MyLibrary.o obj/SingleParticleState.o 
obj/TwoBodyOperator.o : src/TwoBodyOperator.F90 obj/Profiler.o obj/ClassSys.o obj/DefineOperators.o obj/MyLibrary.o obj/OneBodyOperator.o obj/TwoBodyModelSpace.o obj/LinAlgLib.o 
obj/LinAlgLib.o : submodule/LinAlgf90/src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatVecSingle.o obj/MatrixComplex.o obj/MatrixDouble.o obj/MatrixSingle.o obj/VectorComplex.o obj/VectorDouble.o obj/VectorSingle.o obj/SingleDoubleComplex.o obj/LinAlgParameters.o 
obj/LinAlgParameters.o : submodule/LinAlgf90/src/LinAlgParameters.f90 
obj/MatVecComplex.o : submodule/LinAlgf90/src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MatVecDouble.o : submodule/LinAlgf90/src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatVecSingle.o : submodule/LinAlgf90/src/MatVecSingle.f90 obj/MatrixSingle.o obj/VectorSingle.o obj/LinAlgParameters.o 
obj/MatrixComplex.o : submodule/LinAlgf90/src/MatrixComplex.f90 obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MatrixDouble.o : submodule/LinAlgf90/src/MatrixDouble.f90 obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatrixSingle.o : submodule/LinAlgf90/src/MatrixSingle.f90 obj/VectorSingle.o obj/LinAlgParameters.o 
obj/SingleDoubleComplex.o : submodule/LinAlgf90/src/SingleDoubleComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/MatrixDouble.o obj/VectorDouble.o obj/MatrixSingle.o obj/VectorSingle.o 
obj/VectorComplex.o : submodule/LinAlgf90/src/VectorComplex.f90 obj/LinAlgParameters.o 
obj/VectorDouble.o : submodule/LinAlgf90/src/VectorDouble.f90 obj/LinAlgParameters.o 
obj/VectorSingle.o : submodule/LinAlgf90/src/VectorSingle.f90 obj/LinAlgParameters.o 
obj/Atomic.o : main/Atomic.F90 obj/MyLibrary.o obj/HFMBPT.o obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o 
obj/HFInput.o : main/HFInput.F90 obj/ClassSys.o 
obj/HFMain.o : main/HFMain.F90 obj/Atomic.o obj/WriteOperator.o obj/HFMBPT.o obj/HartreeFock.o obj/ThreeBodyMonInteraction.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o obj/Profiler.o 
obj/WriteOperator.o : main/WriteOperator.F90 obj/Profiler.o obj/ClassSys.o obj/Operators.o obj/HFInput.o 
