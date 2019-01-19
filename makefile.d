obj/CommonLibrary.o : src/CommonLibrary.f90 obj/ClassSys.o 
obj/ClassSys.o : src/ClassSys.f90 
obj/Operators.o : src/Operators.F90 obj/Profiler.o obj/DefineOperators.o obj/NOperators.o obj/ModelSpace.o 
obj/NOperators.o : src/NOperators.F90 obj/Profiler.o obj/ClassSys.o obj/DefineOperators.o obj/CommonLibrary.o obj/ModelSpace.o obj/LinAlgLib.o 
obj/MPIFunction.o : src/MPIFunction.F90 
obj/MBPT.o : src/MBPT.F90 obj/CommonLibrary.o obj/Profiler.o obj/Operators.o obj/ModelSpace.o 
obj/Profiler.o : src/Profiler.F90 obj/MPIFunction.o obj/ClassSys.o 
obj/ModelSpace.o : src/ModelSpace.F90 obj/LinAlgLib.o obj/CommonLibrary.o obj/ClassSys.o obj/Profiler.o obj/SingleParticleState.o 
obj/SingleParticleState.o : src/SingleParticleState.F90 
obj/DefineOperators.o : src/DefineOperators.F90 obj/CommonLibrary.o 
obj/HartreeFock.o : src/HartreeFock.F90 obj/MBPT.o obj/ModelSpace.o obj/ClassSys.o obj/CommonLibrary.o obj/Profiler.o obj/Operators.o obj/LinAlgLib.o 
obj/SingleDouble.o : LinAlgf90/src/SingleDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/MatrixSingle.o obj/VectorSingle.o 
obj/LinAlgLib.o : LinAlgf90/src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatVecSingle.o obj/MatrixComplex.o obj/MatrixDouble.o obj/MatrixSingle.o obj/VectorComplex.o obj/VectorDouble.o obj/VectorSingle.o obj/SingleDouble.o obj/LinAlgParameters.o 
obj/MatVecDouble.o : LinAlgf90/src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/LinAlgParameters.o 
obj/VectorSingle.o : LinAlgf90/src/VectorSingle.f90 obj/LinAlgParameters.o 
obj/LinAlgParameters.o : LinAlgf90/src/LinAlgParameters.f90 
obj/VectorComplex.o : LinAlgf90/src/VectorComplex.f90 obj/LinAlgParameters.o 
obj/MatrixSingle.o : LinAlgf90/src/MatrixSingle.f90 obj/VectorSingle.o obj/LinAlgParameters.o 
obj/MatrixDouble.o : LinAlgf90/src/MatrixDouble.f90 obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatVecSingle.o : LinAlgf90/src/MatVecSingle.f90 obj/MatrixSingle.o obj/VectorSingle.o obj/LinAlgParameters.o 
obj/VectorDouble.o : LinAlgf90/src/VectorDouble.f90 obj/LinAlgParameters.o 
obj/MatrixComplex.o : LinAlgf90/src/MatrixComplex.f90 obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MatVecComplex.o : LinAlgf90/src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/LinAlgParameters.o 
