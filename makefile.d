obj/wallclock.o : src/wallclock.c 
obj/RotationGroup.o : src/RotationGroup.f90 
obj/class_sys.o : src/class_sys.f90 
obj/MPIFunction.o : src/MPIFunction.F90 
obj/ModelSpace.o : src/ModelSpace.F90 obj/RotationGroup.o obj/LinAlgLib.o obj/VectorDouble.o obj/MatrixDouble.o obj/MPIFunction.o obj/HFInput.o 
obj/ScalarOperator.o : src/ScalarOperator.F90 obj/class_stopwatch.o obj/RotationGroup.o obj/class_sys.o obj/Read3BME.o obj/LinAlgLib.o obj/VectorDouble.o obj/MatrixDouble.o obj/ModelSpace.o obj/HFInput.o 
obj/HFInput.o : src/HFInput.F90 obj/class_sys.o obj/MPIFunction.o 
obj/class_stopwatch.o : src/class_stopwatch.F90 obj/MPIFunction.o 
obj/Read3BME.o : src/Read3BME.F90 obj/class_sys.o obj/RotationGroup.o obj/ModelSpace.o obj/MPIFunction.o obj/HFInput.o 
obj/HFMain.o : src/HFMain.F90 obj/class_stopwatch.o obj/MPIFunction.o obj/Read3BME.o obj/ModelSpace.o obj/RotationGroup.o obj/HFInput.o 
obj/LinAlgLib.o : LinAlgf90/src/LinAlgLib.f90 obj/Parameters.o obj/MatVecComplex.o obj/MatVecDouble.o obj/MatrixComplex.o obj/MatrixDouble.o obj/VectorComplex.o obj/VectorDouble.o 
obj/MatVecDouble.o : LinAlgf90/src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o 
obj/VectorComplex.o : LinAlgf90/src/VectorComplex.f90 obj/Parameters.o 
obj/MatrixDouble.o : LinAlgf90/src/MatrixDouble.f90 obj/VectorDouble.o 
obj/VectorDouble.o : LinAlgf90/src/VectorDouble.f90 obj/Parameters.o 
obj/MatrixComplex.o : LinAlgf90/src/MatrixComplex.f90 obj/VectorComplex.o 
obj/Parameters.o : LinAlgf90/src/Parameters.f90 
obj/MatVecComplex.o : LinAlgf90/src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o 
