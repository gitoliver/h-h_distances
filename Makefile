#############################################################################
# Makefile for building: h-h.exe
# Generated by qmake (2.01a) (Qt 4.8.7) on: Mo. Jul 10 11:52:13 2017
# Project:  h-h_distances.pro
# Template: app
# Command: /usr/lib/x86_64-linux-gnu/qt4/bin/qmake -o Makefile h-h_distances.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = 
CFLAGS        = -m64 -pipe -O2 -Wall -W $(DEFINES)
CXXFLAGS      = -m64 -pipe -O2 -Wall -W $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++-64 -I. -I../../gems/gmml/bin
LINK          = g++
LFLAGS        = -m64 -Wl,-O1
LIBS          = $(SUBLIBS)   -L/home/oliver/Programs/Cplusplus/h-h_distances/../../gems/gmml/bin/ -lgmml 
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/lib/x86_64-linux-gnu/qt4/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cpp \
		io.cpp 
OBJECTS       = main.o \
		io.o
DIST          = /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		h-h_distances.pro
QMAKE_TARGET  = h-h.exe
DESTDIR       = 
TARGET        = h-h.exe

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: h-h_distances.pro  /usr/share/qt4/mkspecs/linux-g++-64/qmake.conf /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf
	$(QMAKE) -o Makefile h-h_distances.pro
/usr/share/qt4/mkspecs/common/unix.conf:
/usr/share/qt4/mkspecs/common/linux.conf:
/usr/share/qt4/mkspecs/common/gcc-base.conf:
/usr/share/qt4/mkspecs/common/gcc-base-unix.conf:
/usr/share/qt4/mkspecs/common/g++-base.conf:
/usr/share/qt4/mkspecs/common/g++-unix.conf:
/usr/share/qt4/mkspecs/qconfig.pri:
/usr/share/qt4/mkspecs/features/qt_functions.prf:
/usr/share/qt4/mkspecs/features/qt_config.prf:
/usr/share/qt4/mkspecs/features/exclusive_builds.prf:
/usr/share/qt4/mkspecs/features/default_pre.prf:
/usr/share/qt4/mkspecs/features/release.prf:
/usr/share/qt4/mkspecs/features/default_post.prf:
/usr/share/qt4/mkspecs/features/shared.prf:
/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf:
/usr/share/qt4/mkspecs/features/warn_on.prf:
/usr/share/qt4/mkspecs/features/resources.prf:
/usr/share/qt4/mkspecs/features/uic.prf:
/usr/share/qt4/mkspecs/features/yacc.prf:
/usr/share/qt4/mkspecs/features/lex.prf:
/usr/share/qt4/mkspecs/features/include_source_dir.prf:
qmake:  FORCE
	@$(QMAKE) -o Makefile h-h_distances.pro

dist: 
	@$(CHK_DIR_EXISTS) .tmp/h-h.exe1.0.0 || $(MKDIR) .tmp/h-h.exe1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/h-h.exe1.0.0/ && (cd `dirname .tmp/h-h.exe1.0.0` && $(TAR) h-h.exe1.0.0.tar h-h.exe1.0.0 && $(COMPRESS) h-h.exe1.0.0.tar) && $(MOVE) `dirname .tmp/h-h.exe1.0.0`/h-h.exe1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/h-h.exe1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


check: first

compiler_rcc_make_all:
compiler_rcc_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

main.o: main.cpp io.h \
		/home/oliver/Programs/gems/gmml/includes/gmml.hpp \
		../../gems/gmml/includes/common.hpp \
		../../gems/gmml/includes/GeometryTopology/coordinate.hpp \
		../../gems/gmml/includes/Glycan/sugarname.hpp \
		../../gems/gmml/includes/utils.hpp \
		../../gems/gmml/includes/MolecularModeling/atom.hpp \
		../../gems/gmml/includes/MolecularModeling/moleculardynamicatom.hpp \
		../../gems/gmml/includes/MolecularModeling/quantommechanicatom.hpp \
		../../gems/gmml/includes/MolecularModeling/dockingatom.hpp \
		../../gems/gmml/includes/superimposition.hpp \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/Geometry \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/Core \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/DisableStupidWarnings.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/Macros.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/MKL_support.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/misc/blas.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/Constants.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/Meta.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/ForwardDeclarations.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/StaticAssert.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/XprHelper.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/Memory.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/NumTraits.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/GenericPacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/MathFunctionsImpl.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/SSE/PacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AVX/PacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AVX512/PacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AVX512/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/SSE/Complex.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/SSE/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AVX/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AVX/Complex.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AVX/TypeCasting.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/SSE/TypeCasting.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/PacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/Complex.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/NEON/PacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/NEON/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/NEON/Complex.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/PacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/Complex.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/Half.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/PacketMathHalf.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/TypeCasting.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/PacketMath.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/MathFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/Default/Settings.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/functors/TernaryFunctors.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/functors/BinaryFunctors.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/functors/UnaryFunctors.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/functors/NullaryFunctors.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/functors/StlFunctors.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/functors/AssignmentFunctors.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/Complex.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/IO.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/DenseCoeffsBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/DenseBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/plugins/BlockMethods.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/MatrixBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/plugins/CommonCwiseUnaryOps.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/plugins/CommonCwiseBinaryOps.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/plugins/MatrixCwiseUnaryOps.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/plugins/MatrixCwiseBinaryOps.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/EigenBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Product.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CoreEvaluators.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/AssignEvaluator.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Assign.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/ArrayBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/plugins/ArrayCwiseUnaryOps.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/plugins/ArrayCwiseBinaryOps.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/BlasUtil.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/DenseStorage.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/NestByValue.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/ReturnByValue.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/NoAlias.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/PlainObjectBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Matrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Array.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CwiseTernaryOp.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CwiseBinaryOp.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CwiseUnaryOp.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CwiseNullaryOp.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CwiseUnaryView.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/SelfCwiseBinaryOp.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Dot.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/StableNorm.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Stride.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/MapBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Map.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Ref.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Block.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/VectorBlock.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Transpose.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/DiagonalMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Diagonal.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/DiagonalProduct.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Redux.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Visitor.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Fuzzy.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Swap.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CommaInitializer.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/GeneralProduct.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Solve.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Inverse.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/SolverBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/PermutationMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Transpositions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/TriangularMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/SelfAdjointView.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/GeneralBlockPanelKernel.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/Parallelizer.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/ProductEvaluators.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixVector.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/SolveTriangular.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrixTriangular.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixVector.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointProduct.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointRank2Update.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixVector.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverVector.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/BandMatrix.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/CoreIterators.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/ConditionEstimator.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/BooleanRedux.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Select.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/VectorwiseOp.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Random.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Replicate.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Reverse.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/ArrayWrapper.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrix_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixVector_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixVector_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixMatrix_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixVector_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverMatrix_BLAS.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/Assign_MKL.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/GlobalFunctions.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Core/util/ReenableStupidWarnings.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/SVD \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/QR \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/Cholesky \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Cholesky/LLT.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Cholesky/LDLT.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/misc/lapacke.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/misc/lapacke_mangling.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Cholesky/LLT_LAPACKE.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/Jacobi \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Jacobi/Jacobi.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/Householder \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Householder/Householder.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Householder/HouseholderSequence.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Householder/BlockHouseholder.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/QR/HouseholderQR.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/QR/FullPivHouseholderQR.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/QR/ColPivHouseholderQR.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/QR/CompleteOrthogonalDecomposition.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/QR/HouseholderQR_LAPACKE.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/QR/ColPivHouseholderQR_LAPACKE.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/misc/RealSvd2x2.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/SVD/UpperBidiagonalization.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/SVD/SVDBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/SVD/JacobiSVD.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/SVD/BDCSVD.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/SVD/JacobiSVD_LAPACKE.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/LU \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/misc/Kernel.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/misc/Image.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/LU/FullPivLU.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/LU/PartialPivLU.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/LU/PartialPivLU_LAPACKE.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/LU/Determinant.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/LU/InverseImpl.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/LU/arch/Inverse_SSE.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/OrthoMethods.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/EulerAngles.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Homogeneous.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/RotationBase.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Rotation2D.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Quaternion.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/AngleAxis.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Transform.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Translation.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Scaling.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Hyperplane.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/ParametrizedLine.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/AlignedBox.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/Umeyama.h \
		../../gems/gmml/includes/Eigen_Algebra_Template_Library/src/Geometry/arch/Geometry_SSE.h \
		../../gems/gmml/includes/MolecularModeling/assembly.hpp \
		../../gems/gmml/includes/GeometryTopology/plane.hpp \
		../../gems/gmml/includes/Glycan/chemicalcode.hpp \
		../../gems/gmml/includes/Glycan/monosaccharide.hpp \
		../../gems/gmml/includes/Glycan/ontologyvocabulary.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbfile.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		../../gems/gmml/includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		../../gems/gmml/includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		../../gems/gmml/includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		../../gems/gmml/includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		../../gems/gmml/includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		../../gems/gmml/includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbmodel.hpp \
		../../gems/gmml/includes/Glycan/oligosaccharide.hpp \
		../../gems/gmml/includes/Glycan/note.hpp \
		../../gems/gmml/includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		../../gems/gmml/includes/InputSet/CifFileSpace/ciffileatom.hpp \
		../../gems/gmml/includes/InputSet/CifFileSpace/ciffileprocessingexception.hpp \
		../../gems/gmml/includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbatom.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbcompoundcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbconnectcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbdisulfidebondcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbformula.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbformulacard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheadercard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbhelix.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbhelixcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheterogen.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheterogencard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheterogenname.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheterogensynonym.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdblink.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdblinkresidue.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbmatrixn.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbmodeltypecard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdboriginxn.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdboriginxncard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbresidue.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbresiduemodification.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbresiduesequence.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbscalen.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbscalencard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbsheet.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbsheetcard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbsite.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbsitecard.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp \
		../../gems/gmml/includes/InputSet/PdbFileSpace/pdbtitlecard.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp \
		../../gems/gmml/includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologybond.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp \
		../../gems/gmml/includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		../../gems/gmml/includes/GeometryTopology/cell.hpp \
		../../gems/gmml/includes/GeometryTopology/grid.hpp \
		../../gems/gmml/includes/GeometryTopology/InternalCoordinate/angle.hpp \
		../../gems/gmml/includes/GeometryTopology/InternalCoordinate/dihedral.hpp \
		../../gems/gmml/includes/GeometryTopology/InternalCoordinate/distance.hpp \
		../../gems/gmml/includes/MolecularModeling/atomnode.hpp \
		../../gems/gmml/includes/MolecularModeling/element.hpp \
		../../gems/gmml/includes/MolecularModeling/residue.hpp \
		../../gems/gmml/includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		../../gems/gmml/includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp \
		../../gems/gmml/includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		../../gems/gmml/includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		../../gems/gmml/includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		../../gems/gmml/includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		../../gems/gmml/includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		../../gems/gmml/includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		../../gems/gmml/includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp \
		../../gems/gmml/includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessor.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp \
		../../gems/gmml/includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

io.o: io.cpp io.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o io.o io.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

