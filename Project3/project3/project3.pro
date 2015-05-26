TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -larmadillo -lblas -llapack

SOURCES += main.cpp \
    CG.cpp \
    VMCSolver.cpp \
    WaveFunctions.cpp \
    HydrogenOrbitals.cpp \
    lib.cpp \
    GaussianOrbitals.cpp

HEADERS += \
    VMCSolver.h \
    WaveFunctions.h \
    HydrogenOrbitals.h \
    lib.h \
    GaussianOrbitals.h

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
