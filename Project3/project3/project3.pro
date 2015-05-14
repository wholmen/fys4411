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

