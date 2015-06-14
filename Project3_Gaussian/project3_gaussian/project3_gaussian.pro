TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -larmadillo -lblas -llapack

SOURCES += main.cpp \
    GaussianOrbitals.cpp \
    VMCSolver.cpp \
    WaveFunctions.cpp \
    lib.cpp

HEADERS += \
    GaussianOrbitals.h \
    VMCSolver.h \
    WaveFunctions.h \
    lib.h

