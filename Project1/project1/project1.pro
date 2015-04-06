TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    VMCSolver.cpp \
    lib.cpp \
    Beryllium.cpp \
    Helium.cpp

HEADERS += \
    VMCSolver.h \
    lib.h \
    Helium.h \
    Beryllium.h \
    Atom.h

