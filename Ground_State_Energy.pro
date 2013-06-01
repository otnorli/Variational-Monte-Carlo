TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    wavefunction.cpp \
    ../../cppLibrary/lib.cpp \
    slaterimportantsampling.cpp \
    slaterenergy.cpp \
    minimazation.cpp

LIBS += -larmadillo -llapack -lblas

HEADERS += \
    wavefunction.h \
    ../../cppLibrary/lib.h \
    slaterimportantsampling.h \
    slaterenergy.h \
    minimazation.h
