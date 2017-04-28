SET(IONMOVER_SRCS
  IonMoverFactory.cpp
  IonSystem.cpp
  BOSurfaceFactory.cpp
)

ADD_LIBRARY(ionmovers ${IONMOVER_SRCS})
