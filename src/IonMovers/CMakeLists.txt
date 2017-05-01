SET(IONMOVER_SRCS
  IonMoverFactory.cpp
  IonSystem.cpp
  BOSurfaceFactory.cpp
  CEIMC/CEIMCFactory.cpp
  CEIMC/CEIMC.cpp
)

ADD_LIBRARY(ionmovers ${IONMOVER_SRCS})
