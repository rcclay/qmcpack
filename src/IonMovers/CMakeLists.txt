SET(IONMOVER_SRCS
  IonMoverFactory.cpp
  IonSystem.cpp
  BOSurfaceFactory.cpp
  IonDriverBase.cpp
  CEIMC/CEIMC.cpp
  CEIMC/CEIMCUpdateAll.cpp
  BOSurfaces/LJPot.cpp
)

ADD_LIBRARY(ionmovers ${IONMOVER_SRCS})
