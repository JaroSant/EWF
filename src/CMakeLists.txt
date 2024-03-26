# Boost setup (needed for Windows to work)
set(Boost_USE_STATIC_LIBS OFF CACHE BOOL "use static libraries from Boost")
find_package(Boost REQUIRED COMPONENTS)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

if (WIN32)
  # disable autolinking in boost
  add_definitions( -DBOOST_ALL_NO_LIB )

  # force all boost libraries to dynamic link (we already disabled
  # autolinking, so I don't know why we need this, but we do!)
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif()

set(EWF_src myHelpers.cpp Polynomial.cpp PolynomialRootFinder.cpp WrightFisher.cpp)

set(EWF_hdr myHelpers.h Polynomial.h PolynomialRootFinder.h WrightFisher.h)

add_library(EWF STATIC ${EWF_src} ${EWF_hdr})

target_include_directories(
    EWF INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:>
)

set_target_properties(EWF PROPERTIES PUBLIC_HEADER "${EWF_hdr}")

target_link_libraries(EWF PRIVATE Boost::boost)

target_include_directories(EWF PUBLIC ${Boost_INCLUDE_DIRS})

if (PYTHON_BINDINGS)
    set_target_properties(EWF PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    pybind11_add_module(EWF_pybind EWF_pybind.cpp)
    target_link_libraries(EWF_pybind PRIVATE EWF)
endif()


