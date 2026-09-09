# if(TARGET eigen)
#     return()
# endif()

# message(STATUS "Fetch third-party: creating target 'eigen'")

# include(FetchContent) # Load CMake module FetchContent designed for populating content at configure time.
# FetchContent_Declare(
#     eigen
#     GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
#     GIT_TAG tags/3.4.0
#     GIT_SHALLOW TRUE
# )
# FetchContent_MakeAvailable(eigen)

# # Add your target in Project B (Eigen)
# add_library(eigen INTERFACE)
# add_library(Eigen ALIAS eigen)
# target_include_directories(eigen INTERFACE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>)

# message(STATUS "eigen_SOURCE_DIR in CMake-learn now is : ${eigen_SOURCE_DIR}")
# message(STATUS "eigen_BINARY_DIR in CMake-learn now is : ${eigen_BINARY_DIR}")
# # Export the target
# export(TARGETS eigen NAMESPACE geomp:: FILE geompConfig.cmake)