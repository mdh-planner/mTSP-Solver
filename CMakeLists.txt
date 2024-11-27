cmake_minimum_required(VERSION 3.10)

# Project name
project(MR_MT_Solver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")


# Include other source folders if necessary (e.g., 'include' folder in your project)
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}")

# Add the executable target (list all your source files)
add_executable(MR_MT_Solver
    "MR_MT_Solver.cpp"
    "src/chromosome.cpp"
    #"src/CPScheduler.cpp"
    "src/GA.cpp"
    #"src/importProblem.cpp"
    "src/interactionRewards.cpp"
    "src/logger.cpp"
    "src/model.cpp"
    #"src/problem_generator.cpp"
    "src/uva_problems.cpp"
    #"src/vtScheduler.cpp"
)

# Optional: set build type (Debug, Release)
set(CMAKE_BUILD_TYPE Release)  # or Debug

