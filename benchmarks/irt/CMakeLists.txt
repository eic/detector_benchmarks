cmake_minimum_required(VERSION 3.16)

project(irt_analysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the packages
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(EDM4EIC REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist Physics)

# Collect all source files
file(GLOB SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(irt_analysis
    ${SOURCES}
)

target_include_directories(irt_analysis PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Do not let linker drop podioRootIO
target_link_options(irt_analysis PRIVATE
    -Wl,--no-as-needed
)


# Link against the imported targets
target_link_libraries(irt_analysis
    PRIVATE
	podio::podio
	podio::podioIO
        podio::podioRootIO
	EDM4HEP::edm4hep
        EDM4EIC::edm4eic
        ROOT::Core
        ROOT::RIO
        ROOT::Tree
        ROOT::Hist
	ROOT::Physics
)

install(TARGETS irt_analysis DESTINATION bin)
