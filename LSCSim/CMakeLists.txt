set(MODULE_NAME LSCSim)
message(STATUS "Configuring ${MODULE_NAME} ...")

file(GLOB MODULE_HEADERS CONFIGURE_DEPENDS inc/*.hh)
file(GLOB MODULE_SOURCES CONFIGURE_DEPENDS src/*.cc)

add_library(${MODULE_NAME} SHARED
    ${MODULE_SOURCES}
)

target_include_directories(${MODULE_NAME}
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${Geant4_INCLUDE_DIRS}
)

target_link_libraries(${MODULE_NAME}
    PUBLIC
    MCObjs
    GLG4Sim
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ${Geant4_LIBRARIES}
)

install(TARGETS ${MODULE_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY inc/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hh"
)

# ---- Executables ----
file(GLOB MODULE_EXECUTABLES CONFIGURE_DEPENDS test/*.cc)

foreach(testsourcefile ${MODULE_EXECUTABLES})
    get_filename_component(testname ${testsourcefile} NAME_WE)
    add_executable(${testname} ${testsourcefile})
    target_link_libraries(${testname} PRIVATE ${MODULE_NAME})
    install(TARGETS ${testname} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()

install(DIRECTORY mac DESTINATION .)
install(DIRECTORY data DESTINATION .)
