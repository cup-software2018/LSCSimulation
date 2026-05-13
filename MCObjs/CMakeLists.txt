set(MODULE_NAME MCObjs)
message(STATUS "Configuring ${MODULE_NAME} ...")

file(GLOB MODULE_HEADERS CONFIGURE_DEPENDS inc/*.hh)
file(GLOB MODULE_SOURCES CONFIGURE_DEPENDS src/*.cc)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)

set(DICT_HEADERS "")
foreach(_h ${MODULE_HEADERS})
    get_filename_component(_n ${_h} NAME)
    list(APPEND DICT_HEADERS ${_n})
endforeach()

ROOT_GENERATE_DICTIONARY(G__${MODULE_NAME}
    ${DICT_HEADERS}
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h
    OPTIONS
    -I${CMAKE_CURRENT_SOURCE_DIR}/inc
    -Wno-inconsistent-missing-override
    NOINSTALL
)

add_library(${MODULE_NAME} SHARED
    ${MODULE_SOURCES}
    G__${MODULE_NAME}.cxx
)

target_include_directories(${MODULE_NAME}
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(${MODULE_NAME}
    PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
)

install(TARGETS ${MODULE_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY inc/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hh"
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/lib${MODULE_NAME}_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/lib${MODULE_NAME}.rootmap
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
