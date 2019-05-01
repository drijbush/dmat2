find_program(ford NAMES ford)

if (EXISTS ${ford})
  add_custom_target(doc
    COMMAND gen-doc.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
else()
  message(STATUS "Cannot find ford -- disabling docs")
endif()
