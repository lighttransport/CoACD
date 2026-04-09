if(NOT DEFINED DEFAULT_EXEC OR NOT DEFINED CDT_EXEC)
  message(FATAL_ERROR "DEFAULT_EXEC and CDT_EXEC must be set")
endif()

execute_process(
  COMMAND "${DEFAULT_EXEC}"
  RESULT_VARIABLE default_result
  OUTPUT_VARIABLE default_output
  ERROR_VARIABLE default_error
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT default_result EQUAL 0)
  message(FATAL_ERROR
          "Built-in clip runner failed with code ${default_result}\n${default_error}")
endif()

execute_process(
  COMMAND "${CDT_EXEC}"
  RESULT_VARIABLE cdt_result
  OUTPUT_VARIABLE cdt_output
  ERROR_VARIABLE cdt_error
  OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT cdt_result EQUAL 0)
  message(FATAL_ERROR
          "CDT clip runner failed with code ${cdt_result}\n${cdt_error}")
endif()

string(REPLACE "\n" ";" default_lines "${default_output}")
string(REPLACE "\n" ";" cdt_lines "${cdt_output}")

list(LENGTH default_lines default_len)
list(LENGTH cdt_lines cdt_len)
if(NOT default_len EQUAL cdt_len)
  message(FATAL_ERROR
          "Runner outputs have different line counts.\nBuilt-in:\n${default_output}\n\nCDT:\n${cdt_output}")
endif()

math(EXPR last_index "${default_len} - 1")
foreach(index RANGE ${last_index})
  list(GET default_lines ${index} default_line)
  list(GET cdt_lines ${index} cdt_line)

  if(default_line STREQUAL "" AND cdt_line STREQUAL "")
    continue()
  endif()

  if(NOT "${default_line}" MATCHES
         "^([^|]+)\\|success=([0-9]+)\\|plane=(-?[0-9]+)\\|cut=(-?[0-9]+)\\|pos_volume=(-?[0-9]+)\\|neg_volume=(-?[0-9]+)\\|pos_points=([0-9]+)\\|pos_tris=([0-9]+)\\|neg_points=([0-9]+)\\|neg_tris=([0-9]+)$")
    message(FATAL_ERROR "Failed to parse built-in line: ${default_line}")
  endif()
  set(default_name "${CMAKE_MATCH_1}")
  set(default_success "${CMAKE_MATCH_2}")
  set(default_plane "${CMAKE_MATCH_3}")
  set(default_cut "${CMAKE_MATCH_4}")
  set(default_pos_volume "${CMAKE_MATCH_5}")
  set(default_neg_volume "${CMAKE_MATCH_6}")
  set(default_pos_points "${CMAKE_MATCH_7}")
  set(default_neg_points "${CMAKE_MATCH_9}")

  if(NOT "${cdt_line}" MATCHES
         "^([^|]+)\\|success=([0-9]+)\\|plane=(-?[0-9]+)\\|cut=(-?[0-9]+)\\|pos_volume=(-?[0-9]+)\\|neg_volume=(-?[0-9]+)\\|pos_points=([0-9]+)\\|pos_tris=([0-9]+)\\|neg_points=([0-9]+)\\|neg_tris=([0-9]+)$")
    message(FATAL_ERROR "Failed to parse CDT line: ${cdt_line}")
  endif()
  set(cdt_name "${CMAKE_MATCH_1}")
  set(cdt_success "${CMAKE_MATCH_2}")
  set(cdt_plane "${CMAKE_MATCH_3}")
  set(cdt_cut "${CMAKE_MATCH_4}")
  set(cdt_pos_volume "${CMAKE_MATCH_5}")
  set(cdt_neg_volume "${CMAKE_MATCH_6}")
  set(cdt_pos_points "${CMAKE_MATCH_7}")
  set(cdt_neg_points "${CMAKE_MATCH_9}")

  if(NOT default_name STREQUAL cdt_name)
    message(FATAL_ERROR
            "Runner outputs are misaligned.\nBuilt-in line: ${default_line}\nCDT line: ${cdt_line}")
  endif()

  if(cdt_success STREQUAL "1")
    if(NOT default_success STREQUAL "1")
      message(FATAL_ERROR
              "Built-in clip failed where CDT succeeded for ${cdt_name}.\nBuilt-in: ${default_line}\nCDT: ${cdt_line}")
    endif()

    foreach(field IN ITEMS plane cut pos_volume neg_volume pos_points neg_points)
      if(NOT default_${field} STREQUAL cdt_${field})
        message(FATAL_ERROR
                "Built-in clip disagrees with CDT for ${cdt_name} on ${field}.\nBuilt-in: ${default_line}\nCDT: ${cdt_line}")
      endif()
    endforeach()
  endif()
endforeach()
