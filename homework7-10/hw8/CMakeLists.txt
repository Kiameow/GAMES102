add_executable(hw8 hw8.cpp)
target_link_libraries(hw8 PRIVATE hw_common)

add_custom_command(TARGET hw8 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:hw8>
            $<TARGET_RUNTIME_DLLS:hw8>
    COMMAND_EXPAND_LISTS
)
