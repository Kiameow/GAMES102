add_executable(test main.cpp)
target_link_libraries(test PRIVATE hw_common)

add_custom_command(TARGET test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:test>
            $<TARGET_RUNTIME_DLLS:test>
    COMMAND_EXPAND_LISTS
)
