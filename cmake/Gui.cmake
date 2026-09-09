if(EXISTS "${PROJECT_SOURCE_DIR}/.toolchains/gui-sysroot/usr")
  list(PREPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/.toolchains/gui-sysroot/usr")
endif()
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_WAYLAND OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_X11 ON CACHE BOOL "" FORCE)
set(OpenGL_GL_PREFERENCE LEGACY)
FetchContent_Declare(glfw URL https://github.com/glfw/glfw/archive/refs/tags/3.4.tar.gz DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_Declare(imgui URL https://github.com/ocornut/imgui/archive/refs/tags/v1.91.8.tar.gz DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_Declare(implot URL https://github.com/epezent/implot/archive/refs/tags/v0.16.tar.gz DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(glfw imgui implot)
find_package(OpenGL REQUIRED)
add_library(astraflow_ui STATIC
  ${imgui_SOURCE_DIR}/imgui.cpp ${imgui_SOURCE_DIR}/imgui_draw.cpp
  ${imgui_SOURCE_DIR}/imgui_tables.cpp ${imgui_SOURCE_DIR}/imgui_widgets.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
  ${implot_SOURCE_DIR}/implot.cpp ${implot_SOURCE_DIR}/implot_items.cpp)
target_include_directories(astraflow_ui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends ${implot_SOURCE_DIR})
target_link_libraries(astraflow_ui PUBLIC glfw OpenGL::GL)
add_executable(astraflow_gui apps/gui/main.cpp)
target_link_libraries(astraflow_gui PRIVATE astraflow astraflow_ui)
