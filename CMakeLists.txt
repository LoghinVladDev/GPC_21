cmake_minimum_required(VERSION 3.17)
project(GPC)

set(CMAKE_CXX_STANDARD 20)

include_directories()

set(GLUT_LIBS GL GLU glut)
set(GLFW_LIBS GL GLU glfw3 X11 Xxf86vm Xrandr pthread Xi dl Xinerama Xcursor assimp)

include_directories(dependencies/eng/src/GLFW/src/engine/obj/util/shader)

macro(glut_exe)
    add_executable(GLut_${ARGV0} freeglut/${ARGV0}/main.cpp)
    target_link_libraries(GLut_${ARGV0} ${GLUT_LIBS})
endmacro()

macro(glfw_exe)
    add_executable(GLFW_${ARGV0}
            glfw/${ARGV0}/main.cpp
            dependencies/glad/src/glad.c
            dependencies/eng/src/GLFW/src/engine/obj/util/shader/Shader.cpp
    )
    target_link_libraries(GLFW_${ARGV0} ${GLFW_LIBS})
endmacro()

glut_exe(t1)
glut_exe(t2)
glut_exe(t3)
target_link_libraries(GLut_t3 pthread)

glfw_exe(t1)