# CMakeList.txt: проект CMake для polynomials; включите исходный код и определения,
# укажите здесь логику для конкретного проекта.
#

# Добавьте источник в исполняемый файл этого проекта.
add_executable (polynomials "src/polynomials.cpp"  "src/excerpt.cpp" "inc/excerpt.h" "inc/modules.h" "src/modules.cpp" "inc/polynomials.h" "src/main.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET polynomials PROPERTY CXX_STANDARD 20)
endif()

# TODO: Добавьте тесты и целевые объекты, если это необходимо.