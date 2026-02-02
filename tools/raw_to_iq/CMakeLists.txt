cmake_minimum_required(VERSION 3.20)
project(raw_to_iq C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(SOURCES
src/bypass.c
src/packet_decode.c
src/read_file.c
src/reconstruct.c
)

set(HEADERS
src/bypass.h
src/packet_decode.h
src/reconstruct.h
)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4 /EHsc)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -fexceptions)
endif()

# Группируем файлы в Visual Studio
source_group("Source Files" FILES ${SOURCES})
source_group("Header Files" FILES ${HEADERS})

target_include_directories(${PROJECT_NAME} PRIVATE src)

# Линковка
if(WIN32)
    # Аналог ws2_32
    target_link_libraries(${PROJECT_NAME} PRIVATE ws2_32)
endif()