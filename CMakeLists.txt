cmake_minimum_required(VERSION 4.1)
project(Tp3_Jean_Awe)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Tp3_Jean_Awe
        src/main.cpp
        src/Node.cpp
        src/LinkedList.cpp
        src/HashTable.cpp
        src/AlignmentResult.cpp
        src/BestResultsHeap.cpp
        src/BaseLUT.cpp
        src/GenomicPosition.cpp)

#add_executable(test_linkedlist
#        src/test_linkedlist.cpp
#        src/Node.cpp
#        src/LinkedList.cpp
#        src/GenomicPosition.cpp)

#add_executable(test_hashtable
#        src/test_hashtable.cpp
#        src/Node.cpp
#        src/LinkedList.cpp
#        src/HashTable.cpp
#        src/BaseLUT.cpp
#        src/GenomicPosition.cpp)

#add_executable(test_bestresultsheap
#        src/test_bestresultsheap.cpp
#        src/AlignmentResult.cpp
#        src/BestResultsHeap.cpp
#        src/GenomicPosition.cpp)

# === Debug mémoire (désactivé pour la remise) ===
# add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
# add_link_options(-fsanitize=address)