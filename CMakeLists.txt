project( VTKVolumetricSources )

cmake_minimum_required( VERSION 2.8 )

find_package( ParaView REQUIRED )
include( ${PARAVIEW_USE_FILE} )

add_subdirectory( Plugins )