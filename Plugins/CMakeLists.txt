add_paraview_plugin( VolumetricCylinderSource "1.0"
  SERVER_MANAGER_XML VolumetricCylinderSource.xml
  SERVER_MANAGER_SOURCES vtkVolumetricCylinderSource.h vtkVolumetricCylinderSource.cxx
)

add_paraview_plugin( VolumetricHollowCylinderSource "1.0"
  SERVER_MANAGER_XML VolumetricHollowCylinderSource.xml
  SERVER_MANAGER_SOURCES vtkVolumetricHollowCylinderSource.h vtkVolumetricHollowCylinderSource.cxx
)

add_paraview_plugin( VolumetricEllipsoidSource "1.0"
  SERVER_MANAGER_XML VolumetricEllipsoidSource.xml
  SERVER_MANAGER_SOURCES vtkVolumetricEllipsoidSource.h vtkVolumetricEllipsoidSource.cxx
)

add_paraview_plugin( PartialVolumeModeller "1.0"
  SERVER_MANAGER_XML PartialVolumeModeller.xml
  SERVER_MANAGER_SOURCES vtkPartialVolumeModeller.h vtkPartialVolumeModeller.cxx
)
