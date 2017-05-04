#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)

art_make(
          LIBRARY_NAME
                        FlashCharge
          LIB_LIBRARIES
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        art_Framework_Core
                        art_Framework_Services_Registry
                        art_Framework_Principal
                        art_Persistency_Common canvas_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance canvas_Persistency_Provenance       
                        art_Utilities canvas_Utilities canvas_Utilities
                        art_Utilities 
                        art_Persistency_Common 
                        art_Persistency_Provenance 
                        uboone_SpaceCharge
                        uboone_SpaceChargeServices_SpaceChargeServiceMicroBooNE_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CETLIB_EXCEPT}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
                        ${SIMULATIONBASE}

          MODULE_LIBRARIES
          			SpaceCharge
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        art_Framework_Core
                        art_Framework_Services_Registry
                        art_Framework_Principal
                        art_Utilities 
                        art_Persistency_Provenance canvas_Persistency_Provenance canvas_Persistency_Provenance       
                        art_Utilities canvas_Utilities canvas_Utilities
                        art_Utilities 
                        art_Persistency_Common 
                        art_Persistency_Provenance 
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CETLIB_EXCEPT}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
                        ${SIMULATIONBASE}
)

add_subdirectory(job)
add_subdirectory(Spacecharge)

install_headers()
install_source()
install_fhicl()
