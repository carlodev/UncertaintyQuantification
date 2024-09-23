str0 = "
// Simcenter STAR-CCM+ macro: macro_sd7003.java
// Written by Simcenter STAR-CCM+ 18.02.008
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.base.report.*;
import star.turbulence.*;
import star.turbulence.*;
import star.kwturb.*;
import star.cadmodeler.*;
import star.flow.*;
import star.meshing.*;

public class $filename extends StarMacro {

  public void execute() {
    execute0();
    execute1();

  }

  private void execute0() {"

  str1 =   "Simulation simulation_0 = 
      getActiveSimulation();

    CadModel cadModel_0 = 
      ((CadModel) simulation_0.get(SolidModelManager.class).getObject(\"3D-CAD Model 1\"));

    ScalarQuantityDesignParameter scalarQuantityDesignParameter_0 = 
      ((ScalarQuantityDesignParameter) cadModel_0.getDesignParameterManager().getObject(\"AoA\"));

    Units units_0 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"deg\"));

    scalarQuantityDesignParameter_0.getQuantity().setValueAndUnits($AoA, units_0);

    cadModel_0.update();

    SolidModelPart solidModelPart_0 = 
      ((SolidModelPart) simulation_0.get(SimulationPartManager.class).getPart(\"Body 1\"));

    simulation_0.get(SimulationPartManager.class).updateParts(new NeoObjectVector(new Object[] {solidModelPart_0}));

    MeshPart meshPart_0 = 
      ((MeshPart) simulation_0.get(SimulationPartManager.class).getPart(\"Subtract\"));

    simulation_0.get(SimulationPartManager.class).removeParts(new NeoObjectVector(new Object[] {meshPart_0}));

    Region region_0 = 
      simulation_0.getRegionManager().getRegion(\"Region\");

    simulation_0.getRegionManager().removeRegions(Arrays.<Region>asList(region_0));

    MeshActionManager meshActionManager_0 = 
      simulation_0.get(MeshActionManager.class);

    MeshPart meshPart_1 = 
      ((MeshPart) simulation_0.get(SimulationPartManager.class).getPart(\"Intersect\"));

    MeshPart meshPart_2 = 
      meshActionManager_0.subtractParts(new NeoObjectVector(new Object[] {solidModelPart_0, meshPart_1}), meshPart_1, \"Discrete\");

    simulation_0.getRegionManager().newRegionsFromParts(new NeoObjectVector(new Object[] {meshPart_2}), \"OneRegionPerPart\", null, \"OneBoundaryPerPartSurface\", null, RegionManager.CreateInterfaceMode.BOUNDARY, \"OneEdgeBoundaryPerPart\", null);

    Region region_1 = 
      simulation_0.getRegionManager().getRegion(\"Subtract\");

    Boundary boundary_0 = 
      region_1.getBoundaryManager().getBoundary(\"inlet\");

    InletBoundary inletBoundary_0 = 
      ((InletBoundary) simulation_0.get(ConditionTypeManager.class).get(InletBoundary.class));

    boundary_0.setBoundaryType(inletBoundary_0);

    Units units_1 = 
      simulation_0.getUnitsManager().getInternalUnits(new IntVector(new int[] {0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));

    VelocityMagnitudeProfile velocityMagnitudeProfile_0 = 
      boundary_0.getValues().get(VelocityMagnitudeProfile.class);

    velocityMagnitudeProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition(\"\$u0\");

    Boundary boundary_1 = 
      region_1.getBoundaryManager().getBoundary(\"outlet\");

    PressureBoundary pressureBoundary_0 = 
      ((PressureBoundary) simulation_0.get(ConditionTypeManager.class).get(PressureBoundary.class));

    boundary_1.setBoundaryType(pressureBoundary_0);

    Boundary boundary_2 = 
      region_1.getBoundaryManager().getBoundary(\"z+\");

    SymmetryBoundary symmetryBoundary_0 = 
      ((SymmetryBoundary) simulation_0.get(ConditionTypeManager.class).get(SymmetryBoundary.class));

    boundary_2.setBoundaryType(symmetryBoundary_0);

    Boundary boundary_3 = 
      region_1.getBoundaryManager().getBoundary(\"z0\");

    boundary_3.setBoundaryType(symmetryBoundary_0);

    boundary_0.getConditions().get(FlowDirectionOption.class).setSelected(FlowDirectionOption.Type.COMPONENTS);

    TurbulentViscosityRatioProfile turbulentViscosityRatioProfile_0 = 
      boundary_0.getValues().get(TurbulentViscosityRatioProfile.class);

    Units units_2 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"\"));

    turbulentViscosityRatioProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($μr, units_2);

    TurbulenceIntensityProfile turbulenceIntensityProfile_0 = 
      boundary_0.getValues().get(TurbulenceIntensityProfile.class);

    turbulenceIntensityProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($TI, units_2);

    ForceCoefficientReport forceCoefficientReport_0 = 
      ((ForceCoefficientReport) simulation_0.getReportManager().getReport(\"CD\"));

    forceCoefficientReport_0.getParts().setQuery(null);

    Boundary boundary_4 = 
      region_1.getBoundaryManager().getBoundary(\"airfoil\");

    forceCoefficientReport_0.getParts().setObjects(boundary_4);

    ForceCoefficientReport forceCoefficientReport_1 = 
      ((ForceCoefficientReport) simulation_0.getReportManager().getReport(\"CL\"));

    forceCoefficientReport_1.getParts().setQuery(null);

    forceCoefficientReport_1.getParts().setObjects(boundary_4);

    AutoMeshOperation autoMeshOperation_0 = 
      ((AutoMeshOperation) simulation_0.get(MeshOperationManager.class).getObject(\"Automated Mesh\"));

    autoMeshOperation_0.getInputGeometryObjects().setQuery(null);

    autoMeshOperation_0.getInputGeometryObjects().setObjects(meshPart_2);

    SurfaceCustomMeshControl surfaceCustomMeshControl_0 = 
      ((SurfaceCustomMeshControl) autoMeshOperation_0.getCustomMeshControls().getObject(\"NoBL\"));

    surfaceCustomMeshControl_0.getGeometryObjects().setQuery(null);

    PartSurface partSurface_0 = 
      ((PartSurface) meshPart_2.getPartSurfaceManager().getPartSurface(\"inlet\"));

    PartSurface partSurface_1 = 
      ((PartSurface) meshPart_2.getPartSurfaceManager().getPartSurface(\"outlet\"));

    PartSurface partSurface_2 = 
      ((PartSurface) meshPart_2.getPartSurfaceManager().getPartSurface(\"z+\"));

    PartSurface partSurface_3 = 
      ((PartSurface) meshPart_2.getPartSurfaceManager().getPartSurface(\"z0\"));

    surfaceCustomMeshControl_0.getGeometryObjects().setObjects(partSurface_0, partSurface_1, partSurface_2, partSurface_3);

    SurfaceCustomMeshControl surfaceCustomMeshControl_1 = 
      ((SurfaceCustomMeshControl) autoMeshOperation_0.getCustomMeshControls().getObject(\"WT_walls\"));

    surfaceCustomMeshControl_1.getGeometryObjects().setQuery(null);

    PartSurface partSurface_4 = 
      ((PartSurface) meshPart_2.getPartSurfaceManager().getPartSurface(\"walls\"));

    surfaceCustomMeshControl_1.getGeometryObjects().setObjects(partSurface_4);

    MeshPipelineController meshPipelineController_0 = 
      simulation_0.get(MeshPipelineController.class);

    meshPipelineController_0.generateVolumeMesh();
"

str2 = "  ScalarGlobalParameter scalarGlobalParameter_0 = 
((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject(\"u_ts\"));


Units units_3 = 
((Units) simulation_0.getUnitsManager().getObject(\"m/s\"));

scalarGlobalParameter_0.getQuantity().setValueAndUnits($Uinf, units_3);

PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum(\"Physics 1\"));

    GammaReThetaTransitionModel gammaReThetaTransitionModel_0 = 
      physicsContinuum_0.getModelManager().getModel(GammaReThetaTransitionModel.class);

    gammaReThetaTransitionModel_0.setS1($s1);

    gammaReThetaTransitionModel_0.getConset1().setValueAndUnits($C1, units_2);

    SstKwTurbModel sstKwTurbModel_0 = 
      physicsContinuum_0.getModelManager().getModel(SstKwTurbModel.class);

    sstKwTurbModel_0.setSigma_w1($σw1);

    sstKwTurbModel_0.setBetaStar($βstar);

    sstKwTurbModel_0.setA1($α1);


    PlaneSection planeSection_2 = 
    (PlaneSection) simulation_0.getPartManager().createImplicitPart(new NeoObjectVector(new Object[] {}), new DoubleVector(new double[] {0.0, 0.0, 1.0}), new DoubleVector(new double[] {0.0, 0.0, 0.0}), 0, 1, new DoubleVector(new double[] {0.0}));

  planeSection_2.setPresentationName(\"Cf_section\");

  planeSection_2.getInputParts().setQuery(null);

  Region region_2 = 
    simulation_0.getRegionManager().getRegion(\"Region\");

  Boundary boundary_5 = 
    region_2.getBoundaryManager().getBoundary(\"Body 1.airfoil\");

  planeSection_2.getInputParts().setObjects(boundary_5);

  Units units_4 = 
    ((Units) simulation_0.getUnitsManager().getObject(\"m\"));

  planeSection_2.getOriginCoordinate().setCoordinate(units_4, units_4, units_4, new DoubleVector(new double[] {0.0, 0.0, 0.1}));

  XYPlot xYPlot_0 = 
    ((XYPlot) simulation_0.getPlotManager().getPlot(\"Cf_top\"));

  xYPlot_0.getParts().setQuery(null);

  xYPlot_0.getParts().setObjects(planeSection_2);


    simulation_0.getSimulationIterator().run();


    XYPlot xYPlot_1 = 
      ((XYPlot) simulation_0.getPlotManager().getPlot(\"Cf_top\"));

    xYPlot_1.export(resolvePath(\"/home/cbrunelli/STARCCMp/UQ/Results/$CFNAME.csv\"), \",\");

    MonitorPlot monitorPlot_0 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot(\"CL Monitor Plot\"));

    monitorPlot_0.export(resolvePath(\"/home/cbrunelli/STARCCMp/UQ/Results/$CLNAME.csv\"), \",\");

    MonitorPlot monitorPlot_1 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot(\"CD Monitor Plot\"));

    monitorPlot_1.export(resolvePath(\"/home/cbrunelli/STARCCMp/UQ/Results/$CDNAME.csv\"), \",\");
  }

  private void execute1() {
  }
}
"