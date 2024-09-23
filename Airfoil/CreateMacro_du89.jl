using FileIO
using XLSX, DataFrames


function write_macro(Uinf,TI,AoA,μr,σw1,α1,βstar, s1, C1, filename::String)

CLNAME = "CL"*filename
CDNAME =  "CD"*filename
CFNAME =  "CF"*filename

str0 = "// Simcenter STAR-CCM+ macro: macroG.java
// Written by Simcenter STAR-CCM+ 18.02.008
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.base.report.*;
import star.turbulence.*;
import star.cadmodeler.*;
import star.flow.*;
import star.vis.*;
import star.kwturb.*;
import star.meshing.*;
import star.vis.*;

public class $filename extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
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

    MeshActionManager meshActionManager_0 = 
      simulation_0.get(MeshActionManager.class);

    MeshPart meshPart_1 = 
      ((MeshPart) simulation_0.get(SimulationPartManager.class).getPart(\"Intersect\"));

    MeshPart meshPart_2 = 
      meshActionManager_0.subtractParts(new NeoObjectVector(new Object[] {solidModelPart_0, meshPart_1}), meshPart_1, \"Discrete\");

    Region region_0 = 
      simulation_0.getRegionManager().getRegion(\"Region\");

    simulation_0.getRegionManager().removeRegions(Arrays.<Region>asList(region_0));

    AutoMeshOperation autoMeshOperation_0 = 
      ((AutoMeshOperation) simulation_0.get(MeshOperationManager.class).getObject(\"Automated Mesh\"));

    autoMeshOperation_0.getInputGeometryObjects().setQuery(null);

    autoMeshOperation_0.getInputGeometryObjects().setObjects(meshPart_2);

    Region region_1 = 
      simulation_0.getRegionManager().createEmptyRegion();

    region_1.setPresentationName(\"Region\");

    Boundary boundary_0 = 
      region_1.getBoundaryManager().getBoundary(\"Default\");

    region_1.getBoundaryManager().removeBoundaries(new NeoObjectVector(new Object[] {boundary_0}));

    simulation_0.getRegionManager().newRegionsFromParts(new NeoObjectVector(new Object[] {meshPart_2}), \"OneRegion\", region_1, \"OneBoundaryPerPartSurface\", null, RegionManager.CreateInterfaceMode.BOUNDARY, \"OneEdgeBoundaryPerPart\", null);

    Boundary boundary_1 = 
      region_1.getBoundaryManager().getBoundary(\"Subtract.inlet\");

    InletBoundary inletBoundary_0 = 
      ((InletBoundary) simulation_0.get(ConditionTypeManager.class).get(InletBoundary.class));

    boundary_1.setBoundaryType(inletBoundary_0);

    boundary_1.getConditions().get(FlowDirectionOption.class).setSelected(FlowDirectionOption.Type.COMPONENTS);

    TurbulenceIntensityProfile turbulenceIntensityProfile_0 = 
      boundary_1.getValues().get(TurbulenceIntensityProfile.class);

    Units units_1 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"\"));

    turbulenceIntensityProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($TI, units_1);

    TurbulentViscosityRatioProfile turbulentViscosityRatioProfile_0 = 
      boundary_1.getValues().get(TurbulentViscosityRatioProfile.class);

    turbulentViscosityRatioProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($μr, units_1);

    Units units_2 = 
      simulation_0.getUnitsManager().getInternalUnits(new IntVector(new int[] {0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));

    VelocityMagnitudeProfile velocityMagnitudeProfile_0 = 
      boundary_1.getValues().get(VelocityMagnitudeProfile.class);

    velocityMagnitudeProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setDefinition(\"\$u0\");

    Boundary boundary_2 = 
      region_1.getBoundaryManager().getBoundary(\"Subtract.outlet\");

    PressureBoundary pressureBoundary_0 = 
      ((PressureBoundary) simulation_0.get(ConditionTypeManager.class).get(PressureBoundary.class));

    boundary_2.setBoundaryType(pressureBoundary_0);

    Boundary boundary_3 = 
      region_1.getBoundaryManager().getBoundary(\"Subtract.z+\");

    SymmetryBoundary symmetryBoundary_0 = 
      ((SymmetryBoundary) simulation_0.get(ConditionTypeManager.class).get(SymmetryBoundary.class));

    boundary_3.setBoundaryType(symmetryBoundary_0);

    Boundary boundary_4 = 
      region_1.getBoundaryManager().getBoundary(\"Subtract.z0\");

    boundary_4.setBoundaryType(symmetryBoundary_0);

    PlaneSection planeSection_0 = 
      (PlaneSection) simulation_0.getPartManager().createImplicitPart(new NeoObjectVector(new Object[] {}), new DoubleVector(new double[] {0.0, 0.0, 1.0}), new DoubleVector(new double[] {0.0, 0.0, 0.0}), 0, 1, new DoubleVector(new double[] {0.0}));

    Units units_3 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"m\"));

    planeSection_0.getOriginCoordinate().setCoordinate(units_3, units_3, units_3, new DoubleVector(new double[] {0.0, 0.0, 0.0}));

    planeSection_0.getInputParts().setQuery(null);

    Boundary boundary_5 = 
      region_1.getBoundaryManager().getBoundary(\"Subtract.Airfoil\");

    planeSection_0.getInputParts().setObjects(boundary_5);

    planeSection_0.setPresentationName(\"Cf_section\");

    XYPlot xYPlot_0 = 
      ((XYPlot) simulation_0.getPlotManager().getPlot(\"Cf_top\"));

    xYPlot_0.getParts().setQuery(null);

    xYPlot_0.getParts().setObjects(planeSection_0);

    ForceCoefficientReport forceCoefficientReport_0 = 
      ((ForceCoefficientReport) simulation_0.getReportManager().getReport(\"CD\"));

    forceCoefficientReport_0.getParts().setQuery(null);

    forceCoefficientReport_0.getParts().setObjects(boundary_5);

    ForceCoefficientReport forceCoefficientReport_1 = 
      ((ForceCoefficientReport) simulation_0.getReportManager().getReport(\"CL\"));

    forceCoefficientReport_1.getParts().setQuery(null);

    forceCoefficientReport_1.getParts().setObjects(boundary_5);

    ScalarGlobalParameter scalarGlobalParameter_0 = 
      ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject(\"u_ts\"));

    scalarGlobalParameter_0.getQuantity().setValueAndUnits($Uinf, units_2);

    PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum(\"Physics 1\"));

    SstKwTurbModel sstKwTurbModel_0 = 
      physicsContinuum_0.getModelManager().getModel(SstKwTurbModel.class);

    sstKwTurbModel_0.getUpwindOption().setSelected(UpwindOption.Type.SECOND_ORDER);

    sstKwTurbModel_0.setSigma_w1($σw1);

    sstKwTurbModel_0.setA1($α1);

    sstKwTurbModel_0.setBetaStar($βstar);

    GammaReThetaTransitionModel gammaReThetaTransitionModel_0 = 
      physicsContinuum_0.getModelManager().getModel(GammaReThetaTransitionModel.class);

    gammaReThetaTransitionModel_0.setS1($s1);

    gammaReThetaTransitionModel_0.getConset1().setValueAndUnits($C1, units_1);

    autoMeshOperation_0.execute();

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
}

"
macro_text = str0 

write(filename*".java", macro_text);

end

function convert_num2str(n::Int64)
  alphanumeric = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  t1 = Int(ceil(n/26^2))
  t2 = Int(ceil(n/26))
  t3 = mod(n-1,26)+1
  return alphanumeric[t1]*alphanumeric[t2]*alphanumeric[t3]  
end

test_run = DataFrame((XLSX.readtable("Sampling_DU89.xlsx","Sheet1")))
Nsamples = size(test_run)[1]



sbatch_string = open("string_sbatch.txt", "w")
alpha_nums = open("alpha_nums.txt", "w")

for i = 1:1:Nsamples
  Uinf=test_run[i,:].U_freestream
  TI=test_run[i,:].TI
  AoA=test_run[i,:].AoA
  μr=test_run[i,:].Viscosity_Ratio
  σw1=test_run[i,:].σw1
  α1=test_run[i,:].α1
  βstar= test_run[i,:].βstar
  s1= test_run[i,:].s1
  C1= test_run[i,:].C1

  alpha_string = convert_num2str(i)
  
  tstring = "starccm+ -batchsystem slurm -batch Macros/Sim$alpha_string.java DU89_A1_Re500000_3D_ReGammaTheta_Base.sim\n"
  write(sbatch_string, tstring);
  write(alpha_nums,"$(alpha_string)\n");
  filename= "Sim"*alpha_string
  write_macro(Uinf,TI,AoA,μr,σw1,α1,βstar, s1, C1, filename)
end

close(sbatch_string)
close(alpha_nums)

