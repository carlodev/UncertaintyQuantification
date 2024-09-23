// Simcenter STAR-CCM+ macro: TestMacro1.java
// Written by Simcenter STAR-CCM+ 18.06.007
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.turbulence.*;
import star.cadmodeler.*;
import star.kwturb.*;
import star.meshing.*;

public class TestMacro1 extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    CadModel cadModel_0 = 
      ((CadModel) simulation_0.get(SolidModelManager.class).getObject("3D-CAD Model 1"));

    ScalarQuantityDesignParameter scalarQuantityDesignParameter_0 = 
      ((ScalarQuantityDesignParameter) cadModel_0.getDesignParameterManager().getObject("AoA"));

    Units units_2 = 
      ((Units) simulation_0.getUnitsManager().getObject("deg"));

    scalarQuantityDesignParameter_0.getQuantity().setValueAndUnits(5.1, units_2);

    cadModel_0.update();

    SolidModelPart solidModelPart_0 = 
      ((SolidModelPart) simulation_0.get(SimulationPartManager.class).getPart("Body 1"));

    simulation_0.get(SimulationPartManager.class).updateParts(new NeoObjectVector(new Object[] {solidModelPart_0}));

    ScalarGlobalParameter scalarGlobalParameter_0 = 
      ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject("u_ts"));

    Units units_0 = 
      ((Units) simulation_0.getUnitsManager().getObject("m/s"));

    scalarGlobalParameter_0.getQuantity().setValueAndUnits(37.2325, units_0);

    PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum("Physics 1"));

    TurbulenceIntensityProfile turbulenceIntensityProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(TurbulenceIntensityProfile.class);

    Units units_1 = 
      ((Units) simulation_0.getUnitsManager().getObject(""));

    turbulenceIntensityProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(0.01, units_1);

    TurbulentViscosityRatioProfile turbulentViscosityRatioProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(TurbulentViscosityRatioProfile.class);

    turbulentViscosityRatioProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits(9.0, units_1);

    SstKwTurbModel sstKwTurbModel_0 = 
      physicsContinuum_0.getModelManager().getModel(SstKwTurbModel.class);

    sstKwTurbModel_0.setSigma_w1(0.7);

    sstKwTurbModel_0.setA1(0.33);

    sstKwTurbModel_0.setBetaStar(0.07);

    GammaReThetaTransitionModel gammaReThetaTransitionModel_0 = 
      physicsContinuum_0.getModelManager().getModel(GammaReThetaTransitionModel.class);

    gammaReThetaTransitionModel_0.setS1(1.8);

    gammaReThetaTransitionModel_0.getConset1().setValueAndUnits(2.134, units_1);

    MeshPipelineController meshPipelineController_0 = 
      simulation_0.get(MeshPipelineController.class);

    meshPipelineController_0.generateVolumeMesh();

    StepStoppingCriterion stepStoppingCriterion_0 = 
      ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion("Maximum Steps"));

    stepStoppingCriterion_0.setIsUsed(true);

    IntegerValue integerValue_0 = 
      stepStoppingCriterion_0.getMaximumNumberStepsObject();

    integerValue_0.getQuantity().setValue(5000);

    simulation_0.getSimulationIterator().run();


    XYPlot xYPlot_1 = 
      ((XYPlot) simulation_0.getPlotManager().getPlot("Cf_top"));

    xYPlot_1.export(resolvePath("/home/romain.ospital/STARCCMp/UQ/Results/CFSimAAA.csv"), ",");

    MonitorPlot monitorPlot_0 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("CL Monitor Plot"));

    monitorPlot_0.export(resolvePath("/home/romain.ospital/STARCCMp/UQ/Results/CLSimAAA.csv"), ",");

    MonitorPlot monitorPlot_1 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("CD Monitor Plot"));

    monitorPlot_1.export(resolvePath("/home/romain.ospital/STARCCMp/UQ/Results/CDSimAAA.csv"), ",");

  }
}
