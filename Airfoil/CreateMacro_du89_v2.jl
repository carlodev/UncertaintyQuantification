using FileIO
using XLSX, DataFrames


function write_macro(Uinf,TI,AoA,μr,σw1,α1,βstar, s1, C1, filename::String)

CLNAME = "CL"*filename
CDNAME =  "CD"*filename
CFNAME =  "CF"*filename

str0 = "// Simcenter STAR-CCM+ macro: Macro_good.java
// Written by Simcenter STAR-CCM+ 18.06.007
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.turbulence.*;
import star.cadmodeler.*;
import star.kwturb.*;
import star.meshing.*;

public class $filename extends StarMacro {

  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();


  Solution solution_0 = 
      simulation_0.getSolution();

    solution_0.clearSolution();

    
    CadModel cadModel_0 = 
      ((CadModel) simulation_0.get(SolidModelManager.class).getObject(\"3D-CAD Model 1\"));

    ScalarQuantityDesignParameter scalarQuantityDesignParameter_0 = 
      ((ScalarQuantityDesignParameter) cadModel_0.getDesignParameterManager().getObject(\"AoA\"));

    Units units_1 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"deg\"));

    scalarQuantityDesignParameter_0.getQuantity().setValueAndUnits($AoA, units_1);

    cadModel_0.update();

    SolidModelPart solidModelPart_0 = 
      ((SolidModelPart) simulation_0.get(SimulationPartManager.class).getPart(\"Body 1\"));

    simulation_0.get(SimulationPartManager.class).updateParts(new NeoObjectVector(new Object[] {solidModelPart_0}));

    MeshOperationPart meshOperationPart_0 = 
      ((MeshOperationPart) simulation_0.get(SimulationPartManager.class).getPart(\"Subtract\"));

    simulation_0.get(SimulationPartManager.class).updateParts(new NeoObjectVector(new Object[] {meshOperationPart_0}));

    simulation_0.get(MeshOperationManager.class).executeAll();

    ScalarGlobalParameter scalarGlobalParameter_0 = 
      ((ScalarGlobalParameter) simulation_0.get(GlobalParameterManager.class).getObject(\"u_ts\"));

    Units units_2 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"m/s\"));

    scalarGlobalParameter_0.getQuantity().setValueAndUnits($Uinf, units_2);

    PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum(\"Physics 1\"));

    TurbulenceIntensityProfile turbulenceIntensityProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(TurbulenceIntensityProfile.class);

    Units units_3 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"\"));

    turbulenceIntensityProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($TI, units_3);

    TurbulentViscosityRatioProfile turbulentViscosityRatioProfile_0 = 
      physicsContinuum_0.getInitialConditions().get(TurbulentViscosityRatioProfile.class);

    turbulentViscosityRatioProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($μr, units_3);

    GammaReThetaTransitionModel gammaReThetaTransitionModel_0 = 
      physicsContinuum_0.getModelManager().getModel(GammaReThetaTransitionModel.class);

    gammaReThetaTransitionModel_0.getConset1().setValueAndUnits($C1, units_3);

    gammaReThetaTransitionModel_0.setS1($s1);

    SstKwTurbModel sstKwTurbModel_0 = 
      physicsContinuum_0.getModelManager().getModel(SstKwTurbModel.class);

    sstKwTurbModel_0.setSigma_w1($σw1);

    sstKwTurbModel_0.setBetaStar($βstar);

    sstKwTurbModel_0.setA1($α1);

    StepStoppingCriterion stepStoppingCriterion_0 = 
      ((StepStoppingCriterion) simulation_0.getSolverStoppingCriterionManager().getSolverStoppingCriterion(\"Maximum Steps\"));


    stepStoppingCriterion_0.setIsUsed(true);

     IntegerValue integerValue_0 = 
      stepStoppingCriterion_0.getMaximumNumberStepsObject();

    integerValue_0.getQuantity().setValue(10000);

    simulation_0.getSimulationIterator().run();


    XYPlot xYPlot_1 = 
      ((XYPlot) simulation_0.getPlotManager().getPlot(\"Cf_top\"));

    xYPlot_1.export(resolvePath(\"/home/romain.ospital/STARCCMp/UQ/Results/$CFNAME.csv\"), \",\");

    MonitorPlot monitorPlot_0 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot(\"CL Monitor Plot\"));

    monitorPlot_0.export(resolvePath(\"/home/romain.ospital/STARCCMp/UQ/Results/$CLNAME.csv\"), \",\");

    MonitorPlot monitorPlot_1 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot(\"CD Monitor Plot\"));

    monitorPlot_1.export(resolvePath(\"/home/romain.ospital/STARCCMp/UQ/Results/$CDNAME.csv\"), \",\");
  }
}

"
macro_text = str0 

write(filename*".java", macro_text);

end

function convert_num2str(n::Int64)
  if n < 1 || n >= 26^3
    error("Input must be between 1 and 26^3 - 1 (inclusive).")
  end

  # Initialize an array for the letter sequence
  letters = ['A', 'A', 'A']

  # Convert to 0-based index
  n -= 1

  for i in 1:3
      letters[4 - i] = Char('A' + (n % 26))
      n ÷= 26
  end

  return join(letters)
end



test_run = DataFrame((XLSX.readtable("Sampling_DU89_$(Re)_$(AoA_mean).xlsx","Sheet1")))
Nsamples = size(test_run)[1]



sbatch_string = open("string_sbatch.txt", "w")
alpha_nums = open("alpha_nums.txt", "w")

convert_num2str(1)


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
  str_nn = Int(i)
  println(str_nn)
  alpha_string = convert_num2str(str_nn)
  
  tstring = "starccm+ -batchsystem slurm -batch Macros/Sim$alpha_string.java DU89_A1_Re500000_3D_ReGammaTheta_Base.sim\n"
  write(sbatch_string, tstring);
  write(alpha_nums,"$(alpha_string)\n");
  filename= "Sim"*alpha_string
  write_macro(Uinf,TI,AoA,μr,σw1,α1,βstar, s1, C1, filename)
end

close(sbatch_string)
close(alpha_nums)

