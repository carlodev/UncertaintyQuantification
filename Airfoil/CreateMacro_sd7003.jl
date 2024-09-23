using FileIO
using XLSX, DataFrames



function write_macro(TI,μr,σw1,α1,βstar, s1, C1, filename::String)

CLNAME = "CL"*filename
CDNAME =  "CD"*filename
CFNAME =  "CF"*filename

macro_text = "
// Simcenter STAR-CCM+ macro: macro_sd7003.java
// Written by Simcenter STAR-CCM+ 18.02.008
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.turbulence.*;
import star.kwturb.*;

public class $filename extends StarMacro {

  public void execute() {
    execute0();
    execute1();

  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    Region region_0 = 
      simulation_0.getRegionManager().getRegion(\"Region\");

    Boundary boundary_0 = 
      region_0.getBoundaryManager().getBoundary(\"Body 1.inlet\");

    TurbulentViscosityRatioProfile turbulentViscosityRatioProfile_0 = 
      boundary_0.getValues().get(TurbulentViscosityRatioProfile.class);

    Units units_0 = 
      ((Units) simulation_0.getUnitsManager().getObject(\"\"));

    turbulentViscosityRatioProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($μr, units_0);

    TurbulenceIntensityProfile turbulenceIntensityProfile_0 = 
      boundary_0.getValues().get(TurbulenceIntensityProfile.class);

    turbulenceIntensityProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($TI, units_0);

    PhysicsContinuum physicsContinuum_0 = 
      ((PhysicsContinuum) simulation_0.getContinuumManager().getContinuum(\"Physics 1\"));

    GammaReThetaTransitionModel gammaReThetaTransitionModel_0 = 
      physicsContinuum_0.getModelManager().getModel(GammaReThetaTransitionModel.class);

    gammaReThetaTransitionModel_0.setS1($s1);

    gammaReThetaTransitionModel_0.getConset1().setValueAndUnits($C1, units_0);

    SstKwTurbModel sstKwTurbModel_0 = 
      physicsContinuum_0.getModelManager().getModel(SstKwTurbModel.class);

    sstKwTurbModel_0.setSigma_w1($σw1);

    sstKwTurbModel_0.setBetaStar($βstar);

    sstKwTurbModel_0.setA1($α1);

    simulation_0.getSimulationIterator().run();


    XYPlot xYPlot_1 = 
      ((XYPlot) simulation_0.getPlotManager().getPlot(\"Cf_top\"));

    xYPlot_1.export(resolvePath(\"/home/cbrunelli/STARCCMp/sd7003/UQ/Results/$CFNAME.csv\"), \",\");

    MonitorPlot monitorPlot_0 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot(\"CL Monitor Plot\"));

    monitorPlot_0.export(resolvePath(\"/home/cbrunelli/STARCCMp/sd7003/UQ/Results/$CLNAME.csv\"), \",\");

    MonitorPlot monitorPlot_1 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot(\"CD Monitor Plot\"));

    monitorPlot_1.export(resolvePath(\"/home/cbrunelli/STARCCMp/sd7003/UQ/Results/$CDNAME.csv\"), \",\");
  }

  private void execute1() {
  }
}

"

write(filename*".java", macro_text);

end

function convert_num2str(n::Int64)
  alphanumeric = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  t1 = Int(ceil(n/26^2))
  t2 = Int(ceil(n/26))
  t3 = mod(n-1,26)+1
  return alphanumeric[t1]*alphanumeric[t2]*alphanumeric[t3]  
end


test_run = DataFrame((XLSX.readtable("Airfoil/Sampling_sd7003.xlsx","Sheet1")))
Nsamples = size(test_run)[1]

sbatch_string = open("string_sbatch.txt", "w")
alpha_nums = open("alpha_nums.txt", "w")

for i = 1:1:Nsamples
  TI=test_run[i,:].TI
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
  write_macro(TI,μr,σw1,α1,βstar, s1, C1, filename)
end

close(sbatch_string)
close(alpha_nums)


