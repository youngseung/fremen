package System "DMFC systems" 
  model Reference "The reference DMFC system" 
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    inner parameter Real RH_env = 60;
    
    parameter VolumeFlowRate anodeFlow = 30e-6/60 
      "Flow rate for the anodic loop";
    parameter VolumeFlowRate fuelFlow = 30e-7/60 "Flow rate for the fuel";
    parameter VolumeFlowRate cathodeFlow = 30E-5/60 
      "Flow rate for the cathodic loop";
    parameter Temperature T_cooler_out = 350;
    parameter Temperature T_condenser = 320;
    parameter Modelica.SIunits.Current I_cell = 0 "Cell current";
    
    Flow.Mixer mixer annotation (extent=[-16,-40; 4,-20]);
    Flow.Pump pump "The anodic-loop pump" 
              annotation (extent=[-70,-36; -58,-24]);
    annotation (Diagram);
    Flow.Cooler cooler "The solution-loop cooler" 
                  annotation (extent=[-10,-14; 12,6]);
    Flow.PureMethanolSource pureMethanolSource 
      "A substitute for an actual tank"   annotation (extent=[2,-58; 10,-50]);
    Flow.Pump fuelPump "The smaller fuel pump" 
                  annotation (extent=[-10,-58; -2,-50]);
    Flow.Separator degasser "The CO2-degasser" 
                        annotation (extent=[18,-14; 38,6]);
    Flow.SinkPort co2sink "The gas outlet of the degasser" 
                      annotation (extent=[42,-4; 50,4]);
    Flow.ConstantVoltageFuelCell fuelCell 
                      annotation (extent=[-46,-8; -24,14]);
    Flow.EnvironmentPort environment "The air from the environment" 
      annotation (extent=[-84,16; -64,36]);
    Flow.GasFlowController blower "The mass-flow controller" 
      annotation (extent=[-72,4; -64,12]);
    Flow.Cooler cathodeCooler "The cathode-side cooler" 
                  annotation (extent=[-10,4; 12,24]);
    Flow.Separator condenser "The water-recuperating unit" 
                        annotation (extent=[44,4; 64,24]);
    Flow.SinkPort airSink "The gas outlet of the condenser" 
                      annotation (extent=[74,14; 82,22]);
  equation 
    connect(pump.inlet, mixer.outlet) 
      annotation (points=[-64.12,-30; -14,-30], style(pattern=0));
    connect(pureMethanolSource.c, fuelPump.inlet) 
      annotation (points=[6,-54; -6.08,-54], style(pattern=0));
    connect(environment.c, blower.inlet) annotation (points=[-83,21; -89.5,21;
          -89.5,8; -68.08,8], style(pattern=0));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-68,12;
          -58,12; -58,6.3; -46,6.3], style(pattern=0));
    connect(cathodeCooler.inlet, fuelCell.cathode_outlet) annotation (points=
          [-9.34,14; -12,14; -12,6.3; -24,6.3], style(pattern=0));
    connect(cathodeCooler.outlet, condenser.inlet) 
      annotation (points=[11.34,14; 44,14], style(pattern=0));
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (points=[-6,-50; -6,-44; -6,-38; -6,-38], style(pattern=0));
    connect(pump.outlet, fuelCell.anode_inlet) 
      annotation (points=[-64,-24; -64,-0.3; -46,-0.3], style(pattern=0));
    connect(cooler.outlet, degasser.inlet) annotation (points=[11.34,-4; 13.005,
          -4; 13.005,-4; 14.67,-4; 14.67,-4; 18,-4], style(pattern=0));
    pump.V = anodeFlow;
    
    blower.V = cathodeFlow;
    fuelPump.V = fuelFlow;
    cooler.outletTemperature.T = T_cooler_out;
    cathodeCooler.outletTemperature.T = T_condenser;
    fuelCell.I = I_cell;
    
    connect(fuelCell.anode_outlet, cooler.inlet) annotation (points=[-24,-0.08; 
          -17,-0.08; -17,-4; -9.34,-4], style(pattern=0));
    connect(degasser.gasOutlet, co2sink.flowPort) annotation (points=[35,
          6.66134e-16; 38.5,6.66134e-16; 38.5,3.88578e-17; 42.4,3.88578e-17],
        style(pattern=0));
    connect(condenser.gasOutlet, airSink.flowPort) 
      annotation (points=[61,18; 74.4,18], style(pattern=0));
    connect(condenser.liquidOutlet, mixer.waterInlet) 
      annotation (points=[61,10; 62,10; 62,-30; 2,-30], style(pattern=0));
    connect(degasser.liquidOutlet, mixer.loopInlet) annotation (points=[35,-8; 36,
          -8; 36,-14; -6,-14; -6,-22], style(pattern=0));
  end Reference;
  annotation (uses(Modelica(version="2.2.1")));
end System;
