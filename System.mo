package System "DMFC systems" 
  model Reference_NoControl "The reference DMFC system" 
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
    
    Flow.Mixer mixer annotation (extent=[-10,-74; 10,-54]);
    Flow.Pump pump "The anodic-loop pump" 
              annotation (extent=[-74,-70; -62,-58]);
    annotation (Diagram);
    Flow.Cooler cooler "The solution-loop cooler" 
                  annotation (extent=[4,-14; 26,6]);
    Flow.PureMethanolSource pureMethanolSource 
      "A substitute for an actual tank"   annotation (extent=[8,-92; 16,-84]);
    Flow.Pump fuelPump "The smaller fuel pump" 
                  annotation (extent=[-4,-92; 4,-84]);
    Flow.Separator degasser "The CO2-degasser" 
                        annotation (extent=[32,-14; 52,6]);
    Flow.SinkPort co2sink "The gas outlet of the degasser" 
                      annotation (extent=[62,-4; 70,4]);
    Flow.ConstantVoltageFuelCell fuelCell 
                      annotation (extent=[-46,-8; -24,14]);
    Flow.EnvironmentPort environment "The air from the environment" 
      annotation (extent=[-84,16; -64,36]);
    Flow.GasFlowController blower "The mass-flow controller" 
      annotation (extent=[-72,4; -64,12]);
    Flow.Cooler cathodeCooler "The cathode-side cooler" 
                  annotation (extent=[2,66; 24,86]);
    Flow.Separator condenser "The water-recuperating unit" 
                        annotation (extent=[56,66; 76,86]);
    Flow.SinkPort airSink "The gas outlet of the condenser" 
                      annotation (extent=[86,76; 94,84]);
  equation 
    connect(pump.inlet, mixer.outlet) 
      annotation (points=[-68.12,-64; -8,-64], style(color=62, rgbcolor={0,127,
            127}));
    connect(pureMethanolSource.c, fuelPump.inlet) 
      annotation (points=[12,-88; -0.08,-88], style(color=62, rgbcolor={0,127,
            127}));
    connect(environment.c, blower.inlet) annotation (points=[-83,21; -89.5,21;
          -89.5,8; -68.08,8], style(color=62, rgbcolor={0,127,127}));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-68,12;
          -58,12; -58,6.3; -46,6.3], style(color=62, rgbcolor={0,127,127}));
    connect(cathodeCooler.inlet, fuelCell.cathode_outlet) annotation (points=[2.66,76; 
          0,76; 0,6.3; -24,6.3], style(color=62, rgbcolor={0,127,127}));
    connect(cathodeCooler.outlet, condenser.inlet) 
      annotation (points=[23.34,76; 56,76], style(color=62, rgbcolor={0,127,127}));
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (points=[1.05471e-16,-84; 1.05471e-16,-78; 0,-72; 6.10623e-16,
          -72], style(color=62, rgbcolor={0,127,127}));
    connect(cooler.outlet, degasser.inlet) annotation (points=[25.34,-4; 27.005,
          -4; 27.005,-4; 28.67,-4; 28.67,-4; 32,-4], style(color=62, rgbcolor={
            0,127,127}));
    connect(degasser.gasOutlet, co2sink.flowPort) annotation (points=[49,
          6.66134e-16; 38.5,6.66134e-16; 38.5,3.88578e-17; 62.4,3.88578e-17], 
        style(color=62, rgbcolor={0,127,127}));
    connect(condenser.gasOutlet, airSink.flowPort) 
      annotation (points=[73,80; 86.4,80], style(color=62, rgbcolor={0,127,127}));
    connect(condenser.liquidOutlet, mixer.waterInlet) 
      annotation (points=[73,72; 74,72; 74,-64; 8,-64], style(color=62, 
          rgbcolor={0,127,127}));
    connect(degasser.liquidOutlet, mixer.loopInlet) annotation (points=[49,-8; 
          48,-8; 48,-46; 6.10623e-16,-46; 6.10623e-16,-56], style(color=62, 
          rgbcolor={0,127,127}));
    
    pump.V = anodeFlow;
    blower.V = cathodeFlow;
    fuelPump.V = fuelFlow;
    cooler.outletTemperature.T = T_cooler_out;
    cathodeCooler.outletTemperature.T = T_condenser;
    fuelCell.I = I_cell;
    
    connect(fuelCell.anode_outlet, cooler.inlet) annotation (points=[-24,-0.3; 
          -10,-0.3; -10,-4; 4.66,-4], style(
        color=62, 
        rgbcolor={0,127,127}, 
        fillColor=62, 
        rgbfillColor={0,127,127}, 
        fillPattern=1));
    connect(fuelCell.anode_inlet, pump.outlet) annotation (points=[-46,-0.3; 
          -68,-0.3; -68,-58], style(
        color=62, 
        rgbcolor={0,127,127}, 
        fillColor=62, 
        rgbfillColor={0,127,127}, 
        fillPattern=1));
  end Reference_NoControl;
  annotation (uses(Modelica(version="2.2.1")));
end System;
