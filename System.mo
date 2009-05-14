

package System "DMFC systems" 
  partial model Reference "The reference DMFC system, no control applied" 
    
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    inner parameter Real RH_env = 60;
    
    Flow.Mixer mixer annotation (extent=[-10,-70; 10,-50]);
    Flow.Pump pump "The anodic-loop pump" 
              annotation (extent=[-30,-66; -42,-54]);
    annotation (Diagram);
    replaceable Flow.AbstractCooler anodeCooler "The solution-loop cooler" 
                  annotation (extent=[10,-30; 30,-10]);
    Flow.PureMethanolSource pureMethanolSource 
      "A substitute for an actual tank"   annotation (extent=[30,-96; 42,-84]);
    Flow.Pump fuelPump "The smaller fuel pump" 
                  annotation (extent=[20,-96; 8,-84]);
    Flow.Separator degasser "The CO2-degasser" 
                        annotation (extent=[38,-30; 58,-10]);
    Flow.SinkPort co2sink "The gas outlet of the degasser" 
                      annotation (extent=[68,-16; 76,-8]);
    replaceable Flow.FuelCell fuelCell 
                      annotation (extent=[-50,-12; -14,22]);
    Flow.EnvironmentPort environment "The air from the environment" 
      annotation (extent=[-100,0; -80,20]);
    Flow.GasFlowController blower "The mass-flow controller" 
      annotation (extent=[-76,16; -64,4], rotation=270);
    replaceable Flow.AbstractCooler cathodeCooler "The cathode-side cooler" 
                  annotation (extent=[32,30; 52,50]);
    Flow.Separator condenser "The water-recuperating unit" 
                        annotation (extent=[68,30; 86,50]);
    Flow.SinkPort airSink "The gas outlet of the condenser" 
                      annotation (extent=[92,46; 100,54]);
    replaceable Modelica.Electrical.Analog.Interfaces.OnePort load 
      "Load connected to the cell"       annotation (extent=[-52,84; -40,96]);
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (extent=[-8,24; 8,40]);
    Modelica.Electrical.Analog.Sensors.CurrentSensor amperometer 
      "Current in external circuit" annotation (extent=[-32,80; -12,100]);
  equation 
    connect(pump.inlet, mixer.outlet) 
      annotation (points=[-36,-60; -8,-60],    style(color=62, rgbcolor={0,127,
            127}));
    connect(pureMethanolSource.outlet, fuelPump.inlet) 
      annotation (points=[36,-90; 14,-90],    style(color=62, rgbcolor={0,127,
            127}));
    connect(environment.outlet, blower.inlet) 
                                         annotation (points=[-81,10; -70,10],
                              style(color=62, rgbcolor={0,127,127}));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-64,10; 
          -50,10; -50,10.1],         style(color=62, rgbcolor={0,127,127}));
    connect(cathodeCooler.outlet, condenser.inlet) 
      annotation (points=[51.4,40; 68,40],  style(color=62, rgbcolor={0,127,127}));
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (points=[14,-84; 0,-84; 0,-68; 6.10623e-16,-68],
                style(color=62, rgbcolor={0,127,127}));
    connect(anodeCooler.outlet, degasser.inlet) 
                                           annotation (points=[29.4,-20; 38,-20],
                                                     style(color=62, rgbcolor={
            0,127,127}));
    connect(degasser.gasOutlet, co2sink.inlet)    annotation (points=[55,-16;
          62,-16; 62,-12; 68.4,-12],
        style(color=62, rgbcolor={0,127,127}));
    connect(condenser.gasOutlet, airSink.inlet) 
      annotation (points=[83.3,44; 88,44; 88,50; 92.4,50],
                                           style(color=62, rgbcolor={0,127,127}));
    connect(condenser.liquidOutlet, mixer.waterInlet) 
      annotation (points=[83.3,36; 88,36; 88,-60; 8,-60],
                                                        style(color=62,
          rgbcolor={0,127,127}));
    connect(degasser.liquidOutlet, mixer.loopInlet) annotation (points=[55,-24; 
          60,-24; 60,-40; 0,-40; 0,-52; 6.10623e-16,-52],   style(color=62,
          rgbcolor={0,127,127}));
    connect(fuelCell.anode_outlet, anodeCooler.inlet) 
                                                 annotation (points=[-14,-0.1; 
          -14,0; -4,0; -4,-20; 10.6,-20],
                                      style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=62,
        rgbfillColor={0,127,127},
        fillPattern=1));
    connect(fuelCell.anode_inlet, pump.outlet) annotation (points=[-50,-0.1; 
          -50,0; -60,0; -60,-54; -36,-54],
                              style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=62,
        rgbfillColor={0,127,127},
        fillPattern=1));
    connect(cathodeCooler.inlet, fuelCell.cathode_outlet) annotation (points=[32.6,40; 
          12,40; 12,14; -14,14; -14,10.1],    style(color=62, rgbcolor={0,127,
            127}));
    connect(fuelCell.minus, ground.p) annotation (points=[-21.2,15.2; -21.2,40; 
          2.10942e-16,40], style(color=3, rgbcolor={0,0,255}));
    connect(amperometer.p, load.n) 
      annotation (points=[-32,90; -40,90], style(color=3, rgbcolor={0,0,255}));
    connect(fuelCell.minus, amperometer.n) annotation (points=[-21.2,15.2; 
          -21.2,40; 0,40; 0,90; -12,90], style(color=3, rgbcolor={0,0,255}));
    connect(fuelCell.plus, load.p) annotation (points=[-42.8,15.2; -42.8,40; 
          -64,40; -64,90; -52,90], style(color=3, rgbcolor={0,0,255}));
  end Reference;
  
  model Reference_NoControl "Sets manipulable variables with parameters" 
    extends Reference(
      redeclare Modelica.Electrical.Analog.Sources.ConstantCurrent load(I=5),
      redeclare Flow.DiscretisedCooler cathodeCooler(exchanger(A=3.46E-2, U=82)),
      redeclare Flow.DiscretisedCooler anodeCooler,
      redeclare Flow.ConstantVoltageFuelCell fuelCell,
      mixer(
        c(fixed=true),
        T(fixed=true),
        V(fixed=true)));
    
    import Modelica.SIunits.VolumeFlowRate;
    
    parameter VolumeFlowRate V_fuel = 5E-8/60;
    parameter VolumeFlowRate V_anode = 10E-6/60;
    parameter VolumeFlowRate V_cathode = 500E-6/60;
    parameter VolumeFlowRate V_cooler = 10E-3/60;
    parameter VolumeFlowRate V_condenser = 2E-3/60;
    
  equation 
    V_fuel = fuelPump.V;
    V_anode = pump.V;
    V_cathode = blower.V;
    V_cooler = anodeCooler.mfc.V;
    V_condenser = cathodeCooler.mfc.V;
    
  end Reference_NoControl;
  
  annotation (uses(Modelica(version="2.2.1")));
  
  model Reference_Control 
    "The reference DMFC system derived from the one to be presented at ASME FC09" 
    extends Reference(redeclare Flow.ConstantVoltageFuelCell fuelCell,mixer(V(fixed=true),c(fixed=true),T(fixed=true)),
      redeclare Modelica.Electrical.Analog.Sources.ConstantCurrent load(I=5),
      redeclare Flow.DiscretisedCooler cathodeCooler(exchanger(A=3.46E-2, U=82)),
      redeclare Flow.DiscretisedCooler anodeCooler);
    
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Current;
    import Modelica.SIunits.Volume;
    import Modelica.SIunits.Time;
    
    import Flow.MolarFlow;
    
    import Thermo.Molecules.Oxygen;
    import Thermo.Molecules.Water;
    import Thermo.Molecules.Methanol;
    import Thermo.Phases.Liquid;
    import Thermo.rho;
    import Thermo.mw;
    import Thermo.K;
    
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    
  public 
    Control.CathodeLambdaControl K_cath "Cathode lambda controller" 
      annotation (extent=[-74,24; -66,34], rotation=270);
    Control.FuelControl K_fuel annotation (extent=[-16,-94; -4,-86]);
    Control.WaterControl K_cond annotation (extent=[26,8; 38,16], rotation=0);
    Control.AnodeLambdaControl K_an(lambda=5) 
                                    annotation (extent=[-70,-64; -60,-56]);
    Modelica.Blocks.Sources.RealExpression C(y=1000) 
      "Target concentration, mol/m^3" annotation (extent=[-100,52; -80,72]);
    Control.TemperatureControl K annotation (extent=[-10,-38; 2,-26]);
  equation 
    
    annotation (experiment(StopTime=7200), experimentSetupOutput,
      Diagram);
    connect(amperometer.i, K_cath.I) annotation (points=[-22,80; -22,76; -67.6,
          76; -67.6,35],
                   style(
        color=3,
        rgbcolor={0,0,255},
        pattern=3));
    connect(K_fuel.I, amperometer.i) annotation (points=[-17.2,-87.6; -18,-88; 
          -90,-88; -90,40; -68,40; -68,76; -22,76; -22,80], style(
        color=3,
        rgbcolor={0,0,255},
        pattern=3));
    connect(blower.V, K_cath.V) annotation (points=[-70,16; -70,23], style(
        color=2,
        rgbcolor={0,255,0},
        pattern=3));
    connect(cathodeCooler.T_ref, K_cond.T_ref) annotation (points=[42,37; 42,12; 
          39.2,12],          style(
        color=1,
        rgbcolor={255,0,0},
        pattern=3));
    connect(K_cond.p_mix, mixer.p) annotation (points=[24.8,14.4; 18,14; 18,-74; 
          -4.7,-74; -4.7,-66.9], style(
        color=78,
        rgbcolor={127,0,127},
        pattern=3));
    connect(condenser.T, K_cond.T_cond) annotation (points=[86.9,40; 94,40; 94,
          -2; 20,-2; 20,9.6; 24.8,9.6],
                                      style(
        color=1,
        rgbcolor={255,0,0},
        pattern=3));
    connect(degasser.T, K_fuel.T_deg) annotation (points=[59,-20; 66,-20; 66,
          -98; -30,-98; -30,-92.4; -17.2,-92.4], style(
        color=1,
        rgbcolor={255,0,0},
        pattern=3));
    connect(pump.V, K_an.V) annotation (points=[-42,-60; -59,-60],   style(
        color=2,
        rgbcolor={0,255,0},
        pattern=3));
    connect(K_an.I, amperometer.i) annotation (points=[-71,-57.6; -90,-58; -90,
          40; -68,40; -68,76; -22,76; -22,80],     style(
        color=3,
        rgbcolor={0,0,255},
        pattern=3));
    connect(C.y, K_cath.c) annotation (points=[-79,62; -76,62; -76,48; -72.4,48; 
          -72.4,35],
                   style(
        color=62,
        rgbcolor={0,127,127},
        pattern=3));
    connect(C.y, K_an.c) annotation (points=[-79,62; -76,62; -76,-62.4; -71,
          -62.4], style(
        color=62,
        rgbcolor={0,127,127},
        pattern=3));
    connect(C.y, K_fuel.c) annotation (points=[-79,62; -76,62; -76,-90; -17.2,
          -90], style(
        color=62,
        rgbcolor={0,127,127},
        pattern=3));
    connect(K_cath.V, K_cond.V_cath) annotation (points=[-70,23; -70,18; 16,18;
          16,12; 24.8,12], style(
        color=2,
        rgbcolor={0,255,0},
        pattern=3));
    connect(K_fuel.V, fuelPump.V) annotation (points=[-2.8,-90; 8,-90], style(
        color=2,
        rgbcolor={0,255,0},
        pattern=3));
    connect(K.T_m, fuelCell.T) annotation (points=[-11.2,-32; -20,-32; -20,-20; 
          -8,-20; -8,5.34; -12.2,5.34], style(
        color=1, 
        rgbcolor={255,0,0}, 
        pattern=3));
    connect(anodeCooler.T_ref, K.T_deg_ref) annotation (points=[20,-23; 20,-32; 
          3.2,-32], style(
        color=1, 
        rgbcolor={255,0,0}, 
        pattern=3));
  end Reference_Control;
  
  partial model CoolingIntegration 
    "The reference DMFC system, no control applied" 
    
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    inner parameter Real RH_env = 60;
    
    Flow.Mixer mixer annotation (extent=[-10,-42; 10,-22]);
    Flow.Pump pump "The anodic-loop pump" 
              annotation (extent=[-50,-38; -38,-26]);
    annotation (Diagram);
    replaceable Flow.Cooler cooler "The outlet cooler" 
                  annotation (extent=[0,-6; 22,16]);
    Flow.PureMethanolSource pureMethanolSource 
      "A substitute for an actual tank"   annotation (extent=[30,-58; 42,-46]);
    Flow.Pump fuelPump "The smaller fuel pump" 
                  annotation (extent=[8,-58; 20,-46]);
    Flow.Separator degasser "The CO2-degasser" 
                        annotation (extent=[32,-6; 52,16]);
    Flow.SinkPort co2sink "The gas outlet of the degasser" 
                      annotation (extent=[54,20; 62,28]);
    replaceable Flow.FuelCell fuelCell 
                      annotation (extent=[-50,-12; -14,22]);
    Flow.EnvironmentPort environment "The air from the environment" 
      annotation (extent=[-82,22; -62,42]);
    Flow.GasFlowController blower "The mass-flow controller" 
      annotation (extent=[-74,2; -66,10]);
    replaceable Modelica.Electrical.Analog.Interfaces.OnePort load 
      "Load connected to the cell"       annotation (extent=[-52,34; -32,54]);
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (extent=[-32,-30; -18,-16]);
  equation 
    connect(pump.inlet, mixer.outlet) 
      annotation (points=[-44,-32; -8,-32],    style(color=62, rgbcolor={0,127,
            127}));
    connect(pureMethanolSource.outlet, fuelPump.inlet) 
      annotation (points=[36,-52; 14,-52],    style(color=62, rgbcolor={0,127,
            127}));
    connect(environment.outlet, blower.inlet) 
                                         annotation (points=[-63,27; -89.5,27;
          -89.5,6; -70,6],    style(color=62, rgbcolor={0,127,127}));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-70,10;
          -50,10; -50,10.1],         style(color=62, rgbcolor={0,127,127}));
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (points=[14,-46; 0,-46; 0,-40; 6.10623e-16,-40],
                style(color=62, rgbcolor={0,127,127}));
    connect(cooler.outlet, degasser.inlet) annotation (points=[21.34,5; 32,5],
                                                     style(color=62, rgbcolor={
            0,127,127}));
    connect(degasser.gasOutlet, co2sink.inlet)    annotation (points=[49,9.4;
          48,9.4; 48,24; 54.4,24],
        style(color=62, rgbcolor={0,127,127}));
    connect(fuelCell.anode_outlet, cooler.inlet) annotation (points=[-14,-0.1;
          -14,0; -4,0; -4,5; 0.66,5], style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=62,
        rgbfillColor={0,127,127},
        fillPattern=1));
    connect(fuelCell.anode_inlet, pump.outlet) annotation (points=[-50,-0.1;
          -50,0; -76,0; -76,-26; -44,-26],
                              style(
        color=62,
        rgbcolor={0,127,127},
        fillColor=62,
        rgbfillColor={0,127,127},
        fillPattern=1));
    connect(load.n, fuelCell.plus)   annotation (points=[-32,44; -32,29.6; -32,
          15.2; -42.8,15.2],
        style(
        color=3,
        rgbcolor={0,0,255},
        pattern=0,
        fillColor=43,
        rgbfillColor={255,85,85},
        fillPattern=1));
    connect(load.p, fuelCell.minus)   annotation (points=[-52,44; -58,44; -58,
          -16; -21.2,-16; -21.2,15.2],
                                   style(
        color=3,
        rgbcolor={0,0,255},
        pattern=0,
        fillColor=43,
        rgbfillColor={255,85,85},
        fillPattern=1));
    connect(ground.p, fuelCell.minus) annotation (points=[-25,-16; -21.2,-16;
          -21.2,15.2],
                 style(
        color=3,
        rgbcolor={0,0,255},
        pattern=0,
        fillColor=43,
        rgbfillColor={255,85,85},
        fillPattern=1));
    connect(fuelCell.cathode_outlet, cooler.inlet) annotation (points=[-14,10.1;
          -4,10.1; -4,5; 0.66,5], style(color=62, rgbcolor={0,127,127}));
    connect(degasser.liquidOutlet, mixer.waterInlet) annotation (points=[49,0.6;
          48.5,0.6; 48.5,-32; 8,-32], style(color=62, rgbcolor={0,127,127}));
  end CoolingIntegration;
  
end System;
