                                    /**
 * © Federico Zenith, 2008-2009.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


package System "DMFC systems" 
  partial model Reference "The reference DMFC system" 
    
    import Modelica.SIunits.Efficiency;
    import Units.MolarFlow;
    
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Units.Temperature T_env = 298.15;
    inner parameter Units.RelativeHumidity RH_env = 60;
    
    Efficiency eta_to_cell "Fraction of methanol consumed in the cell";
    Efficiency eta_system "Overall system efficiency";
    
    Flow.UnitOperations.Mixer mixer 
                     annotation (extent=[-10,-70; 10,-50]);
    Flow.Measurements.LiquidPump pump "The anodic-loop pump" 
              annotation (extent=[-30,-66; -42,-54]);
    annotation (Diagram, Documentation(info="<html>
<p>This is a generic reference system, with no process integration
whatsoever. Some components, such as the fuel cell, are abstract and
must be specialised in subclasses.</p>
</html>"));
    replaceable Flow.UnitOperations.Coolers.Abstract anodeCooler 
      "The solution-loop cooler" 
                  annotation (extent=[10,-30; 30,-10]);
    Flow.Sources.Methanol pureMethanolSource "A substitute for an actual tank" 
                                          annotation (extent=[30,-96; 42,-84]);
    Flow.Measurements.LiquidPump fuelPump "The smaller fuel pump" 
                  annotation (extent=[20,-96; 8,-84]);
    Flow.UnitOperations.Separator degasser "The CO2-degasser" 
                        annotation (extent=[38,-30; 58,-10]);
    Flow.Sink co2sink "The gas outlet of the degasser" 
                      annotation (extent=[68,-16; 76,-8]);
    replaceable Flow.UnitOperations.FuelCells.Abstract fuelCell 
                      annotation (extent=[-50,-12; -14,22]);
    Flow.Sources.Environment environment "The air from the environment" 
      annotation (extent=[-100,0; -80,20]);
    Flow.Measurements.GasFlowController blower "The mass-flow controller" 
      annotation (extent=[-76,16; -64,4], rotation=270);
    replaceable Flow.UnitOperations.Coolers.Abstract cathodeCooler 
      "The cathode-side cooler" 
                  annotation (extent=[32,30; 52,50]);
    Flow.UnitOperations.Separator condenser "The water-recuperating unit" 
                        annotation (extent=[68,30; 86,50]);
    Flow.Sink airSink "The gas outlet of the condenser" 
                      annotation (extent=[92,46; 100,54]);
    replaceable Modelica.Electrical.Analog.Interfaces.OnePort load 
      "Load connected to the cell"       annotation (extent=[-52,84; -40,96]);
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (extent=[-8,24; 8,40]);
    Modelica.Electrical.Analog.Sensors.CurrentSensor amperometer 
      "Current in external circuit" annotation (extent=[-32,80; -12,100]);
    
  protected 
    MolarFlow inCell = fuelCell.anode_inlet.n[1] - (-fuelCell.anode_outlet.n[1]) 
      "Methanol consumed in the cell";
    MolarFlow inDeg =  -degasser.gasOutlet.n[1] "Methanol lost in the degasser";
    
  equation 
    eta_to_cell = inCell / (inCell + inDeg);
    eta_system = eta_to_cell * fuelCell.eta_total;
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
  
  model Reference_NoControl "The reference system with manual control" 
    extends Reference(
      redeclare Modelica.Electrical.Analog.Sources.ConstantCurrent load(I=5),
      redeclare Flow.UnitOperations.Coolers.Simple cathodeCooler,
      redeclare Flow.UnitOperations.Coolers.Simple anodeCooler,
      redeclare Flow.UnitOperations.FuelCells.ConstantVoltage fuelCell,
                                                                mixer(c(fixed=true),T(fixed=true),V(fixed=true)));
    
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    
    parameter VolumeFlowRate V_fuel = 5E-8/60;
    parameter VolumeFlowRate V_anode = 10E-6/60;
    parameter VolumeFlowRate V_cathode = 500E-6/60;
    parameter Temperature T_cooler = 310;
    parameter Temperature T_condenser = 320;
    
  equation 
    V_fuel = fuelPump.V;
    V_anode = pump.V;
    V_cathode = blower.V;
    T_cooler = anodeCooler.T_ref;
    T_condenser = cathodeCooler.T_ref;
    
    annotation (Documentation(info="<html>
<p>This simple specialisation of the generic reference-system model
allows to set the manipulable variables by hand as parameters, and
see what happens.</p>
</html>"));
  end Reference_NoControl;
  
  annotation (uses(Modelica(version="2.2.1")), Documentation(info="<html>
<p>A collection of complete DMFC systems with and without control algorithms.</p>
</html>"));
  
  model Reference_Control 
    "The reference DMFC system derived from the one to be presented at ASME FC09" 
    extends Reference(redeclare Flow.UnitOperations.FuelCells.Thevenin fuelCell,
      redeclare Modelica.Electrical.Analog.Sources.ConstantCurrent load(I=5),
      redeclare Flow.UnitOperations.Coolers.Simple cathodeCooler,
      redeclare Flow.UnitOperations.Coolers.Simple anodeCooler,
      mixer(T(fixed=true), c(fixed=true)));
    
  public 
    Control.CathodeLambdaControl K_cath "Cathode lambda controller" 
      annotation (extent=[-74,24; -66,34], rotation=270);
    Control.FuelControl K_fuel annotation (extent=[-16,-94; -4,-86]);
    Control.WaterControl K_cond annotation (extent=[26,8; 38,18], rotation=0);
    Control.AnodeLambdaControl K_an annotation (extent=[-70,-64; -60,-56]);
    Control.TemperatureControl K_deg 
                                 annotation (extent=[-16,-34; -4,-22]);
    annotation (experiment(StopTime=7200), experimentSetupOutput,
      Diagram,
      Documentation(info="<html>
<p>This specialisation of the reference system implements a series of
controllers. Note that controller connections are dotted and colour-coded.</p>
</html>"));
  equation 
    connect(amperometer.i, K_cath.I) annotation (points=[-22,80; -22,76; -70,76;
          -70,35], style(
        color=3,
        rgbcolor={0,0,255},
        pattern=3));
    connect(K_fuel.I, amperometer.i) annotation (points=[-17.2,-87.6; -18,-88;
          -90,-88; -90,40; -70,40; -70,76; -22,76; -22,80], style(
        color=3,
        rgbcolor={0,0,255},
        pattern=3));
    connect(blower.V, K_cath.V) annotation (points=[-70,16; -70,23], style(
        color=2,
        rgbcolor={0,255,0},
        pattern=3));
    connect(cathodeCooler.T_ref, K_cond.T_ref) annotation (points=[42,37; 42,13;
          39.2,13],          style(
        color=1,
        rgbcolor={255,0,0},
        pattern=3));
    connect(K_cond.p_mix, mixer.p) annotation (points=[24.8,16; 18,16; 18,8; 4,
          8; 4,-44; 16,-44; 16,-74; -4.7,-74; -4.7,-66.9],
                                 style(
        color=78,
        rgbcolor={127,0,127},
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
    connect(K_an.I, amperometer.i) annotation (points=[-71,-60; -90,-60; -90,40;
          -70,40; -70,76; -22,76; -22,80],         style(
        color=3,
        rgbcolor={0,0,255},
        pattern=3));
    connect(K_cath.V, K_cond.V_cath) annotation (points=[-70,23; -70,18; 16,18;
          16,13; 24.8,13], style(
        color=2,
        rgbcolor={0,255,0},
        pattern=3));
    connect(K_fuel.V, fuelPump.V) annotation (points=[-2.8,-90; 8,-90], style(
        color=2,
        rgbcolor={0,255,0},
        pattern=3));
    connect(K_deg.T_m, fuelCell.T) 
                               annotation (points=[-17.2,-28; -20,-28; -20,-16; 
          -8,-16; -8,5.34; -12.2,5.34], style(
        color=1,
        rgbcolor={255,0,0},
        pattern=3));
    connect(anodeCooler.T_ref, K_deg.T_deg_ref) 
                                            annotation (points=[20,-23; 20,
          -26.8; -2.8,-26.8],
                    style(
        color=1,
        rgbcolor={255,0,0},
        pattern=3));
    connect(cathodeCooler.T_process_out, K_cond.T_cond) annotation (points=[
          51.4,38.6; 60,38; 60,0; 20,0; 20,10; 24.8,10], style(
        color=1,
        rgbcolor={255,0,0},
        pattern=3));
    connect(K_deg.isSaturated, anodeCooler.isSaturated) annotation (points=[-2.8,
          -29.2; 22,-29.2; 22,-24; 23,-23],      style(
        color=5,
        rgbcolor={255,0,255},
        pattern=3));
  end Reference_Control;
  
  partial model Mingled "A DMFC system with outlet mingling" 
    
    import Modelica.SIunits.Efficiency;
    import Units.MolarFlow;
    
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Units.Temperature T_env = 298.15;
    inner parameter Units.RelativeHumidity RH_env = 60;
    
    Efficiency eta_to_cell "Fraction of methanol consumed in the cell";
    Efficiency eta_system "Overall system efficiency";
    
    Flow.UnitOperations.Mixer mixer 
                     annotation (extent=[-10,-70; 10,-50]);
    Flow.Measurements.LiquidPump pump "The anodic-loop pump" 
              annotation (extent=[-30,-66; -42,-54]);
    annotation (Diagram, Documentation(info="<html>
<p>This is a generic reference system, with no process integration
whatsoever. Some components, such as the fuel cell, are abstract and
must be specialised in subclasses.</p>
</html>"));
    replaceable Flow.UnitOperations.Coolers.Abstract cooler "The system cooler"
                  annotation (extent=[14,-4; 32,14]);
    Flow.Sources.Methanol pureMethanolSource "A substitute for an actual tank" 
                                          annotation (extent=[30,-96; 42,-84]);
    Flow.Measurements.LiquidPump fuelPump "The smaller fuel pump" 
                  annotation (extent=[20,-96; 8,-84]);
    Flow.UnitOperations.Separator separator "The loop separator" 
                        annotation (extent=[48,-6; 68,16]);
    Flow.Sink co2sink "The gas outlet of the degasser" 
                      annotation (extent=[80,16; 88,24]);
    replaceable Flow.UnitOperations.FuelCells.Abstract fuelCell 
                      annotation (extent=[-50,-12; -14,22]);
    Flow.Sources.Environment environment "The air from the environment" 
      annotation (extent=[-100,0; -80,20]);
    Flow.Measurements.GasFlowController blower "The mass-flow controller" 
      annotation (extent=[-76,16; -64,4], rotation=270);
    replaceable Modelica.Electrical.Analog.Interfaces.OnePort load 
      "Load connected to the cell"       annotation (extent=[-52,84; -40,96]);
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (extent=[-8,24; 8,40]);
    Modelica.Electrical.Analog.Sensors.CurrentSensor amperometer 
      "Current in external circuit" annotation (extent=[-32,80; -12,100]);
  protected 
    MolarFlow inCell = fuelCell.anode_inlet.n[1] - (-fuelCell.anode_outlet.n[1]) 
      "Methanol consumed in the cell";
    MolarFlow inDeg =  -separator.gasOutlet.n[1] 
      "Methanol lost in the separator";
    
  equation 
    eta_to_cell = inCell / (inCell + inDeg);
    eta_system = eta_to_cell * fuelCell.eta_total;
    
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
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (points=[14,-84; 0,-84; 0,-68; 6.10623e-16,-68],
                style(color=62, rgbcolor={0,127,127}));
    connect(cooler.outlet, separator.inlet) 
                                           annotation (points=[31.46,5; 48,5],
                                                     style(color=62, rgbcolor={
            0,127,127}));
    connect(fuelCell.anode_outlet, cooler.inlet) annotation (points=[-14,-0.1; 
          -14,0; 0,0; 0,5; 14.54,5],  style(
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
    connect(fuelCell.minus, ground.p) annotation (points=[-21.2,15.2; -21.2,40; 
          2.10942e-16,40], style(color=3, rgbcolor={0,0,255}));
    connect(amperometer.p, load.n) 
      annotation (points=[-32,90; -40,90], style(color=3, rgbcolor={0,0,255}));
    connect(fuelCell.minus, amperometer.n) annotation (points=[-21.2,15.2; 
          -21.2,40; 0,40; 0,90; -12,90], style(color=3, rgbcolor={0,0,255}));
    connect(fuelCell.plus, load.p) annotation (points=[-42.8,15.2; -42.8,40; 
          -64,40; -64,90; -52,90], style(color=3, rgbcolor={0,0,255}));
    connect(cooler.inlet, fuelCell.cathode_outlet) annotation (points=[14.54,5; 
          0,5; 0,10.1; -14,10.1], style(color=62, rgbcolor={0,127,127}));
    connect(separator.liquidOutlet, mixer.waterInlet) annotation (points=[65,0.6; 
          65,-60; 8,-60],      style(color=62, rgbcolor={0,127,127}));
    connect(separator.gasOutlet, co2sink.inlet) annotation (points=[65,9.4; 72,
          9.4; 72,20; 80.4,20], style(color=62, rgbcolor={0,127,127}));
  end Mingled;
  
  model Mingled_NoControl "The mingled system with manual control" 
    extends Mingled(
      redeclare Flow.UnitOperations.FuelCells.Thevenin fuelCell,
      redeclare Modelica.Electrical.Analog.Sources.ConstantCurrent load(I=5),
      redeclare Flow.UnitOperations.Coolers.Simple cooler,
      mixer(
        c(fixed=true),
        T(fixed=true),
        V(fixed=true)));
    
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    
    parameter VolumeFlowRate V_fuel = 5E-8/60;
    parameter VolumeFlowRate V_anode = 10E-6/60;
    parameter VolumeFlowRate V_cathode = 500E-6/60;
    parameter Temperature T_cooler = 320;
    
  equation 
    V_fuel = fuelPump.V;
    V_anode = pump.V;
    V_cathode = blower.V;
    T_cooler = cooler.T_ref;
    
    annotation (Diagram, experiment(StopTime=7200));
  end Mingled_NoControl;
end System;
