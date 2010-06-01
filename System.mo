within ;
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
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    inner parameter Units.RelativeHumidity RH_env = 60;

    Efficiency eta_to_cell "Fraction of methanol consumed in the cell";
    Efficiency eta_system "Overall system efficiency";

    Flow.UnitOperations.Mixer mixer 
                     annotation (Placement(transformation(extent={{-10,-70},{10,
              -50}}, rotation=0)));
    Flow.Measurements.LiquidPump pump "The anodic-loop pump" 
              annotation (Placement(transformation(extent={{-30,-66},{-42,-54}},
            rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),      graphics),
                         Documentation(info="<html>
<p>This is a generic reference system, with no process integration
whatsoever. Some components, such as the fuel cell, are abstract and
must be specialised in subclasses.</p>
</html>"));
    replaceable Flow.UnitOperations.Coolers.Abstract anodeCooler
      "The solution-loop cooler" 
                  annotation (Placement(transformation(extent={{10,-30},{30,-10}},
            rotation=0)));
    Flow.Sources.Methanol pureMethanolSource "A substitute for an actual tank" 
                                          annotation (Placement(transformation(
            extent={{30,-96},{42,-84}}, rotation=0)));
    Flow.Measurements.LiquidPump fuelPump "The smaller fuel pump" 
                  annotation (Placement(transformation(extent={{20,-96},{8,-84}},
            rotation=0)));
    Flow.UnitOperations.Separator degasser "The CO2-degasser" 
                        annotation (Placement(transformation(extent={{38,-30},{
              58,-10}}, rotation=0)));
    replaceable Flow.UnitOperations.Stack.Abstract fuelCell 
                      annotation (Placement(transformation(extent={{-50,-12},{
              -14,22}}, rotation=0)));
    Flow.Sources.Environment environment "The air from the environment" 
      annotation (Placement(transformation(extent={{-100,0},{-80,20}}, rotation=
             0)));
    Flow.Measurements.GasFlowController blower "The mass-flow controller" 
      annotation (Placement(transformation(
          origin={-70,10},
          extent={{6,-6},{-6,6}},
          rotation=270)));
    replaceable Flow.UnitOperations.Coolers.Abstract cathodeCooler
      "The cathode-side cooler" 
                  annotation (Placement(transformation(extent={{32,30},{52,50}},
            rotation=0)));
    Flow.UnitOperations.Separator condenser "The water-recuperating unit" 
                        annotation (Placement(transformation(extent={{68,30},{
              86,50}}, rotation=0)));
    Flow.Sink airSink "The gas outlet of the condenser" 
                      annotation (Placement(transformation(extent={{92,64},{100,
              72}}, rotation=0)));
    replaceable Modelica.Electrical.Analog.Interfaces.TwoPin load
      "Load connected to the cell"       annotation (Placement(transformation(
            extent={{-52,84},{-40,96}}, rotation=0)));
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (Placement(transformation(extent={{-8,24},{8,40}}, rotation=0)));
    Modelica.Electrical.Analog.Sensors.CurrentSensor amperometer
      "Current in external circuit" annotation (Placement(transformation(extent=
             {{-32,80},{-12,100}}, rotation=0)));

  protected
    MolarFlow inCell = fuelCell.anode_inlet.n[1] - (-fuelCell.anode_outlet.n[1])
      "Methanol consumed in the cell";
    MolarFlow inDeg =  -degasser.gasOutlet.n[1] "Methanol lost in the degasser";

  public
    Flow.Measurements.MethanolInAir emissions
      "Measurement on methanol emissions"
      annotation (Placement(transformation(extent={{68,60},{88,80}})));
  equation
    eta_to_cell = inCell / (inCell + inDeg);
    eta_system = eta_to_cell * fuelCell.eta_total;
    connect(pump.inlet, mixer.outlet) 
      annotation (Line(points={{-36,-60},{-8,-60}}, color={0,127,127}));
    connect(pureMethanolSource.outlet, fuelPump.inlet) 
      annotation (Line(points={{36,-90},{14,-90}}, color={0,127,127}));
    connect(environment.outlet, blower.inlet) 
                                         annotation (Line(points={{-81,10},{-70,
            10}}, color={0,127,127}));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (Line(points={{-64,10},
            {-50,10},{-50,10.1}},         color={0,127,127}));
    connect(cathodeCooler.outlet, condenser.inlet) 
      annotation (Line(points={{51.4,40},{68,40}}, color={0,127,127}));
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (Line(points={{14,-84},{0,-84},{0,-68},{6.10623e-16,-68}},
          color={0,127,127}));
    connect(anodeCooler.outlet, degasser.inlet) 
                                           annotation (Line(points={{29.4,-20},
            {38,-20}}, color={0,127,127}));
    connect(condenser.liquidOutlet, mixer.waterInlet) 
      annotation (Line(points={{83.3,36},{88,36},{88,-60},{8,-60}}, color={0,
            127,127}));
    connect(degasser.liquidOutlet, mixer.loopInlet) annotation (Line(points={{55,-24},
            {60,-24},{60,-40},{0,-40},{0,-52},{6.10623e-16,-52}},         color=
           {0,127,127}));
    connect(fuelCell.anode_outlet, anodeCooler.inlet) 
                                                 annotation (Line(points={{-14,
            -0.1},{-14,0},{-4,0},{-4,-20},{10.6,-20}}, color={0,127,127}));
    connect(fuelCell.anode_inlet, pump.outlet) annotation (Line(points={{-50,
            -0.1},{-50,0},{-60,0},{-60,-54},{-36,-54}}, color={0,127,127}));
    connect(cathodeCooler.inlet, fuelCell.cathode_outlet) annotation (Line(
          points={{32.6,40},{12,40},{12,14},{-14,14},{-14,10.1}}, color={0,127,
            127}));
    connect(fuelCell.minus, ground.p) annotation (Line(points={{-21.2,15.2},{
            -21.2,40},{1.22125e-16,40}}, color={0,0,255}));
    connect(amperometer.p, load.n) 
      annotation (Line(points={{-32,90},{-40,90}}, color={0,0,255}));
    connect(fuelCell.minus, amperometer.n) annotation (Line(points={{-21.2,15.2},
            {-21.2,40},{0,40},{0,90},{-12,90}}, color={0,0,255}));
    connect(fuelCell.plus, load.p) annotation (Line(points={{-42.8,15.2},{-42.8,
            40},{-64,40},{-64,90},{-52,90}}, color={0,0,255}));
    connect(airSink.inlet, emissions.outlet) annotation (Line(
        points={{92.4,68},{85,68}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(emissions.inlet, degasser.gasOutlet) annotation (Line(
        points={{71,68},{64,68},{64,-16},{55,-16}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(condenser.gasOutlet, emissions.inlet) annotation (Line(
        points={{83.3,44},{83.3,56},{64,56},{64,68},{71,68}},
        color={0,127,127},
        smooth=Smooth.None));
  end Reference;

  model Reference_NoControl "The reference system with manual control"
    extends Reference(
      redeclare Modelica.Electrical.Analog.Sources.ConstantCurrent load(I=5),
      redeclare Flow.UnitOperations.Coolers.Simple cathodeCooler,
      redeclare Flow.UnitOperations.Coolers.Simple anodeCooler,
      redeclare Flow.UnitOperations.Stack.ConstantVoltage fuelCell,
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

  annotation (uses(
      Flow(version="1"),
      Modelica(version="3.1"),
      Control(version="1"),
      Units(version="1")),                     Documentation(info="<html>
<p>A collection of complete DMFC systems with and without control algorithms.</p>
</html>"),
    version="1",
    conversion(noneFromVersion=""));

  model Reference_Control
    "The reference DMFC system derived from the one to be presented at ASME FC09"
    extends Reference(redeclare Flow.UnitOperations.Stack.Thevenin fuelCell(cells=3),
      redeclare Load load(step(I=4, offset=3)),
      redeclare Flow.UnitOperations.Coolers.Simple cathodeCooler,
      redeclare Flow.UnitOperations.Coolers.Simple anodeCooler,
      mixer(T(fixed=true), c(fixed=true)));

      model Load
        extends Modelica.Electrical.Analog.Interfaces.TwoPin;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
              extent={{-100,-20},{100,20}},
              lineColor={0,0,255},
              textString="Load")}),   Diagram(coordinateSystem(preserveAspectRatio=true,
                extent={{-100,-100},{100,100}}), graphics),
        experiment(StopTime=8000));
        Modelica.Electrical.Analog.Sources.StepCurrent step(
        I=2,
        offset=5,
        startTime(displayUnit="h") = 3600) 
          annotation (Placement(transformation(extent={{-20,20},{20,60}})));
      Modelica.Electrical.Analog.Sources.SineCurrent sine(
        I=2,
        freqHz=2E-3,
        startTime(displayUnit="h") = 7200) 
        annotation (Placement(transformation(extent={{-20,-60},{20,-20}})));
      equation
        connect(step.p, p)          annotation (Line(
            points={{-20,40},{-60,40},{-60,5.55112e-16},{-100,5.55112e-16}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(sine.p, p) annotation (Line(
          points={{-20,-40},{-60,-40},{-60,5.55112e-16},{-100,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
        connect(sine.n, n) annotation (Line(
          points={{20,-40},{60,-40},{60,0},{100,0},{100,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
        connect(step.n, n) annotation (Line(
          points={{20,40},{60,40},{60,5.55112e-16},{100,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
      end Load;

  public
    Control.CathodeLambdaControl K_cath(c_est=1200,
      cells=3,
      aA=5.02E-9,
      b=0.24,
      lambda=6) "Cathode lambda controller" 
      annotation (Placement(transformation(
          origin={-70,29},
          extent={{-5,-4},{5,4}},
          rotation=270)));
    Control.ReferenceFuelControl K_fuel(
      cells=3,
      aA=5.02E-9,
      b=0.24)                  annotation (Placement(transformation(extent={{
              -16,-94},{-4,-86}}, rotation=0)));
    Control.WaterControl K_cond(T_0(displayUnit="K") = 305) 
                                annotation (Placement(transformation(extent={{28,4},{
              40,14}},        rotation=0)));
    Control.AnodeLambdaControl K_an(c_est_an=1200, c_est_mix=800,
      cells=3,
      aA=5.02E-9,
      b=0.24)                       annotation (Placement(transformation(extent=
             {{-70,-64},{-60,-56}}, rotation=0)));
    Control.TemperatureControl K_temp(T_FC_ref=340) 
                                 annotation (Placement(transformation(extent={{
              -16,-34},{-4,-22}}, rotation=0)));
    annotation (experiment(StopTime=10800),experimentSetupOutput,
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics),
      Documentation(info="<html>
<p>This specialisation of the reference system implements a series of
controllers. Note that controller connections are dotted and colour-coded.</p>
</html>"));

    output Modelica.SIunits.HeatFlowRate crossover_heat = 725000 * fuelCell.n_x;
    output Modelica.SIunits.HeatFlowRate heat_removal = - (fuelCell.cathode_outlet.H + fuelCell.anode_outlet.H);

  equation
    connect(amperometer.i, K_cath.I) annotation (Line(
        points={{-22,80},{-22,76},{-70,76},{-70,35}},
        color={0,0,255},
        pattern=LinePattern.Dot));
    connect(K_fuel.I, amperometer.i) annotation (Line(
        points={{-17.2,-87.6},{-18,-88},{-90,-88},{-90,40},{-70,40},{-70,76},{
            -22,76},{-22,80}},
        color={0,0,255},
        pattern=LinePattern.Dot));
    connect(blower.V, K_cath.V) annotation (Line(
        points={{-70,16},{-70,23}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(cathodeCooler.T_ref, K_cond.T_ref) annotation (Line(
        points={{42,37},{42,9},{41.2,9}},
        color={255,0,0},
        pattern=LinePattern.Dot));
    connect(K_cond.p_mix, mixer.p) annotation (Line(
        points={{26.8,12},{18,12},{4,12},{4,-44},{16,-44},{16,-74},{-4.7,-74},{
            -4.7,-66.9}},
        color={127,0,127},
        pattern=LinePattern.Dot));
    connect(degasser.T, K_fuel.T_deg) annotation (Line(
        points={{59,-20},{66,-20},{66,-98},{-30,-98},{-30,-92.4},{-17.2,-92.4}},
        color={255,0,0},
        pattern=LinePattern.Dot));
    connect(pump.V, K_an.V) annotation (Line(
        points={{-42,-60},{-59,-60}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(K_an.I, amperometer.i) annotation (Line(
        points={{-71,-60},{-90,-60},{-90,40},{-70,40},{-70,76},{-22,76},{-22,80}},
        color={0,0,255},
        pattern=LinePattern.Dot));
    connect(K_cath.V, K_cond.V_cath) annotation (Line(
        points={{-70,23},{-70,18},{16,18},{16,9},{26.8,9}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(K_fuel.V, fuelPump.V) annotation (Line(
        points={{-2.8,-90},{8,-90}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(K_temp.T_m, fuelCell.T) 
                               annotation (Line(
        points={{-17.2,-28},{-20,-28},{-20,-16},{-8,-16},{-8,5.34},{-12.2,5.34}},
        color={255,0,0},
        pattern=LinePattern.Dot));
    connect(anodeCooler.T_ref, K_temp.T_deg_ref) 
                                            annotation (Line(
        points={{20,-23},{20,-28},{-2.8,-28}},
        color={255,0,0},
        pattern=LinePattern.Dot));
    connect(cathodeCooler.T_process_out, K_cond.T_cond) annotation (Line(
        points={{51.4,38.6},{60,38},{60,0},{20,0},{20,6},{26.8,6}},
        color={255,0,0},
        pattern=LinePattern.Dot));

  end Reference_Control;

  partial model Mingled "A DMFC system with outlet mingling"

    import Modelica.SIunits.Efficiency;
    import Units.MolarFlow;

    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    inner parameter Units.RelativeHumidity RH_env = 60;

    Efficiency eta_to_cell "Fraction of methanol consumed in the cell";
    Efficiency eta_system "Overall system efficiency";

    Flow.UnitOperations.Mixer mixer 
                     annotation (Placement(transformation(extent={{-10,-70},{10,
              -50}}, rotation=0)));
    Flow.Measurements.LiquidPump pump "The anodic-loop pump" 
              annotation (Placement(transformation(extent={{-30,-66},{-42,-54}},
            rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),      graphics),
                         Documentation(info="<html>
<p>This is a generic reference system, with no process integration
whatsoever. Some components, such as the fuel cell, are abstract and
must be specialised in subclasses.</p>
</html>"));
    replaceable Flow.UnitOperations.Coolers.Abstract cooler "The system cooler"
                  annotation (Placement(transformation(extent={{14,-4},{32,14}},
            rotation=0)));
    Flow.Sources.Methanol pureMethanolSource "A substitute for an actual tank" 
                                          annotation (Placement(transformation(
            extent={{30,-96},{42,-84}}, rotation=0)));
    Flow.Measurements.LiquidPump fuelPump "The smaller fuel pump" 
                  annotation (Placement(transformation(extent={{20,-96},{8,-84}},
            rotation=0)));
    Flow.UnitOperations.Separator separator "The loop separator" 
                        annotation (Placement(transformation(extent={{48,-6},{
              68,16}}, rotation=0)));
    Flow.Sink co2sink "The gas outlet of the degasser" 
                      annotation (Placement(transformation(extent={{86,18},{94,
              26}}, rotation=0)));
    replaceable Flow.UnitOperations.Stack.Abstract fuelCell 
                      annotation (Placement(transformation(extent={{-50,-12},{
              -14,22}}, rotation=0)));
    Flow.Sources.Environment environment "The air from the environment" 
      annotation (Placement(transformation(extent={{-100,0},{-80,20}}, rotation=
             0)));
    Flow.Measurements.GasFlowController blower "The mass-flow controller" 
      annotation (Placement(transformation(
          origin={-70,10},
          extent={{6,-6},{-6,6}},
          rotation=270)));
    replaceable Modelica.Electrical.Analog.Interfaces.OnePort load
      "Load connected to the cell"       annotation (Placement(transformation(
            extent={{-52,84},{-40,96}}, rotation=0)));
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (Placement(transformation(extent={{-8,24},{8,40}}, rotation=0)));
    Modelica.Electrical.Analog.Sensors.CurrentSensor amperometer
      "Current in external circuit" annotation (Placement(transformation(extent=
             {{-32,80},{-12,100}}, rotation=0)));
  protected
    MolarFlow inCell = fuelCell.anode_inlet.n[1] - (-fuelCell.anode_outlet.n[1])
      "Methanol consumed in the cell";
    MolarFlow inDeg =  -separator.gasOutlet.n[1]
      "Methanol lost in the separator";

  public
    Flow.Measurements.MethanolInAir emissions
      "Checks whether there is too much methanol in the outlet" 
      annotation (Placement(transformation(extent={{68,16},{84,32}})));
  equation
    eta_to_cell = inCell / (inCell + inDeg);
    eta_system = eta_to_cell * fuelCell.eta_total;

    connect(pump.inlet, mixer.outlet) 
      annotation (Line(points={{-36,-60},{-8,-60}}, color={0,127,127}));
    connect(pureMethanolSource.outlet, fuelPump.inlet) 
      annotation (Line(points={{36,-90},{14,-90}}, color={0,127,127}));
    connect(environment.outlet, blower.inlet) 
                                         annotation (Line(points={{-81,10},{-70,
            10}}, color={0,127,127}));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (Line(points={{-64,10},
            {-50,10},{-50,10.1}},         color={0,127,127}));
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (Line(points={{14,-84},{0,-84},{0,-68},{6.10623e-16,-68}},
          color={0,127,127}));
    connect(cooler.outlet, separator.inlet) 
                                           annotation (Line(points={{31.46,5},{
            48,5}}, color={0,127,127}));
    connect(fuelCell.anode_outlet, cooler.inlet) annotation (Line(points={{-14,
            -0.1},{-14,0},{0,0},{0,5},{14.54,5}}, color={0,127,127}));
    connect(fuelCell.anode_inlet, pump.outlet) annotation (Line(points={{-50,
            -0.1},{-50,0},{-60,0},{-60,-54},{-36,-54}}, color={0,127,127}));
    connect(fuelCell.minus, ground.p) annotation (Line(points={{-21.2,15.2},{
            -21.2,40},{1.22125e-16,40}}, color={0,0,255}));
    connect(amperometer.p, load.n) 
      annotation (Line(points={{-32,90},{-40,90}}, color={0,0,255}));
    connect(fuelCell.minus, amperometer.n) annotation (Line(points={{-21.2,15.2},
            {-21.2,40},{0,40},{0,90},{-12,90}}, color={0,0,255}));
    connect(fuelCell.plus, load.p) annotation (Line(points={{-42.8,15.2},{-42.8,
            40},{-64,40},{-64,90},{-52,90}}, color={0,0,255}));
    connect(cooler.inlet, fuelCell.cathode_outlet) annotation (Line(points={{14.54,5},
            {0,5},{0,10.1},{-14,10.1}},          color={0,127,127}));
    connect(separator.liquidOutlet, mixer.waterInlet) annotation (Line(points={{65,0.6},
            {65,-60},{8,-60}},          color={0,127,127}));
    connect(separator.gasOutlet, emissions.inlet) annotation (Line(
        points={{65,9.4},{65,22.7},{70.4,22.7},{70.4,22.4}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(co2sink.inlet, emissions.outlet) annotation (Line(
        points={{86.4,22},{85.15,22},{85.15,22.4},{81.6,22.4}},
        color={0,127,127},
        smooth=Smooth.None));
  end Mingled;

  model Mingled_NoControl "The mingled system with manual control"
    extends Mingled(
      redeclare Flow.UnitOperations.Stack.Thevenin fuelCell,
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

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),      graphics),
                         experiment(StopTime=7200));
  end Mingled_NoControl;

  model Mingled_Control
    extends Mingled(
      redeclare Flow.UnitOperations.Coolers.Simple cooler,
      redeclare Flow.UnitOperations.Stack.Thevenin fuelCell(cells=3,
        V0=2.1,
        R=.15,
        k_x=2E-6,
        k_m_333=8E-6,
        T(start=303.15)),
      redeclare Load load(step(I=4, offset=3)),
      mixer(T(fixed=true), c(fixed=true)));
    Control.CathodeLambdaControl K_cath(
      cells=3,
      lambda=3,
      c_est=1100,
      aA=4.16E-9,
      b=0.2) 
      annotation (Placement(transformation(
          origin={-70,33},
          extent={{-5,-6},{5,6}},
          rotation=270)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),      graphics),
                          experiment(StopTime=10800));
    Control.WaterControl K_cond annotation (Placement(transformation(extent={{0,
              -32},{12,-20}}, rotation=0)));
    Control.MingledTemperatureControl K_T(
      n=3,
      c=900,
      aA=4.16E-9,
      b=0.2,
      T_r(displayUnit="degC"),
      V_cp(displayUnit="ml"))             annotation (Placement(transformation(
            extent={{-70,-78},{-54,-60}}, rotation=0)));
    Control.MingledFuelControl K_fuel(
      cells=3,
      lambda=3,
      aA=4.16E-9,
      b=0.2)                          annotation (Placement(transformation(
            extent={{-26,-98},{-12,-82}}, rotation=0)));

    model Load
      extends Modelica.Electrical.Analog.Interfaces.TwoPin;
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={Text(
              extent={{-100,-20},{100,20}},
              lineColor={0,0,255},
              textString="Load")}),
        Diagram(coordinateSystem(preserveAspectRatio=true,
                extent={{-100,-100},{100,100}}), graphics),
        experiment(StopTime=8000));
      Modelica.Electrical.Analog.Sources.StepCurrent step(
        I=2,
        offset=5,
        startTime(displayUnit="h") = 3600) 
        annotation (Placement(transformation(extent={{-20,20},{20,60}})));
      Modelica.Electrical.Analog.Sources.SineCurrent sine(
        I=2,
        freqHz=2E-3,
        startTime(displayUnit="h") = 7200) 
        annotation (Placement(transformation(extent={{-20,-60},{20,-20}})));
    equation
      connect(step.p, p) annotation (Line(
          points={{-20,40},{-60,40},{-60,5.55112e-16},{-100,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));

      connect(sine.p, p) annotation (Line(
          points={{-20,-40},{-60,-40},{-60,5.55112e-16},{-100,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(sine.n, n) annotation (Line(
          points={{20,-40},{60,-40},{60,0},{100,0},{100,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(step.n, n) annotation (Line(
          points={{20,40},{60,40},{60,5.55112e-16},{100,5.55112e-16}},
          color={0,0,255},
          smooth=Smooth.None));
    end Load;

  equation
    connect(K_cath.V, blower.V) annotation (Line(
        points={{-70,27},{-70,16}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(K_cath.I, amperometer.i) annotation (Line(
        points={{-70,39},{-70,74},{-22,74},{-22,80}},
        color={0,0,255},
        pattern=LinePattern.Dot));
    connect(K_cath.V, K_cond.V_cath) annotation (Line(
        points={{-70,27},{-70,24},{-56,24},{-56,-26},{-1.2,-26}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(cooler.T_ref, K_cond.T_ref) annotation (Line(
        points={{23,2.3},{23,-26},{13.2,-26}},
        color={255,0,0},
        pattern=LinePattern.Dot));
    connect(K_cond.p_mix, mixer.p) annotation (Line(
        points={{-1.2,-22.4},{-1.2,-22},{-12,-22},{-12,-66.9},{-4.7,-66.9}},
        color={128,0,255},
        pattern=LinePattern.Dot));
    connect(K_cond.T_cond, cooler.T_process_out) annotation (Line(
        points={{-1.2,-29.6},{-10,-30},{-10,-40},{40,-40},{40,4},{31.46,3.74}},
        color={255,0,0},
        pattern=LinePattern.Dot));

    connect(K_T.V, pump.V) annotation (Line(
        points={{-52.4,-69},{-48,-69},{-48,-60},{-42,-60}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(K_T.I, amperometer.i) annotation (Line(
        points={{-71.6,-74.4},{-72,-74},{-80,-74},{-80,74},{-22,74},{-22,80}},
        color={0,0,255},
        pattern=LinePattern.Dot));
    connect(K_T.T_stack, fuelCell.T) annotation (Line(
        points={{-71.6,-63.6},{-76,-64},{-76,-18},{-4,-18},{-4,5.34},{-12.2,
            5.34}},
        color={255,0,0},
        pattern=LinePattern.Dot));
    connect(K_T.T_mix, mixer.T) annotation (Line(
        points={{-71.6,-69},{-78,-69},{-78,-50},{-10,-50},{-9,-53}},
        color={255,0,0},
        pattern=LinePattern.Dot));
    connect(fuelPump.V, K_fuel.V) annotation (Line(
        points={{8,-90},{-10.6,-90}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(K_fuel.V_cath, K_cath.V) annotation (Line(
        points={{-27.4,-85.2},{-46,-86},{-46,-26},{-56,-26},{-56,24},{-70,24},{
            -70,27}},
        color={0,255,0},
        pattern=LinePattern.Dot));
    connect(K_fuel.I, amperometer.i) annotation (Line(
        points={{-27.4,-90},{-80,-90},{-80,74},{-22,74},{-22,80}},
        color={0,0,255},
        pattern=LinePattern.Dot));
    connect(K_fuel.T_sep, cooler.T_process_out) annotation (Line(
        points={{-27.4,-94.8},{-40,-94},{-40,-98},{40,-98},{40,4},{31.46,3.74}},
        color={255,0,0},
        pattern=LinePattern.Dot));

  end Mingled_Control;

  partial model DoubleTank "DMFC system with tanks for water and solution"

    import Modelica.SIunits.Efficiency;
    import Units.MolarFlow;

    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    inner parameter Units.RelativeHumidity RH_env = 60;

    Efficiency eta_to_cell "Fraction of methanol consumed in the cell";
    Efficiency eta_system "Overall system efficiency";

    Flow.UnitOperations.Mixer waterTank(c(start=0))
      "Tank containing the make-up water" 
                     annotation (Placement(transformation(extent={{30,-74},{44,
              -58}}, rotation=0)));
    Flow.Measurements.LiquidPump waterPump "The pure-water pump" 
              annotation (Placement(transformation(extent={{-14,-72},{-26,-60}},
            rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),      graphics),
                         Documentation(info="<html>
<p>This is a generic reference system, with no process integration
whatsoever. Some components, such as the fuel cell, are abstract and
must be specialised in subclasses.</p>
</html>"));
    replaceable Flow.UnitOperations.Coolers.Abstract anodeCooler
      "The solution-loop cooler" 
                  annotation (Placement(transformation(extent={{10,-30},{30,-10}},
            rotation=0)));
    Flow.Sources.Methanol pureMethanolSource "A substitute for an actual tank" 
                                          annotation (Placement(transformation(
            extent={{12,-92},{24,-80}}, rotation=0)));
    Flow.Measurements.LiquidPump fuelPump "The smaller fuel pump" 
                  annotation (Placement(transformation(extent={{-14,-92},{-26,
              -80}},
            rotation=0)));
    Flow.UnitOperations.Separator degasser "The CO2-degasser" 
                        annotation (Placement(transformation(extent={{38,-30},{
              58,-10}}, rotation=0)));
    Flow.Sink co2sink "The gas outlet of the degasser" 
                      annotation (Placement(transformation(extent={{68,-16},{76,
              -8}}, rotation=0)));
    replaceable Flow.UnitOperations.Stack.Abstract fuelCell 
                      annotation (Placement(transformation(extent={{-50,-12},{
              -14,22}}, rotation=0)));
    Flow.Sources.Environment environment "The air from the environment" 
      annotation (Placement(transformation(extent={{-100,0},{-80,20}}, rotation=
             0)));
    Flow.Measurements.GasFlowController blower "The mass-flow controller" 
      annotation (Placement(transformation(
          origin={-70,10},
          extent={{6,-6},{-6,6}},
          rotation=270)));
    replaceable Flow.UnitOperations.Coolers.Abstract cathodeCooler
      "The cathode-side cooler" 
                  annotation (Placement(transformation(extent={{32,30},{52,50}},
            rotation=0)));
    Flow.UnitOperations.Separator condenser "The water-recuperating unit" 
                        annotation (Placement(transformation(extent={{68,30},{
              86,50}}, rotation=0)));
    Flow.Sink airSink "The gas outlet of the condenser" 
                      annotation (Placement(transformation(extent={{92,46},{100,
              54}}, rotation=0)));
    replaceable Modelica.Electrical.Analog.Interfaces.TwoPin load
      "Load connected to the cell"       annotation (Placement(transformation(
            extent={{-52,84},{-40,96}}, rotation=0)));
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (Placement(transformation(extent={{-8,24},{8,40}}, rotation=0)));
    Modelica.Electrical.Analog.Sensors.CurrentSensor amperometer
      "Current in external circuit" annotation (Placement(transformation(extent=
             {{-32,80},{-12,100}}, rotation=0)));

  protected
    MolarFlow inCell = fuelCell.anode_inlet.n[1] - (-fuelCell.anode_outlet.n[1])
      "Methanol consumed in the cell";
    MolarFlow inDeg =  -degasser.gasOutlet.n[1] "Methanol lost in the degasser";

  public
    Flow.UnitOperations.Mixer solutionTank "Tank to gather the outlet solution"
                     annotation (Placement(transformation(extent={{10,-54},{24,
              -38}}, rotation=0)));
    Flow.Measurements.LiquidPump pump "The anodic-loop pump" 
              annotation (Placement(transformation(extent={{-14,-52},{-26,-40}},
            rotation=0)));
    Flow.Measurements.FlowConcentration FC6
      "Measures the concentration fed to the stack" 
      annotation (Placement(transformation(extent={{-64,-68},{-48,-52}})));
  equation
    eta_to_cell = inCell / (inCell + inDeg);
    eta_system = eta_to_cell * fuelCell.eta_total;
    connect(waterPump.inlet, waterTank.outlet) 
      annotation (Line(points={{-20,-66},{-20,-66},{31.4,-66}},
                                                    color={0,127,127}));
    connect(pureMethanolSource.outlet, fuelPump.inlet) 
      annotation (Line(points={{18,-86},{-2,-86},{-20,-86}},
                                                   color={0,127,127}));
    connect(environment.outlet, blower.inlet) 
                                         annotation (Line(points={{-81,10},{-70,
            10}}, color={0,127,127}));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (Line(points={{-64,10},
            {-50,10},{-50,10.1}},         color={0,127,127}));
    connect(cathodeCooler.outlet, condenser.inlet) 
      annotation (Line(points={{51.4,40},{68,40}}, color={0,127,127}));
    connect(anodeCooler.outlet, degasser.inlet) 
                                           annotation (Line(points={{29.4,-20},
            {38,-20}}, color={0,127,127}));
    connect(degasser.gasOutlet, co2sink.inlet)    annotation (Line(points={{55,
            -16},{62,-16},{62,-12},{68.4,-12}}, color={0,127,127}));
    connect(condenser.gasOutlet, airSink.inlet) 
      annotation (Line(points={{83.3,44},{88,44},{88,50},{92.4,50}}, color={0,
            127,127}));
    connect(fuelCell.anode_outlet, anodeCooler.inlet) 
                                                 annotation (Line(points={{-14,
            -0.1},{-14,0},{-4,0},{-4,-20},{10.6,-20}}, color={0,127,127}));
    connect(cathodeCooler.inlet, fuelCell.cathode_outlet) annotation (Line(
          points={{32.6,40},{12,40},{12,14},{-14,14},{-14,10.1}}, color={0,127,
            127}));
    connect(fuelCell.minus, ground.p) annotation (Line(points={{-21.2,15.2},{
            -21.2,40},{1.22125e-16,40}}, color={0,0,255}));
    connect(amperometer.p, load.n) 
      annotation (Line(points={{-32,90},{-40,90}}, color={0,0,255}));
    connect(fuelCell.minus, amperometer.n) annotation (Line(points={{-21.2,15.2},
            {-21.2,40},{0,40},{0,90},{-12,90}}, color={0,0,255}));
    connect(fuelCell.plus, load.p) annotation (Line(points={{-42.8,15.2},{-42.8,
            40},{-64,40},{-64,90},{-52,90}}, color={0,0,255}));
    connect(pump.inlet, solutionTank.outlet) 
      annotation (Line(points={{-20,-46},{-20,-46},{11.4,-46}},
                                                    color={0,127,127}));
    connect(solutionTank.waterInlet, degasser.liquidOutlet) annotation (Line(
        points={{22.6,-46},{56,-46},{56,-24},{55,-24}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(waterTank.waterInlet, condenser.liquidOutlet) annotation (Line(
        points={{42.6,-66},{84,-66},{84,36},{83.3,36}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(pump.outlet, FC6.outlet) annotation (Line(
        points={{-20,-40},{-34,-40},{-34,-60},{-49.6,-60}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(waterPump.outlet, FC6.outlet) annotation (Line(
        points={{-20,-60},{-49.6,-60}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(fuelPump.outlet, FC6.outlet) annotation (Line(
        points={{-20,-80},{-34,-80},{-34,-60},{-49.6,-60}},
        color={0,127,127},
        smooth=Smooth.None));
    connect(FC6.inlet, fuelCell.anode_inlet) annotation (Line(
        points={{-62.4,-60},{-80,-60},{-80,-0.1},{-50,-0.1}},
        color={0,127,127},
        smooth=Smooth.None));
  end DoubleTank;
end System;
