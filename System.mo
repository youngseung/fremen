package System "DMFC systems" 
  partial model Reference "The reference DMFC system, no control applied" 
    
    inner parameter Modelica.SIunits.Pressure p_env = 101325;
    inner parameter Modelica.SIunits.Temperature T_env = 298.15;
    inner parameter Real RH_env = 60;
    
    Flow.Mixer mixer annotation (extent=[-10,-42; 10,-22]);
    Flow.Pump pump "The anodic-loop pump" 
              annotation (extent=[-50,-38; -38,-26]);
    annotation (Diagram);
    Flow.Cooler anodeCooler "The solution-loop cooler" 
                  annotation (extent=[0,-10; 22,10]);
    Flow.PureMethanolSource pureMethanolSource 
      "A substitute for an actual tank"   annotation (extent=[30,-58; 42,-46]);
    Flow.Pump fuelPump "The smaller fuel pump" 
                  annotation (extent=[8,-58; 20,-46]);
    Flow.Separator degasser "The CO2-degasser" 
                        annotation (extent=[32,-10; 52,10]);
    Flow.SinkPort co2sink "The gas outlet of the degasser" 
                      annotation (extent=[54,20; 62,28]);
    replaceable Flow.FuelCell fuelCell 
                      annotation (extent=[-50,-12; -14,22]);
    Flow.EnvironmentPort environment "The air from the environment" 
      annotation (extent=[-82,22; -62,42]);
    Flow.GasFlowController blower "The mass-flow controller" 
      annotation (extent=[-74,2; -66,10]);
    Flow.Cooler cathodeCooler "The cathode-side cooler" 
                  annotation (extent=[0,0; 22,20]);
    Flow.Separator condenser "The water-recuperating unit" 
                        annotation (extent=[54,0; 72,20]);
    Flow.SinkPort airSink "The gas outlet of the condenser" 
                      annotation (extent=[74,28; 82,36]);
    replaceable Modelica.Electrical.Analog.Sources.ConstantCurrent I_cell(I=5) 
      "Current passing through the cell" annotation (extent=[-52,34; -32,54]);
    Modelica.Electrical.Analog.Basic.Ground ground 
      annotation (extent=[-32,-30; -18,-16]);
  equation 
    connect(pump.inlet, mixer.outlet) 
      annotation (points=[-44.12,-32; -8,-32], style(color=62, rgbcolor={0,127,
            127}));
    connect(pureMethanolSource.c, fuelPump.inlet) 
      annotation (points=[36,-52; 13.88,-52], style(color=62, rgbcolor={0,127,
            127}));
    connect(environment.c, blower.inlet) annotation (points=[-81,27; -89.5,27;
          -89.5,6; -70.08,6], style(color=62, rgbcolor={0,127,127}));
    connect(blower.outlet, fuelCell.cathode_inlet) annotation (points=[-70,10; 
          -50,10; -50,10.1],         style(color=62, rgbcolor={0,127,127}));
    connect(cathodeCooler.outlet, condenser.inlet) 
      annotation (points=[21.34,10; 54,10], style(color=62, rgbcolor={0,127,127}));
    connect(fuelPump.outlet, mixer.fuelInlet) 
      annotation (points=[14,-46; 0,-46; 0,-40; 6.10623e-16,-40],
                style(color=62, rgbcolor={0,127,127}));
    connect(anodeCooler.outlet, degasser.inlet) 
                                           annotation (points=[21.34,
          5.32907e-16; 27.005,5.32907e-16; 27.005,6.10623e-16; 32,6.10623e-16],
                                                     style(color=62, rgbcolor={
            0,127,127}));
    connect(degasser.gasOutlet, co2sink.flowPort) annotation (points=[49,4; 48,
          4; 48,24; 54.4,24],
        style(color=62, rgbcolor={0,127,127}));
    connect(condenser.gasOutlet, airSink.flowPort) 
      annotation (points=[69.3,14; 70,14; 70,32; 74.4,32],
                                           style(color=62, rgbcolor={0,127,127}));
    connect(condenser.liquidOutlet, mixer.waterInlet) 
      annotation (points=[69.3,6; 70,6; 70,-32; 8,-32], style(color=62,
          rgbcolor={0,127,127}));
    connect(degasser.liquidOutlet, mixer.loopInlet) annotation (points=[49,-4; 
          48,-4; 48,-14; 0,-14; 0,-24; 6.10623e-16,-24],    style(color=62,
          rgbcolor={0,127,127}));
    connect(fuelCell.anode_outlet, anodeCooler.inlet) 
                                                 annotation (points=[-14,-0.1; 
          -14,0; -4,0; -4,5.32907e-16; 0.66,5.32907e-16],
                                      style(
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
    connect(cathodeCooler.inlet, fuelCell.cathode_outlet) annotation (points=[0.66,10; 
          -14,10; -14,10.1],                  style(color=62, rgbcolor={0,127,
            127}));
    connect(I_cell.n, fuelCell.plus) annotation (points=[-32,44; -32,15.2],
        style(
        color=3,
        rgbcolor={0,0,255},
        pattern=0,
        fillColor=43,
        rgbfillColor={255,85,85},
        fillPattern=1));
    connect(I_cell.p, fuelCell.minus) annotation (points=[-52,44; -58,44; -58,
          -16; -32,-16; -32,-5.2], style(
        color=3,
        rgbcolor={0,0,255},
        pattern=0,
        fillColor=43,
        rgbfillColor={255,85,85},
        fillPattern=1));
    connect(ground.p, fuelCell.minus) annotation (points=[-25,-16; -32,-16; -32,
          -5.2], style(
        color=3,
        rgbcolor={0,0,255},
        pattern=0,
        fillColor=43,
        rgbfillColor={255,85,85},
        fillPattern=1));
  end Reference;
  
  model Reference_NoControl "The reference DMFC system, no control applied" 
    extends Reference(redeclare Flow.OvervoltageFuelCell fuelCell);
    
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    
    parameter VolumeFlowRate anodeFlow = 30e-6/60 
      "Flow rate for the anodic loop";
    parameter VolumeFlowRate fuelFlow = 30e-7/60 "Flow rate for the fuel";
    parameter VolumeFlowRate cathodeFlow = 30E-5/60 
      "Flow rate for the cathodic loop";
    parameter Temperature T_cooler_out = 350;
    parameter Temperature T_condenser = 320;
    
  equation 
    pump.V = anodeFlow;
    blower.V = cathodeFlow;
    fuelPump.V = fuelFlow;
    anodeCooler.outletTemperature.T = T_cooler_out;
    cathodeCooler.outletTemperature.T = T_condenser;
    
  end Reference_NoControl;
  annotation (uses(Modelica(version="2.2.1")));
  model Reference_SimpleControl 
    "The reference DMFC system, some simple control loops" 
    extends Reference(redeclare Flow.OvervoltageFuelCell fuelCell);
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Thermo.Oxygen;
    
    parameter VolumeFlowRate anodeFlow = 30e-6/60 
      "Flow rate for the anodic loop";
    parameter Real lambda = 2 "Air excess ratio";
    parameter Real k_c = 1 "P-control constant for concentration";
    parameter Temperature T_cooler_out = 340;
    parameter Temperature T_condenser = 320;
    
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    
  equation 
    // Air inflow controlled by simple lambda control.
    fuelCell.cathode_inlet.n[Oxygen] = lambda * fuelCell.I/4/F;
    
    // Fuel flow set to maintain initial concentration
    /* Note that this controller is totally bananas, not just because it
   * assumes we can perfectly read the mixer concentration, but also because
   * it case of excess concentration we could pull methanol out of it! */
    fuelPump.V = k_c * (mixer.c_MeOH_0 - fuelCell.cm_in.c);
    
    // Directly specified variables
    pump.V = anodeFlow;
    anodeCooler.outletTemperature.T = T_cooler_out;
    cathodeCooler.outletTemperature.T = T_condenser;
    
  end Reference_SimpleControl;
  
  model Reference_ASME "The reference DMFC system to be presented at ASME FC09" 
    extends Reference(redeclare Flow.ConstantVoltageFuelCell fuelCell);
    import Modelica.SIunits.VolumeFlowRate;
    import Modelica.SIunits.Temperature;
    import Modelica.SIunits.Concentration;
    import Modelica.SIunits.CurrentDensity;
    import Modelica.SIunits.Current;
    import Modelica.SIunits.Time;
    import Thermo.Oxygen;
    import Thermo.Water;
    
    parameter Real lambda_c = 2 "Air excess ratio";
    parameter Real lambda_a = 2 "Methanol solution excess ratio";
    parameter Concentration c_a_ref = 1000 
      "Target anodic methanol concentration";
    
    parameter Real a = 6*F*fuelCell.k_ad*b "Estimated dependence of i_x on c_a";
    parameter Real b = 1/(fuelCell.k_ad*fuelCell.d_M/fuelCell.D_M +1) 
      "Estimated dependence of i_x on i";
    CurrentDensity i_x_est = a*c_a_ref - b*fuelCell.i 
      "Estimate of the cross-over current density";
    CurrentDensity I_x_est = i_x_est*fuelCell.A 
      "Estimate of the cross-over current";
    
    parameter Time tau_w = 600 "Water control design time";
    Real f_M "Methanol degasser loss factor"; // TODO
    Real Kp "Proportionality constant for water content control";
    Real delta_n = sum(mixer.n) - mixer.n_0 "Difference of moles in the mixer";
    
  protected 
    constant Modelica.SIunits.FaradayConstant F = 96485.3415;
    
  equation 
    // Air inflow controlled by simple lambda control.
    fuelCell.cathode_inlet.n[Oxygen] = lambda_c * (fuelCell.I+I_x_est)/4/F;
    
    // Solution inflow, also controlled with lambda control.
    pump.V = lambda_a / mixer.c * (fuelCell.I+I_x_est)/6/F;
    
    // Fuel flow feedforward control TODO define f_M!
    fuelPump.F = (1-b)/6/F + a*fuelCell.A/6/F * c_a_ref + f_M * c_a_ref * fuelCell.I;
    f_M = 10E-11 * 10^((condenser.Tm.T-310)/42);
    
    // Water content feedback control
    condenser.Tm.T = 330 + Kp * delta_n;
    Kp = 1 / (tau_w * blower.F / p_env * Thermo.dp_h2o_dt(condenser.Tm.T, 1));
    
    // Constant temperature in the cell assumed to be maintained by misterious stranger.
    der(fuelCell.Tm.T) = 0;
    
  end Reference_ASME;
end System;
