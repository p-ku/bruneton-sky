using System;
using System.Collections.Generic;
using System.Diagnostics;
using Godot;

using static Definitions;
using static Definitions.DensityProfile;
using static Definitions.DensityProfileLayer;

using static Constants;
using static Func;


public class Main : Control
{

  internal double[] solar_irradiance;
  internal float sun_angular_radius;
  internal float bottom_radius;
  internal float top_radius;
  internal DensityProfile rayleigh_density;
  internal double[] rayleigh_scattering;
  internal DensityProfile mie_density;
  internal double[] mie_scattering;
  internal double[] mie_extinction;
  internal float mie_phase_function_g;
  internal DensityProfile absorption_density;
  internal double[] absorption_extinction;
  // The average albedo of the ground.
  internal double[] ground_albedo;
  internal float mu_s_min;
  // Declare member variables here. Examples:
  // private int a = 2;
  // private string b = "text";

  // Called when the node enters the scene tree for the first time.
  public override void _Ready()
  {
	AtmosphereParameters ATMO = new AtmosphereParameters();
	GD.Print(ATMO);
	//   for (int i = 0; i < TRANSMITTANCE_TEXTURE_HEIGHT; i++)
	//     for (int j = 0; j < TRANSMITTANCE_TEXTURE_WIDTH; j++)
	//     {
	//       double[] transmittance = Func.ComputeTransmittanceToTopAtmosphereBoundary(ATMO, i, j);
	//
	//     }

  }

  //  // Called every frame. 'delta' is the elapsed time since the previous frame.
  //  public override void _Process(float delta)
  //  {
  //      
  //  }
}
