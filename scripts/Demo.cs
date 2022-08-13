using System;
using System.Collections.Generic;
using System.Diagnostics;
using Godot;

using static Definitions;
using static Definitions.DensityProfile;
using static Definitions.DensityProfileLayer;
using static Constants;
using static Func;
using static Model;
using MathNet.Numerics.LinearAlgebra;


public class Demo : Node
{
  private const double kPi = Math.PI;
  private const double kSunAngularRadius = 0.00935 / 2.0;

  private const double kSunSolidAngle = kPi * kSunAngularRadius * kSunAngularRadius;

  private const double kLengthUnitInMeters = 1000.0;
  bool use_constant_solar_spectrum_ = false;
  bool use_ozone_ = true;
  bool use_combined_textures_ = true;
  static bool use_half_precision_ = true;
  //  use_luminance_(NONE),
  bool do_white_balance_ = false;
  bool show_help_ = true;
  int program_ = 0;
  double view_distance_meters_ = 9000.0;
  double view_zenith_angle_radians_ = 1.47;
  double view_azimuth_angle_radians_ = -0.1;
  double sun_zenith_angle_radians_ = 1.3;
  double sun_azimuth_angle_radians_ = 2.9;
  double exposure_ = 10.0;
  private const double kLambdaR = 680.0;
  private const double kLambdaG = 550.0;
  private const double kLambdaB = 440.0;

  private const int kLambdaMin = 360;
  private const int kLambdaMax = 830;

  double[] kSolarIrradiance = {
  1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
  1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
  1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
  1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
  1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
  1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
  };

  double[] kOzoneCrossSection = {
   1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
  8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
   1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
  4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
  2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
  6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
  2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
  };
  static double kDobsonUnit = 2.687e20;
  // Maximum number density of ozone molecules, in m^-3 (computed so at to get
  // 300 Dobson units of ozone - for this we divide 300 DU by the integral of
  // the ozone density profile defined below, which is equal to 15km).
  double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;
  // Wavelength independent solar irradiance "spectrum" (not physically
  // realistic, but was used in the original implementation).
  double kConstantSolarIrradiance = 1.5;
  double kBottomRadius = 6360000.0;
  double kTopRadius = 6420000.0;
  double kRayleigh = 1.24062e-6;
  static double kRayleighScaleHeight = 8000.0;
  static double kMieScaleHeight = 1200.0;
  double kMieAngstromAlpha = 0.0;
  double kMieAngstromBeta = 5.328e-3;
  double kMieSingleScatteringAlbedo = 0.9;
  double kMiePhaseFunctionG = 0.8;
  double kGroundAlbedo = 0.1;
  double max_sun_zenith_angle =
   (use_half_precision_ ? 102.0 / 180.0 * kPi : 120.0 / 180.0 * kPi);
  DensityProfileLayer rayleigh_layer = new DensityProfileLayer(0.0, 1.0, -1.0 / kRayleighScaleHeight, 0.0, 0.0);
  DensityProfileLayer mie_layer = new DensityProfileLayer(0.0, 1.0, -1.0 / kMieScaleHeight, 0.0, 0.0);
  DensityProfile ozone_density = new DensityProfile();
  DensityProfileLayer ozone_density1 = new DensityProfileLayer(25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0);

  DensityProfileLayer ozone_density2 = new DensityProfileLayer(0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0);
  public override void _Ready()
  {
	InitModel();
  }
  void InitModel()
  {
	List<double> wavelengths = new List<double>();
	List<double> solar_irradiance = new List<double>();
	List<double> rayleigh_scattering = new List<double>();
	List<double> mie_scattering = new List<double>();
	List<double> mie_extinction = new List<double>();
	List<double> absorption_extinction = new List<double>();
	List<double> ground_albedo = new List<double>();

	for (int l = kLambdaMin; l <= kLambdaMax; l += 10)
	{
	  double lambda = (double)l * 1.0e-3;  // micro-meters
	  double mie = kMieAngstromBeta / kMieScaleHeight * Math.Pow(lambda, -kMieAngstromAlpha);
	  wavelengths.Add(l);
	  //   wavelengths[count] = l;
	  if (use_constant_solar_spectrum_)
	  {
		//  solar_irradiance[count] = kConstantSolarIrradiance;
		solar_irradiance.Add(kConstantSolarIrradiance);
	  }
	  else
	  {
		//  solar_irradiance[count] = kSolarIrradiance[(l - kLambdaMin) / 10];
		solar_irradiance.Add(kSolarIrradiance[(l - kLambdaMin) / 10]);
	  }
	  rayleigh_scattering.Add(kRayleigh * Math.Pow(lambda, -4.0));
	  mie_scattering.Add(mie * kMieSingleScatteringAlbedo);
	  mie_extinction.Add(mie);
	  absorption_extinction.Add(use_ozone_ ?
	  kMaxOzoneNumberDensity * kOzoneCrossSection[(l - kLambdaMin) / 10] : 0.0);
	  ground_albedo.Add(kGroundAlbedo);
	  //   rayleigh_scattering[count] = kRayleigh * Math.Pow(lambda, -4.0);
	  //   mie_scattering[count] = mie * kMieSingleScatteringAlbedo;
	  //   mie_extinction[count] = mie;
	  //   absorption_extinction[count] = use_ozone_ ?
	  //   kMaxOzoneNumberDensity * kOzoneCrossSection[(l - kLambdaMin) / 10] : 0.0;
	  //   ground_albedo[count] = kGroundAlbedo;
	  //   count++;
	}
	//  double white_point_r = 1.0;
	//  double white_point_g = 1.0;
	//  double white_point_b = 1.0;
	// if (do_white_balance_)
	// {
	//   Model.ConvertSpectrumToLinearSrgb(out white_point_r, out white_point_g, out white_point_b);
	//   double white_point = (white_point_r + white_point_g + white_point_b) / 3.0;
	//   white_point_r /= white_point;
	//   white_point_g /= white_point;
	//   white_point_b /= white_point;
	// }

	ozone_density.layers = new DensityProfileLayer[2] { ozone_density1, ozone_density2 };


	Model myAtmo = new Model();
	myAtmo.Wavelengths = wavelengths;
	myAtmo.SolarIrradiance = Vector<double>.Build.DenseOfArray(solar_irradiance.ToArray());
	myAtmo.SunAngularRadius = kSunAngularRadius;
	myAtmo.BottomRadius = kBottomRadius;
	myAtmo.TopRadius = kTopRadius;
	myAtmo.RayleighDensity = rayleigh_layer;
	myAtmo.RayleighScattering = Vector<double>.Build.DenseOfArray(rayleigh_scattering.ToArray());
	myAtmo.MieDensity = mie_layer;
	myAtmo.MieScattering = Vector<double>.Build.DenseOfArray(mie_scattering.ToArray());
	myAtmo.MieExtinction = Vector<double>.Build.DenseOfArray(mie_extinction.ToArray());
	myAtmo.MiePhaseFunctionG = kMiePhaseFunctionG;
	myAtmo.AbsorptionDensity = ozone_density;
	myAtmo.AbsorptionExtinction = Vector<double>.Build.DenseOfArray(absorption_extinction.ToArray());
	myAtmo.GroundAlbedo = Vector<double>.Build.DenseOfArray(ground_albedo.ToArray());
	myAtmo.MaxSunZenithAngle = sun_zenith_angle_radians_;
	myAtmo.LengthUnitInMeters = kLengthUnitInMeters;
	myAtmo.UseLuminance = LUMINANCE.NONE;
	myAtmo.CombineScatteringTextures = true;
	myAtmo.HalfPrecision = true;
	// myAtmo.NumPrecomputedWavelengths = solar_irradiance.Count;
	myAtmo.Init(2);
	//  new Model(wavelengths, solar_irradiance, kSunAngularRadius,
	//        kBottomRadius, kTopRadius, { rayleigh_layer }, rayleigh_scattering,
	//    { mie_layer}, mie_scattering, mie_extinction, kMiePhaseFunctionG,
	//    ozone_density, absorption_extinction, ground_albedo, max_sun_zenith_angle,
	//    kLengthUnitInMeters, use_luminance_ == PRECOMPUTED ? 15 : 3,
	//    use_combined_textures_, use_half_precision_);
  }
}
