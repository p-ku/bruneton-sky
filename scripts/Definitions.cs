using System;

using System.Collections.Generic;
using static System.Math;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

public static class Definitions
{
  public const double m = 1;
  public const double nm = 1;
  public const double rad = 1;
  public const double sr = 1;
  public const double watt = 1;
  public const double lm = 1;
 // public const double PI = PI;
  public const double km = 1000 * m;
  public const double m2 = m * m;
  public const double m3 = m * m * m;
  public const double pi = PI * rad;
  public const double deg = pi / 180;
  public const double watt_per_square_meter = watt / m2;
  public const double watt_per_square_meter_per_sr = watt / (m2 * sr);
  public const double watt_per_square_meter_per_nm = watt / (m2 * nm);
  public const double watt_per_square_meter_per_sr_per_nm = watt / (m2 * sr * nm);
  public const double watt_per_cubic_meter_per_sr_per_nm = watt / (m3 * sr * nm);
  public const double cd = lm / sr;
  public const double kcd = 1000 * cd;
  public const double cd_per_square_meter = cd / m2;
  public const double kcd_per_square_meter = kcd / m2;

  public struct DensityProfileLayer
  {
    public double Width { get; private set; }
    public double ExpTerm { get; private set; }
    public double ExpScale { get; private set; }
    public double LinearTerm { get; private set; }
    public double ConstantTerm { get; private set; }
    public DensityProfileLayer(double width, double exp_term, double exp_scale, double linear_term, double constant_term)
    {
      Width = width;
      ExpTerm = exp_term;
      ExpScale = exp_scale;
      LinearTerm = linear_term;
      ConstantTerm = constant_term;
    }
  };

  public struct DensityProfile
  {
    public DensityProfileLayer[] layers;
  };
  public struct AtmosphereParameters
  {
    // The solar irradiance at the top of the atmosphere.
    public Vector<double> solar_irradiance;
    // The sun's angular radius. Warning: the implementation uses approximations
    // that are valid only if this angle is smaller than 0.1 radians.
    public double sun_angular_radius;
    // The distance between the planet center and the bottom of the atmosphere.
    public double bottom_radius;
    // The distance between the planet center and the top of the atmosphere.
    public double top_radius;
    // The density profile of air molecules, i.e. a function from altitude to
    // dimensionless values between 0 (null density) and 1 (maximum density).
    public DensityProfile rayleigh_density;
    // The scattering coefficient of air molecules at the altitude where their
    // density is maximum (usually the bottom of the atmosphere), as a function of
    // wavelength. The scattering coefficient at altitude h is equal to
    // 'rayleigh_scattering' times 'rayleigh_density' at this altitude.
    public Vector<double> rayleigh_scattering;
    // The density profile of aerosols, i.e. a function from altitude to
    // dimensionless values between 0 (null density) and 1 (maximum density).
    public DensityProfile mie_density;
    // The scattering coefficient of aerosols at the altitude where their density
    // is maximum (usually the bottom of the atmosphere), as a function of
    // wavelength. The scattering coefficient at altitude h is equal to
    // 'mie_scattering' times 'mie_density' at this altitude.
    public Vector<double> mie_scattering;
    // The extinction coefficient of aerosols at the altitude where their density
    // is maximum (usually the bottom of the atmosphere), as a function of
    // wavelength. The extinction coefficient at altitude h is equal to
    // 'mie_extinction' times 'mie_density' at this altitude.
    public Vector<double> mie_extinction;
    // The asymetry parameter for the Cornette-Shanks phase function for the
    // aerosols.
    public double mie_phase_function_g;
    // The density profile of air molecules that absorb light (e.g. ozone), i.e.
    // a function from altitude to dimensionless values between 0 (null density)
    // and 1 (maximum density).
    public DensityProfile absorption_density;
    // The extinction coefficient of molecules that absorb light (e.g. ozone) at
    // the altitude where their density is maximum, as a function of wavelength.
    // The extinction coefficient at altitude h is equal to
    // 'absorption_extinction' times 'absorption_density' at this altitude.
    public Vector<double> absorption_extinction;
    // The average albedo of the ground.
    public Vector<double> ground_albedo;
    // The cosine of the maximum Sun zenith angle for which atmospheric scattering
    // must be precomputed (for maximum precision, use the smallest Sun zenith
    // angle yielding negligible sky light radiance values. For instance, for the
    // Earth case, 102 degrees is a good choice - yielding mu_s_min = -0.2).
    public double mu_s_min;
  }

}
