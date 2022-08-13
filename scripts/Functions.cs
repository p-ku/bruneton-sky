using System;
using System.Diagnostics;
using Godot;

using static Definitions;
using static Constants;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using static System.Math;
using MathNet.Numerics.LinearAlgebra;


public static class Func
{
  public static double clamp(double num, double min, double max)
  {
	return Math.Min(Math.Max(num, min), max);
  }
  public static double mod(double x, double y)
  {
	return x - y * Math.Floor(x / y);
  }
  public static double ClampCosine(double mu)
  {
	return clamp(mu, -1, 1);
  }
  public static double ClampDistance(double d)
  {
	return Math.Max(d, 0.0 * m);
  }
  public static double ClampRadius(in AtmosphereParameters atmosphere, double r)
  {
	return clamp(r, atmosphere.bottom_radius, atmosphere.top_radius);
  }
  public static double SmoothStep(double edge0, double edge1, double x)
  {
	double t = Math.Min(Math.Max((x - edge0) / (edge1 - edge0), 0.0), 1.0);
	return t * t * (3.0 - 2.0 * t);
  }
  public static double SafeSqrt(double area)
  {
	return Math.Sqrt(Math.Max(area, 0.0 * m2));
  }
  public static double DistanceToTopAtmosphereBoundary(in AtmosphereParameters atmosphere, double r, double mu)
  {
	Debug.Assert(r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	double discriminant = r * r * (mu * mu - 1) + atmosphere.top_radius * atmosphere.top_radius;
	return ClampDistance(-r * mu + SafeSqrt(discriminant));
  }
  public static double DistanceToBottomAtmosphereBoundary(in AtmosphereParameters atmosphere, double r, double mu)
  {
	Debug.Assert(r >= atmosphere.bottom_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	double discriminant = r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius;
	return ClampDistance(-r * mu - SafeSqrt(discriminant));
  }
  public static bool RayIntersectsGround(in AtmosphereParameters atmosphere, double r, double mu)
  {
	Debug.Assert(r >= atmosphere.bottom_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	return mu < 0.0 && r * r * (mu * mu - 1) + atmosphere.bottom_radius * atmosphere.bottom_radius >= 0.0 * m2;
  }
  public static double GetLayerDensity(in DensityProfileLayer layer, double altitude)
  {
	double density = layer.ExpTerm * Math.Exp(layer.ExpScale * altitude) +
	layer.LinearTerm * altitude + layer.ConstantTerm;
	return clamp(density, 0.0, 1);
  }

  public static double GetProfileDensity(in DensityProfile profile, double altitude)
  {
	return altitude < profile.layers[0].Width ?
	  GetLayerDensity(profile.layers[0], altitude) :
	  GetLayerDensity(profile.layers[1], altitude);
  }

  public static double ComputeOpticalLengthToTopAtmosphereBoundary(
   in AtmosphereParameters atmosphere, in DensityProfile profile,
   double r, double mu)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	// Number of intervals for the numerical integration.
	const int SAMPLE_COUNT = 500;
	// The integration step, i.e. the length of each integration interval.
	double dx = DistanceToTopAtmosphereBoundary(atmosphere, r, mu) / SAMPLE_COUNT;
	// Integration loop.
	double result = 0.0 * m;
	for (int i = 0; i <= SAMPLE_COUNT; ++i)
	{
	  double d_i = i * dx;
	  // Distance between the current sample point and the planet center.
	  double r_i = Math.Sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
	  // Number density at the current sample point (divided by the number density
	  // at the bottom of the atmosphere, yielding a dimensionless number).
	  double y_i = GetProfileDensity(profile, r_i - atmosphere.bottom_radius);
	  // Sample weight (from the trapezoidal rule).
	  double weight_i = i == 0.0 || i == SAMPLE_COUNT ? 0.5 : 1.0;
	  result += y_i * weight_i * dx;
	}
	return result;
  }
  public static Vector<double> ComputeTransmittanceToTopAtmosphereBoundary(
  in AtmosphereParameters atmosphere, double r, double mu)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1.0);

	Vector<double> result = Vector<double>.Build.Dense(3);

	result = (-(atmosphere.rayleigh_scattering * ComputeOpticalLengthToTopAtmosphereBoundary(
	 atmosphere, atmosphere.rayleigh_density, r, mu) + atmosphere.mie_extinction *
	 ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.mie_density, r, mu) +
	 atmosphere.absorption_extinction * ComputeOpticalLengthToTopAtmosphereBoundary(
	 atmosphere, atmosphere.absorption_density, r, mu))).PointwiseExp();
	return result;
  }
  public static double GetTextureCoordFromUnitRange(double x, int texture_size)
  {
	return 0.5 / texture_size + x * (1.0 - 1.0 / texture_size);
  }

  static double GetUnitRangeFromTextureCoord(double u, int texture_size)
  {
	return (u - 0.5 / texture_size) / (1.0 - 1.0 / texture_size);
  }
  static Vector<double> GetTransmittanceTextureUvFromRMu(in AtmosphereParameters atmosphere,
  double r, double mu)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1.0);
	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	double H = Math.Sqrt(atmosphere.top_radius * atmosphere.top_radius -
	  atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon.
	double rho =
	  SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
	// and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
	double d = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
	double d_min = atmosphere.top_radius - r;
	double d_max = rho + H;
	double x_mu = (d - d_min) / (d_max - d_min);
	double x_r = rho / H;
	double[] resultArr = { GetTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_TEXTURE_WIDTH),
  GetTextureCoordFromUnitRange(x_r, TRANSMITTANCE_TEXTURE_HEIGHT)};

	return Vector<double>.Build.DenseOfArray(resultArr); ;
  }
  public static void GetRMuFromTransmittanceTextureUv(in AtmosphereParameters atmosphere,
  in Vector<double> uv, out double r, out double mu)
  {
	Debug.Assert(uv[0] >= 0.0 && uv[0] <= 1);
	Debug.Assert(uv[1] >= 0.0 && uv[1] <= 1);
	double x_mu = GetUnitRangeFromTextureCoord(uv[0], TRANSMITTANCE_TEXTURE_WIDTH);
	double x_r = GetUnitRangeFromTextureCoord(uv[1], TRANSMITTANCE_TEXTURE_HEIGHT);
	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	double H = Math.Sqrt(atmosphere.top_radius * atmosphere.top_radius -
	  atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon, from which we can compute r:
	double rho = H * x_r;
	r = Math.Sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
	// and maximum values over all mu - obtained for (r,1) and (r,mu_horizon) -
	// from which we can recover mu:
	double d_min = atmosphere.top_radius - r;
	double d_max = rho + H;
	double d = d_min + x_mu * (d_max - d_min);
	mu = d == 0.0 * m ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
	mu = ClampCosine(mu);
  }
  public static Vector<double> ComputeTransmittanceToTopAtmosphereBoundaryTexture(
   in AtmosphereParameters atmosphere, in double[] frag_coord)
  {
	double r, mu;
	double[] uvArr = { frag_coord[0] / TRANSMITTANCE_TEXTURE_WIDTH, frag_coord[1] / TRANSMITTANCE_TEXTURE_HEIGHT };
	GetRMuFromTransmittanceTextureUv(atmosphere, Vector<double>.Build.DenseOfArray(uvArr), out r, out mu);
	return ComputeTransmittanceToTopAtmosphereBoundary(atmosphere, r, mu);
  }

  public static Vector<double> GetTransmittanceToTopAtmosphereBoundary(
   in AtmosphereParameters atmosphere, in Vector<double>[,] transmittance_texture, double r, double mu)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Vector<double> uv = GetTransmittanceTextureUvFromRMu(atmosphere, r, mu).PointwiseRound();

	return transmittance_texture[(int)uv[0], (int)uv[1]];
  }
  public static Vector<double> GetTransmittance(in AtmosphereParameters atmosphere,
  in Vector<double>[,] transmittance_texture, double r, double mu, double d, bool ray_r_mu_intersects_ground)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1.0);
	Debug.Assert(d >= 0.0 * m);

	double r_d = ClampRadius(atmosphere, Math.Sqrt(d * d + 2.0 * r * mu * d + r * r));
	double mu_d = ClampCosine((r * mu + d) / r_d);

	if (ray_r_mu_intersects_ground)
	{
	  Vector<double> result = GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, r_d, -mu_d) /
	   GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, r, -mu);

	  return result.PointwiseMinimum(Vector<double>.One);
	}
	else
	{
	  Vector<double> result = GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, r, mu) /
	   GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, r_d, mu_d);

	  return result.PointwiseMinimum(Vector<double>.One);
	}
  }
  static Vector<double> GetTransmittanceToSun(
  in AtmosphereParameters atmosphere,
  in Vector<double>[,] transmittance_texture,
  double r, double mu_s)
  {
	double sin_theta_h = atmosphere.bottom_radius / r;
	double cos_theta_h = -Math.Sqrt(Math.Max(1.0 - sin_theta_h * sin_theta_h, 0.0));
	Vector<double> vecIn = GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, r, mu_s);

	return vecIn * SmoothStep(-sin_theta_h * atmosphere.sun_angular_radius / rad,
	 sin_theta_h * atmosphere.sun_angular_radius / rad, mu_s - cos_theta_h);
  }

  public static void ComputeSingleScatteringIntegrand(in AtmosphereParameters atmosphere,
   in Vector<double>[,] transmittance_texture, double r, double mu, double mu_s, double nu, double d,
   bool ray_r_mu_intersects_ground, out Vector<double> rayleigh, out Vector<double> mie)
  {
	double r_d = ClampRadius(atmosphere, Math.Sqrt(d * d + 2.0 * r * mu * d + r * r));
	double mu_s_d = ClampCosine((r * mu_s + d * nu) / r_d);
	Vector<double> sunTran = GetTransmittanceToSun(atmosphere, transmittance_texture, r_d, mu_s_d);
	Vector<double> tran = GetTransmittance(atmosphere, transmittance_texture, r, mu, d, ray_r_mu_intersects_ground);

	Vector<double> transmittance = tran.PointwiseMultiply(sunTran);
	rayleigh = transmittance * GetProfileDensity(atmosphere.rayleigh_density, r_d - atmosphere.bottom_radius);
	mie = transmittance * GetProfileDensity(atmosphere.mie_density, r_d - atmosphere.bottom_radius);
  }

  public static double DistanceToNearestAtmosphereBoundary(in AtmosphereParameters atmosphere,
   double r, double mu, bool ray_r_mu_intersects_ground)
  {
	if (ray_r_mu_intersects_ground)
	{
	  return DistanceToBottomAtmosphereBoundary(atmosphere, r, mu);
	}
	else
	{
	  return DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
	}
  }
  public static void ComputeSingleScattering(in AtmosphereParameters atmosphere,
   in Vector<double>[,] transmittance_texture, double r, double mu, double mu_s, double nu,
   bool ray_r_mu_intersects_ground, out Vector<double> rayleigh, out Vector<double> mie)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	Debug.Assert(mu_s >= -1.0 && mu_s <= 1);
	Debug.Assert(nu >= -1.0 && nu <= 1);

	// Number of intervals for the numerical integration.
	const int SAMPLE_COUNT = 50;
	// The integration step, i.e. the length of each integration interval.
	double dx =
	  DistanceToNearestAtmosphereBoundary(atmosphere, r, mu, ray_r_mu_intersects_ground) / SAMPLE_COUNT;
	// Integration loop.
	Vector<double> rayleigh_sum = Vector<double>.Build.Dense(3);
	Vector<double> mie_sum = Vector<double>.Build.Dense(3);
	for (int i = 0; i <= SAMPLE_COUNT; ++i)
	{
	  double d_i = i * dx;
	  // The Rayleigh and Mie single scattering at the current sample point.
	  Vector<double> rayleigh_i = Vector<double>.Build.Dense(3);
	  Vector<double> mie_i = Vector<double>.Build.Dense(3);
	  ComputeSingleScatteringIntegrand(atmosphere, transmittance_texture,
	  r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, out rayleigh_i, out mie_i);
	  // Sample weight (from the trapezoidal rule).
	  double weight_i = (i == 0.0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;

	  rayleigh_sum = rayleigh_sum + (rayleigh_i * weight_i);
	  mie_sum += mie_i * weight_i;
	}
	rayleigh = rayleigh_sum * dx * atmosphere.solar_irradiance * atmosphere.rayleigh_scattering;
	mie = mie_sum * dx * atmosphere.solar_irradiance * atmosphere.mie_scattering;
  }
  public static double RayleighPhaseFunction(double nu)
  {
	double k = 3.0 / (16.0 * PI * sr);
	return k * (1.0 + nu * nu);
  }

  public static double MiePhaseFunction(double g, double nu)
  {
	double k = 3.0 / (8.0 * PI * sr) * (1.0 - g * g) / (2.0 + g * g);
	return k * (1.0 + nu * nu) / Math.Pow(1.0 + g * g - 2.0 * g * nu, 1.5);
  }
  public static Vector<double> GetScatteringTextureUvwzFromRMuMuSNu(in AtmosphereParameters atmosphere,
   double r, double mu, double mu_s, double nu, bool ray_r_mu_intersects_ground)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	Debug.Assert(mu_s >= -1.0 && mu_s <= 1);
	Debug.Assert(nu >= -1.0 && nu <= 1);

	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	double H = Math.Sqrt(atmosphere.top_radius * atmosphere.top_radius -
	  atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon.
	double rho =
	  SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
	double u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);

	// Discriminant of the quadratic equation for the intersections of the ray
	// (r,mu) with the ground (see RayIntersectsGround).
	double r_mu = r * mu;
	double discriminant = r_mu * r_mu - r * r + atmosphere.bottom_radius * atmosphere.bottom_radius;
	double u_mu;
	if (ray_r_mu_intersects_ground)
	{
	  // Distance to the ground for the ray (r,mu), and its minimum and maximum
	  // values over all mu - obtained for (r,-1) and (r,mu_horizon).
	  double dee = -r_mu - SafeSqrt(discriminant);
	  double dee_min = r - atmosphere.bottom_radius;
	  double dee_max = rho;
	  u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(dee_max == dee_min ? 0.0 :
	  (dee - dee_min) / (dee_max - dee_min), SCATTERING_TEXTURE_MU_SIZE / 2);
	}
	else
	{
	  // Distance to the top atmosphere boundary for the ray (r,mu), and its
	  // minimum and maximum values over all mu - obtained for (r,1) and
	  // (r,mu_horizon).
	  double dee = -r_mu + SafeSqrt(discriminant + H * H);
	  double dee_min = atmosphere.top_radius - r;
	  double dee_max = rho + H;
	  u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
	(dee - dee_min) / (dee_max - dee_min), SCATTERING_TEXTURE_MU_SIZE / 2);
	}

	double d = DistanceToTopAtmosphereBoundary(atmosphere, atmosphere.bottom_radius, mu_s);
	double d_min = atmosphere.top_radius - atmosphere.bottom_radius;
	double d_max = H;
	double a = (d - d_min) / (d_max - d_min);
	double D = DistanceToTopAtmosphereBoundary(atmosphere, atmosphere.bottom_radius, atmosphere.mu_s_min);
	double A = (D - d_min) / (d_max - d_min);
	// An ad-hoc function equal to0.0for mu_s = mu_s_min (because then d = D and
	// thus a = A), equal to 1.0for mu_s = 1.0(because then d = d_min and thus
	// a = 0), and with a large slope around mu_s =0.0, to get more texture 
	// samples near the horizon.
	double u_mu_s = GetTextureCoordFromUnitRange(
	 Math.Max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_SIZE_S);

	double u_nu = (nu + 1.0) / 2.0;

	return Vector<double>.Build.DenseOfArray(new double[] { u_nu, u_mu_s, u_mu, u_r });
  }
  public static void GetRMuMuSNuFromScatteringTextureUvwz(in AtmosphereParameters atmosphere,
  in Vector<double> uvwz, out double r, out double mu, out double mu_s,
  out double nu, out bool ray_r_mu_intersects_ground)
  {
	Debug.Assert(uvwz[0] >= 0.0 && uvwz[0] <= 1.0);
	Debug.Assert(uvwz[1] >= 0.0 && uvwz[1] <= 1.0);
	Debug.Assert(uvwz[2] >= 0.0 && uvwz[2] <= 1.0);
	Debug.Assert(uvwz[3] >= 0.0 && uvwz[3] <= 1.0);

	// Distance to top atmosphere boundary for a horizontal ray at ground level.
	double H = Math.Sqrt(atmosphere.top_radius * atmosphere.top_radius -
	  atmosphere.bottom_radius * atmosphere.bottom_radius);
	// Distance to the horizon.
	double rho = H * GetUnitRangeFromTextureCoord(uvwz[2], SCATTERING_TEXTURE_R_SIZE);
	r = Math.Sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);

	if (uvwz[2] < 0.5)
	{
	  // Distance to the ground for the ray (r,mu), and its minimum and maximum
	  // values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
	  // we can recover mu:
	  double dee_min = r - atmosphere.bottom_radius;
	  double dee_max = rho;
	  double dee = dee_min + (dee_max - dee_min) * GetUnitRangeFromTextureCoord(
	  1.0 - 2.0 * uvwz[3], SCATTERING_TEXTURE_MU_SIZE / 2);
	  mu = dee == 0.0 * m ? -1.0 : ClampCosine(-(rho * rho + dee * dee) / (2.0 * r * dee));
	  ray_r_mu_intersects_ground = true;
	}
	else
	{
	  // Distance to the top atmosphere boundary for the ray (r,mu), and its
	  // minimum and maximum values over all mu - obtained for (r,1) and
	  // (r,mu_horizon) - from which we can recover mu:
	  double dee_min = atmosphere.top_radius - r;
	  double dee_max = rho + H;
	  double dee = dee_min + (dee_max - dee_min) * GetUnitRangeFromTextureCoord(
	  2.0 * uvwz[2] - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
	  mu = dee == 0.0 * m ? 1.0 : ClampCosine((H * H - rho * rho - dee * dee) / (2.0 * r * dee));
	  ray_r_mu_intersects_ground = false;
	}

	double x_mu_s = GetUnitRangeFromTextureCoord(uvwz[3], SCATTERING_TEXTURE_MU_SIZE_S);
	double d_min = atmosphere.top_radius - atmosphere.bottom_radius;
	double d_max = H;
	double D = DistanceToTopAtmosphereBoundary(atmosphere, atmosphere.bottom_radius, atmosphere.mu_s_min);
	double A = (D - d_min) / (d_max - d_min);
	double a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
	double d = d_min + Math.Min(a, A) * (d_max - d_min);
	mu_s = d == 0.0 * m ? 1.0 : ClampCosine((H * H - d * d) / (2.0 * atmosphere.bottom_radius * d));

	nu = ClampCosine(uvwz[2] * 2.0 - 1);
  }
  public static void GetRMuMuSNuFromScatteringTextureFragCoord(
in AtmosphereParameters atmosphere, in double[] frag_coord,
out double r, out double mu, out double mu_s, out double nu, out bool ray_r_mu_intersects_ground)
  {
	double frag_coord_nu = Math.Floor(frag_coord[0] / SCATTERING_TEXTURE_MU_SIZE_S);
	double frag_coord_mu_s = mod(frag_coord[0], SCATTERING_TEXTURE_MU_SIZE_S);
	double[] SCATTERING_TEXTURE_SIZE = {
  SCATTERING_TEXTURE_NU_SIZE - 1,
  SCATTERING_TEXTURE_MU_SIZE_S,
  SCATTERING_TEXTURE_MU_SIZE,
  SCATTERING_TEXTURE_R_SIZE};

	double[] uvwzArr = { frag_coord_nu, frag_coord_mu_s, frag_coord[0], frag_coord[2] };
	Vector<double> uvwz = Vector<double>.Build.DenseOfArray(uvwzArr) /
	 Vector<double>.Build.DenseOfArray(SCATTERING_TEXTURE_SIZE);

	GetRMuMuSNuFromScatteringTextureUvwz(
	  atmosphere, uvwz, out r, out mu, out mu_s, out nu, out ray_r_mu_intersects_ground);
	// Clamp nu to its valid range of values, given mu and mu_s.
	nu = clamp(nu, mu * mu_s - Math.Sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
	  mu * mu_s + Math.Sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
  }
  public static void ComputeSingleScatteringTexture(in AtmosphereParameters atmosphere,
   in Vector<double>[,] transmittance_texture, in double[] frag_coord,
   out Vector<double> rayleigh, out Vector<double> mie)
  {
	double r;
	double mu;
	double mu_s;
	double nu;
	bool ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
	   out r, out mu, out mu_s, out nu, out ray_r_mu_intersects_ground);
	ComputeSingleScattering(atmosphere, transmittance_texture,
	  r, mu, mu_s, nu, ray_r_mu_intersects_ground, out rayleigh, out mie);
  }
  //  TEMPLATE(AbstractSpectrum)
  public static Vector<double> GetScatteringT(in AtmosphereParameters atmosphere,
in Vector<double>[,,] scattering_texture, double r, double mu, double mu_s, double nu,
  bool ray_r_mu_intersects_ground)
  {
	Vector<double> uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
	  atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	double tex_coord_x = uvwz[0] * SCATTERING_TEXTURE_NU_SIZE - 1.0;
	double tex_x = Math.Floor(tex_coord_x);
	double lerp = tex_coord_x - tex_x;
	double[] uvw0 = { (tex_x + uvwz[1]) / SCATTERING_TEXTURE_NU_SIZE, uvwz[3], uvwz[2] };
	double[] uvw1 = { (tex_x + 1.0 + uvwz[1]) / SCATTERING_TEXTURE_NU_SIZE, uvwz[3], uvwz[2] };

	Vector<double> result0 = scattering_texture[
	  (int)Math.Round(uvw0[0]), (int)Math.Round(uvw0[1]), (int)Math.Round(uvw0[2])];
	Vector<double> result1 = scattering_texture[
	  (int)Math.Round(uvw1[0]), (int)Math.Round(uvw1[1]), (int)Math.Round(uvw1[2])];

	return result0 * (1.0 - lerp) + result1 * lerp;
  }

  public static Vector<double> GetScattering(
   in AtmosphereParameters atmosphere,
   in Vector<double>[,,] single_rayleigh_scattering_texture,
   in Vector<double>[,,] single_mie_scattering_texture,
   in Vector<double>[,,] multiple_scattering_texture,
   double r, double mu, double mu_s, double nu,
   bool ray_r_mu_intersects_ground,
   int scattering_order)
  {
	if (scattering_order == 1)
	{
	  Vector<double> rayleigh = GetScatteringT(atmosphere, single_rayleigh_scattering_texture, r, mu, mu_s, nu,
	  ray_r_mu_intersects_ground);
	  Vector<double> mie = GetScatteringT(atmosphere, single_mie_scattering_texture, r, mu, mu_s, nu,
	   ray_r_mu_intersects_ground);

	  return rayleigh * RayleighPhaseFunction(nu) + mie * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
	}
	else
	{
	  return GetScatteringT(atmosphere, multiple_scattering_texture, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
	}
  }

  public static Vector<double> ComputeScatteringDensity(in AtmosphereParameters atmosphere,
  in Vector<double>[,] transmittance_texture, in Vector<double>[,,] single_rayleigh_scattering_texture,
  in Vector<double>[,,] single_mie_scattering_texture, in Vector<double>[,,] multiple_scattering_texture,
  in Vector<double>[,] irradiance_texture, double r, double mu, double mu_s, double nu, int scattering_order)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	Debug.Assert(mu_s >= -1.0 && mu_s <= 1);
	Debug.Assert(nu >= -1.0 && nu <= 1);
	Debug.Assert(scattering_order >= 2);

	// Compute unit direction vectors for the zenith, the view direction omega and
	// and the sun direction omega_s, such that the cosine of the view-zenith
	// angle is mu, the cosine of the sun-zenith angle is mu_s, and the cosine of
	// the view-sun angle is nu. The goal is to simplify computations below.
	Vector<double> zenith_direction = Vector<double>.Build.DenseOfArray(new double[] { 0.0, 0.0, 1.0 });
	Vector<double> omega = Vector<double>.Build.DenseOfArray(new double[] { Math.Sqrt(1.0 - mu * mu), 0.0, mu });
	double sun_dir_x = omega[0] == 0.0 ? 0.0 : (nu - mu * mu_s) / omega[0];
	double sun_dir_y = Math.Sqrt(Math.Max(1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0));
	Vector<double> omega_s = Vector<double>.Build.DenseOfArray(new double[] { sun_dir_x, sun_dir_y, mu_s });
	const int SAMPLE_COUNT = 16;
	const double dphi = pi / SAMPLE_COUNT;
	const double dtheta = pi / SAMPLE_COUNT;
	Vector<double> rayleigh_mie = Vector<double>.Build.Dense(3, 0.0 * watt_per_cubic_meter_per_sr_per_nm);

	// Nested loops for the integral over all the incident directions omega_i.
	for (int l = 0; l < SAMPLE_COUNT; ++l)
	{
	  double theta = (l) + 0.5 * dtheta;
	  double cos_theta = Math.Cos(theta);
	  double sin_theta = Math.Sin(theta);
	  bool ray_r_theta_intersects_ground =
	  RayIntersectsGround(atmosphere, r, cos_theta);

	  // The distance and transmittance to the ground only depend on theta, so we
	  // can compute them in the outer loop for efficiency.
	  double distance_to_ground = 0.0 * m;
	  Vector<double> transmittance_to_ground = Vector<double>.Build.Dense(3);
	  Vector<double> ground_albedo = Vector<double>.Build.Dense(3);
	  if (ray_r_theta_intersects_ground)
	  {
		distance_to_ground = DistanceToBottomAtmosphereBoundary(atmosphere, r, cos_theta);
		transmittance_to_ground = GetTransmittance(atmosphere, transmittance_texture, r, cos_theta,
		  distance_to_ground, true /* ray_intersects_ground */);
		ground_albedo = atmosphere.ground_albedo;
	  }

	  for (int m = 0; m < 2 * SAMPLE_COUNT; ++m)
	  {
		double phi = (m) + 0.5 * dphi;

		Vector<double> omega_i = Vector<double>.Build.DenseOfArray(new[] { Math.Cos(phi) * sin_theta, Math.Sin(phi) * sin_theta, cos_theta });

		double domega_i = (dtheta / rad) * (dphi / rad) * Math.Sin(theta) * sr;

		// The radiance L_i arriving from direction omega_i after n-1.0bounces is
		// the sum of a term given by the precomputed scattering texture for the
		// (n-1)-th order:
		double nu1 = omega_s.DotProduct(omega_i);
		Vector<double> incident_radiance = GetScattering(atmosphere, single_rayleigh_scattering_texture,
		single_mie_scattering_texture, multiple_scattering_texture, r, omega_i[2], mu_s, nu1,
		  ray_r_theta_intersects_ground, scattering_order - 1);

		// and of the contribution from the light paths with n-1.0bounces and whose
		// last bounce is on the ground. This contribution is the product of the
		// transmittance to the ground, the ground albedo, the ground BRDF, and
		// the irradiance received on the ground after n-2 bounces.

		Vector<double> ground_normal = (zenith_direction * r + omega_i * distance_to_ground).Normalize(2);
		Vector<double> ground_irradiance = GetIrradiance(
		  atmosphere, irradiance_texture, atmosphere.bottom_radius, ground_normal.DotProduct(omega_s));

		Vector<double> result = transmittance_to_ground.PointwiseMultiply(ground_albedo);
		incident_radiance += (result * (1.0 / (PI * sr))).PointwiseMultiply(ground_irradiance); ;
		// The radiance finally scattered from direction omega_i towards direction
		// -omega is the product of the incident radiance, the scattering
		// coefficient, and the phase function for directions omega and omega_i
		// (all this summed over all particle types, i.e. Rayleigh and Mie).
		double nu2 = omega.DotProduct(omega_i);
		double rayleigh_density = GetProfileDensity(atmosphere.rayleigh_density, r - atmosphere.bottom_radius);
		double mie_density = GetProfileDensity(atmosphere.mie_density, r - atmosphere.bottom_radius);

		rayleigh_mie += incident_radiance.PointwiseMultiply(atmosphere.rayleigh_scattering * rayleigh_density *
		RayleighPhaseFunction(nu2) + atmosphere.mie_scattering * mie_density *
		MiePhaseFunction(atmosphere.mie_phase_function_g, nu2)) * domega_i;
	  }
	}
	return rayleigh_mie;
  }
  public static Vector<double> ComputeMultipleScattering(in AtmosphereParameters atmosphere,
   in Vector<double>[,] transmittance_texture, in Vector<double>[,,] scattering_density_texture,
   double r, double mu, double mu_s, double nu, bool ray_r_mu_intersects_ground)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu >= -1.0 && mu <= 1);
	Debug.Assert(mu_s >= -1.0 && mu_s <= 1);
	Debug.Assert(nu >= -1.0 && nu <= 1);

	// Number of intervals for the numerical integration.
	const int SAMPLE_COUNT = 50;
	// The integration step, i.e. the length of each integration interval.
	double dx = DistanceToNearestAtmosphereBoundary(atmosphere, r, mu, ray_r_mu_intersects_ground) / SAMPLE_COUNT;
	// Integration loop.
	Vector<double> rayleigh_mie_sum = Vector<double>.Build.Dense(3, 0.0 * watt_per_square_meter_per_sr_per_nm);

	for (int i = 0; i <= SAMPLE_COUNT; ++i)
	{
	  double d_i = (double)i * dx;
	  // The r, mu and mu_s parameters at the current integration point (see the
	  // single scattering section for a detailed explanation).
	  double r_i = ClampRadius(atmosphere, Math.Sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r));
	  double mu_i = ClampCosine((r * mu + d_i) / r_i);
	  double mu_s_i = ClampCosine((r * mu_s + d_i * nu) / r_i);

	  // The Rayleigh and Mie multiple scattering at the current sample point.
	  Vector<double> scat = GetScatteringT(
	  atmosphere, scattering_density_texture, r_i, mu_i, mu_s_i, nu, ray_r_mu_intersects_ground);

	  Vector<double> tran = GetTransmittance(
	atmosphere, transmittance_texture, r, mu, d_i, ray_r_mu_intersects_ground);


	  Vector<double> rayleigh_mie_i = scat.PointwiseMultiply(tran) * dx;

	  // Sample weight (from the trapezoidal rule).
	  double weight_i = (i == 0.0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;

	  rayleigh_mie_sum += rayleigh_mie_i * weight_i;
	}
	return rayleigh_mie_sum;
  }
  public static Vector<double> ComputeScatteringDensityTexture(in AtmosphereParameters atmosphere,
  in Vector<double>[,] transmittance_texture, in Vector<double>[,,] single_rayleigh_scattering_texture,
  in Vector<double>[,,] single_mie_scattering_texture, in Vector<double>[,,] multiple_scattering_texture,
  in Vector<double>[,] irradiance_texture, in double[] frag_coord, int scattering_order)
  {
	double r, mu, mu_s, nu;
	bool ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
	   out r, out mu, out mu_s, out nu, out ray_r_mu_intersects_ground);
	return ComputeScatteringDensity(atmosphere, transmittance_texture,
	  single_rayleigh_scattering_texture, single_mie_scattering_texture,
	  multiple_scattering_texture, irradiance_texture, r, mu, mu_s, nu, scattering_order);
  }

  public static Vector<double> ComputeMultipleScatteringTexture(in AtmosphereParameters atmosphere,
   in Vector<double>[,] transmittance_texture, in Vector<double>[,,] scattering_density_texture,
   in double[] frag_coord, out double nu)
  {
	double r, mu, mu_s;

	bool ray_r_mu_intersects_ground;
	GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord,
	out r, out mu, out mu_s, out nu, out ray_r_mu_intersects_ground);
	return ComputeMultipleScattering(atmosphere, transmittance_texture,
	scattering_density_texture, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
  }
  public static Vector<double> ComputeDirectIrradiance(in AtmosphereParameters atmosphere,
  in Vector<double>[,] transmittance_texture, double r, double mu_s)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu_s >= -1.0 && mu_s <= 1.0);

	double alpha_s = atmosphere.sun_angular_radius / rad;
	// Approximate average of the cosine factor mu_s over the visible fraction of
	// the Sun disc.
	double average_cosine_factor = mu_s < -alpha_s ? 0.0 : (mu_s > alpha_s ? mu_s :
	  (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));
	Vector<double> tran = GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, r, mu_s);
	Vector<double> result = atmosphere.solar_irradiance.PointwiseMultiply(tran) * average_cosine_factor;
	return result;
  }
  public static Vector<double> ComputeIndirectIrradiance(in AtmosphereParameters atmosphere,
  in Vector<double>[,,] single_rayleigh_scattering_texture, in Vector<double>[,,] single_mie_scattering_texture,
  in Vector<double>[,,] multiple_scattering_texture, double r, double mu_s, int scattering_order)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu_s >= -1.0 && mu_s <= 1.0);
	Debug.Assert(scattering_order >= 1.0);

	const int SAMPLE_COUNT = 32;
	const double dphi = pi / SAMPLE_COUNT;
	const double dtheta = pi / SAMPLE_COUNT;

	Vector<double> result = Vector<double>.Build.Dense(3, 0.0 * watt_per_square_meter_per_nm);

	Vector<double> omega_s = Vector<double>.Build.DenseOfArray(new[] { Math.Sqrt(1.0 - mu_s * mu_s), 0.0, mu_s });
	for (int j = 0; j < SAMPLE_COUNT / 2; ++j)
	{
	  double theta = (j) + 0.5 * dtheta;
	  for (int i = 0; i < 2 * SAMPLE_COUNT; ++i)
	  {
		double phi = (i) + 0.5 * dphi;
		double[] omegaArr = { Math.Cos(phi) * Math.Sin(theta), Math.Sin(phi) * Math.Sin(theta), Math.Cos(theta) };
		Vector<double> omega = Vector<double>.Build.DenseOfArray(omegaArr);
		double domega = (dtheta / rad) * (dphi / rad) * Math.Sin(theta) * sr;

		double nu = omega.DotProduct(omega_s);
		Vector<double> scat = GetScattering(atmosphere, single_rayleigh_scattering_texture,
		  single_mie_scattering_texture, multiple_scattering_texture,
		r, omega[2], mu_s, nu, false /* ray_r_theta_intersects_ground */, scattering_order);

		result += scat * omega[2] * domega;
	  }
	}
	return result;
  }
  static Vector<double> GetIrradianceTextureUvFromRMuS(
  in AtmosphereParameters atmosphere, double r, double mu_s)
  {
	Debug.Assert(r >= atmosphere.bottom_radius && r <= atmosphere.top_radius);
	Debug.Assert(mu_s >= -1.0 && mu_s <= 1);
	double x_r = (r - atmosphere.bottom_radius) / (atmosphere.top_radius - atmosphere.bottom_radius);
	double x_mu_s = mu_s * 0.5 + 0.5;

	return Vector<double>.Build.DenseOfArray(new double[] {
  GetTextureCoordFromUnitRange(x_mu_s, IRRADIANCE_TEXTURE_WIDTH),
  GetTextureCoordFromUnitRange(x_r, IRRADIANCE_TEXTURE_HEIGHT)});

  }
  static void GetRMuSFromIrradianceTextureUv(in AtmosphereParameters atmosphere,
   in Vector<double> uv, out double r, out double mu_s)
  {
	Debug.Assert(uv[0] >= 0.0 && uv[0] <= 1);
	Debug.Assert(uv[1] >= 0.0 && uv[1] <= 1);
	double x_mu_s = GetUnitRangeFromTextureCoord(uv[0], IRRADIANCE_TEXTURE_WIDTH);
	double x_r = GetUnitRangeFromTextureCoord(uv[1], IRRADIANCE_TEXTURE_HEIGHT);
	r = atmosphere.bottom_radius + x_r * (atmosphere.top_radius - atmosphere.bottom_radius);
	mu_s = ClampCosine(2.0 * x_mu_s - 1);
  }

  public static Vector<double> ComputeDirectIrradianceTexture(in AtmosphereParameters atmosphere,
  in Vector<double>[,] transmittance_texture, in double[] frag_coord)
  {
	double r, mu_s;
	Vector<double> uv = Vector<double>.Build.DenseOfArray(new[] {frag_coord[0] / IRRADIANCE_TEXTURE_WIDTH,
  frag_coord[1] / IRRADIANCE_TEXTURE_HEIGHT});
	//uv[0] = frag_coord[0] / IRRADIANCE_TEXTURE_WIDTH;
	//uv[1] = frag_coord[1] / IRRADIANCE_TEXTURE_HEIGHT;

	GetRMuSFromIrradianceTextureUv(atmosphere, uv, out r, out mu_s);
	return ComputeDirectIrradiance(atmosphere, transmittance_texture, r, mu_s);
  }
  public static Vector<double> ComputeIndirectIrradianceTexture(in AtmosphereParameters atmosphere,
   in Vector<double>[,,] single_rayleigh_scattering_texture,
   in Vector<double>[,,] single_mie_scattering_texture, in Vector<double>[,,] multiple_scattering_texture,
   in double[] frag_coord, int scattering_order)
  {
	double r, mu_s;
	Vector<double> uv = Vector<double>.Build.DenseOfArray(new[] {frag_coord[0] / IRRADIANCE_TEXTURE_WIDTH,
  frag_coord[1] / IRRADIANCE_TEXTURE_HEIGHT});
	//	uv[0] = frag_coord[0] / IRRADIANCE_TEXTURE_WIDTH;
	//	uv[1] = frag_coord[1] / IRRADIANCE_TEXTURE_HEIGHT;

	GetRMuSFromIrradianceTextureUv(atmosphere, uv, out r, out mu_s);
	return ComputeIndirectIrradiance(atmosphere, single_rayleigh_scattering_texture,
	single_mie_scattering_texture, multiple_scattering_texture, r, mu_s, scattering_order);
  }
  public static Vector<double> GetIrradiance(in AtmosphereParameters atmosphere,
  in Vector<double>[,] irradiance_texture, double r, double mu_s)
  {
	Vector<double> uv = GetIrradianceTextureUvFromRMuS(atmosphere, r, mu_s).PointwiseRound();
	return irradiance_texture[(int)uv[0], (int)uv[1]]; ;
  }
}
